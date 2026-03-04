#!/usr/bin/env python3
"""
Compute base-pair contacts mini-PDBs (from find_base_contacts_NHO.py) as QM candidates.

Filename formats recognised
  New : {ResName1}{Chain1}{Seq1}_{ResName2}{Chain2}{Seq2}_{PDB-ID}.pdb
        e.g.  GD35_CC26_4pco.pdb
  Old : {Seq1}_{Seq2}_{PDB-ID}.pdb
        e.g.  35_26_4pco.pdb

The script is fully self-contained: all geometry and contacts are computed
directly from the mini-PDB coordinates.  No external summary file is needed.

For every mini-PDB the script computes
  d_NO   shortest base-base N···O distance (both directions, no cutoff)
  C1C1   C1'–C1' distance
  tilt   angle between fitted base-plane normals  (degrees)
  H···O  all contacts within --cut  (H in NH-donor set, O in base/ribose)
  ∠      O–H···N angle at H  (degrees)  when donor heavy-atom exists in file

Output categories
─────────────────
GU pairs
  Canonical_Aform      H1–O2 < 2.2 Å  AND  10.3 ≤ C1C1 ≤ 11.2 Å  AND  tilt ≤ 20°
  Bifurcated_contact   H1–O2 present  AND  (H21–O2 or H22–O2 present)
  GU_other             remaining GU pairs with ≥1 contact

GC pairs  (Watson-Crick)
  WatsonCrick_GC       (H21–O2 or H22–O2) < 2.2 Å  AND  (H41–O6 or H42–O6) < 2.2 Å
  GC_other             remaining GC pairs with ≥1 contact

AU pairs  (Watson-Crick)
  WatsonCrick_AU       H61–O4 or H62–O4  < 2.2 Å
  AU_other             remaining AU pairs with ≥1 contact

AG / AC / UC / same-base pairs
  {type}_general       all pairs with ≥1 contact, sorted by d_NO

Usage
─────
  # GU only (backward-compatible):
  python score_base_pairs.py --pairs GU --cut 2.5 contacts/*.pdb

  # GU + GC + AU Watson-Crick:
  python score_base_pairs.py --pairs GU,GC,AU --cut 2.5 contacts/*.pdb

  # All recognised pair types found in directory:
  python score_base_pairs.py --cut 3.5 contacts/*.pdb
"""

import argparse
import glob
import math
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from Bio.PDB import PDBParser

# ══════════════════════════════════════════════════════════════════════════════
# Chemistry constants
# ══════════════════════════════════════════════════════════════════════════════

# Heavy base atoms (no H, no backbone, no ribose)
BASE_ATOMS: Dict[str, Set[str]] = {
    "G": {"N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"},
    "U": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"},
    "C": {"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"},
    "A": {"N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"},
}

# Ribose oxygens included as potential H-bond acceptors
RIBOSE_O: Set[str] = {"O2'", "O3'", "O4'", "O5'"}

# Complete acceptor-O set per residue (base oxygens + ribose)
ACCEPTOR_O: Dict[str, Set[str]] = {
    b: {a for a in atoms if a.startswith("O")} | RIBOSE_O
    for b, atoms in BASE_ATOMS.items()
}

# NH hydrogen name → donor heavy atom on the SAME residue
H_TO_DONOR: Dict[str, str] = {
    "H1":  "N1",   # G N1-H  (canonical wobble donor)
    "H3":  "N3",   # U N3-H  (wobble / WC donor)
    "H21": "N2",   # G N2-H21 (amino)
    "H22": "N2",   # G N2-H22 (amino)
    "H41": "N4",   # C N4-H41 (amino)
    "H42": "N4",   # C N4-H42 (amino)
    "H61": "N6",   # A N6-H61 (amino)
    "H62": "N6",   # A N6-H62 (amino)
    "H2":  "C2",   # A C2-H  (minor-groove C–H; donor listed as C2)
}

ALL_BASES: Tuple[str, ...] = ("G", "U", "C", "A")
ALL_BASES_SET: Set[str] = set(ALL_BASES)

# Canonical pair-code ordering (biochemically meaningful)
_CANON_ORDER: List[str] = ["GC", "GU", "AU", "AG", "AC", "UC", "GG", "AA", "UU", "CC"]
_CANON_SET:   Set[str]  = set(_CANON_ORDER)


def norm_pair(b1: str, b2: str) -> str:
    """Canonical pair code regardless of input order, e.g. norm_pair('U','G') → 'GU'."""
    k = b1 + b2
    if k in _CANON_SET:
        return k
    k2 = b2 + b1
    return k2 if k2 in _CANON_SET else k  # fallback: keep as-is


# ══════════════════════════════════════════════════════════════════════════════
# Data classes
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class Contact:
    """One H···O contact with optional O–H···N angle."""
    h_resname: str
    h_atom:    str
    o_resname: str
    o_atom:    str
    dist_HO:   float
    angle_OHN: Optional[float]   # degrees at H; None when donor heavy-atom absent

    def key(self) -> str:
        """Short identifier, e.g. 'H1-O2'."""
        return f"{self.h_atom}-{self.o_atom}"

    def short(self) -> str:
        """Compact representation for output columns."""
        ang = f",∠{self.angle_OHN:.0f}°" if self.angle_OHN is not None else ""
        return f"{self.h_atom}-{self.o_atom}:{self.dist_HO:.2f}{ang}"

    def verbose(self) -> str:
        """Full representation."""
        ang = (f"  ∠O-H-N={self.angle_OHN:.1f}°"
               if self.angle_OHN is not None else "  ∠=N/A")
        return (f"{self.h_resname}:{self.h_atom}···"
                f"{self.o_resname}:{self.o_atom}  "
                f"d={self.dist_HO:.3f} Å{ang}")


@dataclass
class PairRecord:
    """All geometry and contacts for one mini-PDB."""
    filepath:  str
    pdb_id:    str
    res1name:  str
    res2name:  str
    pair_type: str            # canonical code, e.g. 'GU'
    d_NO:      Optional[float]
    c1c1:      Optional[float]
    tilt:      Optional[float]
    contacts:  List[Contact] = field(default_factory=list)

    # ── Contact queries ────────────────────────────────────────────────────────

    def best_dist(self, h_atom: str, o_atom: str) -> Optional[float]:
        """Shortest H···O distance matching the given h/o atom names."""
        vals = [c.dist_HO for c in self.contacts
                if c.h_atom == h_atom and c.o_atom == o_atom]
        return min(vals) if vals else None

    def has_contact(self, h_atom: str, o_atom: str,
                    max_dist: float = 2.2) -> bool:
        return any(c.h_atom == h_atom and c.o_atom == o_atom
                   and c.dist_HO <= max_dist for c in self.contacts)

    def any_of(self, pairs: List[Tuple[str, str]],
               max_dist: float = 2.2) -> bool:
        """True if any (h_atom, o_atom) pair in list has a contact ≤ max_dist."""
        return any(self.has_contact(h, o, max_dist) for h, o in pairs)

    # ── Geometry predicates ────────────────────────────────────────────────────

    def is_aform(self) -> bool:
        return (self.c1c1 is not None and self.tilt is not None
                and 10.3 <= self.c1c1 <= 11.2 and self.tilt <= 20.0)

    # ── Formatting ─────────────────────────────────────────────────────────────

    def contacts_str(self) -> str:
        return ("  ".join(c.short() for c in
                          sorted(self.contacts, key=lambda x: x.dist_HO))
                or "None")

    def row(self) -> str:
        dno  = f"{self.d_NO:.3f}"  if self.d_NO  is not None else "NA"
        c1c1 = f"{self.c1c1:.3f}" if self.c1c1 is not None else "NA"
        tilt = f"{self.tilt:.1f}" if self.tilt  is not None else "NA"
        return (f"{self.filepath}\t"
                f"d_NO={dno}\tC1C1={c1c1}\ttilt={tilt}\t"
                f"contacts={self.contacts_str()}")


# ══════════════════════════════════════════════════════════════════════════════
# BioPython helpers
# ══════════════════════════════════════════════════════════════════════════════

def get_two_residues(structure):
    """
    Extract exactly two standard residues from a mini-PDB structure.
    Raises ValueError if fewer than two are found.
    """
    reslist = [
        res
        for model in structure
        for chain in model
        for res in chain
        if res.id[0] == " "      # skip HETATM / water
    ]
    if len(reslist) < 2:
        raise ValueError(f"Only {len(reslist)} standard residue(s) found; expected 2.")
    if len(reslist) > 2:
        # Take first 2 (extra residues from unusual files are ignored)
        reslist = reslist[:2]
    return reslist[0], reslist[1]


def atom_coords_array(res, names: Set[str]) -> np.ndarray:
    """Return (N,3) array of coordinates for named atoms present in residue."""
    pts = [res[n].get_coord() for n in names if res.has_id(n)]
    return np.array(pts, dtype=float) if pts else np.zeros((0, 3))


def fit_plane_normal(pts: np.ndarray) -> Optional[np.ndarray]:
    """Least-squares plane normal for a point cloud (returns unit vector or None)."""
    if pts.shape[0] < 3:
        return None
    c = pts.mean(axis=0)
    X = pts - c
    _, v = np.linalg.eigh((X.T @ X) / len(pts))
    n = v[:, 0]
    nrm = np.linalg.norm(n)
    return n / nrm if nrm > 0 else None


# ══════════════════════════════════════════════════════════════════════════════
# Geometry calculations
# ══════════════════════════════════════════════════════════════════════════════

def compute_d_NO(r1, r2, n1: str, n2: str) -> Optional[float]:
    """
    Shortest base-base N···O distance in either direction (no cutoff).
    Only atoms listed in BASE_ATOMS for each residue are considered.
    """
    b1 = BASE_ATOMS.get(n1, set())
    b2 = BASE_ATOMS.get(n2, set())
    best = float("inf")

    for a in r1:
        if a.element.upper() != "N" or a.get_name() not in b1:
            continue
        for b in r2:
            if b.element.upper() != "O" or b.get_name() not in b2:
                continue
            d = float(np.linalg.norm(a.get_coord() - b.get_coord()))
            if d < best:
                best = d

    for a in r2:
        if a.element.upper() != "N" or a.get_name() not in b2:
            continue
        for b in r1:
            if b.element.upper() != "O" or b.get_name() not in b1:
                continue
            d = float(np.linalg.norm(a.get_coord() - b.get_coord()))
            if d < best:
                best = d

    return best if best < float("inf") else None


def compute_c1c1(r1, r2) -> Optional[float]:
    """C1'–C1' distance."""
    if r1.has_id("C1'") and r2.has_id("C1'"):
        return float(np.linalg.norm(r1["C1'"].get_coord() - r2["C1'"].get_coord()))
    return None


def compute_tilt(r1, r2, n1: str, n2: str) -> Optional[float]:
    """Angle (°) between the two fitted base planes."""
    pts1 = atom_coords_array(r1, BASE_ATOMS.get(n1, set()))
    pts2 = atom_coords_array(r2, BASE_ATOMS.get(n2, set()))
    nv1  = fit_plane_normal(pts1)
    nv2  = fit_plane_normal(pts2)
    if nv1 is None or nv2 is None:
        return None
    dot = float(abs(np.dot(nv1, nv2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot))))


# ══════════════════════════════════════════════════════════════════════════════
# H···O contact finder + O–H···N angle
# ══════════════════════════════════════════════════════════════════════════════

def _ohn_angle(h_xyz: np.ndarray,
               o_xyz: np.ndarray,
               donor_xyz: np.ndarray) -> Optional[float]:
    """
    O–H···N angle at H (degrees).
    Vectors from H to O and from H to the donor heavy atom.
    A perfectly linear H-bond gives 180°.
    """
    vo = o_xyz     - h_xyz
    vn = donor_xyz - h_xyz
    do = np.linalg.norm(vo)
    dn = np.linalg.norm(vn)
    if do < 1e-6 or dn < 1e-6:
        return None
    cos_a = float(np.dot(vo, vn) / (do * dn))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_a))))


def find_HO_contacts(
    r1, r2,
    resname1: str, resname2: str,
    cutoff: float,
) -> List[Contact]:
    """
    Find all H···O contacts within *cutoff* Å (both directions).

    For each contact also compute the O–H···N (donor) angle at H.
    H atoms must be present in the PDB file (all-H or at least polar-H models).
    """
    contacts: List[Contact] = []
    acc1 = ACCEPTOR_O.get(resname1, set())
    acc2 = ACCEPTOR_O.get(resname2, set())

    def _search(h_res, o_res, h_rn: str, o_rn: str, o_names: Set[str]) -> None:
        for ah in h_res:
            hname = ah.get_name().strip()
            if ah.element.upper() != "H" or hname not in H_TO_DONOR:
                continue
            hxyz = ah.get_coord()

            for ao in o_res:
                oname = ao.get_name().strip()
                if ao.element.upper() != "O" or oname not in o_names:
                    continue
                d = float(np.linalg.norm(hxyz - ao.get_coord()))
                if d > cutoff:
                    continue

                # O–H···N angle: donor heavy atom is on the SAME residue as H
                dname = H_TO_DONOR[hname]
                angle = None
                if h_res.has_id(dname):
                    angle = _ohn_angle(hxyz, ao.get_coord(),
                                       h_res[dname].get_coord())

                contacts.append(
                    Contact(h_rn, hname, o_rn, oname, d, angle)
                )

    _search(r1, r2, resname1, resname2, acc2)
    _search(r2, r1, resname2, resname1, acc1)
    return contacts


# ══════════════════════════════════════════════════════════════════════════════
# Filename helpers
# ══════════════════════════════════════════════════════════════════════════════

def extract_pdb_id(filepath: str) -> str:
    """
    Extract the originating PDB ID from a mini-PDB filename.
    Works for both old (35_26_4pco.pdb) and new (GD35_CC26_4pco.pdb) formats.
    The PDB ID is always the last underscore-delimited token before '.pdb'.
    """
    stem  = os.path.basename(filepath)
    if stem.lower().endswith(".pdb"):
        stem = stem[:-4]
    parts = stem.split("_")
    return parts[-1] if parts else stem


# ══════════════════════════════════════════════════════════════════════════════
# Output helpers
# ══════════════════════════════════════════════════════════════════════════════

def _print_block(title: str, records: List[PairRecord]) -> None:
    """Print a scored category block.  Silently skips if records is empty."""
    if not records:
        return
    print(f"# {title}")
    for rec in records:
        print(f"{title}\t{rec.row()}")
    print()


# ══════════════════════════════════════════════════════════════════════════════
# Per-pair-type scoring functions
# ══════════════════════════════════════════════════════════════════════════════

def score_GU(records: List[PairRecord]) -> None:
    """
    GU-specific categories:
      Canonical_Aform     H1-O2 < 2.2 Å  + A-form geometry
      Bifurcated_contact  H1-O2 present  + (H21-O2 or H22-O2 present)
      GU_other            everything else with ≥1 contact
    """
    canonical:   List[PairRecord] = []
    bifurcated:  List[PairRecord] = []
    other:       List[PairRecord] = []

    for rec in records:
        is_canon = rec.has_contact("H1", "O2") and rec.is_aform()
        is_bifur = (rec.has_contact("H1", "O2", max_dist=9.9)    # any dist
                    and rec.any_of([("H21", "O2"), ("H22", "O2")],
                                   max_dist=9.9))

        if is_canon:
            canonical.append(rec)
        if is_bifur:
            bifurcated.append(rec)
        if not is_canon and not is_bifur:
            other.append(rec)

    # Sort canonical by closeness to ideal 1.8 Å H···O
    canonical.sort(
        key=lambda r: abs((r.best_dist("H1", "O2") or 9.9) - 1.8))
    # Sort bifurcated by H1-O2 distance
    bifurcated.sort(
        key=lambda r: r.best_dist("H1", "O2") or 9.9)
    other.sort(key=lambda r: r.d_NO or 9.9)

    _print_block("Canonical_Aform",    canonical)
    _print_block("Bifurcated_contact", bifurcated)
    _print_block("GU_other",           other)


def score_GC(records: List[PairRecord]) -> None:
    """
    GC Watson-Crick category:
      WatsonCrick_GC  (H21-O2 or H22-O2) < 2.2 Å  AND  (H41-O6 or H42-O6) < 2.2 Å
      GC_other        anything else with ≥1 contact
    """
    wc:    List[PairRecord] = []
    other: List[PairRecord] = []

    for rec in records:
        g_amino = rec.any_of([("H21", "O2"), ("H22", "O2")])
        c_amino = rec.any_of([("H41", "O6"), ("H42", "O6")])
        if g_amino and c_amino:
            wc.append(rec)
        else:
            other.append(rec)

    wc.sort(key=lambda r: min(
        r.best_dist("H21", "O2") or 9.9,
        r.best_dist("H22", "O2") or 9.9,
    ))
    other.sort(key=lambda r: r.d_NO or 9.9)

    _print_block("WatsonCrick_GC", wc)
    _print_block("GC_other",       other)


def score_AU(records: List[PairRecord]) -> None:
    """
    AU Watson-Crick category:
      WatsonCrick_AU  H61-O4 or H62-O4  < 2.2 Å
      AU_other        anything else with ≥1 contact
    """
    wc:    List[PairRecord] = []
    other: List[PairRecord] = []

    for rec in records:
        if rec.any_of([("H61", "O4"), ("H62", "O4")]):
            wc.append(rec)
        else:
            other.append(rec)

    wc.sort(key=lambda r: min(
        r.best_dist("H61", "O4") or 9.9,
        r.best_dist("H62", "O4") or 9.9,
    ))
    other.sort(key=lambda r: r.d_NO or 9.9)

    _print_block("WatsonCrick_AU", wc)
    _print_block("AU_other",       other)


def score_generic(pair_type: str, records: List[PairRecord]) -> None:
    """Fallback: single block sorted by d_NO."""
    records = sorted(records, key=lambda r: r.d_NO or 9.9)
    _print_block(f"{pair_type}_general", records)


# Dispatch table: pair_type → scorer function
_SCORERS = {
    "GU": score_GU,
    "GC": score_GC,
    "AU": score_AU,
}


# ══════════════════════════════════════════════════════════════════════════════
# --pairs argument parser
# ══════════════════════════════════════════════════════════════════════════════

def parse_pair_arg(raw: Optional[str]) -> Optional[Set[str]]:
    """
    Parse --pairs argument into a set of canonical pair codes.
    Returns None to mean 'all recognised pairs'.
    """
    if raw is None or raw.strip().upper() == "ALL":
        return None
    codes: Set[str] = set()
    for tok in raw.split(","):
        tok = tok.strip().upper()
        if len(tok) != 2:
            print(f"WARNING: ignoring unrecognised pair token '{tok}'",
                  file=sys.stderr)
            continue
        b1, b2 = tok[0], tok[1]
        if b1 not in ALL_BASES_SET or b2 not in ALL_BASES_SET:
            print(f"WARNING: '{tok}' contains an unknown base letter — skipping.",
                  file=sys.stderr)
            continue
        codes.add(norm_pair(b1, b2))
    return codes or None


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Score and rank base-pair mini-PDBs as QM model candidates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # GU pairs only (backward-compatible legacy mode):
  %(prog)s --pairs GU --cut 2.5 contacts/*.pdb > GU_scored.txt

  # GU + GC + AU Watson-Crick:
  %(prog)s --pairs GU,GC,AU --cut 2.5 contacts/*.pdb > WC_scored.txt

  # Every pair type found in directory (all-atom PDBs required for angles):
  %(prog)s --cut 3.5 contacts/*.pdb > all_pairs_scored.txt
""",
    )
    ap.add_argument(
        "pdbs", nargs="+",
        help="Mini-PDB files or glob patterns, e.g. contacts/*.pdb",
    )
    ap.add_argument(
        "--cut", type=float, default=3.5,
        help="H···O distance cutoff in Å for contact search (default: 3.5). "
             "Only pairs with ≥1 contact within this cutoff are reported.",
    )
    ap.add_argument(
        "--pairs", default=None,
        help=(
            "Comma-separated base-pair types to score, e.g. GU,GC,AU. "
            "Use ALL or omit to score every recognised pair found in the files. "
            "Recognised bases: G, U, C, A."
        ),
    )
    ap.add_argument(
        "--wc_cut", type=float, default=2.2,
        help="H···O distance threshold (Å) used for Watson-Crick / canonical "
             "category assignment (default: 2.2). "
             "Does not affect which contacts are searched — use --cut for that.",
    )
    ap.add_argument(
        "--no_angle", action="store_true",
        help="Skip O–H···N angle computation (faster; useful when H atoms absent).",
    )
    ap.add_argument(
        "--verbose", action="store_true",
        help="Print one contact per line with full annotation instead of compact form.",
    )

    args = ap.parse_args()

    # ── Validate arguments ────────────────────────────────────────────────────
    if args.cut <= 0:
        print("ERROR: --cut must be > 0", file=sys.stderr)
        sys.exit(2)
    if args.wc_cut <= 0:
        print("ERROR: --wc_cut must be > 0", file=sys.stderr)
        sys.exit(2)

    requested = parse_pair_arg(args.pairs)

    # ── Banner ────────────────────────────────────────────────────────────────
    print(f"# H···O search cutoff : {args.cut} Å", file=sys.stderr)
    print(f"# WC/canonical cutoff : {args.wc_cut} Å", file=sys.stderr)
    if requested:
        print(f"# Pair types requested: {', '.join(sorted(requested))}",
              file=sys.stderr)
    else:
        print("# Scoring all recognised pair types.", file=sys.stderr)

    # ── Collect PDB files ─────────────────────────────────────────────────────
    pdb_files: List[str] = []
    for pat in args.pdbs:
        matched = glob.glob(pat)
        pdb_files.extend(matched if matched else [pat])
    pdb_files = sorted(set(pdb_files))

    if not pdb_files:
        print("ERROR: no PDB files found.", file=sys.stderr)
        sys.exit(1)
    print(f"# Mini-PDB files      : {len(pdb_files)}", file=sys.stderr)

    # ── Parse mini-PDBs ───────────────────────────────────────────────────────
    bparser = PDBParser(QUIET=True)
    by_pair: Dict[str, List[PairRecord]] = {}

    for filepath in pdb_files:
        if not os.path.isfile(filepath):
            print(f"# SKIP (not found): {filepath}", file=sys.stderr)
            continue

        pdb_id = extract_pdb_id(filepath)

        try:
            structure = bparser.get_structure("pair", filepath)
            r1, r2   = get_two_residues(structure)
        except Exception as exc:
            print(f"# SKIP {filepath}: {exc}", file=sys.stderr)
            continue

        n1 = r1.get_resname().strip()
        n2 = r2.get_resname().strip()

        if n1 not in ALL_BASES_SET or n2 not in ALL_BASES_SET:
            print(f"# SKIP {filepath}: unrecognised resname(s) '{n1}'/'{n2}'",
                  file=sys.stderr)
            continue

        pt = norm_pair(n1, n2)
        if requested is not None and pt not in requested:
            continue    # filtered out by --pairs

        # ── Geometry ──────────────────────────────────────────────────────────
        d_NO = compute_d_NO(r1, r2, n1, n2)
        c1c1 = compute_c1c1(r1, r2)
        tilt = compute_tilt(r1, r2, n1, n2)

        # ── H···O contacts + angles ───────────────────────────────────────────
        contacts = find_HO_contacts(r1, r2, n1, n2, args.cut)
        if args.no_angle:
            # Zero out angles so downstream code still works
            contacts = [
                Contact(c.h_resname, c.h_atom, c.o_resname,
                        c.o_atom, c.dist_HO, None)
                for c in contacts
            ]

        if not contacts:
            continue    # no H···O contact within cutoff → not a candidate

        rec = PairRecord(
            filepath=filepath,
            pdb_id=pdb_id,
            res1name=n1,
            res2name=n2,
            pair_type=pt,
            d_NO=d_NO,
            c1c1=c1c1,
            tilt=tilt,
            contacts=contacts,
        )
        by_pair.setdefault(pt, []).append(rec)

    # ── Summary ───────────────────────────────────────────────────────────────
    total = sum(len(v) for v in by_pair.values())
    print(f"# Qualifying pairs    : {total}", file=sys.stderr)

    if not by_pair:
        print("# No qualifying pairs found — check --cut and --pairs.", file=sys.stderr)
        sys.exit(0)

    # ── Header ────────────────────────────────────────────────────────────────
    print("# " + "═" * 72)
    print("# Base-pair contact scoring")
    print(f"#   H···O search cutoff  = {args.cut} Å")
    print(f"#   WC/canonical cutoff  = {args.wc_cut} Å")
    print(f"#   Pair types found     = {', '.join(sorted(by_pair.keys()))}")
    print(f"#   Total pairs scored   = {total}")
    print("# " + "═" * 72)
    print()
    print("# Column format:")
    print("#   category  filepath  d_NO=<Å>  C1C1=<Å>  tilt=<°>  "
          "contacts=<H-O:dist,∠angle°  …>")
    print()

    # ── Score each pair type in canonical order ────────────────────────────────
    for pt in _CANON_ORDER:
        recs = by_pair.get(pt)
        if not recs:
            continue

        print("═" * 72)
        print(f"# {pt} pairs  (n = {len(recs)})")
        print("═" * 72)

        scorer = _SCORERS.get(pt)
        if scorer:
            scorer(recs)
        else:
            score_generic(pt, recs)

        print()

    # ── Optional verbose dump ─────────────────────────────────────────────────
    if args.verbose:
        print()
        print("# " + "─" * 72)
        print("# VERBOSE: per-contact detail for every qualifying pair")
        print("# " + "─" * 72)
        for pt in _CANON_ORDER:
            recs = by_pair.get(pt)
            if not recs:
                continue
            for rec in sorted(recs, key=lambda r: r.d_NO or 9.9):
                print(f"\n## {pt}  {rec.filepath}  "
                      f"d_NO={rec.d_NO:.3f if rec.d_NO else 'NA'}  "
                      f"C1C1={rec.c1c1:.3f if rec.c1c1 else 'NA'}  "
                      f"tilt={rec.tilt:.1f if rec.tilt else 'NA'}")
                for c in sorted(rec.contacts, key=lambda x: x.dist_HO):
                    print(f"    {c.verbose()}")


if __name__ == "__main__":
    main()
