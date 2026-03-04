#!/usr/bin/env python3
"""
Find residue pairs that have BASE–BASE (and sugar) H···O contacts within a cutoff.

Supported base types : G, U, C, A
Supported pair types : any combination, e.g. GU, GC, AU, AG, GG, AA …

Behavior:
- The distance cutoff (--cut) is applied ONLY to H···O contacts (NH···O / NH2···O).
- Pairs are selected based on the presence of at least one H···O contact within
  the cutoff (H attached to N: H1, H3, H21, H22, H41, H42, H61, H62, etc.).
- For each selected pair:
    * ALL N···O contacts (base–base and base–sugar) are reported with NO distance
      cutoff, in both directions (res1 N -> res2 O and res2 N -> res1 O).
    * H···O contacts within the cutoff are reported.
    * One mini-PDB per pair is written.

Outputs:
- --summary      : N···O contacts (unfiltered distances) for H···O-selected pairs.
- --summary_HO   : H···O (NH···O) contacts within the cutoff.
- mini-PDB files : two-residue PDBs for each selected pair.

Usage:
    # GU pairs only (backward-compatible default):
    ./find_base_contacts_NHO.py --pdbin 1msy.pdb --cut 3.5 \\
        --summary contacts/NO_summary.txt --summary_HO contacts/HO_summary.txt

    # GU + GC pairs:
    ./find_base_contacts_NHO.py --pdbin 1msy.pdb --cut 3.5 --pairs GU,GC \\
        --summary contacts/NO_summary.txt --summary_HO contacts/HO_summary.txt

    # All Watson-Crick + wobble pairs:
    ./find_base_contacts_NHO.py --pdbin 1msy.pdb --cut 3.5 --pairs GU,GC,AU,AG \\
        --summary contacts/NO_summary.txt --summary_HO contacts/HO_summary.txt

    # Every possible combination among G, U, C, A:
    ./find_base_contacts_NHO.py --pdbin 1msy.pdb --cut 3.5 --pairs ALL \\
        --summary contacts/NO_summary.txt --summary_HO contacts/HO_summary.txt
"""

import argparse
import itertools
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Base-atom sets (heavy atoms + ribose oxygens considered for H-bonding)
# ---------------------------------------------------------------------------
BASE_ATOMS: Dict[str, Set[str]] = {
    "G": {
        "N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6",
        "N7", "C8", "N9",
        "O2'", "O3'", "O4'", "O5'",
    },
    "U": {
        "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6",
        "O2'", "O3'", "O4'", "O5'",
    },
    "C": {
        "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6",
        "O2'", "O3'", "O4'", "O5'",
    },
    "A": {
        "N1", "C2", "N3", "C4", "C5", "C6", "N6",
        "N7", "C8", "N9",
        "O2'", "O3'", "O4'", "O5'",
    },
}

# Default NH hydrogen names per base (atoms attached to N that can donate H-bonds)
DEFAULT_NH_HNAMES: Dict[str, Set[str]] = {
    "G": {"H1", "H21", "H22"},          # N1-H, N2-H21/H22
    "U": {"H3"},                          # N3-H
    "C": {"H41", "H42"},                  # N4-H41/H42
    "A": {"H2", "H61", "H62"},            # C2-H2 (minor groove), N6-H61/H62
}

# All four recognized base names
ALL_BASES = ("G", "U", "C", "A")

# All unordered pairs (with repetition) among ALL_BASES
ALL_PAIR_CODES = [
    a + b
    for i, a in enumerate(ALL_BASES)
    for b in ALL_BASES[i:]
]


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def norm_icode(ic: str) -> str:
    if ic is None:
        return " "
    s = ic.strip()
    if s == "" or s.upper() == "NA":
        return " "
    return s[0]


def infer_element_from_atomname(atomname: str) -> str:
    s = atomname.strip()
    i = 0
    while i < len(s) and s[i].isdigit():
        i += 1
    s = s[i:]
    if not s:
        return ""
    return s[0].upper() if s[0].isalpha() else ""


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ResidKey:
    chain: str
    resseq: int
    icode: str
    resname: str

    def label(self) -> str:
        ic = self.icode.strip()
        return f"{self.resseq}{ic}" if ic else f"{self.resseq}"


@dataclass
class AtomRec:
    atomname: str
    element: str
    x: float
    y: float
    z: float
    resid: ResidKey
    line: str


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def parse_pdb_atoms(pdb_path: str) -> List[AtomRec]:
    atoms: List[AtomRec] = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not (line.startswith("ATOM ") or line.startswith("HETATM")):
                continue
            if len(line) < 54:
                continue
            try:
                atomname = line[12:16].strip()
                resname  = line[17:20].strip()
                chain    = line[21].strip() or " "
                resseq   = int(line[22:26])
                icode    = norm_icode(line[26:27])
                x        = float(line[30:38])
                y        = float(line[38:46])
                z        = float(line[46:54])
            except ValueError:
                continue
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                element = infer_element_from_atomname(atomname)
            resid = ResidKey(chain, resseq, icode, resname)
            atoms.append(AtomRec(atomname, element, x, y, z, resid, line))
    return atoms


def group_atoms_by_residue(atoms: Iterable[AtomRec]) -> Dict[ResidKey, List[AtomRec]]:
    d: Dict[ResidKey, List[AtomRec]] = {}
    for a in atoms:
        d.setdefault(a.resid, []).append(a)
    return d


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def dist2(a: AtomRec, b: AtomRec) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return dx * dx + dy * dy + dz * dz


def filter_base_atoms(atoms: List[AtomRec], resname: str) -> List[AtomRec]:
    allowed = BASE_ATOMS.get(resname, set())
    return [a for a in atoms if a.atomname in allowed]


def is_sequence_neighbor(r1: ResidKey, r2: ResidKey) -> bool:
    if r1.chain != r2.chain:
        return False
    if r1.icode.strip() or r2.icode.strip():
        return False
    return abs(r1.resseq - r2.resseq) == 1


def base_heavy_coords(base_atoms: List[AtomRec]) -> np.ndarray:
    pts = [(a.x, a.y, a.z) for a in base_atoms if a.element in {"C", "N", "O"}]
    if not pts:
        return np.zeros((0, 3), dtype=float)
    return np.array(pts, dtype=float)


def fit_plane_normal_and_centroid(
    pts: np.ndarray,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    if pts.shape[0] < 3:
        return None
    c = pts.mean(axis=0)
    X = pts - c
    C = (X.T @ X) / float(X.shape[0])
    w, v = np.linalg.eigh(C)
    n = v[:, 0]
    norm = np.linalg.norm(n)
    if norm == 0.0:
        return None
    return n / norm, c


def stacking_veto(
    r1_base_atoms: List[AtomRec],
    r2_base_atoms: List[AtomRec],
    stack_theta_deg: float,
    stack_dz: float,
) -> bool:
    r1_pts = base_heavy_coords(r1_base_atoms)
    r2_pts = base_heavy_coords(r2_base_atoms)
    r1_fit = fit_plane_normal_and_centroid(r1_pts)
    r2_fit = fit_plane_normal_and_centroid(r2_pts)
    if r1_fit is None or r2_fit is None:
        return False
    n1, c1 = r1_fit
    n2, c2 = r2_fit

    dot   = float(abs(np.dot(n1, n2)))
    dot   = max(-1.0, min(1.0, dot))
    theta = math.degrees(math.acos(dot))
    dz    = float(abs(np.dot((c2 - c1), n1)))

    return (theta < stack_theta_deg) and (dz > stack_dz)


# ---------------------------------------------------------------------------
# Contact finders
# ---------------------------------------------------------------------------

def find_all_NO_contacts_unfiltered(
    r1_atoms: List[AtomRec],
    r2_atoms: List[AtomRec],
) -> List[Tuple[AtomRec, AtomRec, float]]:
    """
    All N···O base-base (and sugar) contacts between r1 and r2.
    No distance cutoff is applied.  Both directions are checked.
    Returns (N_atom, O_atom, distance).
    """
    contacts: List[Tuple[AtomRec, AtomRec, float]] = []
    for a1 in r1_atoms:
        if a1.element != "N":
            continue
        for a2 in r2_atoms:
            if a2.element != "O":
                continue
            d = math.dist((a1.x, a1.y, a1.z), (a2.x, a2.y, a2.z))
            contacts.append((a1, a2, d))

    for a2 in r2_atoms:
        if a2.element != "N":
            continue
        for a1 in r1_atoms:
            if a1.element != "O":
                continue
            d = math.dist((a2.x, a2.y, a2.z), (a1.x, a1.y, a1.z))
            contacts.append((a2, a1, d))

    return contacts


def find_all_HO_contacts(
    r1_atoms: List[AtomRec],
    r2_atoms: List[AtomRec],
    cutoff2: float,
    nh_hnames: Set[str],
    resname1: str,
    resname2: str,
) -> List[Tuple[AtomRec, AtomRec, float]]:
    """
    H···O contacts where:
      - H is a hydrogen with name in nh_hnames (NH donors),
      - O is any base/sugar oxygen listed in BASE_ATOMS for the partner residue.
    Both directions (H on r1 -> O on r2, and H on r2 -> O on r1) are checked.
    Returns (H_atom, O_atom, distance), filtered by cutoff2.
    """
    contacts: List[Tuple[AtomRec, AtomRec, float]] = []

    def is_target_H(a: AtomRec) -> bool:
        return a.element == "H" and a.atomname in nh_hnames

    def is_base_O(a: AtomRec, resname: str) -> bool:
        return a.element == "O" and a.atomname in BASE_ATOMS.get(resname, set())

    # H on r1 -> O on r2
    for ah in r1_atoms:
        if not is_target_H(ah):
            continue
        for ao in r2_atoms:
            if not is_base_O(ao, resname2):
                continue
            d2 = dist2(ah, ao)
            if d2 <= cutoff2:
                contacts.append((ah, ao, math.sqrt(d2)))

    # H on r2 -> O on r1
    for ah in r2_atoms:
        if not is_target_H(ah):
            continue
        for ao in r1_atoms:
            if not is_base_O(ao, resname1):
                continue
            d2 = dist2(ah, ao)
            if d2 <= cutoff2:
                contacts.append((ah, ao, math.sqrt(d2)))

    return contacts


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def make_outname(pdbin: str, r1: ResidKey, r2: ResidKey) -> str:
    base = os.path.basename(pdbin)
    d    = os.path.dirname(os.path.abspath(pdbin))
    tag  = f"{r1.resname}{r1.chain}{r1.label()}_{r2.resname}{r2.chain}{r2.label()}"
    return os.path.join(d, f"{tag}_{base}")


def write_pair_pdb(
    outpath: str, r1_atoms: List[AtomRec], r2_atoms: List[AtomRec]
) -> None:
    with open(outpath, "w", encoding="utf-8") as f:
        for a in r1_atoms:
            f.write(a.line)
        for a in r2_atoms:
            f.write(a.line)
        f.write("END\n")


# ---------------------------------------------------------------------------
# Pair-code parsing
# ---------------------------------------------------------------------------

def parse_pair_codes(raw: str) -> List[Tuple[str, str]]:
    """
    Parse user input like "GU,GC,AU" (or "ALL") into a list of (base1, base2) tuples.
    The returned list contains CANONICAL pairs: base1 <= base2 alphabetically,
    unless both are the same (e.g. GG).  Duplicates are removed.
    Ordering within a pair is preserved as typed; only duplicates removed.
    """
    if raw.strip().upper() == "ALL":
        codes = ALL_PAIR_CODES
    else:
        codes = [c.strip().upper() for c in raw.split(",") if c.strip()]

    seen: Set[str] = set()
    result: List[Tuple[str, str]] = []
    for code in codes:
        if len(code) != 2:
            print(f"WARNING: '{code}' is not a 2-character base pair code — skipping.",
                  file=sys.stderr)
            continue
        b1, b2 = code[0], code[1]
        if b1 not in ALL_BASES:
            print(f"WARNING: '{b1}' is not a recognized base (G/U/C/A) — skipping '{code}'.",
                  file=sys.stderr)
            continue
        if b2 not in ALL_BASES:
            print(f"WARNING: '{b2}' is not a recognized base (G/U/C/A) — skipping '{code}'.",
                  file=sys.stderr)
            continue
        key = b1 + b2
        if key not in seen:
            seen.add(key)
            result.append((b1, b2))
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Report BASE–BASE and sugar NH···O contacts between user-specified base pairs; "
            "select pairs based on H···O cutoff; stacking filter; write mini-PDBs; "
            "also report ALL N···O contacts for selected pairs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # GU pairs only (backward-compatible):
  %(prog)s --pdbin 1msy.pdb --cut 3.5 --pairs GU

  # GU and GC pairs:
  %(prog)s --pdbin 1msy.pdb --cut 3.5 --pairs GU,GC

  # Watson-Crick + wobble:
  %(prog)s --pdbin 1msy.pdb --cut 3.5 --pairs GU,GC,AU,AG

  # Every combination among G, U, C, A:
  %(prog)s --pdbin 1msy.pdb --cut 3.5 --pairs ALL

  # Symmetric same-base pairs (e.g. all GG contacts):
  %(prog)s --pdbin 1msy.pdb --cut 3.5 --pairs GG
""",
    )

    ap.add_argument("--pdbin", required=True, help="Input PDB file.")
    ap.add_argument(
        "--cut",
        type=float,
        required=True,
        help="Distance cutoff for H···O contacts (Å). N···O contacts are always unfiltered.",
    )
    ap.add_argument(
        "--pairs",
        default="GU",
        help=(
            "Comma-separated base-pair codes to search for, e.g. GU,GC,AU,AG. "
            "Use ALL to search every combination among G, U, C, A. "
            "Default: GU (backward-compatible)."
        ),
    )
    ap.add_argument("--dryrun", action="store_true", help="Do not write mini-PDB files.")

    ap.add_argument(
        "--stack_theta",
        type=float,
        default=25.0,
        help="Stacking veto: reject if base-plane angle theta < this (degrees). Default: 25.",
    )
    ap.add_argument(
        "--stack_dz",
        type=float,
        default=2.0,
        help="Stacking veto: reject if inter-plane separation dz > this (Å). Default: 2.0.",
    )
    ap.add_argument(
        "--no_stack_filter",
        action="store_true",
        help="Disable stacking veto entirely.",
    )

    ap.add_argument(
        "--summary",
        default=None,
        help="Path to text file where all N···O contacts (unfiltered) will be appended.",
    )
    ap.add_argument(
        "--summary_HO",
        default=None,
        help="Path to text file where all H···O (NH···O) contacts within cutoff will be appended.",
    )
    ap.add_argument(
        "--nh_names",
        default=None,
        help=(
            "Comma-separated list of H atom names treated as N–H donors. "
            "If omitted, a sensible default is derived from the selected base types: "
            "G→H1,H21,H22  U→H3  C→H41,H42  A→H2,H61,H62."
        ),
    )

    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    if args.cut <= 0:
        print("ERROR: --cut must be > 0", file=sys.stderr)
        sys.exit(2)
    if args.stack_theta < 0 or args.stack_theta > 90:
        print("ERROR: --stack_theta should be in [0, 90]", file=sys.stderr)
        sys.exit(2)
    if args.stack_dz < 0:
        print("ERROR: --stack_dz must be >= 0", file=sys.stderr)
        sys.exit(2)

    # ------------------------------------------------------------------
    # Parse pair codes
    # ------------------------------------------------------------------
    pair_codes = parse_pair_codes(args.pairs)
    if not pair_codes:
        print("ERROR: No valid base-pair codes provided. Use e.g. --pairs GU,GC", file=sys.stderr)
        sys.exit(2)

    involved_bases: Set[str] = set()
    for b1, b2 in pair_codes:
        involved_bases.add(b1)
        involved_bases.add(b2)

    print(f"# Searching for pair type(s): {', '.join(b1+b2 for b1,b2 in pair_codes)}",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # NH hydrogen names
    # ------------------------------------------------------------------
    if args.nh_names:
        nh_hnames: Set[str] = {h.strip().upper() for h in args.nh_names.split(",") if h.strip()}
    else:
        nh_hnames = set()
        for base in involved_bases:
            nh_hnames |= DEFAULT_NH_HNAMES.get(base, set())

    print(f"# NH hydrogen names used : {', '.join(sorted(nh_hnames))}", file=sys.stderr)

    # ------------------------------------------------------------------
    # Load PDB
    # ------------------------------------------------------------------
    atoms  = parse_pdb_atoms(args.pdbin)
    by_res = group_atoms_by_residue(atoms)

    # Index residues by base name for quick lookup
    resids_by_base: Dict[str, List[ResidKey]] = {b: [] for b in ALL_BASES}
    for r in by_res:
        if r.resname in resids_by_base:
            resids_by_base[r.resname].append(r)
    for b in resids_by_base:
        resids_by_base[b].sort(key=lambda r: (r.chain, r.resseq, r.icode))

    cutoff2 = args.cut * args.cut
    written = 0
    summary_NO_lines:  List[str] = []
    summary_HO_lines:  List[str] = []

    pdb_base = os.path.basename(args.pdbin)
    pdb_id   = os.path.splitext(pdb_base)[0]

    # ------------------------------------------------------------------
    # Iterate over requested pair types
    # ------------------------------------------------------------------
    for (b1, b2) in pair_codes:
        r1_list = resids_by_base.get(b1, [])
        r2_list = resids_by_base.get(b2, [])

        if not r1_list:
            print(f"# WARNING: No '{b1}' residues found in {args.pdbin}", file=sys.stderr)
            continue
        if not r2_list:
            print(f"# WARNING: No '{b2}' residues found in {args.pdbin}", file=sys.stderr)
            continue

        same_base = (b1 == b2)

        # Generate candidate pairs, avoiding self-pairing and double-counting
        # for symmetric (same-base) cases.
        for idx1, r1 in enumerate(r1_list):
            r1_all      = by_res[r1]
            r1_base_all = filter_base_atoms(r1_all, b1)
            if not r1_base_all:
                continue

            for idx2, r2 in enumerate(r2_list):
                # For same-base pairs: skip lower-triangle to avoid duplicates
                if same_base and idx2 <= idx1:
                    continue
                # Same residue is never a partner of itself
                if r1 == r2:
                    continue
                if is_sequence_neighbor(r1, r2):
                    continue

                r2_all      = by_res[r2]
                r2_base_all = filter_base_atoms(r2_all, b2)
                if not r2_base_all:
                    continue

                if (not args.no_stack_filter) and stacking_veto(
                    r1_base_all, r2_base_all, args.stack_theta, args.stack_dz
                ):
                    continue

                # H···O contacts (cutoff-filtered) — selects the pair
                ho_contacts = find_all_HO_contacts(
                    r1_all, r2_all, cutoff2, nh_hnames, b1, b2
                )
                if not ho_contacts:
                    continue

                # N···O contacts (no cutoff) for the selected pair
                no_contacts = find_all_NO_contacts_unfiltered(r1_base_all, r2_base_all)

                # Print N···O contacts to stdout
                for aN, aO, d in sorted(no_contacts, key=lambda x: x[2]):
                    print(
                        f"{aN.resid.resname} {aN.resid.chain} {aN.resid.resseq} "
                        f"{aN.resid.icode:>2} {aN.atomname:<4} <-> "
                        f"{aO.resid.resname} {aO.resid.chain} {aO.resid.resseq} "
                        f"{aO.resid.icode:>2} {aO.atomname:<4} {d:6.3f}"
                    )
                    summary_NO_lines.append(
                        "\t".join([
                            pdb_id,
                            aN.resid.chain, str(aN.resid.resseq),
                            aN.resid.icode.strip() or ".",
                            aN.resid.resname, aN.atomname,
                            aO.resid.chain, str(aO.resid.resseq),
                            aO.resid.icode.strip() or ".",
                            aO.resid.resname, aO.atomname,
                            f"{d:.3f}",
                        ])
                    )

                # H···O contacts
                for aH, aO, d in sorted(ho_contacts, key=lambda x: x[2]):
                    summary_HO_lines.append(
                        "\t".join([
                            pdb_id,
                            aH.resid.chain, str(aH.resid.resseq),
                            aH.resid.icode.strip() or ".",
                            aH.resid.resname, aH.atomname,
                            aO.resid.chain, str(aO.resid.resseq),
                            aO.resid.icode.strip() or ".",
                            aO.resid.resname, aO.atomname,
                            f"{d:.3f}",
                        ])
                    )

                # Write mini-PDB
                if not args.dryrun:
                    outpath = make_outname(args.pdbin, r1, r2)
                    write_pair_pdb(outpath, r1_all, r2_all)
                    written += 1

    # ------------------------------------------------------------------
    # Write summary files
    # ------------------------------------------------------------------
    if args.summary and summary_NO_lines:
        header = (
            "#PDB_ID\tN_chain\tN_resSeq\tN_iCode\tN_resName\tN_atom\t"
            "O_chain\tO_resSeq\tO_iCode\tO_resName\tO_atom\tdist\n"
        )
        file_exists = os.path.exists(args.summary)
        with open(args.summary, "a", encoding="utf-8") as fout:
            if not file_exists:
                fout.write(header)
            for line in summary_NO_lines:
                fout.write(line + "\n")
            fout.write(f"#SUMMARY\t{pdb_id}\t{len(summary_NO_lines)}\n")

    if args.summary_HO and summary_HO_lines:
        header_ho = (
            "#PDB_ID\tH_chain\tH_resSeq\tH_iCode\tH_resName\tH_atom\t"
            "O_chain\tO_resSeq\tO_iCode\tO_resName\tO_atom\tdist\n"
        )
        file_exists = os.path.exists(args.summary_HO)
        with open(args.summary_HO, "a", encoding="utf-8") as fout:
            if not file_exists:
                fout.write(header_ho)
            for line in summary_HO_lines:
                fout.write(line + "\n")
            fout.write(f"#SUMMARY\t{pdb_id}\t{len(summary_HO_lines)}\n")

    if not args.dryrun:
        print(f"# pairs written: {written}", file=sys.stderr)


if __name__ == "__main__":
    main()
