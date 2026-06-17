#!/usr/bin/env python3
"""
Find residue pairs that have BASE–BASE (and sugar) H···O contacts within a cutoff.

Supported base types : G, U, C, A
Supported pair types : any combination, e.g. GU, GC, AU, AG, GG, AA …

Speed improvements over v1
──────────────────────────
1. Accepts MULTIPLE PDB files in one invocation (--pdbin / --pdblist / --pdbdir).
   No more shell-loop overhead of starting Python 5000 times.
2. Parallel processing via --workers (default: all CPU cores).
3. Centroid pre-filter: residue pairs whose centroids are > --centroid_cut apart
   are skipped before any per-atom work is done.
4. Optional scipy KDTree for the H···O radius search (auto-detected at import).
5. Live progress bar on stderr: N/total  elapsed  ETA  rate.

Behavior (unchanged from v1)
────────────────────────────
- --cut is applied ONLY to H···O contacts.
- Pairs are selected based on >= 1 H···O contact within --cut.
- For each selected pair: ALL N···O contacts (no cutoff) + H···O within cut
  are reported; one mini-PDB is written.

Outputs
───────
--summary      N···O contacts (unfiltered) for H···O-selected pairs.
--summary_HO   H···O contacts within cutoff.
mini-PDB files two-residue PDBs for each selected pair.

Usage
─────
  # Single PDB (backward-compatible):
  python find_rna_base_contacts_NHO.py --pdbin 1msy.pdb --cut 3.5 --pairs GU

  # 5000 PDBs via list file, 8 cores:

  python3 find_rna_base_contacts_NHO.py \
        --pdblist pdblist_reduce.txt \
        --cut 3.0 \
	    --pairs GU \
		--workers 8 \
        --summary    NO_summary.txt \
        --summary_HO HO_summary.txt
        
  # Whole directory:
  python find_rna_base_contacts_NHO.py --pdbdir /data/rna/ --cut 3.5 \\
      --pairs GU,GC --workers 8 \\
      --summary results/NO_summary.txt --summary_HO results/HO_summary.txt
"""

import argparse
import concurrent.futures
import math
import os
import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np

# Optional scipy KDTree — significant speedup for large structures
try:
    from scipy.spatial import KDTree as _KDTree
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ══════════════════════════════════════════════════════════════════════════════
# Chemistry constants
# ══════════════════════════════════════════════════════════════════════════════

BASE_ATOMS: Dict[str, Set[str]] = {
    "G": {"N1","C2","N2","N3","C4","C5","C6","O6","N7","C8","N9",
          "O2'","O3'","O4'","O5'"},
    "U": {"N1","C2","O2","N3","C4","O4","C5","C6",
          "O2'","O3'","O4'","O5'"},
    "C": {"N1","C2","O2","N3","C4","N4","C5","C6",
          "O2'","O3'","O4'","O5'"},
    "A": {"N1","C2","N3","C4","C5","C6","N6","N7","C8","N9",
          "O2'","O3'","O4'","O5'"},
}

DEFAULT_NH_HNAMES: Dict[str, Set[str]] = {
    "G": {"H1", "H21", "H22"},
    "U": {"H3"},
    "C": {"H41", "H42"},
    "A": {"H2", "H61", "H62"},
}

# Donor heavy atom (D) covalently bonded to each H donor atom.
# Used to compute the D–H···A angle.  A geometric fallback (nearest
# heavy atom within bond distance) is used if a name is not listed.
DONOR_OF_H: Dict[str, Dict[str, str]] = {
    "G": {"H1": "N1", "H21": "N2", "H22": "N2"},
    "U": {"H3": "N3"},
    "C": {"H41": "N4", "H42": "N4"},
    "A": {"H2": "C2", "H61": "N6", "H62": "N6"},
}

ALL_BASES = ("G", "U", "C", "A")
ALL_PAIR_CODES = [
    a + b
    for i, a in enumerate(ALL_BASES)
    for b in ALL_BASES[i:]
]

# ══════════════════════════════════════════════════════════════════════════════
# Small helpers
# ══════════════════════════════════════════════════════════════════════════════

def norm_icode(ic: str) -> str:
    if ic is None:
        return " "
    s = ic.strip()
    return " " if (s == "" or s.upper() == "NA") else s[0]


def infer_element_from_atomname(atomname: str) -> str:
    s = atomname.strip().lstrip("0123456789")
    return s[0].upper() if s and s[0].isalpha() else ""


# ══════════════════════════════════════════════════════════════════════════════
# Data structures
# ══════════════════════════════════════════════════════════════════════════════

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


# ══════════════════════════════════════════════════════════════════════════════
# PDB parsing
# ══════════════════════════════════════════════════════════════════════════════

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


def group_atoms_by_residue(
    atoms: Iterable[AtomRec],
) -> Dict[ResidKey, List[AtomRec]]:
    d: Dict[ResidKey, List[AtomRec]] = {}
    for a in atoms:
        d.setdefault(a.resid, []).append(a)
    return d


# ══════════════════════════════════════════════════════════════════════════════
# Geometry helpers
# ══════════════════════════════════════════════════════════════════════════════

def dist2(a: AtomRec, b: AtomRec) -> float:
    dx, dy, dz = a.x - b.x, a.y - b.y, a.z - b.z
    return dx*dx + dy*dy + dz*dz


def find_donor_atom(
    h_atom: AtomRec,
    residue_atoms: List[AtomRec],
) -> Optional[AtomRec]:
    """
    Heavy atom (D) covalently bonded to a hydrogen.
    Tries the DONOR_OF_H name map first, then falls back to the nearest
    heavy atom within bond distance (1.3 A).
    """
    dname = DONOR_OF_H.get(h_atom.resid.resname, {}).get(h_atom.atomname)
    if dname:
        for a in residue_atoms:
            if a.atomname == dname:
                return a
    best: Optional[AtomRec] = None
    best_d2 = 1.69  # (1.3 A)^2
    for a in residue_atoms:
        if a is h_atom or a.element == "H":
            continue
        d2 = dist2(h_atom, a)
        if d2 < best_d2:
            best_d2 = d2
            best = a
    return best


def angle_at_h(donor: AtomRec, h: AtomRec, acceptor: AtomRec) -> Optional[float]:
    """D–H···A angle (degrees), measured at the H atom."""
    v1 = (donor.x - h.x, donor.y - h.y, donor.z - h.z)
    v2 = (acceptor.x - h.x, acceptor.y - h.y, acceptor.z - h.z)
    n1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    n2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    if n1 == 0.0 or n2 == 0.0:
        return None
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    c = max(-1.0, min(1.0, dot / (n1 * n2)))
    return math.degrees(math.acos(c))


def filter_base_atoms(atoms: List[AtomRec], resname: str) -> List[AtomRec]:
    allowed = BASE_ATOMS.get(resname, set())
    return [a for a in atoms if a.atomname in allowed]


def is_sequence_neighbor(r1: ResidKey, r2: ResidKey) -> bool:
    if r1.chain != r2.chain:
        return False
    if r1.icode.strip() or r2.icode.strip():
        return False
    return abs(r1.resseq - r2.resseq) == 1


def c1c1_distance(
    r1_atoms: List[AtomRec],
    r2_atoms: List[AtomRec],
) -> Optional[float]:
    """
    C1'–C1' distance between two residues.
    Returns None if either residue lacks a C1' atom.
    """
    c1_r1 = next((a for a in r1_atoms if a.atomname == "C1'"), None)
    c1_r2 = next((a for a in r2_atoms if a.atomname == "C1'"), None)
    if c1_r1 is None or c1_r2 is None:
        return None
    return math.sqrt(dist2(c1_r1, c1_r2))


def centroid_xyz(atoms: List[AtomRec]) -> Optional[np.ndarray]:
    if not atoms:
        return None
    return np.mean([[a.x, a.y, a.z] for a in atoms], axis=0)


def base_heavy_coords(base_atoms: List[AtomRec]) -> np.ndarray:
    pts = [(a.x, a.y, a.z) for a in base_atoms if a.element in {"C","N","O"}]
    return np.array(pts, dtype=float) if pts else np.zeros((0, 3))


def fit_plane_normal_and_centroid(
    pts: np.ndarray,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    if pts.shape[0] < 3:
        return None
    c = pts.mean(axis=0)
    X = pts - c
    C = (X.T @ X) / float(X.shape[0])
    _, v = np.linalg.eigh(C)
    n = v[:, 0]
    nrm = np.linalg.norm(n)
    return (n / nrm, c) if nrm > 0 else None


def stacking_veto(
    r1_base: List[AtomRec],
    r2_base: List[AtomRec],
    theta_deg: float,
    dz_thresh: float,
) -> bool:
    r1_pts = base_heavy_coords(r1_base)
    r2_pts = base_heavy_coords(r2_base)
    r1_fit = fit_plane_normal_and_centroid(r1_pts)
    r2_fit = fit_plane_normal_and_centroid(r2_pts)
    if r1_fit is None or r2_fit is None:
        return False
    n1, c1 = r1_fit
    n2, c2 = r2_fit
    dot   = float(abs(np.dot(n1, n2)))
    theta = math.degrees(math.acos(max(-1.0, min(1.0, dot))))
    dz    = float(abs(np.dot((c2 - c1), n1)))
    return (theta < theta_deg) and (dz > dz_thresh)


# ══════════════════════════════════════════════════════════════════════════════
# Contact finders
# ══════════════════════════════════════════════════════════════════════════════

def find_all_NO_contacts_unfiltered(
    r1_atoms: List[AtomRec],
    r2_atoms: List[AtomRec],
) -> List[Tuple[AtomRec, AtomRec, float]]:
    contacts = []
    for a1 in r1_atoms:
        if a1.element != "N":
            continue
        for a2 in r2_atoms:
            if a2.element != "O":
                continue
            contacts.append((a1, a2, math.sqrt(dist2(a1, a2))))
    for a2 in r2_atoms:
        if a2.element != "N":
            continue
        for a1 in r1_atoms:
            if a1.element != "O":
                continue
            contacts.append((a2, a1, math.sqrt(dist2(a2, a1))))
    return contacts


def find_all_HO_contacts(
    r1_atoms: List[AtomRec],
    r2_atoms: List[AtomRec],
    cutoff2: float,
    nh_hnames: Set[str],
    resname1: str,
    resname2: str,
) -> List[Tuple[AtomRec, AtomRec, float, Optional[float], Optional[AtomRec]]]:
    """H···O contacts within cutoff. Uses scipy KDTree when available.

    Each contact is
        (H_atom, O_atom, H···O_distance, D-H···A_angle_deg, donor_atom).
    The angle / donor are None if the donor heavy atom cannot be located.
    The donor atom lets callers also report the N···O distance and the same
    N–H···O angle keyed on the donor nitrogen (see NHO summary).
    """
    acc1 = ({a for a in BASE_ATOMS.get(resname1, set()) if a.startswith("O")}
            | {"O2'","O3'","O4'","O5'"})
    acc2 = ({a for a in BASE_ATOMS.get(resname2, set()) if a.startswith("O")}
            | {"O2'","O3'","O4'","O5'"})
    contacts: List[Tuple[AtomRec, AtomRec, float, Optional[float],
                         Optional[AtomRec]]] = []

    def _kdtree_search(h_list, o_list, o_acc):
        o_filt = [a for a in o_list if a.element == "O" and a.atomname in o_acc]
        if not o_filt:
            return
        coords = np.array([[a.x, a.y, a.z] for a in o_filt])
        tree   = _KDTree(coords)
        cut    = math.sqrt(cutoff2)
        for ah in h_list:
            if ah.element != "H" or ah.atomname not in nh_hnames:
                continue
            don = find_donor_atom(ah, h_list)
            for i in tree.query_ball_point([ah.x, ah.y, ah.z], cut):
                ao = o_filt[i]
                ang = angle_at_h(don, ah, ao) if don is not None else None
                contacts.append((ah, ao, math.sqrt(dist2(ah, ao)), ang, don))

    def _naive_search(h_list, o_list, o_acc):
        for ah in h_list:
            if ah.element != "H" or ah.atomname not in nh_hnames:
                continue
            don = find_donor_atom(ah, h_list)
            for ao in o_list:
                if ao.element != "O" or ao.atomname not in o_acc:
                    continue
                d2 = dist2(ah, ao)
                if d2 <= cutoff2:
                    ang = angle_at_h(don, ah, ao) if don is not None else None
                    contacts.append((ah, ao, math.sqrt(d2), ang, don))

    _search = _kdtree_search if HAS_SCIPY else _naive_search
    _search(r1_atoms, r2_atoms, acc2)
    _search(r2_atoms, r1_atoms, acc1)
    return contacts


# ══════════════════════════════════════════════════════════════════════════════
# Output helpers
# ══════════════════════════════════════════════════════════════════════════════

def make_outname(pdbin: str, r1: ResidKey, r2: ResidKey) -> str:
    base = os.path.basename(pdbin)
    d    = os.path.dirname(os.path.abspath(pdbin))
    tag  = f"{r1.resname}{r1.chain}{r1.label()}_{r2.resname}{r2.chain}{r2.label()}"
    return os.path.join(d, f"{tag}_{base}")


def write_pair_pdb(outpath, r1_atoms, r2_atoms):
    with open(outpath, "w", encoding="utf-8") as f:
        for a in r1_atoms:
            f.write(a.line)
        for a in r2_atoms:
            f.write(a.line)
        f.write("END\n")


# ══════════════════════════════════════════════════════════════════════════════
# Pair-code parsing
# ══════════════════════════════════════════════════════════════════════════════

def parse_pair_codes(raw: str) -> List[Tuple[str, str]]:
    all_bases_set = set(ALL_BASES)
    codes = ALL_PAIR_CODES if raw.strip().upper() == "ALL" \
            else [c.strip().upper() for c in raw.split(",") if c.strip()]
    seen: Set[str] = set()
    result: List[Tuple[str, str]] = []
    for code in codes:
        if len(code) != 2 or code[0] not in all_bases_set or code[1] not in all_bases_set:
            print(f"WARNING: skipping unrecognised pair code '{code}'", file=sys.stderr)
            continue
        if code not in seen:
            seen.add(code)
            result.append((code[0], code[1]))
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Worker function — processes ONE PDB file
# Must be module-level for ProcessPoolExecutor pickling
# ══════════════════════════════════════════════════════════════════════════════

def _process_one_pdb(args_tuple):
    """Returns (pdb_id, NO_lines, HO_lines, NHO_lines, n_written, warnings)."""
    (pdb_path, pair_codes, cutoff2, nh_hnames,
     dryrun, stack_theta, stack_dz, no_stack,
     centroid_cut, c1c1_min, c1c1_max) = args_tuple

    NO_lines: List[str] = []
    HO_lines: List[str] = []
    NHO_lines: List[str] = []
    warnings: List[str] = []
    written = 0
    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]

    try:
        atoms = parse_pdb_atoms(pdb_path)
    except Exception as exc:
        return pdb_id, [], [], [], 0, [f"SKIP {pdb_path}: {exc}"]

    by_res = group_atoms_by_residue(atoms)

    resids_by_base: Dict[str, List[ResidKey]] = {b: [] for b in ALL_BASES}
    for r in by_res:
        if r.resname in resids_by_base:
            resids_by_base[r.resname].append(r)
    for b in resids_by_base:
        resids_by_base[b].sort(key=lambda r: (r.chain, r.resseq, r.icode))

    # Precompute residue centroids for coarse distance pre-filter
    centroids: Dict[ResidKey, np.ndarray] = {}
    for r, ats in by_res.items():
        c = centroid_xyz(ats)
        if c is not None:
            centroids[r] = c

    centroid_cut2 = centroid_cut * centroid_cut

    for (b1, b2) in pair_codes:
        r1_list = resids_by_base.get(b1, [])
        r2_list = resids_by_base.get(b2, [])
        if not r1_list or not r2_list:
            continue

        same_base = (b1 == b2)

        # Build KDTree of r2 centroids for O(log N) coarse filtering
        if HAS_SCIPY and len(r2_list) > 4:
            r2_with_c = [(r, centroids[r]) for r in r2_list if r in centroids]
            r2_keys   = [x[0] for x in r2_with_c]
            r2_coords = np.array([x[1] for x in r2_with_c])
            r2_tree   = _KDTree(r2_coords) if len(r2_coords) else None
        else:
            r2_tree = None
            r2_keys = r2_list

        for idx1, r1 in enumerate(r1_list):
            c1 = centroids.get(r1)
            if c1 is None:
                continue

            r1_all      = by_res[r1]
            r1_base_all = filter_base_atoms(r1_all, b1)
            if not r1_base_all:
                continue

            # Coarse centroid filter
            if r2_tree is not None:
                nearby = [r2_keys[i]
                          for i in r2_tree.query_ball_point(c1, centroid_cut)]
            else:
                nearby = [r for r in r2_keys
                          if r not in centroids or
                          float(np.sum((centroids[r] - c1)**2)) <= centroid_cut2]

            for r2 in nearby:
                if same_base:
                    idx2 = r2_list.index(r2)
                    if idx2 <= idx1:
                        continue
                if r1 == r2 or is_sequence_neighbor(r1, r2):
                    continue

                r2_all      = by_res[r2]
                r2_base_all = filter_base_atoms(r2_all, b2)
                if not r2_base_all:
                    continue

                if (not no_stack) and stacking_veto(
                    r1_base_all, r2_base_all, stack_theta, stack_dz
                ):
                    continue

                # ── C1'–C1' distance filter ───────────────────────────────
                # Genuine base pairs have C1'-C1' in a tight range.
                # This rejects residues that are close in 3D space by chance
                # (crystal contacts, backbone proximity) but are not paired.
                if c1c1_min is not None or c1c1_max is not None:
                    d_c1c1 = c1c1_distance(r1_all, r2_all)
                    if d_c1c1 is not None:
                        if c1c1_min is not None and d_c1c1 < c1c1_min:
                            continue
                        if c1c1_max is not None and d_c1c1 > c1c1_max:
                            continue
                    # If C1' is absent (rare), let the pair through

                ho = find_all_HO_contacts(
                    r1_all, r2_all, cutoff2, nh_hnames, b1, b2)
                if not ho:
                    continue

                no = find_all_NO_contacts_unfiltered(r1_base_all, r2_base_all)

                for aN, aO, d in sorted(no, key=lambda x: x[2]):
                    NO_lines.append("\t".join([
                        pdb_id,
                        aN.resid.chain, str(aN.resid.resseq),
                        aN.resid.icode.strip() or ".",
                        aN.resid.resname, aN.atomname,
                        aO.resid.chain, str(aO.resid.resseq),
                        aO.resid.icode.strip() or ".",
                        aO.resid.resname, aO.atomname,
                        f"{d:.3f}",
                    ]))

                for aH, aO, d, ang, _ in sorted(ho, key=lambda x: x[2]):
                    HO_lines.append("\t".join([
                        pdb_id,
                        aH.resid.chain, str(aH.resid.resseq),
                        aH.resid.icode.strip() or ".",
                        aH.resid.resname, aH.atomname,
                        aO.resid.chain, str(aO.resid.resseq),
                        aO.resid.icode.strip() or ".",
                        aO.resid.resname, aO.atomname,
                        f"{d:.3f}",
                        f"{ang:.1f}" if ang is not None else ".",
                    ]))

                # NHO summary: same contacts keyed on the donor nitrogen.
                # Report the N···O heavy-atom distance and the N–H···O angle
                # (identical DHA angle as in the H···O summary).
                nho = [
                    (aN, aO, math.sqrt(dist2(aN, aO)), ang)
                    for aH, aO, d, ang, aN in ho
                    if aN is not None
                ]
                for aN, aO, dNO, ang in sorted(nho, key=lambda x: x[2]):
                    NHO_lines.append("\t".join([
                        pdb_id,
                        aN.resid.chain, str(aN.resid.resseq),
                        aN.resid.icode.strip() or ".",
                        aN.resid.resname, aN.atomname,
                        aO.resid.chain, str(aO.resid.resseq),
                        aO.resid.icode.strip() or ".",
                        aO.resid.resname, aO.atomname,
                        f"{dNO:.3f}",
                        f"{ang:.1f}" if ang is not None else ".",
                    ]))

                if not dryrun:
                    write_pair_pdb(make_outname(pdb_path, r1, r2), r1_all, r2_all)
                    written += 1

    return pdb_id, NO_lines, HO_lines, NHO_lines, written, warnings


# ══════════════════════════════════════════════════════════════════════════════
# Progress bar
# ══════════════════════════════════════════════════════════════════════════════

def _progress(done: int, total: int, t0: float) -> None:
    elapsed = time.time() - t0
    rate    = done / elapsed if elapsed > 0 else 0.0
    eta     = (total - done) / rate if rate > 0 else float("inf")
    pct     = 100.0 * done / total if total else 0.0
    bar_w   = 28
    filled  = int(bar_w * done / total) if total else 0
    bar     = "█" * filled + "░" * (bar_w - filled)
    eta_str = (f"{eta/60:.1f} min" if eta < 3600
               else ("done" if eta == 0 else f"{eta/3600:.1f} hr"))
    print(
        f"\r  [{bar}] {done}/{total} ({pct:.0f}%)  "
        f"{elapsed/60:.1f} min elapsed  ETA {eta_str}  "
        f"{rate:.1f} PDB/s   ",
        end="", file=sys.stderr, flush=True,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Find RNA base-base NH···O contacts. "
                    "Accepts single or many PDB files; parallelises with --workers.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Input (mutually exclusive)
    inp = ap.add_mutually_exclusive_group(required=True)
    inp.add_argument("--pdbin",   nargs="+", metavar="PDB",
                     help="One or more PDB files.")
    inp.add_argument("--pdblist", metavar="FILE",
                     help="Text file: one PDB path per line.")
    inp.add_argument("--pdbdir",  metavar="DIR",
                     help="Process every *.pdb in this directory.")

    # Core
    ap.add_argument("--cut",    type=float, required=True,
                    help="H···O distance cutoff (A).")
    ap.add_argument("--pairs",  default="GU",
                    help="Pair codes, e.g. GU,GC,AU or ALL. Default: GU.")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1,
                    help="Parallel worker processes (default: all CPU cores).")
    ap.add_argument("--centroid_cut", type=float, default=15.0,
                    help="Residue centroid pre-filter distance A (default: 15). "
                         "Pairs further apart are skipped before per-atom work.")

    # Output
    ap.add_argument("--summary",    default=None,
                    help="Append N···O contacts here.")
    ap.add_argument("--summary_HO", default=None,
                    help="Append H···O contacts within --cut here.")
    ap.add_argument("--summary_NHO", default=None,
                    help="Append the same --cut-selected contacts here, but "
                         "keyed on the donor N: N···O heavy-atom distance and "
                         "the N–H···O (DHA) angle.")
    ap.add_argument("--dryrun", action="store_true",
                    help="Skip mini-PDB writing.")

    # Stacking filter
    ap.add_argument("--stack_theta",     type=float, default=25.0)
    ap.add_argument("--stack_dz",        type=float, default=2.0)
    ap.add_argument("--no_stack_filter", action="store_true")

    # C1'-C1' distance filter  (the key guard against non-base-pair false positives)
    ap.add_argument("--c1c1_min", type=float, default=9.0,
                    help="Minimum C1'–C1' distance (Å) for a valid base pair. "
                         "Default: 9.0. Set to 0 to disable lower bound.")
    ap.add_argument("--c1c1_max", type=float, default=12.5,
                    help="Maximum C1'–C1' distance (Å) for a valid base pair. "
                         "Default: 12.5. Set to 0 to disable upper bound.")

    # NH donors
    ap.add_argument("--nh_names", default=None,
                    help="Comma-separated H donor atom names. Default: auto from --pairs.")

    args = ap.parse_args()

    # ── Collect PDB files ──────────────────────────────────────────────────────
    if args.pdbin:
        pdb_files = list(args.pdbin)
    elif args.pdblist:
        with open(args.pdblist, encoding="utf-8") as f:
            pdb_files = [ln.strip() for ln in f
                         if ln.strip() and not ln.startswith("#")]
    else:
        pdb_files = sorted(
            os.path.join(args.pdbdir, fn)
            for fn in os.listdir(args.pdbdir)
            if fn.lower().endswith(".pdb")
        )
    pdb_files = [p for p in pdb_files if os.path.isfile(p)]
    if not pdb_files:
        print("ERROR: no PDB files found.", file=sys.stderr); sys.exit(1)

    # ── Pairs & NH donors ─────────────────────────────────────────────────────
    pair_codes = parse_pair_codes(args.pairs)
    if not pair_codes:
        print("ERROR: no valid pair codes.", file=sys.stderr); sys.exit(2)

    involved: Set[str] = set()
    for b1, b2 in pair_codes:
        involved.update([b1, b2])

    if args.nh_names:
        nh_hnames: Set[str] = {h.strip().upper()
                               for h in args.nh_names.split(",") if h.strip()}
    else:
        nh_hnames = set()
        for b in involved:
            nh_hnames |= DEFAULT_NH_HNAMES.get(b, set())

    cutoff2 = args.cut * args.cut
    c1c1_min = args.c1c1_min if args.c1c1_min > 0 else None
    c1c1_max = args.c1c1_max if args.c1c1_max > 0 else None

    # ── Banner ────────────────────────────────────────────────────────────────
    print(f"# PDB files         : {len(pdb_files)}", file=sys.stderr)
    print(f"# Pair type(s)      : {', '.join(b1+b2 for b1,b2 in pair_codes)}", file=sys.stderr)
    print(f"# H···O cutoff      : {args.cut} A", file=sys.stderr)
    print(f"# C1'-C1' range     : {c1c1_min or 'none'} – {c1c1_max or 'none'} A", file=sys.stderr)
    print(f"# Centroid pre-cut  : {args.centroid_cut} A", file=sys.stderr)
    print(f"# Workers           : {args.workers}", file=sys.stderr)
    print(f"# scipy KDTree      : {'yes' if HAS_SCIPY else 'no  (pip install scipy)'}", file=sys.stderr)
    print(f"# NH donors         : {', '.join(sorted(nh_hnames))}", file=sys.stderr)

    # ── Dispatch ──────────────────────────────────────────────────────────────
    worker_args = [
        (p, pair_codes, cutoff2, nh_hnames,
         args.dryrun, args.stack_theta, args.stack_dz,
         args.no_stack_filter, args.centroid_cut, c1c1_min, c1c1_max)
        for p in pdb_files
    ]

    all_NO:  List[str] = []
    all_HO:  List[str] = []
    all_NHO: List[str] = []
    total_written = 0
    total = len(pdb_files)
    done  = 0
    t0    = time.time()
    _progress(0, total, t0)

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_process_one_pdb, a): a[0] for a in worker_args}
        for fut in concurrent.futures.as_completed(futs):
            done += 1
            try:
                _, NO, HO, NHO, n, _ = fut.result()
                all_NO.extend(NO)
                all_HO.extend(HO)
                all_NHO.extend(NHO)
                total_written += n
            except Exception as exc:
                print(f"\n# ERROR {futs[fut]}: {exc}", file=sys.stderr)
            _progress(done, total, t0)

    elapsed = time.time() - t0
    print(f"\n# Finished in {elapsed/60:.1f} min  |  pairs written: {total_written}",
          file=sys.stderr)

    # ── Write summary files ───────────────────────────────────────────────────
    if args.summary and all_NO:
        header = ("#PDB_ID\tN_chain\tN_resSeq\tN_iCode\tN_resName\tN_atom\t"
                  "O_chain\tO_resSeq\tO_iCode\tO_resName\tO_atom\tdist\n")
        mode = "a" if os.path.exists(args.summary) else "w"
        with open(args.summary, mode, encoding="utf-8") as fout:
            if mode == "w":
                fout.write(header)
            fout.writelines(ln + "\n" for ln in all_NO)
            fout.write(f"#SUMMARY\t{total}_pdbs\t{len(all_NO)}_contacts\n")
        print(f"# N···O summary  -> {args.summary}  ({len(all_NO)} lines)", file=sys.stderr)

    if args.summary_HO and all_HO:
        header_ho = ("#PDB_ID\tH_chain\tH_resSeq\tH_iCode\tH_resName\tH_atom\t"
                     "O_chain\tO_resSeq\tO_iCode\tO_resName\tO_atom\tdist\t"
                     "angle_DHA\n")
        mode = "a" if os.path.exists(args.summary_HO) else "w"
        with open(args.summary_HO, mode, encoding="utf-8") as fout:
            if mode == "w":
                fout.write(header_ho)
            fout.writelines(ln + "\n" for ln in all_HO)
            fout.write(f"#SUMMARY\t{total}_pdbs\t{len(all_HO)}_contacts\n")
        print(f"# H···O summary  -> {args.summary_HO}  ({len(all_HO)} lines)", file=sys.stderr)

    if args.summary_NHO and all_NHO:
        header_nho = ("#PDB_ID\tN_chain\tN_resSeq\tN_iCode\tN_resName\tN_atom\t"
                      "O_chain\tO_resSeq\tO_iCode\tO_resName\tO_atom\tdist\t"
                      "angle_DHA\n")
        mode = "a" if os.path.exists(args.summary_NHO) else "w"
        with open(args.summary_NHO, mode, encoding="utf-8") as fout:
            if mode == "w":
                fout.write(header_nho)
            fout.writelines(ln + "\n" for ln in all_NHO)
            fout.write(f"#SUMMARY\t{total}_pdbs\t{len(all_NHO)}_contacts\n")
        print(f"# N···O+angle summary -> {args.summary_NHO}  ({len(all_NHO)} lines)", file=sys.stderr)


if __name__ == "__main__":
    main()
