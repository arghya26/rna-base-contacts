#!/usr/bin/env python3
"""
Rank base-pair mini-PDBs from compute_rna_base_miniPDBs.py output and select top N QM models.

Input
-----
scored_pairs.txt : output from compute_rna_base_miniPDBs.py
Lines look like (tab-separated):
    Category  filepath  d_NO=2.95  C1C1=10.8  tilt=12.3  contacts=H1-O2:1.82,∠165°  H3-O6:1.78,∠172°

Pair-specific scoring functions
--------------------------------
GU  (wobble)
    2.0 * |H1–O2  – 1.8|    (absent → +4.0)
  + 2.0 * |H3–O6  – 1.8|    (absent → +2.0)
  + 0.05 * max(0, 170 – ∠H1O2)  (absent → +1.0)   ← angle term, GU only
  + 0.05 * max(0, 170 – ∠H3O6)  (absent → +0.5)   ← angle term, GU only
  + 1.0  * |d_NO  – 2.9|
  + 0.2  * |C1C1  – 10.8|   (absent → +1.0)        ← lowered weight
  + 0.5  * max(0, tilt – 15) (absent → +2.0)

GC  (Watson-Crick)
    2.0 * best(H21/H22–O2) from 1.8   (absent → +4.0)
  + 2.0 * best(H41/H42–O6) from 1.8   (absent → +4.0)
  + 1.0  * |d_NO  – 2.9|
  + 0.2  * |C1C1  – 10.4|   (absent → +1.0)
  + 0.5  * max(0, tilt – 15) (absent → +2.0)

AU  (Watson-Crick)
    2.0 * best(H61/H62–O4) from 1.8   (absent → +4.0)
  + 1.0  * |d_NO  – 2.9|
  + 0.2  * |C1C1  – 10.4|   (absent → +1.0)
  + 0.5  * max(0, tilt – 15) (absent → +2.0)

All other pair types  (AG, AC, UC, GG, …)
    2.0 * best H···O contact from 1.8  (absent → +4.0)
  + 1.0  * |d_NO  – 2.9|
  + 0.2  * |C1C1  – 10.6|   (absent → +1.0)
  + 0.5  * max(0, tilt – 15) (absent → +2.0)

Lower score = better QM candidate.

Usage
-----
    # GU only (backward-compatible):
    python rank_base_pairs.py scored_pairs.txt 10 --pairs GU

    # GU + GC + AU:
    python rank_base_pairs.py scored_pairs.txt 10 --pairs GU,GC,AU

    # Everything in the file:
    python rank_base_pairs.py scored_pairs.txt 10
"""

import argparse
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

# ══════════════════════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════════════════════

ALL_BASES: Set[str] = {"G", "U", "C", "A"}

# Canonical pair code ordering for output
_CANON_ORDER: List[str] = [
    "GC", "GU", "AU", "AG", "AC", "UC", "GG", "AA", "UU", "CC"
]

# Ideal C1'–C1' distance per pair type (Å)
_IDEAL_C1C1: Dict[str, float] = {
    "GU": 10.8,
    "GC": 10.4,
    "AU": 10.4,
}
_IDEAL_C1C1_DEFAULT = 10.6

# Category name → pair type mapping (covers score_base_pairs.py output)
_CAT_TO_PAIR: Dict[str, str] = {
    "Canonical_Aform":    "GU",
    "Bifurcated_contact": "GU",
    "GU_other":           "GU",
    "WatsonCrick_GC":     "GC",
    "GC_other":           "GC",
    "WatsonCrick_AU":     "AU",
    "AU_other":           "AU",
}

# ══════════════════════════════════════════════════════════════════════════════
# Data class
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class Entry:
    category:  str
    pair_type: str
    filepath:  str
    d_NO:      Optional[float]
    c1c1:      Optional[float]
    tilt:      Optional[float]
    # contacts: {key: (dist, angle_or_None)}  e.g. {"H1-O2": (1.82, 165.0)}
    contacts:  Dict[str, Tuple[float, Optional[float]]] = field(default_factory=dict)
    score:     float = 0.0

    # ── Contact queries ──────────────────────────────────────────────────────

    def dist(self, key: str) -> Optional[float]:
        """H···O distance for a contact key like 'H1-O2', or None."""
        return self.contacts[key][0] if key in self.contacts else None

    def angle(self, key: str) -> Optional[float]:
        """O–H···N angle (°) for a contact key, or None."""
        return self.contacts[key][1] if key in self.contacts else None

    def best_dist(self, *keys: str) -> Optional[float]:
        """Shortest H···O distance among the given contact keys."""
        vals = [self.contacts[k][0] for k in keys if k in self.contacts]
        return min(vals) if vals else None

    # ── Output ───────────────────────────────────────────────────────────────

    def contacts_str(self) -> str:
        if not self.contacts:
            return "None"
        parts = []
        for k, (d, a) in sorted(self.contacts.items(),
                                 key=lambda x: x[1][0]):
            ang = f",∠{a:.0f}°" if a is not None else ""
            parts.append(f"{k}:{d:.2f}{ang}")
        return "  ".join(parts)

    def row(self, top_n_contacts: int = 4) -> str:
        dno  = f"{self.d_NO:.3f}"  if self.d_NO  is not None else "NA"
        c1c1 = f"{self.c1c1:.3f}" if self.c1c1 is not None else "NA"
        tilt = f"{self.tilt:.1f}" if self.tilt  is not None else "NA"
        return (
            f"{self.filepath}\t"
            f"pair={self.pair_type}\t"
            f"category={self.category}\t"
            f"score={self.score:.3f}\t"
            f"d_NO={dno}\t"
            f"C1C1={c1c1}\t"
            f"tilt={tilt}\t"
            f"contacts={self.contacts_str()}"
        )


# ══════════════════════════════════════════════════════════════════════════════
# Pair-type inference
# ══════════════════════════════════════════════════════════════════════════════

def _norm_pair(b1: str, b2: str) -> str:
    """Canonical two-letter pair code; canonical order from _CANON_ORDER."""
    k  = b1 + b2
    k2 = b2 + b1
    for c in _CANON_ORDER:
        if c == k:
            return k
        if c == k2:
            return k2
    return k   # fallback


def infer_pair_type(filepath: str, category: str) -> Optional[str]:
    """
    Try (in order):
      1. New filename format  GD35_CC26_4pco.pdb  → first letter of each prefix
      2. Category name lookup (hardcoded table)
      3. Category prefix  e.g. 'AG_general' → 'AG'
    Returns canonical pair code or None if undetermined.
    """
    # ── 1. Filename ───────────────────────────────────────────────────────────
    stem  = os.path.basename(filepath)
    if stem.lower().endswith(".pdb"):
        stem = stem[:-4]
    parts = stem.split("_")
    # New format: ≥3 parts, first two start with a base letter
    if len(parts) >= 3:
        b1 = parts[0][0].upper() if parts[0] else ""
        b2 = parts[1][0].upper() if parts[1] else ""
        if b1 in ALL_BASES and b2 in ALL_BASES:
            return _norm_pair(b1, b2)

    # ── 2. Category table ─────────────────────────────────────────────────────
    if category in _CAT_TO_PAIR:
        return _CAT_TO_PAIR[category]

    # ── 3. Category prefix pattern  e.g. 'AG_general' ────────────────────────
    m = re.match(r"^([A-Z]{2})_", category)
    if m:
        tok = m.group(1)
        b1, b2 = tok[0], tok[1]
        if b1 in ALL_BASES and b2 in ALL_BASES:
            return _norm_pair(b1, b2)

    return None


# ══════════════════════════════════════════════════════════════════════════════
# Contact string parser
# ══════════════════════════════════════════════════════════════════════════════

# Matches:  H1-O2:1.82,∠165°   or   H1-O2:1.82
_CONTACT_RE = re.compile(
    r"([A-Za-z0-9']+)-([A-Za-z0-9']+)"   # H_atom - O_atom
    r":([0-9]+\.[0-9]+)"                   # :dist
    r"(?:,∠([0-9]+(?:\.[0-9]+)?)°)?"      # optional ,∠angle°
)

def parse_contacts(s: str) -> Dict[str, Tuple[float, Optional[float]]]:
    """
    Parse the contacts field from score_base_pairs.py output.
    Handles both old format (comma-separated, no angles) and
    new format (double-space-separated, optional ∠ tokens).

    Returns {key: (dist, angle_or_None)}  e.g. {"H1-O2": (1.82, 165.0)}
    """
    out: Dict[str, Tuple[float, Optional[float]]] = {}
    if not s or s.strip() == "None":
        return out
    for m in _CONTACT_RE.finditer(s):
        h_atom  = m.group(1)
        o_atom  = m.group(2)
        dist    = float(m.group(3))
        angle   = float(m.group(4)) if m.group(4) is not None else None
        key     = f"{h_atom}-{o_atom}"
        # Keep shortest distance if duplicate keys exist
        if key not in out or dist < out[key][0]:
            out[key] = (dist, angle)
    return out


# ══════════════════════════════════════════════════════════════════════════════
# Line parser
# ══════════════════════════════════════════════════════════════════════════════

def parse_line(line: str) -> Optional[Entry]:
    """
    Parse one scored-pair line into an Entry (score not yet computed).
    Returns None if the line cannot be parsed.
    """
    fields = line.strip().split("\t")
    if len(fields) < 5:
        return None

    category = fields[0].strip()
    filepath  = fields[1].strip()

    d_NO: Optional[float] = None
    c1c1: Optional[float] = None
    tilt: Optional[float] = None
    contacts_raw = ""

    for tok in fields[2:]:
        tok = tok.strip()
        if tok.startswith("d_NO="):
            try:
                d_NO = float(tok.split("=", 1)[1])
            except ValueError:
                pass
        elif tok.startswith("C1C1="):
            v = tok.split("=", 1)[1]
            if v != "NA":
                try:
                    c1c1 = float(v)
                except ValueError:
                    pass
        elif tok.startswith("tilt="):
            v = tok.split("=", 1)[1]
            if v != "NA":
                try:
                    tilt = float(v)
                except ValueError:
                    pass
        elif tok.startswith("contacts="):
            contacts_raw = tok.split("=", 1)[1]

    if d_NO is None:
        return None     # not a data line

    pair_type = infer_pair_type(filepath, category)
    if pair_type is None:
        return None     # cannot determine pair type

    contacts = parse_contacts(contacts_raw)

    return Entry(
        category=category,
        pair_type=pair_type,
        filepath=filepath,
        d_NO=d_NO,
        c1c1=c1c1,
        tilt=tilt,
        contacts=contacts,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Scoring functions  (lower = better)
# ══════════════════════════════════════════════════════════════════════════════

def _c1c1_term(c1c1: Optional[float], ideal: float) -> float:
    """C1'–C1' penalty with lowered weight (0.2)."""
    if c1c1 is not None:
        return 0.2 * abs(c1c1 - ideal)
    return 1.0    # unknown penalty (lowered from original 2.0)


def _tilt_term(tilt: Optional[float]) -> float:
    if tilt is not None:
        return 0.5 * max(0.0, tilt - 15.0)
    return 2.0


def _dNO_term(d_NO: Optional[float], ideal: float = 2.9) -> float:
    return 1.0 * abs(d_NO - ideal) if d_NO is not None else 2.0


def score_GU(e: Entry) -> float:
    """
    GU wobble scorer — includes O–H···N angle terms (GU only).
    Ideal angle for a linear H-bond = 180°; penalise below 170°.
    """
    s = 0.0

    # H1–O2 (canonical wobble, G H1 → U O2)
    d_H1O2 = e.dist("H1-O2")
    if d_H1O2 is not None:
        s += 2.0 * abs(d_H1O2 - 1.8)
        ang = e.angle("H1-O2")
        s += 0.05 * max(0.0, 170.0 - ang) if ang is not None else 1.0
    else:
        s += 4.0 + 1.0    # distance + angle penalty

    # H3–O6 (U H3 → G O6)
    d_H3O6 = e.dist("H3-O6")
    if d_H3O6 is not None:
        s += 2.0 * abs(d_H3O6 - 1.8)
        ang = e.angle("H3-O6")
        s += 0.05 * max(0.0, 170.0 - ang) if ang is not None else 0.5
    else:
        s += 2.0 + 0.5

    s += _dNO_term(e.d_NO)
    s += _c1c1_term(e.c1c1, _IDEAL_C1C1["GU"])
    s += _tilt_term(e.tilt)
    return s


def score_GC(e: Entry) -> float:
    """GC Watson-Crick scorer — no angle terms."""
    s = 0.0

    # G N2(H)···C O2 — use whichever amino H is closer
    best_g_amino = e.best_dist("H21-O2", "H22-O2")
    s += 2.0 * abs(best_g_amino - 1.8) if best_g_amino is not None else 4.0

    # C N4(H)···G O6
    best_c_amino = e.best_dist("H41-O6", "H42-O6")
    s += 2.0 * abs(best_c_amino - 1.8) if best_c_amino is not None else 4.0

    s += _dNO_term(e.d_NO)
    s += _c1c1_term(e.c1c1, _IDEAL_C1C1["GC"])
    s += _tilt_term(e.tilt)
    return s


def score_AU(e: Entry) -> float:
    """AU Watson-Crick scorer — no angle terms."""
    s = 0.0

    # A N6(H)···U O4
    best_a_amino = e.best_dist("H61-O4", "H62-O4")
    s += 2.0 * abs(best_a_amino - 1.8) if best_a_amino is not None else 4.0

    s += _dNO_term(e.d_NO)
    s += _c1c1_term(e.c1c1, _IDEAL_C1C1["AU"])
    s += _tilt_term(e.tilt)
    return s


def score_generic(e: Entry) -> float:
    """Generic scorer for pair types without a dedicated function — no angle terms."""
    s = 0.0
    best_HO = (min(v for v, _ in e.contacts.values())
               if e.contacts else None)
    s += 2.0 * abs(best_HO - 1.8) if best_HO is not None else 4.0

    s += _dNO_term(e.d_NO)
    s += _c1c1_term(e.c1c1, _IDEAL_C1C1_DEFAULT)
    s += _tilt_term(e.tilt)
    return s


_SCORERS = {
    "GU": score_GU,
    "GC": score_GC,
    "AU": score_AU,
}


def compute_score(e: Entry) -> float:
    scorer = _SCORERS.get(e.pair_type, score_generic)
    return scorer(e)


# ══════════════════════════════════════════════════════════════════════════════
# I/O
# ══════════════════════════════════════════════════════════════════════════════

def read_candidates(path: str,
                    keep_pairs: Optional[Set[str]]) -> List[Entry]:
    entries: List[Entry] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("═"):
                continue
            e = parse_line(line)
            if e is None:
                continue
            if keep_pairs is not None and e.pair_type not in keep_pairs:
                continue
            e.score = compute_score(e)
            entries.append(e)
    return entries


def parse_pairs_arg(raw: Optional[str]) -> Optional[Set[str]]:
    """Returns None (= all) or a set of canonical pair codes."""
    if raw is None or raw.strip().upper() == "ALL":
        return None
    codes: Set[str] = set()
    for tok in raw.split(","):
        tok = tok.strip().upper()
        if len(tok) != 2:
            print(f"WARNING: ignoring unrecognised token '{tok}'", file=sys.stderr)
            continue
        b1, b2 = tok[0], tok[1]
        if b1 not in ALL_BASES or b2 not in ALL_BASES:
            print(f"WARNING: '{tok}' contains unknown base letter — skipping.",
                  file=sys.stderr)
            continue
        codes.add(_norm_pair(b1, b2))
    return codes or None


# ══════════════════════════════════════════════════════════════════════════════
# Output helpers
# ══════════════════════════════════════════════════════════════════════════════

def _print_section(title: str, entries: List[Entry], n: int) -> None:
    if not entries:
        return
    print(f"# {title}")
    for e in sorted(entries, key=lambda x: x.score)[:n]:
        print(e.row())
    print()


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Rank base-pair mini-PDBs from score_base_pairs.py output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # GU only:
  %(prog)s scored.txt 10 --pairs GU

  # GU + GC + AU:
  %(prog)s scored.txt 10 --pairs GU,GC,AU

  # All pair types:
  %(prog)s scored.txt 10
""",
    )
    ap.add_argument("input",  help="Input file from score_base_pairs.py")
    ap.add_argument("n_top",  type=int, help="Number of top models to print per section")
    ap.add_argument(
        "--pairs", default=None,
        help=(
            "Comma-separated pair types to rank, e.g. GU,GC,AU. "
            "Omit or use ALL to rank every pair type found in the file."
        ),
    )

    args = ap.parse_args()

    if args.n_top < 1:
        print("ERROR: n_top must be >= 1", file=sys.stderr)
        sys.exit(2)

    keep_pairs = parse_pairs_arg(args.pairs)

    entries = read_candidates(args.input, keep_pairs)
    if not entries:
        print("# No candidates found — check --pairs and input file.", file=sys.stderr)
        sys.exit(0)

    # ── Score function legend ─────────────────────────────────────────────────
    print("# " + "═" * 72)
    print("# Base-pair QM candidate ranking")
    print(f"#   Input         : {args.input}")
    pair_label = ", ".join(sorted(keep_pairs)) if keep_pairs else "ALL"
    print(f"#   Pair types    : {pair_label}")
    print(f"#   Top N         : {args.n_top}")
    print("#")
    print("#   Scoring summary (lower = better):")
    print("#     GU  : 2×|H1-O2−1.8| + 2×|H3-O6−1.8| + angle penalties")
    print("#           + 1×|d_NO−2.9| + 0.2×|C1C1−10.8| + tilt penalty")
    print("#     GC  : 2×best(H21/H22-O2) + 2×best(H41/H42-O6)")
    print("#           + 1×|d_NO−2.9| + 0.2×|C1C1−10.4| + tilt penalty")
    print("#     AU  : 2×best(H61/H62-O4)")
    print("#           + 1×|d_NO−2.9| + 0.2×|C1C1−10.4| + tilt penalty")
    print("#     other: 2×best_HO + 1×|d_NO−2.9| + 0.2×|C1C1−10.6| + tilt")
    print("#")
    print("#   C1C1 weight lowered to 0.2 (was 0.5); unknown C1C1 penalty = 1.0")
    print("#   Angle term (GU only): 0.05×max(0, 170−∠)°")
    print("# " + "═" * 72)
    print()

    # ── Overall top N ─────────────────────────────────────────────────────────
    _print_section(f"Top {args.n_top} overall (all pair types, lowest score)",
                   entries, args.n_top)

    # ── Per pair type ─────────────────────────────────────────────────────────
    by_pair: Dict[str, List[Entry]] = defaultdict(list)
    for e in entries:
        by_pair[e.pair_type].append(e)

    for pt in _CANON_ORDER:
        recs = by_pair.get(pt)
        if not recs:
            continue
        print("─" * 72)
        print(f"# {pt} pairs  (n = {len(recs)})")
        print("─" * 72)
        print()

        # Top N overall within this pair type
        _print_section(f"Top {args.n_top} {pt} overall", recs, args.n_top)

        # Top N within each category
        by_cat: Dict[str, List[Entry]] = defaultdict(list)
        for e in recs:
            by_cat[e.category].append(e)

        for cat, cat_recs in by_cat.items():
            _print_section(f"Top {args.n_top} [{cat}]", cat_recs, args.n_top)


if __name__ == "__main__":
    main()