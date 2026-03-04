# rna-base-contacts

A four-script pipeline to **find**, **score**, **rank**, and **plot** RNA base–base hydrogen-bond contacts in PDB structures, with a focus on selecting high-quality QM model candidates.

---

## Pipeline overview

```
PDB file(s)
    │
    ▼
find_rna_base_contacts_NHO.py
    │  Finds NH···O contacts; writes mini-PDBs + summary files
    │
    ├──► rna_base_NO_summary.txt    (all N···O contacts, unfiltered)
    ├──► rna_base_HO_summary.txt    (H···O contacts within cutoff)
    └──► mini-PDB files             (one 2-residue PDB per pair)
              │
              ▼
    compute_rna_base_miniPDBs.py
              │  Computes geometry + contacts directly from mini-PDBs
              │
              └──► scored_pairs.txt
                        │
                        ▼
             rank_rna_base_contact_miniPDBs.py
                        │  Ranks by composite score; selects top N
                        │
                        └──► ranked_candidates.txt

    rna_base_HO_summary.txt ──► plot_gu_HO.py
                                    │  Plots contact-type frequencies
                                    │  and distance distributions
                                    └──► *.png figures
```

All scripts support **any combination** of the four RNA bases (G, U, C, A) via a `--pairs` argument. GU wobble pairs, GC/AU Watson-Crick pairs, and non-canonical pairs (AG, AC, UC, …) are all handled.

---

## Requirements

| Package | Purpose |
|---|---|
| [NumPy](https://numpy.org/) | Geometry (all scripts) |
| [Biopython](https://biopython.org/) | PDB parsing (`compute_*`, `rank_*`) |
| [Matplotlib](https://matplotlib.org/) | Plotting (`plot_gu_HO.py`) |

Install all dependencies:

```bash
pip install -r requirements.txt
```

> **Note:** `find_rna_base_contacts_NHO.py` only requires NumPy.

---

## Installation

```bash
git clone https://github.com/arghya26/rna-base-contacts.git
cd rna-base-contacts
pip install -r requirements.txt
```

---

## Quick start

### GU wobble pairs only

```bash
# Step 1 — find contacts, write mini-PDBs
python scripts/find_rna_base_contacts_NHO.py \
    --pdbin  1msy.pdb \
    --pairs  GU \
    --cut    3.5 \
    --summary    results/NO_summary.txt \
    --summary_HO results/HO_summary.txt

# Step 2 — score every mini-PDB
python scripts/compute_rna_base_miniPDBs.py \
    --pairs GU --cut 3.5 \
    results/*GU*.pdb > results/scored.txt

# Step 3 — rank and select top 10
python scripts/rank_rna_base_contact_miniPDBs.py \
    results/scored.txt 10 --pairs GU > results/top10_GU.txt

# Step 4 — plot H···O contact statistics
python scripts/plot_gu_HO.py results/HO_summary.txt
```

### Multiple pair types

```bash
python scripts/find_rna_base_contacts_NHO.py \
    --pdbin 1msy.pdb --pairs GU,GC,AU --cut 3.5 \
    --summary results/NO_summary.txt --summary_HO results/HO_summary.txt

python scripts/compute_rna_base_miniPDBs.py \
    --pairs GU,GC,AU --cut 3.5 results/*.pdb > results/scored.txt

python scripts/rank_rna_base_contact_miniPDBs.py \
    results/scored.txt 10 --pairs GU,GC,AU > results/top10.txt

python scripts/plot_gu_HO.py results/HO_summary.txt
```

### All pair types

```bash
python scripts/find_rna_base_contacts_NHO.py \
    --pdbin 1msy.pdb --pairs ALL --cut 3.5 \
    --summary results/NO_summary.txt --summary_HO results/HO_summary.txt

python scripts/compute_rna_base_miniPDBs.py --cut 3.5 results/*.pdb > results/scored.txt

python scripts/rank_rna_base_contact_miniPDBs.py results/scored.txt 10
```

---

## Script reference

### 1. `find_rna_base_contacts_NHO.py`

Scans a PDB file for base pairs with at least one **NH···O** contact within the cutoff. For each qualifying pair it:

- reports all **N···O** contacts (no distance cutoff) to stdout and `--summary`
- reports **H···O** contacts within `--cut` to `--summary_HO`
- writes a **mini-PDB** containing just the two residues

| Argument | Default | Description |
|---|---|---|
| `--pdbin` | *(required)* | Input PDB file |
| `--cut` | *(required)* | H···O distance cutoff (Å) |
| `--pairs` | `GU` | Pair codes: e.g. `GU,GC,AU` or `ALL` |
| `--summary` | None | Append N···O contacts here |
| `--summary_HO` | None | Append H···O contacts here |
| `--nh_names` | auto | Override NH hydrogen atom names |
| `--stack_theta` | 25° | Stacking veto: reject if base-plane angle < this |
| `--stack_dz` | 2.0 Å | Stacking veto: reject if inter-plane separation > this |
| `--no_stack_filter` | off | Disable stacking filter entirely |
| `--dryrun` | off | Skip mini-PDB writing |

**Mini-PDB filename format:**
```
{ResName1}{Chain1}{SeqNum1}_{ResName2}{Chain2}{SeqNum2}_{PDBID}.pdb
# e.g.  GD35_UC26_4pco.pdb
```

---

### 2. `compute_rna_base_miniPDBs.py`

Reads mini-PDB files and computes geometry and contacts entirely from coordinates — no external summary file needed.

**Computed quantities:**

| Field | Description |
|---|---|
| `d_NO` | Shortest base–base N···O distance (both directions, no cutoff) |
| `C1C1` | C1'–C1' distance (Å) |
| `tilt` | Angle between fitted base-plane normals (°) |
| `H···O` | All contacts within `--cut`; H must be present in the file |
| `∠O–H···N` | Linearity angle at H (°); NA when donor heavy atom absent |

**Output categories:**

| Category | Criteria |
|---|---|
| `Canonical_Aform` | GU: H1–O2 < 2.2 Å AND 10.3 ≤ C1C1 ≤ 11.2 AND tilt ≤ 20° |
| `Bifurcated_contact` | GU: H1–O2 present AND (H21–O2 or H22–O2 present) |
| `GU_other` | Remaining GU pairs with ≥ 1 contact |
| `WatsonCrick_GC` | GC: G-amino···O2 < 2.2 Å AND C-amino···O6 < 2.2 Å |
| `WatsonCrick_AU` | AU: A-amino···O4 < 2.2 Å |
| `{type}_other` | Remaining WC pairs with ≥ 1 contact |
| `{type}_general` | Non-canonical pair types (AG, AC, UC, …) |

| Argument | Default | Description |
|---|---|---|
| `pdbs` | *(required)* | Mini-PDB files (positional, glob-friendly) |
| `--cut` | 3.5 Å | H···O search cutoff |
| `--pairs` | ALL | Pair types to process |
| `--wc_cut` | 2.2 Å | Threshold for WC/canonical category assignment |
| `--no_angle` | off | Skip O–H···N angle computation |
| `--verbose` | off | One contact per line with full annotation |

---

### 3. `rank_rna_base_contact_miniPDBs.py`

Reads the output of `compute_rna_base_miniPDBs.py` and ranks entries by a composite score (lower = better QM candidate).

**Scoring functions:**

| Pair | H···O terms | Angle | d_NO | C1C1 | Tilt |
|---|---|---|---|---|---|
| **GU** | 2×\|H1-O2−1.8\| + 2×\|H3-O6−1.8\| | 0.05×max(0,170−∠)° *(GU only)* | 1× | 0.2× | 0.5× |
| **GC** | 2×best(H21/H22-O2) + 2×best(H41/H42-O6) | — | 1× | 0.2× | 0.5× |
| **AU** | 2×best(H61/H62-O4) | — | 1× | 0.2× | 0.5× |
| **other** | 2×best H···O | — | 1× | 0.2× | 0.5× |

Output: **Top N overall** → **Top N per pair type** → **Top N per category**.

| Argument | Default | Description |
|---|---|---|
| `input` | *(required)* | Output file from `compute_rna_base_miniPDBs.py` |
| `n_top` | *(required)* | Number of top models per section |
| `--pairs` | ALL | Filter to specific pair types |

---

### 4. `plot_gu_HO.py`

Reads `rna_base_HO_summary.txt` and produces three figures:

| Figure | Description |
|---|---|
| `HO_type_frequency.png` | Bar chart: count per H–O contact type |
| `HO_distance_histograms_stepstyle.png` | Transparent overlapping histograms per type |
| `HO_distance_histograms_linestyle.png` | Line-style distance distributions per type |

Colors are consistent across all three figures.

```bash
python scripts/plot_gu_HO.py results/HO_summary.txt
```

---

## Supported pair types

| Code | Description |
|---|---|
| `GU` | G–U wobble |
| `GC` | G–C Watson-Crick |
| `AU` | A–U Watson-Crick |
| `AG` | A–G non-canonical |
| `AC` | A–C non-canonical |
| `UC` | U–C non-canonical |
| `GG`, `AA`, `UU`, `CC` | Same-base pairs |
| `ALL` | All of the above |

---

## Notes on hydrogen atoms

H···O distances and O–H···N angles require **hydrogen atoms to be present** in the PDB. Most deposited RNA structures are heavy-atom only. Recommended tools for adding H:

- [MolProbity](http://molprobity.biochem.duke.edu/)
- [PHENIX](https://phenix-online.org/) (`phenix.reduce`)
- [REFMAC](https://www.ccp4.ac.uk/)

When working with heavy-atom-only structures, use `--no_angle` in `compute_rna_base_miniPDBs.py`. N···O distances and geometry (d_NO, C1C1, tilt) are always computed regardless of H availability.

---

## License

MIT License — see [LICENSE](LICENSE).
