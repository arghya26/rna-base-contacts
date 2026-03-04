# Changelog

## [1.0.0] — 2026-03

Initial release. Four-script pipeline for finding, scoring, ranking, and
plotting RNA base–base NH···O contacts in PDB structures.

### Scripts

- **`find_rna_base_contacts_NHO.py`**
  - Finds H···O contacts for any combination of RNA bases (G, U, C, A)
    via `--pairs` (e.g. `GU`, `GU,GC,AU`, or `ALL`).
  - Stacking filter (base-plane angle + inter-plane separation) removes
    stacked pairs that are not genuine H-bond partners.
  - Writes per-pair mini-PDB files and two summary tables:
    N···O contacts (unfiltered) and H···O contacts (within `--cut`).
  - Mini-PDB filename format: `{Res1}{Chain1}{Seq1}_{Res2}{Chain2}{Seq2}_{PDBID}.pdb`

- **`compute_rna_base_miniPDBs.py`**
  - Self-contained: all geometry and contacts computed directly from
    mini-PDB coordinates — no external summary file required.
  - Computes d_NO, C1C1, tilt, H···O distances, and O–H···N angles.
  - Categorises pairs: `Canonical_Aform`, `Bifurcated_contact` (GU);
    `WatsonCrick_GC`, `WatsonCrick_AU`; `{type}_other`; `{type}_general`.
  - Supports `--pairs`, `--cut`, `--wc_cut`, `--no_angle`, `--verbose`.

- **`rank_rna_base_contact_miniPDBs.py`**
  - Reads output of `compute_rna_base_miniPDBs.py`.
  - Pair-specific composite scoring functions:
    - GU: distance + angle terms (O–H···N linearity), d_NO, C1C1, tilt.
    - GC: two Watson-Crick amino contacts, d_NO, C1C1, tilt.
    - AU: one Watson-Crick amino contact, d_NO, C1C1, tilt.
    - Other: best H···O, d_NO, C1C1, tilt.
  - C1C1 weight 0.2 (lower than distance terms); unknown C1C1 penalty 1.0.
  - Outputs: Top N overall → Top N per pair type → Top N per category.

- **`plot_gu_HO.py`**
  - Reads `rna_base_HO_summary.txt` produced by `find_rna_base_contacts_NHO.py`.
  - Produces three figures with consistent per-type colors:
    1. `HO_type_frequency.png` — bar chart of contact-type counts.
    2. `HO_distance_histograms_stepstyle.png` — overlapping transparent histograms.
    3. `HO_distance_histograms_linestyle.png` — line-style distance distributions.
