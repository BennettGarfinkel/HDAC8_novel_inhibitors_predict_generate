# hHDAC8 Non-Hydroxamate Inhibitor Discovery Pipeline

Computational drug discovery pipeline for identifying potent, drug-like, isoform-selective
non-hydroxamate inhibitors of human HDAC8 (hHDAC8), using ML-driven QSAR/docking-score
prediction and genetic-algorithm-based generative chemistry with built-in docking-score
selectivity for HDAC8 over HDAC1/HDAC6. 

## Project constraints and goals

- **Hard constraint: no hydroxamic acid zinc-binding groups (ZBGs).** All other
  non-hydroxamate ZBG chemotypes are open for exploration.
- **Goal:** identify candidate molecules with real potency (target ≤150-500 nM depending
  on triage stage), acceptable drug-likeness (Lipinski/Veber), synthetic accessibility,
  no reactive/toxic liabilities, and evidence of HDAC8 selectivity over HDAC1/HDAC6 —
  both by predicted IC50 and by predicted docking score.
- **Therapeutic rationale:** HDAC8-selective inhibition avoids the toxicity profile of
  pan-HDAC inhibitors; relevant to pediatric T-cell ALL, neuroblastoma, Cornelia de Lange
  syndrome (SMC3/cohesin deacetylation), and antiparasitic applications.

## Architecture

Cap–linker–ZBG (L-shaped) design paradigm. Key PDB
structures: 1T64 (HDAC8), 4BKX (HDAC1), 5EEM (HDAC6 CD2).

## Pipeline components

| File | Purpose |
|---|---|
| `scripts/pipeline.py` | ZBG whitelist + Tier-1 hard gates, liability engine, data cleaning, IC50/docking model training (`train_model_from_df` generic trainer), single/batch SMILES prediction, ZBG-precedent tracking, potency tiering, 2D L-shape proxy |
| `scripts/ga.py` | ZBG-locked mutation/crossover, NSGA-II Pareto GA with docking-selectivity objective + threshold gate, IC50-selectivity flag, lineage-protected elite selection |
| `scripts/data_prep.py` | Consolidates the 73k-row merged dataset into clean per-target training sets |
| `scripts/train_all_models.py` | **Optional retrain entry point**: retrains all six models from `data/` and rebuilds `models/run_state.pkl` |
| `scripts/run_docking_selective.py` | Generation run: reuses trained models, applies docking gate + IC50 selectivity, exports candidates |
| `scripts/make_figures.py` | Generates all nine figures from the candidate CSVs + model bundles |
| `notebooks/hHDAC8_docking_selective.ipynb` | Full runnable notebook: (optional retrain) → load models → (optional regenerate) GA → assembly → figures → `predict_from_smiles()` |
| `models/run_state.pkl` | Trained model bundles (HDAC8 IC50/docking, HDAC1/HDAC6 IC50/docking) + `ic50_df` |
| `models/docking_selective_run_state.pkl` | Full state from the generation run (hits_df + all model bundles) |
| `candidates/*.csv` | GA output: raw hits, balanced primary deliverable, cluster representatives, strict clean screen |
| `figures/*.png` | Nine diagnostic visualizations (fig1–fig9) |

## Docking + IC50 selectivity gate

The GA uses a **docking-score-based isoform selectivity** system:

1. **Docking gate (hard, defines a qualifying hit):** HDAC8 predicted docking must be
   **≤ -7** (strong on-target binding, and not below a -13 applicability-domain floor);
   HDAC1 and HDAC6 predicted docking must each be **> -7** (weak off-target binding). This
   is a clean directional threshold per isoform, replacing the earlier ±1.5 tolerance band.

2. **IC50 selectivity (reported + preferred, not a hard hit gate):** off-target predicted
   IC50 ≥ 1000 nM and HDAC8 ≥ 0.7 log more potent than each off-target sets a
   `passes_ic50_selectivity` flag. It is **not** a hard gate — as a hard gate it collapsed
   yield to ~1 hit/lineage, because HDAC8's own median training IC50 is ~920 nM. Instead the
   existing `selectivity` Pareto objective pushes the search toward IC50-selective molecules,
   and the flag is enforced only where a fully-selective subset is wanted.

3. **Pareto objective:** `docking_selectivity = min(dock_HDAC1, dock_HDAC6) - dock_HDAC8`
   is maximized alongside pIC50, HDAC8 docking, liability score, and IC50-based selectivity.

The docking gate is deliberately not a hard rejection *inside* the GA loop — that would
freeze lineages whose seeds start outside the target region (confirmed bug, found and fixed).
It gates only what counts as a returned hit; the Pareto objective drives the population toward
selectivity across generations.

## ZBG whitelist (safety-audited)

| ZBG | Status | Rationale |
|---|---|---|
| Hydroxamate | Hard-excluded | Project constraint |
| Hydrazide | Removed | N-N bond, oxidative cleavage, hydrazine-release precedent |
| Thiol (free) | Removed | Oxidizes to reactive disulfide |
| Benzisothiazolone | Removed | Ring-opens to reactive electrophile |
| 8-Hydroxyquinoline | Removed | Promiscuous non-Zn metal chelation; neurotoxicity precedent |
| Salicylamide | Kept | Established HDAC8 ZBG |
| Ortho-aminoanilide | Kept | Clinically precedented (entinostat/mocetinostat class) |
| Acylurea | Kept | Current FragBreed lead chemotype |
| 3-HPT | Kept, flagged | Pyrithione-adjacent risk; open SAR area |
| Carboxylate | Added | Metabolically inert, weak/reversible chelator |
| Trifluoromethyl-ketone | Added | Reversible hydrate-forming ZBG |
| Cyclic-thione | Added | Generalizes 3-HPT; 72+ real examples |
| Triazolopyridine | Added | Published non-hydroxamate HDAC8 ZBG |
| Triazolo[4,3-a]quinoline | Added | Published non-hydroxamate HDAC8 ZBG |

## Data

- `data/HDAC_Docking_Inhibition.csv` — merged 73,467-row dataset
- `data/hdac8_ic50_clean_MERGED.csv` — 6,082 unique HDAC8 IC50 compounds
- `data/hdac1_ic50_clean.csv` / `data/hdac6_ic50_clean.csv` — 10,340 / 9,035 off-target IC50
- `data/hdac1_dock_clean.csv` / `data/hdac6_dock_clean.csv` — 946 / 1,003 off-target docking

Model performance (scaffold-grouped CV R²): HDAC8 IC50 0.54, HDAC8 docking 0.43, HDAC1
IC50 0.51, HDAC6 IC50 0.50. HDAC1/HDAC6 docking models are modest (0.13/0.15) due to
small training sets with high scaffold diversity — treat as directional, not precise.

## Results (current run)

291 total docking-gated hits across 8 of 9 ZBG lineages (Triazoloquinoline has zero
seeds). The strict all-constraints screen (no PAINS/Brenk, Lipinski/Veber, SA ≤ 3.5,
HDAC8 dock ≤ -7, HDAC1/HDAC6 dock > -7, IC50 ≤ 500 nM) yields **6 clean candidates**
(5 Carboxylate, 1 Triazolopyridine). Each output CSV tracks MW, cLogP, HBD, HBA, TPSA,
RotB, both off-target IC50/docking predictions, and the `passes_ic50_selectivity` flag.

## Output columns

Every candidate CSV includes: `SMILES`, `zbg_tag`, `pIC50_pred`, `IC50_nM_est`,
`IC50_tier`, `docking_pred`, `docking_hdac1_pred`, `docking_hdac6_pred`,
`docking_selectivity`, `pIC50_hdac1_pred`, `pIC50_hdac6_pred`, `selectivity_vs_hdac1/6`,
`passes_ic50_selectivity`, `MW`, `cLogP`, `HBD`, `HBA`, `TPSA`, `RotB`, `sa_score`,
`pains_brenk_flagged`, plus applicability-domain and precedent flags.

## Running the notebook

```bash
cd notebooks
jupyter notebook hHDAC8_docking_selective.ipynb
# Run All
```

The notebook has two optional switches near the top:

- **`RETRAIN`** (default `False`): retrain all six models from `../data/` and rebuild
  `../models/run_state.pkl` (a few minutes). Off by default — the shipped models are used.
- **`REGENERATE`** (default `False`): re-run the full two-stage island GA (~25 min). Off by
  default — the cached hit pool in `../models/docking_selective_run_state.pkl` is loaded so
  "Run All" completes in seconds and reproduces the shipped results.

With both `False`, the notebook loads cached models + hits, rebuilds all deliverables and
all nine figures, and defines `predict_from_smiles()`. It reads `../models/`, imports from
`../scripts/`, and writes to `../candidates/` and `../figures/`.

Retraining or regenerating from the command line instead:

```bash
cd scripts
python3 train_all_models.py        # rebuild ../models/run_state.pkl from ../data
python3 run_docking_selective.py   # re-run the GA and export candidates
python3 make_figures.py            # regenerate ../figures from the candidate CSVs
```

## Known limitations

- Docking scores are **ML-predicted** (not live Glide/Maestro). Cross-docking against
  4BKX (HDAC1) / 5EEM (HDAC6) in Maestro and MM-GBSA rescoring recommended before
  wet-lab prioritization.
- HDAC1/HDAC6 docking models have modest scaffold-CV R² (0.13/0.15): small training sets
  (946/1,003 compounds) with ~70% singleton scaffolds and a compressed score range.
  Treat off-target docking numbers as directional, not precise.
- None of the 6 clean candidates also clear the strict IC50-selectivity margin; that is
  reported honestly in the `passes_ic50_selectivity` column rather than forced.
- Triazoloquinoline has zero seeds and produced zero hits.
