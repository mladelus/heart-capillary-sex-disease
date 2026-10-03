# Endothelial state defined without EDN1 (analysis plan, section M)

Python (numpy, pandas, scipy). Reads the per-donor count files saved by the R scripts; run R scripts 00-19 first.

```
python cache_cells.py      # one-time: converts dl/*.rds to fast .npz files
python m1_signature.py     # M1: cross-fitted signature of EDN1-positive capillary cells; per-cell state scores
python m2_state_tests.py   # M2-M4: abundance and intensity by sex, validation, continuum, strict capillary definitions
python m5_genetics.py      # M5: GWAS Catalog test of the state genes (needs gene_positions.csv, see the script header)
```

Data folders are `~/Documents/heart_capillary_project` and `~/Documents/heart_capillary_disease` (set `HCAP_DATA` and `HCAP2_DATA` to change). `common.py` re-implements the donor-level model and the meta-analysis of the R pipeline; it reproduces the R estimates exactly (checked on the share of EDN1-positive cells). Results are written to `state/out/`.
