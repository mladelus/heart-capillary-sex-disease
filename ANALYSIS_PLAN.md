# Analysis plan: do female and male heart capillaries differ in disease?

**Author:** Maria Adelus · **Draft, to be time-stamped on GitHub before any expression data are analyzed.** Later changes are listed under *Deviations*.

## Background and question

In healthy adult human hearts, capillary-function programs did not differ substantially between females and males (Adelus, paper 1). Sex differences in coronary microvascular disease may therefore arise in disease rather than at baseline.

**Primary question:** in diseased human ventricles, do female and male capillary endothelial cells differ in capillary-function programs?

**Secondary question:** is any sex difference in disease larger than in health (sex × disease)?

## Sex bias of the available diseases, and expectations (fixed before analysis)

All four diseases available in open single-cell data are more common or more severe in males. The female-predominant microvascular conditions (coronary microvascular dysfunction, HFpEF, INOCA, takotsubo) have no open single-nucleus heart data. This shapes the hypotheses; citations are to be verified before the manuscript is written.

| Disease | Sex pattern (clinical literature) | Expectation for capillaries in females vs males |
|---|---|---|
| DCM | More common in males; females often have better outcomes | Females less affected: a smaller capillary disease response |
| ACM | Males more often and more severely affected | Same as DCM (descriptive only; 3 females) |
| Myocarditis | Strong male predominance, especially in young adults | Females show less inflammatory/interferon activation in capillaries (descriptive only; 3 females) |
| MI | More common in males; females have worse microvascular outcomes (more microvascular obstruction and MINOCA) | Females show more capillary dysfunction (descriptive only; 4 females) |

**Directional secondary hypothesis (DCM, the only well-powered stratum):** in the same-consortium sex × disease model, the DCM-related change in capillary programs is smaller in females than in males. This is tested as the sex × disease interaction, two-sided, and reported with its direction.

## Data (all open; CELLxGENE Census release 2025-11-08)

- **Donors:** adults (18 years or older) with recorded sex and a non-normal disease label.
- **Tissue:** ventricular tissue only (same labels as paper 1).
- **Minimum cells:** at least 100 endothelial cells per donor. Only primary data are used.
- **Analysis strata:** one stratum per dataset × disease with at least 3 eligible females and 3 eligible males. The metadata inventory (run before this plan) found four:

| Stratum | Source | Females | Males |
|---|---|---|---|
| Dilated cardiomyopathy (DCM) | Reichart et al. 2022 | 12 | 34 |
| Arrhythmogenic cardiomyopathy (ACM) | Reichart et al. 2022 | 3 | 5 |
| Myocardial infarction (MI) | Kuppe et al. 2022 | 4 | 12 |
| Myocarditis | Census myocarditis dataset | 3 | 20 |
| **Total** | | **22** | **71** |

- **Healthy comparison:** the 45 healthy donors of paper 1 (23 females, 22 males), plus any healthy donors in the disease datasets.

## Methods (identical to paper 1 unless stated)

- **Cell definitions:** same broad cell classes; the same capillary, arterial and venous marker rule, applied within each dataset; pseudobulk minimums of at least 30 cells.
- **Outcome programs:** the six paper 1 programs (NO/eNOS, endothelin/ACE, prostacyclin, barrier/junction, fatty-acid transport, angiogenic tip), with the same genes. Scores are mean z-scores, standardized within each stratum.
- **Positive control:** the X-escape program.

## Primary analysis

- **Model:** within each stratum, score ~ sex + age (+ assay when it varies). The female-minus-male estimates are pooled across strata by random-effects meta-analysis (REML), since 4 strata are available.
- **Multiplicity:** Holm correction across the scorable programs.
- **Equivalence:** bounds of ±0.8 SD, as in paper 1. Whether the 90% CI also lies within ±0.5 SD is reported.
- **Verdicts:** *difference*, *equivalent* or *inconclusive*, by the paper 1 rules.
- **Calibration:** 1,000 sex-label permutations within stratum × assay.
- **Sensitivity analyses:**
  - adding the stress score and the ambient heart-muscle RNA score as covariates;
  - DCM alone;
  - leaving each stratum out in turn.

## Secondary analyses

1. **Sex × disease, same consortium.** This uses Reichart DCM together with the Heart Cell Atlas healthy donors (Litviňuková 2020), which come from the same consortium and snRNA-seq protocol.
   - **Model:** score ~ sex × disease + age + assay.
   - **Estimand:** the interaction, i.e., the change in the female-minus-male difference from health to DCM.
   - This analysis is reported as supportive, because disease and study are partly confounded.
2. **Sex × disease, all data.** The pooled disease estimate (primary) is compared with the pooled healthy estimate from paper 1 by a z-test of the difference. It is labeled cross-study.
3. **Pericytes:** the pericyte-to-capillary ratio and the five paper 1 pericyte programs, female minus male within disease.
4. **X–Y dosage:** KDM6A/UTY and KDM5C/KDM5D (X alone and X + Y) in capillary endothelium, female minus male within disease, compared with health.
5. **Disease positive control:** DCM vs healthy capillaries (same-consortium model). It must show clear disease effects, such as a significant shift in at least one program or many genes at FDR < 0.05. This shows that the pipeline detects biology in these data.

## What counts as a result

- **Positive:** the pooled Holm p < 0.05 with the same direction in most strata, and calibration passes. If the secondary sex × disease interaction also has p < 0.05, this is described as "sex differences emerge in disease".
- **Informative null:** equivalence (90% CI within ±0.8 SD).
- **Inconclusive:** anything else.

## Known limitations (stated in advance)

- Few females: 22 in total, 12 with DCM.
- All available diseases are male-biased, so results cannot be generalized to the female-predominant microvascular diseases, for which no open single-nucleus data exist.
- Disease severity, genetic cause and treatment vary and are not recorded.
- Health vs disease comparisons partly cross studies.
- Ages are partly recorded by decade.
- Menopause and hormone data are not available.

## Deviations

1. **Mural cells instead of pericytes (made before any pericyte or sex result was seen).** After download, the cell labels showed that two of the three datasets (DCM/ACM atlas; myocarditis) label pericytes and smooth muscle cells together as "mural cell", and only the myocardial infarction dataset labels pericytes separately. Using the "pericyte" label alone would have given most donors zero pericytes. For aim 3 we therefore use one rule for every dataset: mural cells = pericyte + mural cell + smooth muscle labels (cardiac mural cells are mostly pericytes). These cells are downloaded by a new script, `02b_download_mural.R`. The ratio becomes mural cells per capillary EC, and the five pericyte programs are scored in the mural-cell pseudobulk. The primary analysis (aims 1 and 2) is unchanged.

2. **Bug fix (no change to any estimate):** in script 04 the leave-one-stratum-out rows were labeled with the number of strata ("without 3") instead of the stratum left out, because the loop variable shared its name with a column. The loop variable was renamed.


## Stage 3 — where do sex differences lie? (added after the stage 1 and stage 2 results; committed before any stage 3 analysis was run)

Stages 1 (healthy, paper 1) and 2 (disease) found no large sex differences in the six pre-chosen capillary programs, but a reproducible X–Y dosage pattern. Stage 3 asks where sex differences in the heart are, without restricting to pre-chosen programs. Everything below is exploratory; the healthy cohort (45 donors) is the discovery set and the disease cohort (93 donors) the replication set, and vice versa. No new data are downloaded.

**A. Discovery and replication by cell type (script 08).** Cell types: capillary, arterial and venous endothelium; mural cells (pericyte label in stage 1, mural label in stage 2); fibroblasts; cardiomyocytes; myeloid and lymphoid cells. In each cohort and cell type, donor pseudobulks are analysed with limma-voom (expression ~ sex + age [+ assay]) within each dataset (stage 1) or disease (stage 2), and combined by inverse-variance meta-analysis. Genes are classed as autosomal, X or Y (EnsDb.Hsapiens.v86). Measures of autosomal sex signal per cell type, reported side by side: (i) the Spearman correlation of gene-level sex z-scores between the two cohorts; (ii) the number of autosomal genes at FDR < 0.05 in one cohort that replicate (p < 0.05, same direction) in the other, with 2.5% expected by chance. Pathways: Hallmark gene sets (MSigDB) tested on the combined (Stouffer) autosomal z-scores with a rank-based gene-set test, Benjamini–Hochberg across sets, per cell type. Expectation stated in advance: the autosomal sex signal in capillary endothelium is no larger than in other heart cell types.

**B. Sex-hormone receptor map (script 09).** ESR1, ESR2, GPER1, AR, PGR and CYP19A1: expression (log2 CPM) by cell type and sex in both cohorts; in capillaries, female − male by age band (< 50, ≥ 50 years); and Hallmark estrogen-response (early, late) and androgen-response scores in capillaries, female − male, combined across both cohorts.

**C. Capillary cell states (script 10).** Every capillary EC is scored for stress, interferon response, inflammatory activation, angiogenic tip and proliferation (mean log(1 + CP10k) of fixed gene lists). A cell is "high" for a state if its score is in the top 10% of capillary cells of its dataset; the outcome is each donor's logit fraction of high cells. The capillary-to-arterial/venous position of each capillary cell (z_venous − z_arterial, from script 03) is summarized per donor as its mean. Female − male with the primary methods (per dataset or disease, random-effects meta-analysis), in both cohorts and combined.

**D. Donor-to-donor variability (script 10).** For each program and the X-escape control, log(SD in females / SD in males) of the age-adjusted score within each dataset or disease, combined by meta-analysis; and the mean transcriptome distance (1 − Pearson correlation, 2,000 most variable autosomal genes) between donors of the same sex, females vs males, tested by sex-label permutation within dataset.

**Combined health + disease summaries (script 11).** Female − male estimates of the same outcome from both cohorts are combined by random-effects meta-analysis across all six strata (two healthy datasets, four disease strata), with disease state as a moderator. Figures for stage 3.

**E. Immune cells (script 12; added after stage 3 A–D results, committed before this analysis was run).** Stage 3 showed higher estrogen-response gene expression in females in myeloid and lymphoid cells and fibroblasts, the cell types with the highest ESR1 expression, but no gene-level sex differences replicated in immune cells. Question: do immune cells differ by sex in pre-chosen immune programs? Same methods as the primary analysis (donor pseudobulk, score standardized within dataset or disease, female − male adjusted for age and assay, random-effects meta-analysis per cohort and across all six strata, Holm across programs within a cell type, X-escape positive control).
Myeloid programs: type I interferon response (IFIT1, IFIT3, ISG15, MX1, IFI6, IFI44L, OAS1, STAT1, IRF7, XAF1); MHC class II antigen presentation (HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, CD74, CIITA); inflammatory (IL1B, TNF, NLRP3, CCL3, CCL4, CXCL8, NFKBIA, PTGS2); resident macrophage (LYVE1, F13A1, MRC1, CD163, FOLR2, STAB1); recruited monocyte-derived (CCR2, PLAC8, S100A8, S100A9, VCAN, FCN1); complement/efferocytosis (C1QA, C1QB, C1QC, MERTK, TREM2, APOE).
Lymphoid programs: cytotoxic (GZMB, GZMK, GZMA, PRF1, NKG7, GNLY, CCL5); type I interferon response (as above); B cell (MS4A1, CD79A, CD79B, BANK1).
Also: (i) immune composition per donor (myeloid and lymphoid share of all cells, logit); (ii) immune X-linked genes reported to escape X inactivation (TLR7, TLR8, CXCR3, CD40LG, BTK, IL2RG, CYBB), female − male from script 08; (iii) donor-level Hallmark estrogen-response (early) and oxidative-phosphorylation scores in every cell type, female − male, with and without a heart-muscle ambient RNA score as covariate, to test whether the estrogen signal tracks ESR1 expression and whether the male-higher oxidative-phosphorylation signal is contamination. Limitation known in advance: the myocardial infarction dataset labels many lymphoid cells with terms outside our lymphoid class, so lymphoid results there are incomplete.
