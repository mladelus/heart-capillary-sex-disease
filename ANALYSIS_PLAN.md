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

*(none yet)*
