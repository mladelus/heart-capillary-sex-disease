# =============================================================================
# 10  STAGE 3C + 3D (exploratory)
# C. Capillary cell states: each capillary EC is scored for stress, interferon response, inflammatory activation,
#    angiogenic tip and proliferation (mean log(1 + CP10k)). "High" = top 10% of capillary cells in its dataset
#    (proliferation: any MKI67/TOP2A/CENPF count, because most cells have none). Outcome per donor: logit fraction
#    of high cells. Also the capillary's position on the arterial-venous axis (z_venous - z_arterial, script 03).
# D. Donor-to-donor variability: log(SD females / SD males) of age-adjusted program scores; and the mean
#    transcriptome distance between donors of the same sex (females minus males), permutation test.
# Female - male within each dataset/disease, random-effects meta per cohort and across both cohorts.
# Makes: tables/S3C_capillary_states.csv, tables/S3D_variability.csv, stage3_states_var.rds
# =============================================================================
set.seed(2030)
STATES <- list(
  stress        = STRESS,
  interferon    = c("IFIT1", "IFIT3", "ISG15", "MX1", "IFI6", "IFI44L", "OAS1", "STAT1", "IFITM1", "XAF1"),
  inflammatory  = c("SELE", "VCAM1", "ICAM1", "CXCL2", "CCL2", "IL6"),
  angiogenic    = programs$Angiogenic_tip,
  proliferating = c("MKI67", "TOP2A", "CENPF"))
C <- load_cohorts()
dirs <- c(healthy = PAPER1_DIR, disease = getwd())

# ---- C. per-cell state scores (one donor file at a time) ----
cells <- bind_rows(lapply(names(dirs), function(co) {
  ec <- readRDS(file.path(dirs[[co]], "ec_subtypes.rds"))
  S  <- C[[co]]$samples
  bind_rows(lapply(list.files(file.path(dirs[[co]], "dl"), full.names = TRUE), function(f) {
    x <- readRDS(f); s <- paste(substr(x$dataset_id, 1, 8), x$donor_id, sep = "|")
    if (!s %in% S$sample) return(NULL)
    a <- ec |> filter(dataset_id == x$dataset_id, donor_id == x$donor_id, subtype == "capillary")
    if (!nrow(a)) return(NULL)
    m <- x$ec_counts; lib <- Matrix::colSums(m); m <- m[, a$soma_joinid, drop = FALSE]; lib <- lib[a$soma_joinid]
    sc <- sapply(STATES, function(g) { g <- intersect(g, rownames(m)); if (!length(g)) return(rep(0, ncol(m)))
      Matrix::colMeans(log1p(t(t(m[g, , drop = FALSE]) / lib) * 1e4)) })
    if (is.null(dim(sc))) sc <- matrix(sc, nrow = 1, dimnames = list(NULL, names(STATES)))
    tibble(cohort = co, sample = s, dataset_id = x$dataset_id, axis = a$z_venous - a$z_arterial,
           log_umi = log(as.numeric(lib))) |> bind_cols(as_tibble(sc))
  }))
}))
cells <- cells |> group_by(cohort, dataset_id) |>
  mutate(across(c(stress, interferon, inflammatory, angiogenic),           # random tie-breaking, so exactly ~10% are "high"
                ~ percent_rank(.x + runif(length(.x), 0, 1e-9)) > 0.9, .names = "hi_{.col}"),
         hi_proliferating = proliferating > 0) |> ungroup()
don <- cells |> group_by(cohort, sample) |>
  summarise(n = n(), across(starts_with("hi_"), ~ log((sum(.x) + 0.5) / (n() - sum(.x) + 0.5))),
            axis = mean(axis), log_umi = median(log_umi), .groups = "drop") |>
  left_join(bind_rows(core_cols(C$healthy$samples), core_cols(C$disease$samples)) |> select(sample, stratum, sex, age, assay),
            by = "sample") |>
  filter(n >= MIN_CAP)
outs_c <- c(paste0("hi_", names(STATES)), "axis")
meta_cohorts <- function(dat, outs, extra = NULL) bind_rows(lapply(outs, function(v) {
  per <- bind_rows(lapply(split(dat, dat$stratum), function(d) {
    r <- sex_fit(d, v, extra); if (is.null(r)) NULL else mutate(r, stratum = d$stratum[1], cohort = d$cohort[1]) }))
  if (!nrow(per)) return(NULL)
  bind_rows(per |> group_by(cohort) |> group_modify(~ meta_one(.x)) |> ungroup(),
            meta_one(per) |> mutate(cohort = "both")) |> mutate(outcome = v)
}))
st <- bind_rows(meta_cohorts(don, outs_c) |> mutate(adjusted = "none"),
                meta_cohorts(don, outs_c, "log_umi") |> mutate(adjusted = "+ median log UMI per cell"))
cat("\nRealised % of capillary cells flagged 'high', by cohort (should be about 10, except proliferating):\n")
cells |> group_by(cohort) |> summarise(across(starts_with("hi_"), ~ round(100 * mean(.x), 1))) |> print(width = Inf)
cat("\n== C. Capillary cell states: female - male (SD units of logit fraction; axis = venous-minus-arterial position) ==\n")
st |> select(adjusted, outcome, cohort, k, n_w, n_m, est, lo, hi, p, I2) |> r3() |> arrange(adjusted, outcome, cohort) |> print(n = Inf, width = Inf)
cat("\nMedian % of capillary cells 'high' per donor, by sex:\n")
cells |> left_join(don |> select(sample, sex), by = "sample") |> filter(!is.na(sex)) |>
  group_by(sample, sex) |> summarise(across(starts_with("hi_"), ~ 100 * mean(.x)), .groups = "drop") |>
  group_by(sex) |> summarise(across(starts_with("hi_"), median)) |> r3() |> print(width = Inf)

# ---- D1. variability of program scores ----
sc_h <- readRDS(file.path(PAPER1_DIR, "aim1.rds"))$scores |> mutate(stratum = dataset_id, cohort = "healthy")
sc_d <- readRDS("p2_primary.rds")$scores |> mutate(cohort = "disease")
vprog <- c(PRIMARY, "X_escape")
lvr <- bind_rows(lapply(vprog, function(v) bind_rows(lapply(list(sc_h, sc_d), function(sc) bind_rows(lapply(split(sc, sc$stratum), function(d) {
  d <- d[is.finite(d[[v]]) & !is.na(d$age), ]
  nf <- sum(d$sex == "female"); nm <- sum(d$sex == "male"); if (nf < MIN_PER_SEX || nm < MIN_PER_SEX) return(NULL)
  r <- resid(lm(as.formula(paste(v, "~ age")), data = d))
  tibble(program = v, stratum = d$stratum[1], cohort = d$cohort[1], n_w = nf, n_m = nm,
         est = log(sd(r[d$sex == "female"]) / sd(r[d$sex == "male"])), se = sqrt(1 / (2 * (nf - 1)) + 1 / (2 * (nm - 1))))
}))))))
var_res <- bind_rows(lvr |> group_by(program, cohort) |> group_modify(~ meta_one(.x)) |> ungroup(),
                     lvr |> group_by(program) |> group_modify(~ meta_one(.x)) |> ungroup() |> mutate(cohort = "both")) |>
  mutate(SD_ratio = exp(est), SD_ratio_lo = exp(lo), SD_ratio_hi = exp(hi))
cat("\n== D. Variability: SD in females / SD in males of age-adjusted program scores ==\n")
var_res |> select(program, cohort, k, SD_ratio, SD_ratio_lo, SD_ratio_hi, p, I2) |> r3() |> arrange(program, cohort) |>
  print(n = Inf, width = Inf)

# ---- D2. transcriptome-wide distance between donors of the same sex ----
ch <- chrom_map()
dist_test <- function(co, nperm = 1000) {
  m <- C[[co]]$mats$capillary; S <- core_cols(C[[co]]$samples) |> filter(sample %in% colnames(m), n_capillary >= MIN_CAP, !is.na(age))
  y <- DGEList(m[, S$sample]); y <- calcNormFactors(y[filterByExpr(y, group = S$stratum), , keep.lib.sizes = FALSE])
  lc <- cpm(y, log = TRUE, prior.count = 1)
  lc <- removeBatchEffect(lc, batch = S$stratum, batch2 = if (n_distinct(S$assay) > 1) S$assay else NULL,
                          covariates = cbind(S$age, log(S$n_capillary)), design = model.matrix(~ sex, S))
  au <- intersect(rownames(lc), ch$gene[ch$chr_class == "autosome"])
  top <- au[order(apply(lc[au, ], 1, var), decreasing = TRUE)][1:min(2000, length(au))]
  D <- 1 - cor(lc[top, ])
  stat <- function(sex) {                     # mean within-sex distance, same stratum pairs only
    f <- c(); mm <- c()
    for (s in unique(S$stratum)) { i <- which(S$stratum == s)
      fi <- i[sex[i] == "female"]; mi <- i[sex[i] == "male"]
      if (length(fi) > 1) f <- c(f, D[fi, fi][upper.tri(D[fi, fi])]); if (length(mi) > 1) mm <- c(mm, D[mi, mi][upper.tri(D[mi, mi])]) }
    mean(f) - mean(mm)
  }
  obs <- stat(as.character(S$sex))
  perm <- replicate(nperm, { sx <- as.character(S$sex); for (s in unique(S$stratum)) { i <- which(S$stratum == s); sx[i] <- sample(sx[i]) }; stat(sx) })
  tibble(cohort = co, genes = length(top), median_cap_F = median(S$n_capillary[S$sex == "female"]),
         median_cap_M = median(S$n_capillary[S$sex == "male"]), females_minus_males = obs, p_perm = (1 + sum(abs(perm) >= abs(obs))) / (1 + nperm))
}
dist_res <- bind_rows(dist_test("healthy"), dist_test("disease"))
cat("\n== D. Mean distance between donors of the same sex, females minus males (positive = females more variable) ==\n")
dist_res |> r3() |> print(width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(st, "tables/S3C_capillary_states.csv", row.names = FALSE)
write.csv(bind_rows(var_res, dist_res |> mutate(program = "transcriptome distance")), "tables/S3D_variability.csv", row.names = FALSE)
saveRDS(list(states = st, donors = don, variability = var_res, distance = dist_res), "stage3_states_var.rds")
