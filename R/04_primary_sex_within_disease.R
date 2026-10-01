# =============================================================================
# 04  PRIMARY - female vs male capillary ECs within diseased hearts
# Per disease stratum: capillary pseudobulk -> TMM log-CPM -> program score (mean z) standardized within stratum
#   -> lm(score ~ sex + age [+ assay]) -> female - male in SD units.
# Strata pooled by random-effects meta-analysis (REML). Holm across scorable programs; equivalence +/- 0.8 SD.
# Checks: X-escape positive control; 1,000 sex-label permutations within stratum x assay;
#         + stress + ambient heart-muscle RNA covariates; DCM alone; leave one stratum out.
# Makes: p2_primary.rds, tables/P2_Table_primary.csv, tables/P2_TableS_*.csv
# =============================================================================
set.seed(2029)
P   <- readRDS("pb.rds")
cap <- P$mats$capillary
S   <- P$samples |> filter(n_capillary >= MIN_CAP, sample %in% colnames(cap), !grepl("^normal", stratum))
all_sets <- c(programs, pos_control, list(Stress = STRESS, CM_ambient = CM_AMBIENT))

by_st <- split(S, S$stratum)
by_st <- by_st[vapply(by_st, function(d) sum(d$sex == "female") >= MIN_PER_SEX && sum(d$sex == "male") >= MIN_PER_SEX, logical(1))]
if (!length(by_st)) stop("No disease stratum has enough females and males")
cat("\n== Strata entering the primary analysis ==\n")
bind_rows(lapply(by_st, function(d) tibble(stratum = d$stratum[1], females = sum(d$sex == "female"),
  males = sum(d$sex == "male"), median_cap = median(d$n_capillary), assays = paste(unique(d$assay), collapse = " | ")))) |>
  print(width = Inf)

lc_of <- function(d) {
  y <- DGEList(cap[, d$sample, drop = FALSE])
  y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE])
  cpm(y, log = TRUE, prior.count = 1)
}
lc_st  <- lapply(by_st, lc_of)
scores <- bind_rows(lapply(names(by_st), function(k)
  bind_cols(by_st[[k]], as_tibble(sapply(all_sets, function(g) prog_score(lc_st[[k]], g))))))

cat("\n== Program genes passing expression filter, per stratum ==\n")
print(sapply(lc_st, function(lc) sapply(all_sets, function(g) paste0(sum(g %in% rownames(lc)), "/", length(g)))))

fit_all <- function(sc, sets, extra = NULL) bind_rows(lapply(sets, function(prog) bind_rows(lapply(
  split(sc, sc$stratum), function(d) {
    r <- sex_fit(d, prog, extra); if (is.null(r)) NULL else mutate(r, stratum = d$stratum[1], program = prog)
  }))))
run_meta <- function(per) per |> group_by(program) |> group_modify(~ meta_one(.x)) |> ungroup()

per_st <- fit_all(scores, names(all_sets))
scorable <- intersect(PRIMARY, unique(per_st$program))
main <- run_meta(per_st) |>
  mutate(p_holm = ifelse(program %in% scorable, p.adjust(ifelse(program %in% scorable, p, NA), "holm"), NA),
         MDE80 = (qnorm(0.975) + qnorm(0.8)) * (hi - lo) / (2 * qnorm(0.975)),
         verdict = ifelse(program %in% scorable, verdict(p_holm, lo90, hi90, same_dir, k), NA))
xe <- main |> filter(program == "X_escape")
if (!nrow(xe) || !(xe$est > 0 && xe$p < 0.05)) warning("POSITIVE CONTROL FAILED: X-escape not higher in females")

r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))
cat("\n== PRIMARY: female - male in diseased hearts, SD units (random-effects meta across strata) ==\n")
main |> select(program, k, n_w, n_m, est, lo, hi, p, p_holm, lo90, hi90, within_0.5, pi_lo, pi_hi, I2, same_dir,
               MDE80, verdict) |> r3() |> print(n = Inf, width = Inf)
cat("\n== Per-stratum estimates ==\n")
per_st |> select(program, stratum, n_w, n_m, est, se, p) |> r3() |> arrange(program, stratum) |> print(n = Inf, width = Inf)

# ---- sensitivity ----
sens <- bind_rows(
  run_meta(fit_all(scores, scorable, c("Stress", "CM_ambient"))) |> mutate(analysis = "+ stress + ambient RNA"),
  if (any(grepl("^dilated", names(by_st))))
    run_meta(fit_all(scores |> filter(grepl("^dilated", stratum)), scorable)) |> mutate(analysis = "DCM only"),
  bind_rows(lapply(names(by_st), function(left_out) run_meta(per_st |> filter(stratum != left_out, program %in% scorable)) |>
    mutate(analysis = paste("without", left_out)))))
cat("\n== Sensitivity analyses ==\n")
sens |> select(analysis, program, k, est, lo, hi, p) |> r3() |> print(n = Inf, width = Inf)

# ---- calibration ----
perm <- bind_rows(lapply(1:1000, function(i) {
  if (i %% 100 == 0) message("permutation ", i)
  sc <- scores |> group_by(stratum, assay) |> mutate(sex = sample(sex)) |> ungroup()
  run_meta(fit_all(sc, c(scorable, "X_escape"))) |> select(program, est, p)
}))
calib <- perm |> group_by(program) |> summarise(false_pos_rate_at_0.05 = mean(p < 0.05), pp = list(p)) |>
  left_join(main |> select(program, p_obs = p), by = "program") |>
  rowwise() |> mutate(p_perm = (1 + sum(unlist(pp) <= p_obs)) / (1 + length(unlist(pp)))) |> ungroup() |> select(-pp)
cat("\n== Calibration (1,000 permutations) ==\n")
calib |> r3() |> print(width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(main, "tables/P2_Table_primary.csv", row.names = FALSE)
write.csv(per_st, "tables/P2_TableS_primary_per_stratum.csv", row.names = FALSE)
write.csv(sens, "tables/P2_TableS_primary_sensitivity.csv", row.names = FALSE)
write.csv(calib, "tables/P2_TableS_primary_calibration.csv", row.names = FALSE)
saveRDS(list(scores = scores, per_st = per_st, main = main, sens = sens, calib = calib, perm = perm), "p2_primary.rds")
