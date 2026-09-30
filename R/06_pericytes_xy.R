# =============================================================================
# 06  SECONDARY - pericytes and X-Y dosage in diseased hearts
# (a) pericytes per capillary EC (log ratio) and five pericyte programs: female - male within disease strata
# (b) X-Y paralogs in capillary ECs: X copy alone and X + Y (log2 CPM), female - male within disease strata,
#     compared with healthy capillaries from paper 1 (cross-study difference)
# Makes: p2_pericytes_xy.rds, tables/P2_TableS_pericytes.csv, tables/P2_TableS_xy.csv
# =============================================================================
P <- readRDS("pb.rds")
S <- P$samples |> filter(!grepl("^normal", stratum), n_EC >= MIN_EC) |>
  mutate(peri_cap_ratio = log((n_pericyte + 0.5) / (n_capillary + 0.5)))
r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))
fit_strata <- function(sc, vars) bind_rows(lapply(vars, function(v) bind_rows(lapply(split(sc, sc$stratum), function(d) {
  r <- sex_fit(d, v); if (is.null(r)) NULL else mutate(r, stratum = d$stratum[1], program = v) }))))
meta_all <- function(per) per |> group_by(program) |> group_modify(~ meta_one(.x)) |> ungroup()

# ---- (a) pericytes ----
ratio <- meta_all(fit_strata(S, "peri_cap_ratio"))
pe <- P$mats$pericyte
Sp <- S |> filter(n_pericyte >= MIN_CT, sample %in% colnames(pe))
peri_sc <- bind_rows(lapply(split(Sp, Sp$stratum), function(d) {
  if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX) return(NULL)
  y <- DGEList(pe[, d$sample, drop = FALSE]); y <- calcNormFactors(y[filterByExpr(y, group = d$sex), , keep.lib.sizes = FALSE])
  lc <- cpm(y, log = TRUE, prior.count = 1)
  bind_cols(d, as_tibble(sapply(c(peri_programs, pos_control), function(g) prog_score(lc, g))))
}))
peri <- meta_all(fit_strata(peri_sc, names(c(peri_programs, pos_control))))
cat("\n== (a) Pericytes in disease: female - male (SD units) ==\n")
bind_rows(ratio, peri) |> select(program, k, n_w, n_m, est, lo, hi, p, lo90, hi90, I2) |> r3() |> print(width = Inf)
cat("\nRaw pericytes per 100 capillary ECs (median) by stratum and sex:\n")
S |> group_by(stratum, sex) |> summarise(donors = n(), median = median(100 * n_pericyte / pmax(n_capillary, 1)),
                                         .groups = "drop") |> r3() |> print(n = Inf, width = Inf)

# ---- (b) X-Y dosage in capillary ECs ----
cap <- P$mats$capillary
Sc <- S |> filter(n_capillary >= MIN_CAP, sample %in% colnames(cap))
xy <- bind_rows(lapply(split(Sc, Sc$stratum), function(d) {
  cp <- edgeR::cpm(calcNormFactors(DGEList(cap[, d$sample, drop = FALSE])))
  bind_rows(lapply(seq_len(nrow(xy_pairs)), function(i) {
    X <- if (xy_pairs$X[i] %in% rownames(cp)) cp[xy_pairs$X[i], ] else rep(0, ncol(cp))
    Y <- if (xy_pairs$Y[i] %in% rownames(cp)) cp[xy_pairs$Y[i], ] else rep(0, ncol(cp))
    tibble(d, pair = paste(xy_pairs$X[i], xy_pairs$Y[i], sep = "/"), X_only = log2(X + 1), X_plus_Y = log2(X + Y + 1),
           Y_share = Y / (X + Y))
  }))
}))
lfit <- function(d, y) {
  d <- droplevels(d)
  if (sum(d$sex == "female") < MIN_PER_SEX || sum(d$sex == "male") < MIN_PER_SEX) return(tibble())
  f <- as.formula(paste(y, if (n_distinct(d$assay) > 1) "~ sex + age + assay" else "~ sex + age"))
  co <- summary(lm(f, data = d))$coefficients
  if (!"sexfemale" %in% rownames(co)) return(tibble())
  tibble(est = co["sexfemale", 1], se = co["sexfemale", 2], p = co["sexfemale", 4],
         n_w = sum(d$sex == "female"), n_m = sum(d$sex == "male"))
}
xy_res <- xy |> pivot_longer(c(X_only, X_plus_Y), names_to = "measure", values_to = "value") |>
  group_by(pair, measure, stratum) |> group_modify(~ lfit(.x, "value")) |> ungroup() |>
  group_by(pair, measure) |> group_modify(~ meta_one(.x)) |> ungroup()
healthy <- tryCatch(readRDS(file.path(PAPER1_DIR, "aim4.rds"))$by_group |> filter(group == "capillary") |>
  transmute(pair, measure, healthy = est, se_h = (hi - lo) / (2 * qnorm(0.975))), error = function(e) NULL)
if (!is.null(healthy)) xy_res <- xy_res |> left_join(healthy, by = c("pair", "measure")) |>
  mutate(se_d = (hi - lo) / (2 * qnorm(0.975)), change = est - healthy, change_p = 2 * pnorm(-abs(change / sqrt(se_d^2 + se_h^2))))
cat("\n== (b) X-Y paralogs in capillary ECs, female - male (log2), disease vs health (paper 1) ==\n")
xy_res |> select(any_of(c("pair", "measure", "k", "est", "lo", "hi", "p", "healthy", "change", "change_p"))) |>
  r3() |> arrange(measure, pair) |> print(n = Inf, width = Inf)
cat("\nY share in males' capillary ECs (median %), disease:\n")
xy |> filter(sex == "male") |> group_by(pair) |> summarise(Y_share_pct = round(100 * median(Y_share, na.rm = TRUE))) |> print()

dir.create("tables", showWarnings = FALSE)
write.csv(bind_rows(ratio, peri), "tables/P2_TableS_pericytes.csv", row.names = FALSE)
write.csv(xy_res, "tables/P2_TableS_xy.csv", row.names = FALSE)
saveRDS(list(ratio = ratio, peri = peri, xy_res = xy_res), "p2_pericytes_xy.rds")
