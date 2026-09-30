# =============================================================================
# 05  SECONDARY - does the female-male difference change from health to disease?
# (a) Same consortium: Heart Cell Atlas healthy donors (paper 1) + DCM atlas healthy and DCM donors.
#     score (standardized) ~ sex * disease + age + assay; estimand = sex:disease interaction.
#     Disease positive control: the DCM main effect (programs) and genome-wide DCM vs healthy genes.
#     Sensitivity: + study term; + stress + ambient RNA.
# (b) All data, cross-study: pooled disease estimate (script 04) minus pooled healthy estimate (paper 1).
# Makes: p2_interaction.rds, tables/P2_Table_interaction.csv, tables/P2_TableS_cross_study.csv
# =============================================================================
P2 <- readRDS("pb.rds")
P1 <- readRDS(file.path(PAPER1_DIR, "pb.rds"))
merge_mats <- function(a, b) {
  g <- union(rownames(a), rownames(b))
  fill <- function(m) { out <- matrix(0, length(g), ncol(m), dimnames = list(g, colnames(m))); out[rownames(m), ] <- m; out }
  cbind(fill(a), fill(b))
}

# paper 1 used two datasets; the Heart Cell Atlas is the one not titled "Heart" (as in paper 1, script 10)
s1 <- P1$samples |> filter(!grepl("^Heart$", dataset_title), n_capillary >= MIN_CAP) |>
  mutate(disease = "normal", study = "Heart Cell Atlas")
dcm_ds <- unique(as.character(P2$samples$dataset_id[grepl("^dilated", P2$samples$stratum)]))
if (length(dcm_ds) != 1) stop("Expected exactly one DCM dataset; found ", length(dcm_ds))
s2 <- P2$samples |> filter(dataset_id == dcm_ds, disease %in% c("dilated cardiomyopathy", "normal"),
                           n_capillary >= MIN_CAP) |> mutate(study = "DCM atlas")
if (!nrow(s1)) stop("No Heart Cell Atlas donors found in paper 1 pb.rds")
if (!any(s2$disease == "normal")) warning("No healthy donors in the DCM atlas: the + study sensitivity cannot be estimated")
cat("\n== Same-consortium model: donors ==\n")
bind_rows(s1, s2) |> count(study, disease, sex) |> pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print()

cap <- merge_mats(P1$mats$capillary[, s1$sample, drop = FALSE], P2$mats$capillary[, s2$sample, drop = FALSE])
S <- bind_rows(s1 |> mutate(across(c(sex, assay), as.character)), s2 |> mutate(across(c(sex, assay), as.character))) |>
  mutate(sex = relevel(factor(sex), ref = "male"),
         disease = relevel(factor(ifelse(disease == "normal", "healthy", "DCM")), ref = "healthy"),
         assay = factor(assay), study = factor(study)) |> filter(!is.na(age))
cap <- cap[, S$sample]

y  <- DGEList(cap); y <- calcNormFactors(y[filterByExpr(y, group = paste(S$sex, S$disease)), , keep.lib.sizes = FALSE])
lc <- cpm(y, log = TRUE, prior.count = 1)
sets <- c(programs, pos_control, list(Stress = STRESS, CM_ambient = CM_AMBIENT))
sc <- bind_cols(S, as_tibble(sapply(sets, function(g) prog_score(lc, g))))

int_fit <- function(v, extra = NULL) {
  d <- sc[is.finite(sc[[v]]), ]; d$yy <- as.numeric(scale(d[[v]]))
  f <- as.formula(paste("yy ~ sex * disease + age", if (n_distinct(d$assay) > 1) "+ assay" else "",
                        if (length(extra)) paste("+", paste(extra, collapse = " + ")) else ""))
  co <- summary(lm(f, data = d))$coefficients
  g <- function(t) if (t %in% rownames(co)) co[t, c(1, 2, 4)] else c(NA, NA, NA)
  tibble(program = v,
         sex_in_healthy = g("sexfemale")[1], sex_in_healthy_p = g("sexfemale")[3],
         DCM_effect_in_males = g("diseaseDCM")[1], DCM_effect_p = g("diseaseDCM")[3],
         interaction = g("sexfemale:diseaseDCM")[1], int_se = g("sexfemale:diseaseDCM")[2],
         int_p = g("sexfemale:diseaseDCM")[3]) |>
    mutate(int_lo = interaction - 1.96 * int_se, int_hi = interaction + 1.96 * int_se,
           sex_in_DCM = sex_in_healthy + interaction)
}
scorable <- names(sets)[sapply(names(sets), function(v) sum(is.finite(sc[[v]])) > 0)]
inter <- bind_rows(lapply(scorable, int_fit)) |>
  mutate(int_p_holm = ifelse(program %in% PRIMARY, p.adjust(ifelse(program %in% PRIMARY, int_p, NA), "holm"), NA))
r3 <- function(d) d |> mutate(across(where(is.double), ~ signif(.x, 3)))
cat("\n== (a) Sex x DCM, same consortium (SD units; interaction = change in female - male from health to DCM) ==\n")
inter |> select(program, sex_in_healthy, sex_in_DCM, interaction, int_lo, int_hi, int_p, int_p_holm,
                DCM_effect_in_males, DCM_effect_p) |> r3() |> print(n = Inf, width = Inf)

sens <- bind_rows(
  bind_rows(lapply(intersect(PRIMARY, scorable), int_fit, extra = "study")) |> mutate(analysis = "+ study"),
  bind_rows(lapply(intersect(PRIMARY, scorable), int_fit, extra = c("Stress", "CM_ambient"))) |>
    mutate(analysis = "+ stress + ambient RNA"))
cat("\n== (a) Sensitivity ==\n")
sens |> select(analysis, program, interaction, int_lo, int_hi, int_p) |> r3() |> print(n = Inf, width = Inf)

# genome-wide disease positive control (and genome-wide interaction count)
X  <- model.matrix(if (n_distinct(S$assay) > 1) ~ sex * disease + age + assay else ~ sex * disease + age, S)
X  <- X[, colSums(abs(X)) > 0 & !duplicated(t(X)), drop = FALSE]
fit <- eBayes(lmFit(voom(y, X), X))
cat("\n== Genome-wide (capillary ECs): genes at FDR < 0.05 ==\n")
for (cf in intersect(c("diseaseDCM", "sexfemale", "sexfemale:diseaseDCM"), colnames(X)))
  cat(sprintf("  %-22s %5d\n", cf, sum(topTable(fit, coef = cf, number = Inf)$adj.P.Val < 0.05)))
gw_int <- topTable(fit, coef = "sexfemale:diseaseDCM", number = Inf) |> tibble::rownames_to_column("gene")
cat("\nTop 15 genes for the sex x DCM interaction:\n")
gw_int |> head(15) |> select(gene, logFC, P.Value, adj.P.Val) |> r3() |> print()

# ---- (b) cross-study: pooled disease vs pooled healthy (paper 1) ----
A1 <- readRDS(file.path(PAPER1_DIR, "aim1.rds"))$main
A2 <- readRDS("p2_primary.rds")$main
se_of <- function(lo, hi) (hi - lo) / (2 * qnorm(0.975))
cross <- inner_join(A1 |> transmute(program, healthy = est, se_h = se_of(lo, hi)),
                    A2 |> transmute(program, disease = est, se_d = se_of(lo, hi)), by = "program") |>
  filter(program %in% c(PRIMARY, "X_escape")) |>
  mutate(difference = disease - healthy, se = sqrt(se_h^2 + se_d^2),
         lo = difference - 1.96 * se, hi = difference + 1.96 * se, p = 2 * pnorm(-abs(difference / se)))
cat("\n== (b) Cross-study: female - male in disease minus female - male in health ==\n")
cross |> r3() |> print(width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(inter, "tables/P2_Table_interaction.csv", row.names = FALSE)
write.csv(sens, "tables/P2_TableS_interaction_sensitivity.csv", row.names = FALSE)
write.csv(gw_int, "tables/P2_TableS_genomewide_interaction.csv", row.names = FALSE)
write.csv(cross, "tables/P2_TableS_cross_study.csv", row.names = FALSE)
saveRDS(list(scores = sc, inter = inter, sens = sens, gw_int = gw_int, cross = cross), "p2_interaction.rds")
