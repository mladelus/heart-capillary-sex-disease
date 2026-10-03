# =============================================================================
# 19  STAGE 3K (exploratory) - which capillary cells make EDN1, and the endothelin variant in GTEx
# (1) EDN1-positive vs depth-matched EDN1-negative capillary ECs within each donor (donor = unit of analysis)
#     i   paired differential expression of pseudobulks   ii  state scores, arterial-venous axis, EDNRB co-expression
#     iii EDN1 transcripts per positive cell, females vs males
# (2) GTEx: rs9349379 as an eQTL for EDN1 / PHACTR1; sex-biased eQTL list (if the open file downloads)
# Makes: tables/S3K_*.csv, stage3_edn1_cells.rds
# =============================================================================
set.seed(2031)
C <- load_cohorts()
SETS <- list(stress = STRESS, interferon = c("IFIT1", "IFIT3", "ISG15", "MX1", "IFI6", "IFI44L", "OAS1", "STAT1", "IFITM1", "XAF1"),
             inflammatory = c("SELE", "VCAM1", "ICAM1", "CXCL2", "CCL2", "IL6"), angiogenic = programs$Angiogenic_tip,
             hypoxia = c("VEGFA", "ADM", "SLC2A1", "PDK1", "BNIP3", "EGLN3", "ANKRD37"), shear_KLF2 = c("KLF2", "KLF4", "NOS3", "THBD"))
MIN_POS <- 5
dirs <- c(healthy = PAPER1_DIR, disease = getwd())
pos_pb <- list(); neg_pb <- list(); don <- list()
for (co in names(dirs)) {
  ec <- readRDS(file.path(dirs[[co]], "ec_subtypes.rds")); S <- C[[co]]$samples
  for (f in list.files(file.path(dirs[[co]], "dl"), full.names = TRUE)) {
    x <- readRDS(f); s <- paste(substr(x$dataset_id, 1, 8), x$donor_id, sep = "|"); if (!s %in% S$sample) next
    a <- ec |> filter(dataset_id == x$dataset_id, donor_id == x$donor_id, subtype == "capillary")
    m <- x$ec_counts; if (!"EDN1" %in% rownames(m) || nrow(a) < MIN_CAP) next
    m <- m[, a$soma_joinid, drop = FALSE]; umi <- Matrix::colSums(m); e <- m["EDN1", ]
    ip <- which(e > 0); if (length(ip) < MIN_POS) next
    dec <- cut(rank(umi, ties.method = "first"), 10, labels = FALSE)            # within-donor depth deciles
    ineg <- unlist(lapply(unique(dec[ip]), function(dd) { cand <- which(e == 0 & dec == dd); k <- min(length(cand), 3 * sum(dec[ip] == dd))
      if (k == 0) integer(0) else cand[sample.int(length(cand), k)] }))
    if (length(ineg) < MIN_POS) next
    cp10k <- function(idx, genes) { g <- intersect(genes, rownames(m)); if (!length(g)) return(NA_real_)
      mean(Matrix::colMeans(log1p(t(t(m[g, idx, drop = FALSE]) / umi[idx]) * 1e4))) }
    sc <- sapply(SETS, function(g) cp10k(ip, g) - cp10k(ineg, g))
    don[[s]] <- tibble(sample = s, n_cap = ncol(m), n_pos = length(ip), n_neg = length(ineg),
                       depth_pos = median(umi[ip]), depth_neg = median(umi[ineg]),
                       axis_diff = mean(a$z_venous[ip] - a$z_arterial[ip]) - mean(a$z_venous[ineg] - a$z_arterial[ineg]),
                       capscore_diff = mean(a$z_capillary[ip]) - mean(a$z_capillary[ineg]),
                       EDNRB_pct_pos = if ("EDNRB" %in% rownames(m)) 100 * mean(m["EDNRB", ip] > 0) else NA_real_,
                       EDNRB_pct_neg = if ("EDNRB" %in% rownames(m)) 100 * mean(m["EDNRB", ineg] > 0) else NA_real_,
                       EDN1_per_pos_cell = mean(log1p(e[ip] / umi[ip] * 1e4))) |> bind_cols(as_tibble(t(sc)) |> rename_with(~ paste0(.x, "_diff")))
    pos_pb[[s]] <- Matrix::rowSums(m[, ip, drop = FALSE]); neg_pb[[s]] <- Matrix::rowSums(m[, ineg, drop = FALSE])
  }
}
don <- bind_rows(don) |> left_join(bind_rows(core_cols(C$healthy$samples), core_cols(C$disease$samples)) |>
                                    select(sample, stratum, cohort, sex, age, assay), by = "sample")
cat("\n== Donors with >= ", MIN_POS, " EDN1-positive capillary cells ==\n", sep = "")
don |> group_by(cohort, sex) |> summarise(donors = n(), median_positive_cells = median(n_pos), median_matched_negatives = median(n_neg),
                                         depth_pos = median(depth_pos), depth_neg = median(depth_neg), .groups = "drop") |> r3() |> print()

# ---- (1 ii) how do EDN1-positive cells differ from matched neighbours? one value per donor, one-sample t test ----
dv <- c(paste0(names(SETS), "_diff"), "axis_diff", "capscore_diff")
st <- bind_rows(lapply(dv, function(v) { z <- don[[v]]; z <- z[is.finite(z)]; if (length(z) < 10) return(NULL); tt <- t.test(z)
  tibble(feature = sub("_diff$", "", v), donors = length(z), mean_difference = mean(z), lo = tt$conf.int[1], hi = tt$conf.int[2], p = tt$p.value) }))
eb <- don |> filter(is.finite(EDNRB_pct_pos)); tb <- t.test(eb$EDNRB_pct_pos - eb$EDNRB_pct_neg)
st <- bind_rows(st, tibble(feature = "% of cells with EDNRB (percentage points)", donors = nrow(eb), mean_difference = mean(eb$EDNRB_pct_pos - eb$EDNRB_pct_neg),
                           lo = tb$conf.int[1], hi = tb$conf.int[2], p = tb$p.value)) |> mutate(FDR = p.adjust(p, "BH"))
cat("\n== EDN1-positive minus matched EDN1-negative capillary cells (per-donor differences; axis: positive = more venous-like) ==\n")
st |> r3() |> arrange(p) |> print(width = Inf)
cat(sprintf("\nMedian %% of cells with EDNRB: EDN1-positive %.1f, matched negative %.1f\n", median(eb$EDNRB_pct_pos), median(eb$EDNRB_pct_neg)))

# ---- (1 i) paired differential expression, donor as blocking factor ----
genes <- sort(unique(unlist(lapply(c(pos_pb, neg_pb), names))))
stk <- function(L) sapply(L, function(v) { z <- setNames(numeric(length(genes)), genes); z[names(v)] <- v; z })
M <- cbind(stk(pos_pb), stk(neg_pb)); ids <- names(pos_pb)
colnames(M) <- c(paste0(ids, "_pos"), paste0(ids, "_neg"))
meta <- tibble(donor = factor(rep(ids, 2)), status = factor(rep(c("pos", "neg"), each = length(ids)), levels = c("neg", "pos")))
y <- DGEList(M[rownames(M) != "EDN1", ]); y <- calcNormFactors(y[filterByExpr(y, group = meta$status), , keep.lib.sizes = FALSE])
X <- model.matrix(~ donor + status, meta)
tt <- topTable(eBayes(lmFit(voom(y, X), X)), coef = "statuspos", number = Inf) |> tibble::rownames_to_column("gene") |> as_tibble()
cat("\n== Paired DE, EDN1-positive vs matched negative capillary cells: ", sum(tt$adj.P.Val < 0.05), " genes at FDR < 0.05 of ", nrow(tt),
    " tested (", length(ids), " donors) ==\n", sep = "")
cat("Top 30 (positive logFC = higher in EDN1-positive cells):\n")
tt |> select(gene, logFC, P.Value, adj.P.Val) |> head(30) |> r3() |> print(n = 30)
cat("\nGenes of interest:\n")
tt |> filter(gene %in% c("EDNRB", "ECE1", "KLF2", "KLF4", "NOS3", "VCAM1", "ICAM1", "VEGFA", "ANGPT2", "ESM1", "APLN", "HEY1", "GJA5", "ACKR1", "CA4")) |>
  select(gene, logFC, P.Value, adj.P.Val) |> r3() |> print(n = Inf)

# ---- (1 iii) EDN1 per positive cell, females vs males ----
per_cell <- {
  per <- bind_rows(lapply(split(don, don$stratum), function(d) { r <- sex_fit(d, "EDN1_per_pos_cell")
    if (is.null(r)) NULL else mutate(r, stratum = d$stratum[1], cohort = d$cohort[1]) }))
  if (nrow(per)) bind_rows(per |> group_by(cohort) |> group_modify(~ meta_one(.x)) |> ungroup(), meta_one(per) |> mutate(cohort = "both")) else tibble()
}
cat("\n== EDN1 transcripts per EDN1-positive capillary cell, female - male (SD units) ==\n")
if (nrow(per_cell)) per_cell |> select(cohort, k, n_w, n_m, est, lo, hi, p) |> r3() |> print() else cat("Too few donors per group to estimate.\n")

# ---- (2) GTEx look-ups ----
for (p in c("jsonlite", "data.table")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
api <- "https://gtexportal.org/api/v2/"
getj <- function(u) tryCatch(jsonlite::fromJSON(u), error = function(e) { message("GTEx request failed: ", conditionMessage(e)); NULL })
eq <- tibble()
v <- getj(paste0(api, "dataset/variant?snpId=rs9349379&datasetId=gtex_v8"))
if (!is.null(v) && length(v$data) && nrow(as.data.frame(v$data))) {
  vid <- as.data.frame(v$data)$variantId[1]
  eq <- bind_rows(lapply(c("EDN1", "PHACTR1"), function(gn) {
    r <- getj(paste0(api, "reference/gene?geneId=", gn, "&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38")); if (is.null(r) || !length(r$data)) return(NULL)
    gid <- as.data.frame(r$data)$gencodeId[1]
    bind_rows(lapply(c("Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Heart_Left_Ventricle", "Heart_Atrial_Appendage"), function(ti) {
      d <- getj(paste0(api, "association/dyneqtl?tissueSiteDetailId=", ti, "&gencodeId=", gid, "&variantId=", vid, "&datasetId=gtex_v8"))
      if (is.null(d) || is.null(d$pValue)) return(NULL)
      tibble(variant = "rs9349379", gene = gn, tissue = ti, p = as.numeric(d$pValue), effect = as.numeric(d$nes)) })) }))
}
cat("\n== GTEx: rs9349379 as an eQTL (effect = change per alternative allele; both sexes together) ==\n")
if (nrow(eq)) eq |> r3() |> arrange(gene, p) |> print(n = Inf) else cat("Could not be retrieved.\n")

sb <- NULL; dir.create("gtex", showWarnings = FALSE); sbf <- list.files("gtex", pattern = "sbeQTL.*\\.txt(\\.gz)?$", full.names = TRUE)
if (!length(sbf)) { options(timeout = 1800)
  for (u in c("https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/sex-biased-eqtl/GTEx_Analysis_v8_sbeQTLs.tar.gz",
              "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/sb-eqtl/GTEx_Analysis_v8_sbeQTLs.tar.gz")) {
    dest <- file.path("gtex", "sbeQTLs.tar.gz")
    ok <- tryCatch({ download.file(u, dest, mode = "wb", quiet = TRUE); file.size(dest) > 1e5 }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) { untar(dest, exdir = "gtex"); break } else unlink(dest) }
  sbf <- list.files("gtex", pattern = "sbeQTL.*\\.txt(\\.gz)?$", full.names = TRUE, recursive = TRUE) }
if (length(sbf)) {
  sb <- data.table::fread(sbf[1], data.table = FALSE); gcol <- grep("hugo|gene_name|symbol", names(sb), ignore.case = TRUE, value = TRUE)[1]
  cm <- tryCatch(readRDS("stage3_cm_pathways.rds")$genes$gene, error = function(e) character(0))
  want <- unique(c("EDN1", "PHACTR1", "EDNRA", "EDNRB", "ECE1", xy_pairs$X, cm))
  cat("\n== GTEx sex-biased eQTL list: rows for our genes in heart and artery tissues ==\n")
  hit <- sb[sb[[gcol]] %in% want, ]; tcol <- grep("tissue", names(hit), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(tcol)) hit <- hit[grepl("Heart|Artery", hit[[tcol]]), ]
  qcol <- grep("^qval|q_val", names(hit), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(qcol)) hit <- hit[order(hit[[qcol]]), ]
  print(as_tibble(hit) |> head(60), n = 60, width = Inf)
  if (!is.na(qcol)) cat("\nRows with q < 0.25 (the threshold used by GTEx for sex-biased eQTLs):", sum(hit[[qcol]] < 0.25, na.rm = TRUE), "\n")
  sb <- hit
} else cat("\nThe GTEx sex-biased eQTL file could not be downloaded. If you want this look-up: gtexportal.org -> Downloads -> QTL ->\n",
           "'Sex-biased eQTLs' (v8), save the .tar.gz or .txt into ", normalizePath("gtex"), " and run this script again.\n", sep = "")

dir.create("tables", showWarnings = FALSE)
write.csv(st, "tables/S3K_EDN1_positive_cell_features.csv", row.names = FALSE)
write.csv(tt, "tables/S3K_EDN1_positive_vs_negative_DE.csv", row.names = FALSE)
if (nrow(eq)) write.csv(eq, "tables/S3K_rs9349379_eqtl.csv", row.names = FALSE)
saveRDS(list(donors = don, features = st, de = tt, per_cell = per_cell, eqtl = eq, sex_biased = sb), "stage3_edn1_cells.rds")
