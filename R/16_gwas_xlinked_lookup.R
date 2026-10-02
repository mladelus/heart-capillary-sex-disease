# =============================================================================
# 16  STAGE 3H (exploratory) - human genetics near the X-linked paralogs
# Step 1  GWAS Catalog associations within 500 kb of KDM6A, KDM5C, USP9X, DDX3X, EIF1AX, ZFX, RPS4X, NLGN4X (GRCh38)
# Step 2  GTEx v8 eQTLs for the same genes in heart, artery and blood; overlap with catalogued GWAS variants
# Needs internet. The GWAS Catalog file is large (several hundred MB) and is downloaded once into <data folder>/gwas.
# If the download fails: gwas catalog website -> Downloads -> "All associations ... with added ontology annotations",
# save the .tsv (or .zip) into that folder and run again.
# Makes: tables/S3H_gwas_near_xy_genes.csv, S3H_gwas_summary.csv, S3H_gtex_eqtl_xy_genes.csv, stage3_gwas.rds
# =============================================================================
for (p in c("data.table", "jsonlite")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
WINDOW <- 5e5
XG <- xy_pairs$X
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)
gg <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, return.type = "data.frame",
                       columns = c("gene_name", "seq_name", "gene_seq_start", "gene_seq_end", "gene_biotype"))
loc <- as_tibble(gg) |> filter(gene_name %in% XG, seq_name == "X", gene_biotype == "protein_coding") |>
  group_by(gene = gene_name) |> summarise(start = min(gene_seq_start), end = max(gene_seq_end), .groups = "drop")
cat("\n== X-linked genes (GRCh38) ==\n"); print(loc)

# ---- Step 1: GWAS Catalog ----
dir.create("gwas", showWarnings = FALSE); options(timeout = 3600)
have <- list.files("gwas", pattern = "\\.(tsv|zip)$", full.names = TRUE)
if (!length(have)) {
  urls <- c("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv",
            "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
            "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip")
  for (u in urls) {
    dest <- file.path("gwas", if (grepl("zip$", u)) "gwas_catalog.zip" else "gwas_catalog.tsv")
    ok <- tryCatch({ download.file(u, dest, mode = "wb"); file.size(dest) > 1e7 }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) break else unlink(dest)
  }
  have <- list.files("gwas", pattern = "\\.(tsv|zip)$", full.names = TRUE)
  if (!length(have)) stop("Could not download the GWAS Catalog. Download 'All associations with added ontology annotations' ",
                          "by hand from the GWAS Catalog Downloads page and put the file in ", normalizePath("gwas"))
}
f <- have[1]
if (grepl("zip$", f)) { unzip(f, exdir = "gwas"); f <- list.files("gwas", pattern = "\\.tsv$", full.names = TRUE)[1] }
G <- data.table::fread(f, sep = "\t", quote = "", data.table = FALSE, fill = TRUE)
need <- c("DISEASE/TRAIT", "CHR_ID", "CHR_POS", "SNPS", "P-VALUE", "PUBMEDID", "FIRST AUTHOR", "MAPPED_GENE", "MAPPED_TRAIT")
if (!all(need %in% names(G))) stop("Unexpected GWAS Catalog columns. Found: ", paste(head(names(G), 40), collapse = ", "))
X <- as_tibble(G[G$CHR_ID == "X", ]) |>
  transmute(trait = `DISEASE/TRAIT`, mapped_trait = MAPPED_TRAIT, snp = SNPS, pos = suppressWarnings(as.numeric(CHR_POS)),
            p = suppressWarnings(as.numeric(`P-VALUE`)), mapped_gene = MAPPED_GENE, author = `FIRST AUTHOR`, pubmed = PUBMEDID,
            sample = if ("INITIAL SAMPLE SIZE" %in% names(G)) `INITIAL SAMPLE SIZE` else NA_character_) |>
  filter(!is.na(pos))
cat("\nGWAS Catalog:", nrow(G), "associations in total;", nrow(X), "on the X chromosome with a position.\n")
txt <- tolower(paste(X$trait, X$mapped_trait))
X$category <- case_when(
  grepl("coronary|myocard|heart|cardi|angina|blood pressure|hypertens|stroke|atrial|aort|vascular|arter|venous|thrombo|pulse|qt interval|electrocardio", txt) ~ "cardiovascular",
  grepl("lipid|cholesterol|ldl|hdl|triglycer|body mass|bmi|obes|waist|adipos|fat |diabet|glucose|hba1c|insulin|metabolic", txt) ~ "cardiometabolic",
  grepl("testosterone|estradiol|oestradiol|shbg|sex hormone|menopause|menarche|puberty|androgen", txt) ~ "sex hormone or reproductive",
  grepl("blood cell|platelet|hemoglobin|haemoglobin|lymphocyte|neutrophil|monocyte|eosinophil|immun|autoimmun|lupus|arthritis", txt) ~ "blood or immune",
  TRUE ~ "other")
X$tier <- ifelse(!is.na(X$p) & X$p < 5e-8, "genome-wide (p < 5e-8)", "suggestive (5e-8 to 1e-5)")   # fixed before running
X$sex_specific <- grepl("\\b(women|men|female|male|females|males)\\b", tolower(paste(X$trait, X$sample)))
near_fun <- function(w) bind_rows(lapply(seq_len(nrow(loc)), function(i) X |> filter(pos >= loc$start[i] - w, pos <= loc$end[i] + w) |>
  mutate(gene = loc$gene[i], distance_kb = round(pmax(0, loc$start[i] - pos, pos - loc$end[i]) / 1000))))
near <- near_fun(WINDOW)
X$in_window <- X$pos %in% near$pos
cat("\n== Associations within 500 kb of each gene, by category ==\n")
summ <- near |> distinct(gene, snp, trait, .keep_all = TRUE) |> count(gene, category) |>
  pivot_wider(names_from = category, values_from = n, values_fill = 0) |> right_join(loc |> select(gene), by = "gene") |>
  mutate(across(-gene, ~ tidyr::replace_na(.x, 0)))
print(summ, width = Inf)
cat("\n== Cardiovascular, cardiometabolic and sex-hormone associations near the genes (best p first) ==\n")
near |> filter(category != "other", category != "blood or immune") |> arrange(gene, p) |>
  select(gene, distance_kb, snp, trait, p, tier, category, sex_specific, mapped_gene, author, pubmed) |>
  mutate(trait = substr(trait, 1, 60)) |> print(n = 150, width = Inf)
cat("\n== The same counts by evidence tier (cardiovascular + cardiometabolic only) ==\n")
near |> filter(category %in% c("cardiovascular", "cardiometabolic")) |> distinct(gene, snp, trait, .keep_all = TRUE) |>
  count(gene, tier) |> pivot_wider(names_from = tier, values_from = n, values_fill = 0) |> print(width = Inf)
cat("\n== Sensitivity: 1 Mb windows, cardiovascular + cardiometabolic associations per gene ==\n")
near_fun(1e6) |> filter(category %in% c("cardiovascular", "cardiometabolic")) |> distinct(gene, snp, trait, .keep_all = TRUE) |>
  count(gene, tier) |> pivot_wider(names_from = tier, values_from = n, values_fill = 0) |> print(width = Inf)
cm <- X$category %in% c("cardiovascular", "cardiometabolic")
ft <- fisher.test(table(in_window = X$in_window, cardiometabolic = cm))
cat(sprintf("\nShare of X-chromosome associations that are cardiovascular or cardiometabolic: %.1f%% inside the windows (%d of %d) vs %.1f%% elsewhere on X (%d of %d); odds ratio %.2f, p = %.3g.\n",
            100 * mean(cm[X$in_window]), sum(cm & X$in_window), sum(X$in_window), 100 * mean(cm[!X$in_window]), sum(cm & !X$in_window), sum(!X$in_window),
            ft$estimate, ft$p.value))

# ---- Step 2: GTEx v8 eQTLs for the same genes (GTEx Portal API) ----
api <- "https://gtexportal.org/api/v2/"
getj <- function(u) tryCatch(jsonlite::fromJSON(u), error = function(e) { message("GTEx request failed: ", conditionMessage(e)); NULL })
tissues <- c("Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Artery_Coronary", "Artery_Aorta", "Artery_Tibial", "Whole_Blood")
eq <- bind_rows(lapply(XG, function(gn) {
  r <- getj(paste0(api, "reference/gene?geneId=", gn, "&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38"))
  if (is.null(r) || !length(r$data) || !nrow(as.data.frame(r$data))) return(NULL)
  gid <- as.data.frame(r$data)$gencodeId[1]
  bind_rows(lapply(tissues, function(ti) {
    e <- getj(paste0(api, "association/singleTissueEqtl?gencodeId=", gid, "&tissueSiteDetailId=", ti, "&datasetId=gtex_v8&itemsPerPage=250"))
    if (is.null(e) || !length(e$data)) return(NULL)
    d <- as.data.frame(e$data); if (!nrow(d)) return(NULL)
    tibble(gene = gn, tissue = ti, snp = d$snpId, variant = d$variantId, p = d$pValue, effect = d$nes)
  }))
}))
if (!is.null(eq) && nrow(eq)) {
  cat("\n== GTEx v8: significant eQTL variants per gene and tissue (0 = none reported) ==\n")
  eq |> count(gene, tissue) |> pivot_wider(names_from = tissue, values_from = n, values_fill = 0) |> print(width = Inf)
  cat("\nStrongest eQTL per gene and tissue (effect = change in expression per alternative allele):\n")
  eq |> group_by(gene, tissue) |> slice_min(p, n = 1, with_ties = FALSE) |> ungroup() |> r3() |> print(n = Inf, width = Inf)
  ov <- eq |> inner_join(X |> select(snp, trait, gwas_p = p, category, author, pubmed), by = "snp")
  cat("\n== eQTL variants that are themselves catalogued GWAS variants ==\n")
  if (nrow(ov)) ov |> mutate(trait = substr(trait, 1, 60)) |> r3() |> print(n = Inf, width = Inf) else cat("None.\n")
} else { cat("\nNo GTEx eQTLs were returned for these genes (none significant, or the GTEx service could not be reached).\n"); ov <- tibble() }

# ---- Step 2b: candidate lookup. Each catalogued cardiovascular / cardiometabolic variant near a gene is tested as an
# eQTL for that gene (GTEx "dynamic eQTL" calculation). Because only these few variants are tested, the threshold is
# Bonferroni across the tests made here, not the genome-wide eQTL threshold. ----
cand <- near |> filter(category %in% c("cardiovascular", "cardiometabolic"), grepl("^rs[0-9]+$", snp)) |> distinct(gene, snp)
dyn <- tibble()
if (nrow(cand)) dyn <- bind_rows(lapply(seq_len(nrow(cand)), function(i) {
  v <- getj(paste0(api, "dataset/variant?snpId=", cand$snp[i], "&datasetId=gtex_v8"))
  if (is.null(v) || !length(v$data) || !nrow(as.data.frame(v$data))) return(NULL)
  vid <- as.data.frame(v$data)$variantId[1]
  r <- getj(paste0(api, "reference/gene?geneId=", cand$gene[i], "&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38"))
  if (is.null(r) || !length(r$data)) return(NULL)
  gid <- as.data.frame(r$data)$gencodeId[1]
  bind_rows(lapply(c("Heart_Left_Ventricle", "Artery_Coronary", "Artery_Aorta"), function(ti) {
    d <- getj(paste0(api, "association/dyneqtl?tissueSiteDetailId=", ti, "&gencodeId=", gid, "&variantId=", vid, "&datasetId=gtex_v8"))
    if (is.null(d) || is.null(d$pValue)) return(NULL)
    tibble(gene = cand$gene[i], snp = cand$snp[i], tissue = ti, p = as.numeric(d$pValue), effect = as.numeric(d$nes))
  }))
}))
if (nrow(dyn)) { dyn <- dyn |> mutate(p_bonferroni = pmin(1, p * n()))
  cat("\n== Candidate lookup: catalogued cardiovascular / cardiometabolic variants as eQTLs for the nearby gene (", nrow(dyn), " tests) ==\n", sep = "")
  dyn |> arrange(p) |> r3() |> print(n = Inf, width = Inf)
} else cat("\nCandidate eQTL lookup: no variants could be tested (none found in GTEx, or the service could not be reached).\n")

dir.create("tables", showWarnings = FALSE)
if (nrow(dyn)) write.csv(dyn, "tables/S3H_candidate_eqtl.csv", row.names = FALSE)
write.csv(near, "tables/S3H_gwas_near_xy_genes.csv", row.names = FALSE)
write.csv(summ, "tables/S3H_gwas_summary.csv", row.names = FALSE)
if (!is.null(eq) && nrow(eq)) write.csv(eq, "tables/S3H_gtex_eqtl_xy_genes.csv", row.names = FALSE)
saveRDS(list(genes = loc, near = near, summary = summ, fisher = ft, eqtl = eq, overlap = ov, candidate_eqtl = dyn), "stage3_gwas.rds")
