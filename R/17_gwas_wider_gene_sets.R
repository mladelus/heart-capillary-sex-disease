# =============================================================================
# 17  STAGE 3I (exploratory, widened after seeing script 16) - GWAS Catalog, wider gene sets
# Are our sex-biased genes mapped to cardiovascular GWAS associations more often than comparable genes?
# Gene sets: cardiomyocyte sex-biased genes; X-linked genes higher in females in both cohorts; X members of X-Y pairs;
# capillary program genes; endothelin genes. Also lists every association at the EDN1 / PHACTR1 locus.
# Needs the GWAS Catalog file downloaded by script 16 (<data folder>/gwas) and stage3_de.rds, stage3_cm_pathways.rds.
# Makes: tables/S3I_gene_set_enrichment.csv, S3I_gene_set_hits.csv, S3I_EDN1_locus.csv, stage3_gwas_wide.rds
# =============================================================================
f <- list.files("gwas", pattern = "\\.tsv$", full.names = TRUE)[1]
if (is.na(f)) stop("GWAS Catalog file not found in ", normalizePath("gwas"), " - run script 16 first")
G <- data.table::fread(f, sep = "\t", quote = "", data.table = FALSE, fill = TRUE)
A <- tibble(trait = G$`DISEASE/TRAIT`, mapped_trait = G$MAPPED_TRAIT, snp = G$SNPS, chr = as.character(G$CHR_ID),
            pos = suppressWarnings(as.numeric(G$CHR_POS)), p = suppressWarnings(as.numeric(G$`P-VALUE`)),
            mapped_gene = G$MAPPED_GENE, author = G$`FIRST AUTHOR`, pubmed = G$PUBMEDID)
rm(G); gc(verbose = FALSE)
txt <- tolower(paste(A$trait, A$mapped_trait))
A$category <- case_when(
  grepl("testosterone|estradiol|oestradiol|shbg|sex hormone|menopause|menarche|puberty|androgen", txt) ~ "sex hormone or reproductive",
  grepl("coronary|myocard|heart|cardi|angina|blood pressure|hypertens|stroke|atrial|aort|vascular|arter|venous|thrombo|pulse|qt interval|electrocardio|pr interval|qrs", txt) ~ "cardiovascular",
  grepl("lipid|cholesterol|ldl|hdl|triglycer|body mass|bmi|obes|waist|adipos|fat |diabet|glucose|hba1c|insulin|metabolic", txt) ~ "cardiometabolic",
  TRUE ~ "other")
gw <- A |> filter(!is.na(p), p < 5e-8, !is.na(mapped_gene), mapped_gene != "")
long <- gw |> mutate(gene = strsplit(mapped_gene, " - |, |; |;|,")) |> tidyr::unnest(gene) |> mutate(gene = trimws(gene)) |>
  filter(gene != "") |> distinct(gene, snp, trait, .keep_all = TRUE)
cat("\nGenome-wide significant associations with a mapped gene:", nrow(gw), "; distinct mapped genes:", n_distinct(long$gene), "\n")
g_cv  <- unique(long$gene[long$category == "cardiovascular"])
g_cm  <- unique(long$gene[long$category == "cardiometabolic"])
g_any <- unique(long$gene)

# ---- gene sets and their backgrounds ----
S3 <- readRDS("stage3_de.rds")$de
cmg <- readRDS("stage3_cm_pathways.rds")$genes
bg_cm  <- S3$cardiomyocyte |> filter(chr_class == "autosome") |> pull(gene)
bg_cap <- S3$capillary |> filter(chr_class == "autosome") |> pull(gene)
xall <- bind_rows(lapply(names(S3), function(ct) S3[[ct]] |> filter(chr_class == "X") |> mutate(cell_type = ct)))
x_ours <- xall |> filter(FDR_h < 0.05, FDR_d < 0.05, z_h > 0, z_d > 0) |> pull(gene) |> unique()
XY18 <- c("KDM6A", "KDM5C", "USP9X", "DDX3X", "EIF1AX", "ZFX", "RPS4X", "NLGN4X", "TXLNG", "TBL1X", "PCDH11X", "TMSB4X",
          "PRKX", "AMELX", "SOX3", "RBMX", "HSFX1", "TSPYL2")
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)
gg <- as_tibble(ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, return.type = "data.frame",
        columns = c("gene_name", "seq_name", "gene_seq_start", "gene_seq_end", "gene_biotype")))
bg_x <- gg |> filter(seq_name == "X", gene_biotype == "protein_coding") |> pull(gene_name) |> unique()
sets <- list(
  "Cardiomyocyte sex-biased (55 consistent)" = list(genes = cmg$gene, bg = bg_cm),
  "Cardiomyocyte sex-biased (31 replicated)" = list(genes = cmg$gene[cmg$replicated], bg = bg_cm),
  "X-linked, higher in females in both cohorts" = list(genes = x_ours, bg = unique(xall$gene)),
  "X members of X-Y pairs (18)" = list(genes = XY18, bg = bg_x),
  "Capillary program genes" = list(genes = unique(unlist(programs)), bg = bg_cap),
  "Endothelin genes" = list(genes = c("EDN1", "EDNRA", "EDNRB", "ECE1"), bg = bg_cap))
one <- function(set_genes, bg, hits, what) {
  bg <- union(bg, set_genes); inset <- bg %in% set_genes; hit <- bg %in% hits
  ft <- fisher.test(table(factor(inset, c(FALSE, TRUE)), factor(hit, c(FALSE, TRUE))))
  tibble(trait_class = what, set_genes = sum(inset), set_with_hit = sum(inset & hit), pct_set = 100 * mean(hit[inset]),
         pct_background = 100 * mean(hit[!inset]), odds_ratio = unname(ft$estimate), p = ft$p.value)
}
enr <- bind_rows(lapply(names(sets), function(nm) bind_rows(one(sets[[nm]]$genes, sets[[nm]]$bg, g_cv, "cardiovascular"),
  one(sets[[nm]]$genes, sets[[nm]]$bg, g_cm, "cardiometabolic"), one(sets[[nm]]$genes, sets[[nm]]$bg, g_any, "any trait (control)")) |>
  mutate(gene_set = nm, .before = 1)))
cat("\n== Are the gene sets mapped to GWAS associations more often than comparable genes? (genome-wide tier) ==\n")
enr |> r3() |> print(n = Inf, width = Inf)

hits <- bind_rows(lapply(names(sets), function(nm) long |> filter(gene %in% sets[[nm]]$genes, category %in% c("cardiovascular", "cardiometabolic")) |>
  mutate(gene_set = nm))) |> select(gene_set, gene, category, trait, p, snp, author, pubmed)
cat("\n== Cardiovascular associations mapped to genes in each set (best per gene and trait) ==\n")
hits |> filter(category == "cardiovascular") |> group_by(gene_set, gene, trait) |> slice_min(p, n = 1, with_ties = FALSE) |> ungroup() |>
  arrange(gene_set, gene, p) |> mutate(trait = substr(trait, 1, 55)) |> select(-category) |> print(n = 200, width = Inf)
cat("\nGenes in each set with at least one cardiovascular association:\n")
hits |> filter(category == "cardiovascular") |> group_by(gene_set) |> summarise(genes = paste(sort(unique(gene)), collapse = ", ")) |> print(width = Inf)

# ---- the EDN1 / PHACTR1 locus ----
ed <- gg |> filter(gene_name %in% c("EDN1", "PHACTR1"), gene_biotype == "protein_coding")
lo <- min(ed$gene_seq_start) - 2e5; hi <- max(ed$gene_seq_end) + 2e5; chr_ed <- as.character(ed$seq_name[1])
loc <- A |> filter(chr == chr_ed, pos >= lo, pos <= hi, !is.na(p)) |>
  mutate(women_predominant = grepl("migraine|fibromuscular|dissection|microvascular|angina", tolower(paste(trait, mapped_trait))),
         tier = ifelse(p < 5e-8, "genome-wide", "suggestive"))
cat("\n== EDN1 / PHACTR1 locus (chr", chr_ed, ":", lo, "-", hi, "): associations by trait ==\n", sep = "")
loc |> group_by(trait = substr(trait, 1, 60)) |> summarise(associations = n(), best_p = min(p), category = category[1],
                                                           women_predominant = any(women_predominant), .groups = "drop") |>
  arrange(best_p) |> print(n = 60, width = Inf)
cat("\nThe regulatory variant rs9349379 (reported to control EDN1 expression):\n")
A |> filter(grepl("rs9349379", snp)) |> transmute(trait = substr(trait, 1, 60), p, author, pubmed) |> arrange(p) |> print(n = 40, width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(enr, "tables/S3I_gene_set_enrichment.csv", row.names = FALSE)
write.csv(hits, "tables/S3I_gene_set_hits.csv", row.names = FALSE)
write.csv(loc, "tables/S3I_EDN1_locus.csv", row.names = FALSE)
saveRDS(list(enrichment = enr, hits = hits, edn1_locus = loc, x_ours = x_ours), "stage3_gwas_wide.rds")
