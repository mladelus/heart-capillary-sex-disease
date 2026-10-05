# =============================================================================
# 22  PLAN SECTION U (exploratory) - are genes controlled by the X-Y chromatin regulators shifted between women and men?
# Run:  source("~/Desktop/heart-capillary-sex-disease/R/00_setup.R"); source("~/Desktop/heart-capillary-sex-disease/R/22_xy_target_marks.R")
# Needs stage3_de.rds, stage3_gtex.rds, kidney/tables/kidney_sex_z_by_celltype.csv and internet access (gene sets).
# Makes: tables/S3U_target_marks.csv, tables/S3U_output.txt
# =============================================================================
dir.create("genesets", showWarnings = FALSE)
while (sink.number() > 0) sink()
sink("tables/S3U_output.txt", split = TRUE)
gg <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, columns = c("gene_name", "seq_name", "gene_seq_start", "gene_seq_end"), return.type = "data.frame")
gg <- gg[gg$seq_name %in% as.character(1:22) & !duplicated(gg$gene_name), ]; glen <- setNames(log10(gg$gene_seq_end - gg$gene_seq_start + 1), gg$gene_name)

# ---- sex z-scores per cell type -------------------------------------------------------------------------------------------
Z <- list()
D <- readRDS("stage3_de.rds")$de
for (ct in names(D)) { a <- D[[ct]] |> filter(chr_class == "autosome"); Z[[paste("heart", ct)]] <- setNames(a$z_comb, a$gene) }
gt <- readRDS("stage3_gtex.rds")$sex_all_genes; Z[["GTEx bulk left ventricle"]] <- setNames(gt$t_gtex, gt$gene)
kf <- "kidney/tables/kidney_sex_z_by_celltype.csv"
if (file.exists(kf)) { k <- read.csv(kf, row.names = 1, check.names = FALSE); for (ct in colnames(k)) { v <- setNames(k[[ct]], rownames(k)); Z[[paste("kidney", ct)]] <- v[!is.na(v)] } } else cat("Kidney z-scores not found:", kf, "\n")
Z <- lapply(Z, function(v) { v <- v[names(v) %in% names(glen) & is.finite(v)]; r <- residuals(lm(v ~ glen[names(v)])); setNames(as.numeric(r), names(v)) })   # autosomal genes; gene length removed
cat("Cell types:", length(Z), "\n")

# ---- gene sets --------------------------------------------------------------------------------------------------------------
SETS <- list(); kind <- c()
cgp <- tryCatch({ m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP"), error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")); split(m$gene_symbol, m$gs_name) }, error = function(e) NULL)
bp <- c("BENPORATH_PRC2_TARGETS", "BENPORATH_ES_WITH_H3K27ME3", "BENPORATH_SUZ12_TARGETS", "BENPORATH_EED_TARGETS")
if (!is.null(cgp) && all(bp %in% names(cgp))) { SETS[["PRIMARY polycomb/H3K27me3 targets (union of four Ben-Porath sets)"]] <- unique(unlist(cgp[bp])); for (b in bp) SETS[[b]] <- cgp[[b]] } else cat("Ben-Porath sets not available from msigdbr.\n")
kind[names(SETS)] <- "1 polycomb"
enrichr <- function(name) {
  f <- file.path("genesets", paste0(name, ".gmt"))
  if (!file.exists(f) || file.size(f) < 1e4) tryCatch(download.file(paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", name), f, mode = "wb", quiet = TRUE), error = function(e) NULL)
  if (!file.exists(f) || file.size(f) < 1e4) { cat("Library not available:", name, "\n"); return(list()) }
  L <- strsplit(readLines(f, warn = FALSE), "\t"); setNames(lapply(L, function(x) unique(sub(",.*$", "", x[-(1:2)][nzchar(x[-(1:2)])]))), vapply(L, `[`, "", 1))
}
rm_ <- enrichr("Epigenomics_Roadmap_HM_ChIP-seq"); keep <- grep("(H3K27me3|H3K4me3).*(left ventricle|kidney)|(left ventricle|kidney).*(H3K27me3|H3K4me3)", names(rm_), ignore.case = TRUE, value = TRUE)
for (s in keep) { SETS[[paste("Roadmap:", s)]] <- rm_[[s]]; kind[paste("Roadmap:", s)] <- "2 tissue histone marks" }
tp <- enrichr("TF_Perturbations_Followed_by_Expression"); keep <- grep("^(KDM6A|UTX|UTY|KDM5C|KDM5D|ZFX)[ _]", names(tp), ignore.case = TRUE, value = TRUE); keep <- keep[grepl("human", keep, ignore.case = TRUE) | !grepl("mouse", keep, ignore.case = TRUE)]
for (s in keep) { SETS[[paste("Perturbation:", s)]] <- tp[[s]]; kind[paste("Perturbation:", s)] <- "3 regulator perturbation" }
cat("Gene sets:", length(SETS), "\n"); print(data.frame(set = substr(names(SETS), 1, 90), genes = lengths(SETS), row.names = NULL))

# ---- test ---------------------------------------------------------------------------------------------------------------------
res <- bind_rows(lapply(names(Z), function(ct) bind_rows(lapply(names(SETS), function(s) {
  v <- Z[[ct]]; idx <- which(names(v) %in% SETS[[s]]); if (length(idx) < 25) return(NULL)
  r <- limma::cameraPR(v, list(s = idx), inter.gene.cor = 0.01)
  tibble(cell_type = ct, set = s, kind = kind[[s]], genes = length(idx), mean_z_set = mean(v[idx]), mean_z_rest = mean(v[-idx]),
         direction = ifelse(r$Direction == "Up", "higher in women", "higher in men"), p = r$PValue)
}))))
write.csv(res, "tables/S3U_target_marks.csv", row.names = FALSE)
prim <- res |> filter(grepl("^PRIMARY", set)) |> mutate(organ = sub(" .*$", "", cell_type))
cat("\n== U. PRIMARY set: polycomb/H3K27me3 targets, by cell type (women minus men; gene length removed) ==\n"); prim |> select(cell_type, genes, mean_z_set, mean_z_rest, direction, p) |> r3() |> print(n = Inf, width = Inf)
cat("\nCell types with p < 0.05, by organ and direction:\n"); print(prim |> group_by(organ, direction) |> summarise(cell_types = n(), p_below_0.05 = sum(p < 0.05)), n = Inf)
cat("\n== U. Secondary sets: number of cell types by direction (p < 0.05) ==\n")
res |> filter(!grepl("^PRIMARY", set)) |> group_by(kind, set = substr(set, 1, 80)) |>
  summarise(cell_types = n(), women_higher_sig = sum(p < 0.05 & direction == "higher in women"), men_higher_sig = sum(p < 0.05 & direction == "higher in men"), median_p = median(p)) |> r3() |> print(n = Inf, width = Inf)
sink()
message("Finished. Output saved in tables/S3U_output.txt")
