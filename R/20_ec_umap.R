# =============================================================================
# 20  STAGE 3L (descriptive) - UMAP of the blood endothelial cells of all 138 donors
# A picture, not a test: no p values come from this script. It shows where the capillary, arterial and venous
# cells sit, that female and male cells intermix, and where the EDN1-expressing cells lie.
# Up to 300 endothelial cells per donor (random, fixed seed) so that no donor dominates the map.
# X, Y and mitochondrial genes are left out of the variable genes, so sex is not built into the map.
# Datasets are aligned with Harmony (batch = dataset); if Harmony is not available the plain PCA is used and this is printed.
# Run 00_setup.R first. Makes: stage3_umap.rds, tables/S3L_ec_umap_cells.csv.gz
# =============================================================================
set.seed(2032)
MAX_PER_DONOR <- 300
C  <- load_cohorts()
ch <- chrom_map()
dirs <- c(healthy = PAPER1_DIR, disease = getwd())
S_all <- bind_rows(core_cols(C$healthy$samples), core_cols(C$disease$samples))

mats <- list(); cells <- list()
for (co in names(dirs)) {
  ec <- readRDS(file.path(dirs[[co]], "ec_subtypes.rds")); S <- C[[co]]$samples
  for (f in list.files(file.path(dirs[[co]], "dl"), full.names = TRUE)) {
    x <- readRDS(f); s <- paste(substr(x$dataset_id, 1, 8), x$donor_id, sep = "|"); if (!s %in% S$sample) next
    a <- ec |> filter(dataset_id == x$dataset_id, donor_id == x$donor_id, soma_joinid %in% colnames(x$ec_counts))
    if (nrow(a) < 10) next
    if (nrow(a) > MAX_PER_DONOR) a <- a[sort(sample.int(nrow(a), MAX_PER_DONOR)), ]
    m <- x$ec_counts[, a$soma_joinid, drop = FALSE]
    id <- paste0(co, "_", s, "_", a$soma_joinid); colnames(m) <- id
    mats[[paste(co, s)]] <- m
    cells[[paste(co, s)]] <- tibble(cell = id, cohort = co, sample = s, dataset_id = x$dataset_id, subtype = a$subtype,
                                    z_capillary = a$z_capillary, z_arterial = a$z_arterial, z_venous = a$z_venous,
                                    umis = as.numeric(Matrix::colSums(m)))
  }
}
cells <- bind_rows(cells) |> left_join(S_all |> select(sample, stratum, sex, age), by = "sample")
message(nrow(cells), " endothelial cells from ", n_distinct(cells$sample), " donors")

# one matrix with the union of genes (a gene missing from a donor's file had zero counts there)
genes <- sort(unique(unlist(lapply(mats, rownames))))
align <- function(m) { t <- as(m, "TsparseMatrix")
  Matrix::sparseMatrix(i = match(rownames(m), genes)[t@i + 1], j = t@j + 1, x = t@x,
                       dims = c(length(genes), ncol(m)), dimnames = list(genes, colnames(m))) }
M <- do.call(cbind, unname(lapply(mats, align))); rm(mats); gc(verbose = FALSE)
M <- M[, cells$cell]

seu <- CreateSeuratObject(counts = M, meta.data = as.data.frame(cells |> tibble::column_to_rownames("cell")))
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
drop <- c(ch$gene[ch$chr_class != "autosome"], grep("^MT-", rownames(seu), value = TRUE))
hvg  <- head(setdiff(VariableFeatures(seu), drop), 2000); VariableFeatures(seu) <- hvg
seu <- ScaleData(seu, features = hvg, verbose = FALSE)
seu <- RunPCA(seu, features = hvg, npcs = 30, verbose = FALSE)
red <- "pca"
if (!requireNamespace("harmony", quietly = TRUE)) try(install.packages("harmony"), silent = TRUE)
if (requireNamespace("harmony", quietly = TRUE)) {
  ok <- tryCatch({ seu <- harmony::RunHarmony(seu, group.by.vars = "dataset_id", verbose = FALSE); TRUE },
                 error = function(e) { message("Harmony failed: ", conditionMessage(e)); FALSE })
  if (ok) red <- "harmony"
}
message("UMAP computed on: ", red, if (red == "pca") " (datasets NOT aligned)" else " (datasets aligned)")
seu <- RunUMAP(seu, reduction = red, dims = 1:30, seed.use = 1, verbose = FALSE)
seu <- FindNeighbors(seu, reduction = red, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.4, random.seed = 1, verbose = FALSE)

show <- intersect(c("EDN1", "EDNRB", "ECE1", "KLF2", "SMAD6", "CXCL12", "NEBL", "AQP1", "CA4", "RGCC", "HEY1", "SEMA3G", "GJA5",
                    "ACKR1", "NR2F2", "PLVAP", "VWF", "PECAM1", "XIST", "UTY"), rownames(seu))
nd <- GetAssayData(seu, layer = "data")[show, , drop = FALSE]
out <- cells |> mutate(UMAP1 = Embeddings(seu, "umap")[cell, 1], UMAP2 = Embeddings(seu, "umap")[cell, 2],
                       cluster = as.character(seu$seurat_clusters[cell]),
                       EDN1_pos = if ("EDN1" %in% rownames(M)) as.numeric(M["EDN1", cell] > 0) else NA_real_) |>
  bind_cols(as_tibble(as.matrix(Matrix::t(nd))[cells$cell, , drop = FALSE]))

cat("\n== Cells on the map ==\n")
out |> count(cohort, sex, subtype) |> pivot_wider(names_from = subtype, values_from = n, values_fill = 0) |> print()
cat("\n== Clusters: size, % capillary / arterial / venous, % female, % EDN1-positive ==\n")
out |> group_by(cluster) |> summarise(cells = n(), donors = n_distinct(sample), pct_capillary = 100 * mean(subtype == "capillary"),
    pct_arterial = 100 * mean(subtype == "arterial"), pct_venous = 100 * mean(subtype == "venous"),
    pct_female = 100 * mean(sex == "female"), pct_EDN1_pos = 100 * mean(EDN1_pos), .groups = "drop") |>
  arrange(desc(cells)) |> r3() |> print(n = Inf, width = Inf)
cat(sprintf("\nFemale cells overall: %.1f%%. EDN1-positive cells overall: %.2f%%\n", 100 * mean(out$sex == "female"), 100 * mean(out$EDN1_pos)))

dir.create("tables", showWarnings = FALSE)
gz <- gzfile("tables/S3L_ec_umap_cells.csv.gz", "w"); write.csv(out, gz, row.names = FALSE); close(gz)
saveRDS(list(cells = out, reduction = red, hvg = hvg), "stage3_umap.rds")
message("Saved stage3_umap.rds and tables/S3L_ec_umap_cells.csv.gz")
