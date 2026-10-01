# =============================================================================
# 13  STAGE 3F (exploratory) - cardiomyocyte sex-biased genes: pathways and protein network
# Needs stage3_de.rds (script 08) and internet access (MSigDB gene sets via msigdbr; STRING via its web API).
# Makes: tables/S3F_cm_gene_lists.csv, S3F_cm_ORA.csv, S3F_cm_camera.csv, S3F_cm_string_edges.csv,
#        S3F_endothelin_genes.csv, figures/S3_Fig_cm_network.pdf/.png, stage3_cm_pathways.rds
# =============================================================================
D  <- readRDS("stage3_de.rds")$de
ch <- chrom_map()
cm <- D$cardiomyocyte |> filter(chr_class == "autosome") |>
  mutate(p_comb = 2 * pnorm(-abs(z_comb)), FDR_comb = p.adjust(p_comb, "BH"),
         replicated = rep_h2d | rep_d2h,
         consistent = p_h < 0.05 & p_d < 0.05 & sign(z_h) == sign(z_d) & FDR_comb < 0.05,
         direction = ifelse(z_comb > 0, "higher in females", "higher in males"))
cat("\n== Cardiomyocyte gene lists (autosomal; background =", nrow(cm), "genes tested in both cohorts) ==\n")
cm |> summarise(replicated = sum(replicated), replicated_female = sum(replicated & z_comb > 0), replicated_male = sum(replicated & z_comb < 0),
                consistent = sum(consistent), consistent_female = sum(consistent & z_comb > 0), consistent_male = sum(consistent & z_comb < 0)) |> print()
genes_tab <- cm |> filter(replicated | consistent) |>
  select(gene, chr, direction, logFC_healthy = logFC_h, FDR_healthy = FDR_h, logFC_disease = logFC_d, FDR_disease = FDR_d,
         z_comb, FDR_comb, replicated, consistent) |> arrange(desc(abs(z_comb)))
cat("\n== All replicated or consistent cardiomyocyte genes ==\n")
genes_tab |> r3() |> print(n = Inf, width = Inf)

# ---- gene sets ----
msig <- function(cat, sub = NULL) {
  m <- tryCatch(if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", collection = cat)
                else msigdbr::msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub),
                error = function(e) if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", category = cat)
                                    else msigdbr::msigdbr(species = "Homo sapiens", category = cat, subcategory = sub))
  split(m$gene_symbol, m$gs_name)
}
collections <- list(Hallmark = tryCatch(msig("H"), error = function(e) NULL),
                    Reactome = tryCatch(msig("C2", "CP:REACTOME"), error = function(e) NULL),
                    GO_BP    = tryCatch(msig("C5", "GO:BP"), error = function(e) NULL))
collections <- collections[!sapply(collections, is.null)]
cat("\nGene-set collections loaded:", paste(names(collections), sapply(collections, length), collapse = "; "), "\n")
bg <- cm$gene
sets <- lapply(collections, function(L) { L <- lapply(L, intersect, bg); L[lengths(L) >= 10 & lengths(L) <= 500] })

# ---- over-representation, with the tested cardiomyocyte genes as background ----
lists <- list(replicated = cm$gene[cm$replicated], consistent = cm$gene[cm$consistent],
              consistent_higher_in_males = cm$gene[cm$consistent & cm$z_comb < 0],
              consistent_higher_in_females = cm$gene[cm$consistent & cm$z_comb > 0])
ora <- bind_rows(lapply(names(lists), function(ln) { g <- lists[[ln]]; if (length(g) < 5) return(NULL)
  bind_rows(lapply(names(sets), function(cn) {
    x <- bind_rows(lapply(names(sets[[cn]]), function(s) {
      S <- sets[[cn]][[s]]; ov <- intersect(g, S)
      tibble(list = ln, n_list = length(g), collection = cn, set = s, set_size = length(S), overlap = length(ov),
             expected = length(g) * length(S) / length(bg),
             p = phyper(length(ov) - 1, length(S), length(bg) - length(S), length(g), lower.tail = FALSE),
             genes = paste(sort(ov), collapse = ", "))
    }))
    if (!nrow(x)) return(NULL)
    x |> mutate(FDR = p.adjust(p, "BH")) |> filter(overlap >= 2)   # BH across every set in the collection
  }))
}))
cat("\n== Over-representation: sets with FDR < 0.10 (or the 5 best per list and collection if none) ==\n")
if (nrow(ora)) ora |> group_by(list, collection) |> arrange(p, .by_group = TRUE) |>
  filter(FDR < 0.10 | row_number() <= 5) |> ungroup() |>
  select(list, collection, set, set_size, overlap, expected, p, FDR, genes) |> r3() |> print(n = Inf, width = Inf)

# ---- competitive rank-based test allowing for inter-gene correlation (replaces the script 08 pathway test) ----
stat <- setNames(cm$z_comb, cm$gene)
cam <- bind_rows(lapply(names(sets), function(cn) {
  r <- limma::cameraPR(stat, lapply(sets[[cn]], function(g) which(names(stat) %in% g)), inter.gene.cor = 0.01)
  tibble(collection = cn, set = rownames(r), genes = r$NGenes,
         direction = ifelse(r$Direction == "Up", "higher in females", "higher in males"), p = r$PValue, FDR = r$FDR)
}))
cat("\n== cameraPR on all tested cardiomyocyte genes: sets with FDR < 0.05 (or 5 best per collection if none) ==\n")
cam |> group_by(collection) |> arrange(p, .by_group = TRUE) |> filter(FDR < 0.05 | row_number() <= 5) |> ungroup() |>
  r3() |> print(n = Inf, width = Inf)

# ---- STRING protein network ----
string_get <- function(method, genes, extra = "") {
  u <- paste0("https://string-db.org/api/tsv/", method, "?identifiers=", paste(genes, collapse = "%0d"),
              "&species=9606&caller_identity=heart_capillary_sex_disease", extra)
  tryCatch(read.delim(url(u), stringsAsFactors = FALSE), error = function(e) { message("STRING request failed: ", conditionMessage(e)); NULL })
}
net_genes <- unique(c(lists$replicated, lists$consistent))
edges <- NULL; ppi <- NULL; edges_plus <- NULL
if (length(net_genes) >= 3) {
  edges      <- string_get("network", net_genes, "&required_score=400")
  ppi        <- string_get("ppi_enrichment", net_genes, "&required_score=400")
  edges_plus <- string_get("network", net_genes, "&required_score=400&add_nodes=10")
}
if (!is.null(ppi)) { cat("\n== STRING: do these genes interact more than expected? ==\n"); print(ppi) }
if (!is.null(edges) && nrow(edges)) {
  edges <- edges |> distinct(preferredName_A, preferredName_B, .keep_all = TRUE)
  cat("\n== STRING interactions among the cardiomyocyte genes (confidence >= 0.4) ==\n")
  edges |> select(preferredName_A, preferredName_B, score) |> arrange(desc(score)) |> print()
} else cat("\nNo STRING interactions among the listed genes at confidence >= 0.4.\n")

if (!is.null(edges_plus) && nrow(edges_plus)) {
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  e <- edges_plus |> distinct(preferredName_A, preferredName_B, .keep_all = TRUE)
  g <- igraph::graph_from_data_frame(e[, c("preferredName_A", "preferredName_B", "score")], directed = FALSE)
  nm <- igraph::V(g)$name; z <- stat[nm]
  col <- ifelse(is.na(z) | !nm %in% net_genes, "grey85", ifelse(z > 0, "#D6604D", "#4393C3"))
  draw <- function() { set.seed(1)
    plot(g, layout = igraph::layout_with_fr(g), vertex.color = col, vertex.size = ifelse(nm %in% net_genes, 14, 9),
         vertex.label.cex = 0.6, vertex.label.color = "black", vertex.frame.color = "grey40",
         edge.width = 1 + 3 * igraph::E(g)$score, edge.color = "grey60")
    legend("bottomleft", c("higher in females", "higher in males", "added STRING partner"), pt.bg = c("#D6604D", "#4393C3", "grey85"),
           pch = 21, bty = "n", cex = 0.7) }
  pdf("figures/S3_Fig_cm_network.pdf", width = 6.7, height = 6.7); draw(); dev.off()
  png("figures/S3_Fig_cm_network.png", width = 170, height = 170, units = "mm", res = 300); draw(); dev.off()
  cat("\nNetwork figure saved (", igraph::vcount(g), " proteins, ", igraph::ecount(g), " interactions, including added partners).\n", sep = "")
  cat("Proteins in the network that are in our gene list:", paste(intersect(nm, net_genes), collapse = ", "), "\n")
}

# ---- endothelin / ACE program genes in capillary ECs, gene level ----
ET <- c("EDN1", "ECE1", "EDNRB", "ACE", "EDNRA", "ECE2")
et <- D$capillary |> filter(gene %in% ET) |>
  transmute(gene, in_program = gene %in% programs$Endothelin_ACE, logFC_healthy = logFC_h, p_healthy = p_h,
            logFC_disease = logFC_d, p_disease = p_d, z_comb, p_comb = 2 * pnorm(-abs(z_comb)))
cat("\n== Endothelin/ACE genes in capillary ECs, female - male (log2); genes not listed were too lowly expressed ==\n")
et |> r3() |> arrange(p_comb) |> print(width = Inf)

dir.create("tables", showWarnings = FALSE)
write.csv(genes_tab, "tables/S3F_cm_gene_lists.csv", row.names = FALSE)
write.csv(ora, "tables/S3F_cm_ORA.csv", row.names = FALSE)
write.csv(cam, "tables/S3F_cm_camera.csv", row.names = FALSE)
if (!is.null(edges)) write.csv(edges, "tables/S3F_cm_string_edges.csv", row.names = FALSE)
write.csv(et, "tables/S3F_endothelin_genes.csv", row.names = FALSE)
saveRDS(list(genes = genes_tab, ora = ora, camera = cam, string_edges = edges, string_ppi = ppi, endothelin = et), "stage3_cm_pathways.rds")
