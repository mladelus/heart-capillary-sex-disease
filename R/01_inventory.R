# =============================================================================
# 01  INVENTORY - diseased adult human ventricles in the Census (metadata only)
# Strata = dataset x disease with >= MIN_PER_SEX eligible females and males (>= MIN_EC ventricular ECs).
# Healthy donors from the same datasets are also kept (for the same-consortium sex x disease model).
# Makes: inv_cells.rds, inv_donors.rds, inv_datasets.rds
# =============================================================================
ds <- as.data.frame(census$get("census_info")$get("datasets")$read()$concat()) |>
  select(dataset_id, dataset_title, collection_name, collection_doi)

obs <- census$get("census_data")$get("homo_sapiens")$obs$read(
  value_filter = "tissue_general == 'heart' & is_primary_data == TRUE",
  column_names = c("soma_joinid", "dataset_id", "donor_id", "sex", "development_stage", "tissue", "disease",
                   "cell_type", "assay", "suspension_type", "raw_sum"))$concat() |>
  as.data.frame() |> as_tibble() |>
  mutate(across(where(is.factor), as.character),
         age = parse_age(development_stage), age_decade = development_stage %in% names(dec_mid),
         class = cell_class(cell_type)) |>
  left_join(ds, by = "dataset_id")

vent <- obs |> filter(age >= 18, str_detect(tissue, VENT), sex %in% c("female", "male"))
donors <- vent |> group_by(dataset_id, dataset_title, collection_name, donor_id, sex, age, age_decade, disease) |>
  summarise(cells = n(), ECs = sum(class == "blood EC"),
            assay = paste(sort(unique(assay)), collapse = "; "),
            suspension = paste(sort(unique(suspension_type)), collapse = "; "), .groups = "drop") |>
  filter(ECs >= MIN_EC)
dup <- donors |> count(dataset_id, donor_id) |> filter(n > 1)
if (nrow(dup)) message(nrow(dup), " donor(s) carry more than one disease label; keeping the label with most cells")
donors <- donors |> arrange(desc(cells)) |> distinct(dataset_id, donor_id, .keep_all = TRUE)

strata <- donors |> filter(disease != "normal") |> count(dataset_id, dataset_title, disease, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |>
  filter(female >= MIN_PER_SEX, male >= MIN_PER_SEX) |>
  mutate(stratum = paste0(disease, " (", substr(dataset_title, 1, 30), ")"))
if (anyDuplicated(strata$stratum)) strata <- strata |> mutate(stratum = paste0(disease, " (", substr(dataset_title, 1, 30), ", ", substr(dataset_id, 1, 6), ")"))
cat("\n== Disease strata (>= ", MIN_PER_SEX, " females and males) ==\n", sep = "")
print(strata |> select(stratum, female, male), width = Inf)

donors <- donors |>
  left_join(strata |> select(dataset_id, disease, stratum), by = c("dataset_id", "disease")) |>
  mutate(stratum = ifelse(disease == "normal" & dataset_id %in% strata$dataset_id,
                          paste0("normal (", substr(dataset_title, 1, 30), ")"), stratum),
         included = !is.na(stratum))
cat("\n== Included donors by stratum and sex ==\n")
donors |> filter(included) |> count(stratum, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print(n = Inf, width = Inf)
cat("\n== Age by stratum and sex (median, range) ==\n")
donors |> filter(included) |> group_by(stratum, sex) |>
  summarise(median_age = median(age), min = min(age), max = max(age), .groups = "drop") |> print(n = Inf, width = Inf)
cat("\n== Chemistry by stratum and sex ==\n")
donors |> filter(included) |> count(stratum, assay, sex) |>
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> print(n = Inf, width = Inf)

inc <- donors |> filter(included)
cat("\nTotal: ", nrow(inc), " donors; ", sum(vent$donor_id %in% inc$donor_id & vent$dataset_id %in% inc$dataset_id),
    " cells to download.\n", sep = "")
saveRDS(vent |> semi_join(inc, by = c("dataset_id", "donor_id")), "inv_cells.rds")
saveRDS(donors, "inv_donors.rds")
saveRDS(ds |> filter(dataset_id %in% inc$dataset_id), "inv_datasets.rds")
