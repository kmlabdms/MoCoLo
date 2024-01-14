# Supplement Figure S1-Table S2
# The distribution for the expected and the observed

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read MoCoLo results ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder <- "v3"
folder_path <- sprintf("output/%s", folder)

# 8-oxo-G as pivot
map_shuffle_tbl <- read_tsv(sprintf("output/%s/mapping_oxi_All_simuOXO.tsv", folder))
map_raw <- map_shuffle_tbl %>%
  mutate(group = ifelse(group == "oxidative", "observed", "expected")) %>%
  mutate(motif_type = factor(motif_type, levels = str_sort(unique(motif_type)))) %>%
  arrange(motif_type) %>%
  select(-ho_rate)
sig_8oxoG <- PlotVlineBatch(
  all_nonB = map_raw, metric = "total_overlapped_region",
  ncol = 2, gtitle = "H01: Number of 8-oxo-dG motifs in non-B "
)

# Non-B as pivot
map_shuffle_tbl_simuNB <- read_tsv(sprintf("output/%s/mapping_oxi_All_simuNB.tsv", folder))
map_raw_simuNB <- map_shuffle_tbl_simuNB %>%
  mutate(group = ifelse(group == "oxidative", "observed", "expected")) %>%
  mutate(motif_type = factor(motif_type, levels = str_sort(unique(motif_type)))) %>%
  arrange(motif_type) %>%
  select(-ho_rate)
sig_nonB_simuNB <- PlotVlineBatch(
  all_nonB = map_raw_simuNB, metric = "total_overlapped_motif",
  ncol = 2, gtitle = "H02: Number of non-B motifs in 8-oxo-dG"
)

# Combine the figures
sig_8oxoG | sig_nonB_simuNB




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 2: Summary table  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GetMappingSummary <- function(map_raw, pivot = NULL) {
  mapping <- map_raw %>%
    left_join(total_nonb) %>%
    dplyr::rename(
      colo_nonB = total_overlapped_motif,
      colo_8oxoG = total_overlapped_region
    ) %>%
    mutate(motif_type = str_replace_all(motif_type, "_", " ") %>% str_wrap(10))
  # 4 Get P-value----
  pvalue_dt <- CalcPvalueTable_noEmbed(map_raw)
  pvalue_long <- pvalue_dt %>%
    purrr::set_names(c("motif_type", "8oxoG", "nonB")) %>%
    gather(group, pvalue, - motif_type)

  # output summary table
  mapping_summary <- map_raw %>%
    group_by(motif_type, group) %>%
    summarise(
      mean_colo_nonB = mean(total_overlapped_motif),
      mean_colo_8oxoG = mean(total_overlapped_region)
    ) %>%
    gather("metrics", "value", -motif_type, -group) %>%
    spread(group, value) %>%
    arrange(metrics) %>%
    mutate(group = str_remove(metrics, "mean_colo_")) %>%
    # attach more info
    left_join(total_nonb) %>%
    left_join(pvalue_long) %>%
    ungroup()

  if(!is.null(pivot)) {
    mapping_summary = filter(mapping_summary, group %in% pivot)
  }

  return(mapping_summary)
}
summary_observed <- bind_rows(
  GetMappingSummary(map_raw , pivot = "8oxoG"),
  GetMappingSummary(map_raw_simuNB, pivot = "nonB")
)






