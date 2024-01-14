#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Colocalization Testing ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# input for running  ----
folder <- "v3"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0 data preparation ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8oxoG
oxidative_all <- readRDS("data/first_case/oxoG/simuG_listed_by_permutation_50027x100.rds")
oxog_real_list = list(oxog = oxidative_all[[1]])
# NB
nonbDB

if(!exists("nonb_real_simu_list")) {
nonb_simu <- SimuPosistion(all_nonb_property[1:4], hg19.size, 100, add_ng = F)
nonb_real_simu_list =
  lapply(names(nonbDB),
         \(nonb_type) {
           c(list(nonbDB[[nonb_type]]),
             map(nonb_simu, ~ filter(.x, group == nonb_type)))
         })
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1 Colocal: Simulation 8-oxoG (Not by chr) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1.1 Running the mapping ----
map_shuffle_tbl <- MapGenomewide(
  feature_simulation = oxidative_all,
  feature_condition_list = nonbDB,
  n_core = 36
)
## 1.2 Get P-value ----
pvalue_dt <- CalcPvalueTable(map_shuffle_tbl)
## 1.3 Save output ----
if (F) {
  write_tsv(map_shuffle_tbl, sprintf("output/%s/mapping_oxi_All.tsv", folder))
  write_csv(pvalue_dt, sprintf("output/%s/mapping_oxi_pvalue_All.csv", folder))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3 Colocal: Simulation NonB (Not by Chr) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1.1 Running the mapping ----

#### single nonb
map_shuffle_tbl <- MapGenomewide(
  feature_simulation = nonb_real_simu_list[[7]],
  feature_condition_list = oxog_real_list,
  n_core = 36L
)

#### seven nonb
map_shuffle_tbl <-
  map(
    nonb_real_simu_list,
    ~MapGenomewide(
      feature_simulation = .x,
      feature_condition_list = oxog_real_list,
      n_core = 36L
    )
  )
map_shuffle_tbl_nonb =
  map_shuffle_tbl %>%
  set_names(names(nonbDB)) %>%
  map(~select(.x, -motif_type)) %>%
  bind_rows(.id = "motif_type")


## 1.2 Get P-value ----
pvalue_dt_nonb <- CalcPvalueTable(map_shuffle_tbl_nonb)
## 1.3 Save output ----
if (F) {
  write_tsv(map_shuffle_tbl_nonb, sprintf("output/%s/mapping_oxi_All_simuNB.tsv", folder))
  write_csv(pvalue_dt_nonb, sprintf("output/%s/mapping_oxi_pvalue_All_simuNB.csv", folder))
}



