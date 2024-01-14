# Figure 5 / Figure S3
library(tidyverse)
library(ggpubr)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1 Test One (within non-B) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Are there more non-Bs in 8-oxoG than non-8-oxoG regions?
# Comparision: (8-oxoG vs random)
# for each NonB

if (F) {
  metric = "metric1"
  result_path <- sprintf("output/v3/mapping_oxi_All_simuOXO.tsv")
} else {
  metric = "metric2"
  result_path <- sprintf("output/v3/mapping_oxi_All_simuNB.tsv")
}
map_shuffle_tbl <- read_tsv(result_path)
if (metric == "metric1") {
  pvalue_monte_carlo <- CalcPvalueTable(map_shuffle_tbl) %>%
    select(nonb_type, pvalue_oxidative) %>%
    mutate(monte_carlo = case_when(
      pvalue_oxidative <= 0.001 ~ "***",
      pvalue_oxidative <= 0.05 & pvalue_oxidative > 0.001 ~ "*",
      pvalue_oxidative > 0.05 ~ "ns"
    )) %>%
    mutate(y.position = 0.9)
}
if (metric == "metric2") {
  pvalue_monte_carlo <- CalcPvalueTable(map_shuffle_tbl) %>%
    select(motif_type = nonb_type, pvalue_motif) %>%
    mutate(monte_carlo = case_when(
      pvalue_motif <= 0.001 ~ "***",
      pvalue_motif <= 0.05 & pvalue_motif > 0.001 ~ "*",
      pvalue_motif > 0.05 ~ "ns"
    ))

}

map_raw <-
  map_shuffle_tbl %>%
  mutate(group = ifelse(group == "oxidative", "observed", "expected")) %>%
  mutate(group = factor(group, levels = c("observed", "expected"))) %>%
  mutate(motif_type = factor(motif_type, levels = str_sort(unique(motif_type)))) %>%
  arrange(motif_type) %>% select(-ho_rate) %>%
  mutate(
    motif_short = recode(
      motif_type,
      A_Phased_Repeat = "APR",
      Direct_Repeat = "DR",
      G_Quadruplex_Motif = "G4",
      Inverted_Repeat = "IR",
      Mirror_Repeat = "MR",
      Short_Tandem_Repeat = "STR",
      Z_DNA_Motif = "ZDNA"))

# when simulate oxoG
if (metric == "metric1") {
  # metric1
  oxog_test1 <- map_raw %>%
    select(motif_type, count = total_overlapped_region, group, motif_short) %>%
    mutate(pct = round( count / 50027, 5)) %>%
    mutate(motif_wrap = sprintf("8-oxo-dG\nin %s", motif_short))

  pvalue_monte_carlo_named =
    distinct(oxog_test1, motif_type, motif_wrap) %>%
    left_join(pvalue_monte_carlo, by = c("motif_type" = "nonb_type"))

} else if (metric == "metric2") {
  # metric2
  oxog_test1 <- map_raw %>%
    select(motif_type, count = total_overlapped_motif, group, motif_short) %>%
    left_join(total_nonb) %>%
    mutate(pct = round( count / total_nonB, 5)) %>%
    mutate(motif_wrap = sprintf("%s\nin 8-oxo-dG", motif_short))
  pvalue_monte_carlo_named =
    distinct(oxog_test1, motif_type, motif_wrap) %>%
    left_join(pvalue_monte_carlo)
}




# mean value
pct_mean <- oxog_test1 %>% group_by(motif_wrap, group) %>% dplyr::summarise(pct = mean(pct))
# p value
# this p-value is from compare_mean is not the real monte_carlo pvalue,
# it is just to formatt the table for visualization

stat_test <- ggpubr::compare_means(
  pct ~ group, data = oxog_test1, group.by = "motif_wrap"
) %>% mutate(y.position = 0.9)
if (metric == "metric1") {
  # metric1
  stat_test <- ggpubr::compare_means(
    pct ~ group, data = oxog_test1, group.by = "motif_wrap") %>%
    # adjust y postion to 90%
    mutate(y.position = 0.9) %>%
    left_join(pvalue_monte_carlo_named)
} else if  (metric == "metric2") {
  stat_test <- ggpubr::compare_means(
    pct ~ group, data = oxog_test1, group.by = "motif_wrap") %>%
    # adjust y postion to 9%
    mutate(y.position = 0.09) %>%
    left_join(pvalue_monte_carlo_named)

}




p <- ggplot(pct_mean, aes(group, pct)) +
  geom_bar(stat = "identity",
           aes(fill = group, color = group), alpha = 0.5) +
  #geom_jitter(data=oxog_test1, aes(color = group), size = 2, alpha = 0.5) +
  facet_wrap("motif_wrap", nrow = 1, scale = "free_x") +
  geom_text(aes(label = scales::percent(pct, accuracy = 0.01)), vjust = -0.5, fontface = "bold.italic", size = 3.5)
final_p <-
  p +
  theme_linedraw() +
  #rremove("x.grid")+
  theme(
    axis.text = element_text(size = 10,  face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 10, face = "bold", colour = "black"),

  )+
  theme(
    panel.spacing = unit(1, "lines"),
    panel.grid = element_blank()
  )

if (metric == "metric1") {
  # metric1
  final_p +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggpubr::stat_pvalue_manual(stat_test, label = "monte_carlo", label.size = 5) +
    scale_fill_manual(values = c("red", "grey50")) +
    scale_color_manual(values = c("red", "grey50")) +
    labs(#title = "Overlapped 8-oxo-dG regions",
         #subtitle = "Comparision: Observed (MCF10A, DIP-seq) vs Expected (Simulation)",
         y = "8-oxo-dG regions \n(Overlapped / Total)", x= "" )
} else if (metric == "metric2") {
  # metric2
  final_p +
    scale_y_continuous(labels = scales::percent, limits = c(0,0.1)) +
    ggpubr::stat_pvalue_manual(stat_test, label = "monte_carlo", label.size = 5) +
    scale_fill_manual(values = c("dodgerblue", "grey50")) +
    scale_color_manual(values = c("dodgerblue", "grey50")) +
    labs(#title = "Overlapped non-B motifs",
         #subtitle = "Comparision: Observed (MCF10A, DIP-seq) vs Expected (Simulation)",
         y = "Non-B motifs\n(Overlapped / Total)", x= ""  )
}






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 Visualize the testing results comparing across non-B ----
## computing script: script/case2_5-compare-nonb-permutation.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Are there more 8-oxoG region overalpped one type of non-Bs than another?
# Comparison: (paired non-B)
metric = "metric1"
result_path <- sprintf("output/v3/mapping_oxi_All_simuOXO.tsv")
map_shuffle_tbl <- read_tsv(result_path)

## CALCULATE P-VALUE ----

file_paths <- list.files("output/pct_diff", full.names = TRUE)
file_names <- str_remove(file_paths, "output/pct_diff/") %>%
  str_remove(".tsv") %>%
  str_replace(" - ", "__")
tbls <- lapply(file_paths, read_tsv) %>% set_names(file_names)
diff <- tbls %>% bind_rows(.id = "comparison") %>%
  separate("comparison", c("motif_A", "motif_B"), sep = "__", remove = FALSE) %>%
  mutate(across(where(is.numeric), abs))
diff_list <- diff %>% split(.$comparison)


# CALCULATE P-VALUE
mr = diff_list[[8]]
CalP <-  function(mr, metric="pct_oxog_diff") {
  z <-  mr %>% split(.$group) %>% lapply(pull, metric)
  z0=z[["expected"]]
  z1=z[["observed"]]
  p_value <- mean(z0 > z1)
  return(p_value)
}
diff_p_oxog <- sapply(diff_list, CalP, metric="pct_oxog_diff") %>%
  tibble(comparison = names(.), pvalue_oxog = .)
diff_p_nonB <- sapply(diff_list, CalP, metric="pct_nonb_diff") %>%
  tibble(comparison = names(.), pvlaue_nonb = .)

diff_all <- diff %>%  left_join(diff_p_oxog) %>%  left_join(diff_p_nonB)
diff_obs <- diff_all %>% filter(group == "observed") %>%
  mutate(pct_nonb_diff = (100 * pct_nonb_diff) %>% round(1),
         pct_oxog_diff = (100 * pct_oxog_diff) %>% round(0)
  )


# Using the "Overallap 8-oxo-G" as an example -----



if (metric == "metric1") {
  pvalue_res <- CalcPvalueTable(map_shuffle_tbl) %>%
    select(nonb_type, pvalue_oxidative) %>%
    mutate(monte_carlo = case_when(
      pvalue_oxidative <= 0.001 ~ "***",
      pvalue_oxidative <= 0.05 & pvalue_oxidative > 0.001 ~ "*",
      pvalue_oxidative > 0.05 ~ "ns"
    )) %>%
    mutate(y.position = 0.9)
}

map_raw <-
  map_shuffle_tbl %>%
  mutate(group = ifelse(group == "oxidative", "observed", "expected")) %>%
  mutate(group = factor(group, levels = c("observed", "expected"))) %>%
  mutate(motif_type = factor(motif_type, levels = str_sort(unique(motif_type)))) %>%
  arrange(motif_type) %>% select(-ho_rate) %>%
  mutate(
    motif_short = recode(
      motif_type,
      A_Phased_Repeat = "APR",
      Direct_Repeat = "DR",
      G_Quadruplex_Motif = "G4",
      Inverted_Repeat = "IR",
      Mirror_Repeat = "MR",
      Short_Tandem_Repeat = "STR",
      Z_DNA_Motif = "ZDNA"))

# when simulate oxoG
if (metric == "metric1") {
  # metric1
  oxog_test1 <- map_raw %>%
    select(motif_type, count = total_overlapped_region, group, motif_short) %>%
    mutate(pct = round( count / 50027, 5)) %>%
    mutate(motif_wrap = sprintf("8-oxo-dG\nin %s", motif_short))

  pvalue_monte_carlo_named =
    distinct(oxog_test1, motif_type, motif_wrap) %>%
    left_join(pvalue_monte_carlo)

}

## read mapping results from previous
oxog_seq <- oxog_test1 %>%  filter(group == "observed") %>%
  mutate(motif_wrap = fct_reorder(motif_wrap, pct,.desc = TRUE))

TidyPvalue <- function(diff_obs, col_tag) {
  st_diff <- diff_obs %>% select(group1 = motif_A, group2 = motif_B, pvalue_oxog, pvlaue_nonb) %>%
    mutate(
      group1 = recode(
        group1,
        A_Phased_Repeat = "APR",
        Direct_Repeat = "DR",
        G_Quadruplex_Motif = "G4",
        Inverted_Repeat = "IR",
        Mirror_Repeat = "MR",
        Short_Tandem_Repeat = "STR",
        Z_DNA_Motif = "ZDNA")) %>%
    mutate(group1 = sprintf(col_tag, group1)) %>%
    mutate(
      group2 = recode(
        group2,
        A_Phased_Repeat = "APR",
        Direct_Repeat = "DR",
        G_Quadruplex_Motif = "G4",
        Inverted_Repeat = "IR",
        Mirror_Repeat = "MR",
        Short_Tandem_Repeat = "STR",
        Z_DNA_Motif = "ZDNA")) %>%
    mutate(group2 = sprintf(col_tag, group2)) %>%
    mutate(pvalue_oxog_sig = case_when(
      pvalue_oxog < 0.001 ~ "***",
      pvalue_oxog < 0.05 ~ "*",
      TRUE ~ "ns"
    )) %>%
    mutate(pvlaue_nonb_sig = case_when(
      pvlaue_nonb < 0.001 ~ "***",
      pvlaue_nonb < 0.05 ~ "*",
      TRUE ~ "ns"
    ))

  return(st_diff)

}
if (metric == "metric1") {col_tag <- "8-oxo-dG\nin %s"}
st_diff_sig <- TidyPvalue(diff_obs, col_tag) %>%
  filter(pvalue_oxog_sig < 0.05)

## remove APR
st_diff_sig <- st_diff_sig %>%
  filter(!str_detect(group1, "APR|ZDNA") ) %>%
  filter(!str_detect(group2, "APR|ZDNA"))
oxog_seq <- oxog_seq %>% filter(motif_type != "A_Phased_Repeat")

if (metric == "metric1") {
  q <- ggplot(oxog_seq, aes(motif_wrap, pct)) +
    geom_bar(stat = "identity", aes(fill = motif_wrap), color = 'black') +
    geom_text(aes(label = scales::percent(pct, accuracy = 0.1)), vjust = -0.7) +
    ggpubr::stat_pvalue_manual(st_diff_sig, y.position = 0.7,
                       label = "pvalue_oxog_sig", step.increase  = 0.08 )+
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_fill_brewer(palette = "Set2")+
    theme_linedraw() +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9,  ),
      axis.text.y = element_text(size = 12,  face = "bold")
    )
  q + labs(#title = "Overalapped 8-oxoG in ",
           subtitle = "Comparision across non-B types (permuation)",
           y = "8-oxoG regions % \n(Overlapped / Total)", x= "" )
}







