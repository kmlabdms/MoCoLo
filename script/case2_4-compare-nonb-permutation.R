# 0 Libraries ----
library(tidyverse)
library(patchwork)
library(readr)
library(GenomicRanges)
library(vroom)
library(ggprism)
library(liftOver)
library(scales)
library(tictoc)
library(valr)
library(furrr)
library(parallel)
library(ComplexHeatmap)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename



## functions
Convert_Metric <- function(bt_coverage_output){
  colnames_bt_coverage <- c(
    "chr", "start", "end", "n_motif",
    "length_overlap", "length_region", "coverage_overlap")

  map_motif_region <-
    bt_coverage_output %>%
    set_names(colnames_bt_coverage) %>% as_tibble() %>%
    filter(n_motif>0)

  total_Ol_motif <- sum(map_motif_region$n_motif)
  total_Ol_region <- nrow(map_motif_region)
  ho_rate <- total_Ol_motif / total_Ol_region

  metrics_output <- tibble(total_overlapped_motif=total_Ol_motif,
                           total_overlapped_region=total_Ol_region,
                           ho_rate=ho_rate)
  metrics_output
}
Overlap_Metric <- function(host_feature, nest_feature ){
  bt_coverage_output <- bed_coverage(host_feature, nest_feature)
  metrics_output <- Convert_Metric(bt_coverage_output)
  metrics_output
}
## mapping function
Feature_Mapping_vector <- function(host_feature_list, nest_feature) {
  total_output <-
    map(host_feature_list, Overlap_Metric,
        nest_feature = nest_feature) %>%
    bind_rows() %>%
    mutate(group = c("oxidative", rep("control", length(host_feature_list)-1)))
  total_output
}
Feature_Mapping_loop <- function(host_feature_list, nest_feature) {

  mapping_output <- list()

  for (i in 1:length(host_feature_list)) {
    host_feature <- host_feature_list[[i]]
    mapping_output[[i]] <- Overlap_Metric(host_feature, nest_feature)

  }

  total_output <-
    mapping_output %>%
    bind_rows() %>%
    mutate(group = c("oxidative", rep("control", length(host_feature_list)-1)))
  total_output
}

## pvalue
CalcPvalueTable <-  function(map_shuffle_tbl) {
  CalcPvalue <- function(map_output, metric) {
    values <- map_output %>% dplyr::pull(metric)
    mean( values[2:length(values)] > values[1])
  }

  map_shuffle <- split(map_shuffle_tbl, map_shuffle_tbl$motif_type)

  tibble(
    nonb_type = names(map_shuffle),
    pvalue_oxidative = sapply(map_shuffle, CalcPvalue, metric="total_overlapped_region"),
    pvalue_motif =  sapply(map_shuffle, CalcPvalue, metric="total_overlapped_motif"),
    pvalue_embed = sapply(map_shuffle, CalcPvalue, metric="ho_rate")
  )
}
CalcPctDiff <- function(twoShuffle, oxidative_one) {

  mapping <- MapGenomewide(
    oxidative_one,
    twoShuffle,
    n_core = 36L
  )


  pct_nonb <- mapping$total_overlapped_motif / unlist(lapply(twoShuffle, nrow))
  pct_oxog <- mapping$total_overlapped_region / nrow(oxidative_one[[1]])
  diff <- tibble(
    pct_nonb_diff = pct_nonb[1] - pct_nonb[2],
    pct_oxog_diff = pct_oxog[1] - pct_oxog[2]
  )

  mapping_pct <- mapping %>%
    mutate(total_nonb = unlist(lapply(twoShuffle, nrow))) %>%
    mutate(total_8oxoG = nrow(oxidative_one[[1]])) %>%
    mutate(pct_nonb = sprintf("%1.3f%%", pct_nonb*100),
           pct_8oxoG = sprintf("%1.3f%%", pct_oxog*100),
    ) %>%
    mutate(`8oxoG (overlap/total)` = sprintf("%i / %i", total_overlapped_region, total_8oxoG)) %>%
    mutate(`Non-B (overlap/total)` = sprintf("%i / %i", total_overlapped_motif, total_nonb)) %>%
    select(motif_type, `8oxoG (overlap/total)`, pct_8oxoG, `Non-B (overlap/total)`, pct_nonb)

  return(list(mapping=mapping_pct, diff=diff))
}
PlotVlineNonB <- function(mr, metric="pct_oxog_diff", title = "STR and ZDNA"){

  z <-  mr %>% split(.$group) %>% lapply(pull, metric)
  z0=z[["expected"]]
  z1=z[["observed"]]

  p_value <- mean(z0 > z1)

  if (p_value > 0.05) {
    vline_col <- "blue"
    vline_size = 0.5
    vline_type = "dashed"
  } else  {
    vline_col <- "red"
    vline_size <- 1
    vline_type = "solid"
  }

  p <- ggplot()+
    geom_vline(xintercept = z0, color="grey80")+
    geom_density(aes(x = z0, y=..scaled..), color="black")+
    geom_vline(xintercept = z1,color=vline_col, size = vline_size, linetype = vline_type)+
    theme_linedraw()+
    theme(axis.text  = element_text(size = 10, face="bold"),
          axis.title  = element_text(size = 10, face="bold"),
          plot.title = element_text(size=12, face = "bold")) +
    theme(panel.grid.major = element_line(linetype = "blank"),
          panel.grid.minor = element_line(linetype = "blank"),
          #axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()
    )+
    labs( y = "scaled") +
    #xlim(0,1)+
    ggtitle(label = waiver(), subtitle = title)

  if (p_value < 0.001) {
    pwithsig <- p+
      geom_text(aes(x = (quantile(z0, 0.5) + z1)/2, y = 0.1,
                    label = "p < 0.001, ***"),
                size = 5
      )
  } else if (p_value >= 0.05) {
    pwithsig <- p+
      geom_text(aes(x = (quantile(z0, 0.5) + z1)/2, y = 0.1,
                    label = sprintf("NS")),
                size = 5
      )
  } else {
    pwithsig <- p+
      geom_text(aes(x = (quantile(z0, 0.5) + z1)/2, y = 0.1,
                    label = "p < 0.05, *"),
                size = 5
      )
  }

  return(pwithsig)

}


# 1 Load Data ----
oxidative_one <- oxidative_all[1]
nonbDB


# Select Non-B Pairs for comparison ----
combx <- combn(1:7, 2) %>% as.data.frame()
twoMotifs_list <- lapply(combx, function(comb){nonbDB[comb]})

i = 20

perumutation_times = 100
for (i in 1:length(twoMotifs_list)) {
  twoMotifs <- twoMotifs_list[[i]]
  title_table <- names(twoMotifs) %>% paste(collapse = " - ")


  # RUNNING MAPPING ----
  twoShuffle_list <- lapply(1:perumutation_times, function(x) {
    bind_rows(twoMotifs) %>%
      # permutation (resample non-B type between the two)
      mutate(Type = sample(Type, replace = F)) %>%
      split(.$Type)
  })
  mapping_output <- CalcPctDiff(twoMotifs, oxidative_one)
  observed <- mapping_output$diff %>% mutate(group = "observed")
  expected <- lapply(twoShuffle_list, function(tbl) {
    CalcPctDiff(tbl, oxidative_one)[["diff"]]
  }) %>%
    bind_rows() %>%
    mutate(group = "expected")

  pct_diff <- bind_rows(observed, expected)
  plot_oxog_diff <- PlotVlineNonB(pct_diff, "pct_oxog_diff", title = sprintf("%s (8oxoG)", title_plot))
  plot_nonb_diff <- PlotVlineNonB(pct_diff, "pct_nonb_diff", title = sprintf("%s (Non-B)", title_plot))
  all_p <-
    patchwork::wrap_plots(plot_oxog_diff, plot_nonb_diff) /
    gridExtra::tableGrob(mapping_output$mapping)

  write_tsv(pct_diff, sprintf("output/pct_diff/%s.tsv", title_table))
}


