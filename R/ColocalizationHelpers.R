## functions
Convert_Metric <- function(bt_coverage_output, min_overlapped_counts = 1){
  colnames_bt_coverage <- c(
    "chr" = "chrom",
    "start",
    "end",
    "n_motif" = ".ints",
    "length_overlap" = ".cov",
    "length_region" = ".len",
    "coverage_overlap" = ".frac"
    )

  map_motif_region <-
    bt_coverage_output %>%
    select(colnames_bt_coverage) %>%
    filter(n_motif>=min_overlapped_counts)

  # region is the host_feature (being simulated)
  total_Ol_region <- nrow(map_motif_region)
  # motif is the nest_feature (not being simulated)
  # TODO: Same motif maybe calculated for twice
  # if we just simple sum it up
  total_Ol_motif <- sum(map_motif_region$n_motif)

  ho_rate <- total_Ol_motif / total_Ol_region

  metrics_output <- tibble(total_overlapped_motif=total_Ol_motif,
                           total_overlapped_region=total_Ol_region,
                           ho_rate=ho_rate)
  metrics_output
}
Overlap_Metric <- function(host_feature, nest_feature ){
  # perform overlap
  bt_coverage_output <- bed_coverage(host_feature, nest_feature)
  # covert to readable output
  metrics_output <- Convert_Metric(bt_coverage_output)
  metrics_output
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Feature_Mapping_loop ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Feature_Mapping_loop
#'
#'
#'
#' @param host_feature_list
#' @param nest_feature
#'
#' @return
#' @export
#'
#' @examples
Feature_Mapping_loop <- function(feature_simulation, feature_condition) {

  mapping_output <- list()

  for (i in 1:length(feature_simulation)) {
    mapping_output[[i]] <- Overlap_Metric(
      host_feature = feature_simulation[[i]],
      nest_feature = feature_condition)
  }

  total_output <-
    mapping_output %>%
    bind_rows() %>%
    # TODO: change to observed and expected
    mutate(group = c("oxidative", rep("control",
                                      length(feature_simulation)-1)))
  return(total_output)
}


#' MapGenomewide
#'
#' The main function to perform co-localization
#'
#' @param region_list
#' @param nonB_list
#' @param n_core
#' @param save
#'
#' @return
#' @export
#'
#' @examples
MapGenomewide <- function(feature_simulation, # a list
                          feature_condition_list,  # a table
                          n_core = 12L, save = FALSE){


  feature_simulation <- map(feature_simulation,
                            ~select(.x, chrom, start, end))

  mapping_out <-
    parallel::mclapply(
      X = feature_condition_list,
      FUN = Feature_Mapping_loop,
      feature_simulation = feature_simulation,
      mc.cores = n_core
    )


  ## TODO change to "motif_type" to "condition_feature_type"
  mapping_out_tbl <- mapping_out %>%  bind_rows(.id = "motif_type")

  return(mapping_out_tbl)
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
CalcPvalueTable_noEmbed <-  function(map_shuffle_tbl) {

  CalcPvalue <- function(map_output, metric) {
    values <- map_output %>% dplyr::pull(metric)
    mean(values[2:length(values)] > values[1])
  }

  map_shuffle <- split(map_shuffle_tbl, map_shuffle_tbl$motif_type)

  tibble(
    nonb_type = names(map_shuffle),
    pvalue_oxidative = sapply(map_shuffle, CalcPvalue, metric="total_overlapped_region"),
    pvalue_motif =  sapply(map_shuffle, CalcPvalue, metric="total_overlapped_motif"),
    #TODO
    #pvalue_embed = sapply(map_shuffle, CalcPvalue, metric="ho_rate")

  )
}
