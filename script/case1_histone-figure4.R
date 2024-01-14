# Case 1
# This is the first case for histone data
# Figure 4B


#' TestCoLocalHistone
#'
#' Test Co-Localization for Histone data by genomic regions
#' histone: H3K9me3, h4k20me3
#'
#'
#' @param loci_annotation  "Intergenic", "TTS", "exon",  "intron", "promoter-TSS"
#' @param output_counts TRUE or FALSE  whether to output the overlapped counts
#' @param pivot_name  H3K9me3 or h4k20me3
#' @param n_simulation default=100
#'
#' @return
#' @export
#'
#' @examples
TestCoLocalHistone <- function (
    loci_annotation = "intron",
    output_counts = FALSE,
    pivot_name = "H3K9me3",
    n_simulation = 100
    ) {

  # Build hg19
  genome_build <- hg19.size

  # peak annotations
  h3k9Anno <- read_tsv("data/second_case/peaks/H3K9me3_R1_peaks.annotatePeaks.txt") %>% clean_names() %>%
    dplyr::rename(name = peak_id_cmd_annotate_peaks_pl_h3k9me3_r1_peaks_broad_peak_genome_fa_gid_gtf_genes_gtf_cpu_6) %>%
    select(name, annotation) %>%
    mutate(annotation = str_remove(annotation, " \\(.*"))
  h4k20Anno <-
    read_tsv("data/second_case/peaks/H4K20me3_R1_peaks.annotatePeaks.txt") %>% clean_names() %>%
    dplyr::rename(name = peak_id_cmd_annotate_peaks_pl_h4k20me3_r1_peaks_broad_peak_genome_fa_gid_gtf_genes_gtf_cpu_6) %>%
    select(name, annotation) %>%
    mutate(annotation = str_remove(annotation, " \\(.*"))


  h3k9me3 <-
    read_csv("data/second_case/peaks/H3K9me3_R1_peaks.csv") %>%
    clean_names() %>%
    left_join(h3k9Anno) %>% dplyr::rename(chrom = chr) %>%
    mutate(chrom = paste0("chr", chrom))  %>%
    filter(annotation == loci_annotation)
  h4k20me3 <- read_csv("data/second_case/peaks/H4K20me3_R1_peaks.csv") %>%
    clean_names() %>%
    left_join(h4k20Anno) %>% dplyr::rename(chrom = chr) %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    filter(annotation == loci_annotation)

  #  2 Simulation  ----------------------------------------------------------

  if (pivot_name == "H3K9me3"){
    pivot <- h3k9me3
    query <- h4k20me3
  } else if (pivot_name == "h4k20me3") {
    pivot <- h4k20me3
    query <- h3k9me3
  }

  # Prep for simulation by shuffle Coordinate
  n <- n_simulation
  pivot_simu <- list()
  pivot_simu[[1]] <- pivot
  for (i in 2:n) {
    pivot_simu[[i]] <- valr::bed_shuffle(pivot, genome_build)
  }


  # Check overlapps by regions ----------------------------------------------
  list_overlap <- list()
  # Evaluate the number of overlapped regions by setting different coverage cutoffs.
  for (fractionCut in seq(0,1,by = 0.1)) {
    colname <- glue::glue("{as.character(fractionCut*100)}%")
    list_overlap[[colname]] <-
      # bed_coverage outputs the same number of lines as
      # the "x" input table
      # in the case, the "pivot_simu" below,
      # which is the simulated histone marker.
      # simulate and evalute on the same marker.
      lapply(pivot_simu, bed_coverage, y = query) %>%
      lapply(filter, .frac >= fractionCut) %>%
      sapply(nrow)
  }
  tbl_overlap <- list_overlap %>% bind_rows() %>%
    mutate(group = c("observed", rep("expected", nrow(.) - 1))) %>%
    relocate(group)
  # output table
  if (output_counts) {
    write_tsv(tbl_overlap, sprintf("output/output_histone/tbl_overlap_%s.tsv", pivot_name))
  }

  # visualization
  data <- tbl_overlap %>%
    gather("coverage", "overlapped_pivot_peaks", -group) %>%
    mutate(coverage = factor(coverage, levels = sprintf("%i%%",0:100))) %>%
    mutate(group = factor(group, levels = c("observed", "expected")))
  data_nozero <- data %>%
    filter(coverage != "0%")

  plot_dot <- ggplot(data_nozero, aes(x = coverage, y = overlapped_pivot_peaks,
                                 label = overlapped_pivot_peaks,
                                 color = group)) +
    geom_point(size = 4)+
    geom_point(size = 4, shape = 1, color = "black")+
    geom_hline(yintercept = data$overlapped_pivot_peaks[1])+
    scale_color_manual(values = c("#EF8887", "#96BBFA")) +
    theme_linedraw() +
    theme(legend.position =  c(0.8, 0.7))+
    ggtitle(sprintf("Pivot: %s at %s", pivot_name, loci_annotation))

  return(plot_dot)
}


# Run the test and visualization ------------------------------------------
# 1 test single locus
TestCoLocalHistone(loci_annotation = "intron", pivot_name = "H3K9me3")
# 2 test all loci
loci_list <-c("Intergenic",   "TTS", "exon" , "intron", "promoter-TSS")
plog_list <- lapply(loci_list, TestCoLocalHistone, pivot_name = "H3K9me3")
patchwork::wrap_plots(plog_list, nrow = 2)



