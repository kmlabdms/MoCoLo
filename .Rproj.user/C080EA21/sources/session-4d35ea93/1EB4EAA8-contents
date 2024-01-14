options(bedtools.path = "/usr/bin/")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0 Libraries ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(GenomicRanges)
library(scales)
library(tictoc)
library(valr)
library(furrr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(ggpubr)
library(vroom)
library(ggprism)
library(parallel)
library(patchwork)
library(ComplexHeatmap)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename


theme_set(
  theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "bottom") +
    theme(
      legend.position = c(0.99, 0.99),
      legend.justification = c('right', "top"),
      legend.background = element_blank(),
      legend.key.size = unit(2, "mm")
    )
)



