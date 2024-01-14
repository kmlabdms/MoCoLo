options(bedtools.path = "/opt/homebrew/bin/")
library(bedtoolsr)
library(GenomicRanges)
library(tidyverse)
library(scales)
library(tictoc)
library(valr)
library(furrr)
library(liftOver)
library(BSgenome.Hsapiens.UCSC.hg19)
select <- dplyr::select
filter <- dplyr::filter


# data
# read 8-oxoG region ----
# convert hg18 to hg19
hg19.size <- read.delim("data-raw/hg19.size", header=FALSE)[1:24,] %>%
  set_names(c('chrom', 'size'))
Hg18_To_Hg19 <- function (file_bed) {

  Bed_To_GRange <-function(i){
    GRanges(seqnames = i$chr,
            ranges = IRanges(start = i$start,end = i$end))
  }
  file_GRange <- Bed_To_GRange(file_bed)

  chromsome <- rtracklayer::import.chain('data-raw/hg18ToHg19.over.chain')
  gr_list <- liftOver(file_GRange, chromsome)
  file_GRange_hg19 <-
    liftOver(file_GRange, chromsome) %>%
    as_tibble() %>%
    select(seqnames,start,end)
}
MCF10A_raw <-
  read_delim("data/first_case/GSE100234_DIP_8oxodG_peaks_MCF10A.bed",
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(`-log10(pvalue)`>5, fold_enrichment>7)
# Perform filter and convert from hg18 to hg19
MCF10A_filtered <-  Hg18_To_Hg19(file_bed = MCF10A_raw)

# host feature (Clean)
oxidative <-  MCF10A_filtered %>% set_names("chrom", "start", "end")



# 1  The observed ----
oxidative_nuc <- readRDS("data/first_case/MCFT.nuc.rds") %>%
  set_names("chrom", "start", "end", "nG", "length", "pG")  %>%
  mutate(length = length+1)
regions_real_raw <- oxidative_nuc %>% arrange(desc(length), desc(nG),chrom, start, end)
# Two features was removed here. 52,489 -> 52,487
# Feature (chr17:4756097-4756097) has length = 0, Skipping.
# Feature (chr9:141122440-141122440) has length = 0, Skipping.


# 2 Run simulation -----
## 2.1 Build peak pools
region_real <- regions_real_raw %>%
  mutate(peak_id = paste0("peak", 1:nrow(.)))

# Simulation 8-oxo-G with G%
region_real_simu <-
  SimulateRegions(
    region_real,
    simu_fuc = SimulatePool,
    n_test = 100, cores = 50,
  )
saveRDS(region_real_simu,
        sprintf("data/first_case/oxog_realG%%_simu.%sx%s.rds",
                nrow(region_real_simu[[1]]), length(region_real_simu)))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Simulation (no restrict on G%) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
real <- oxidative_all[[1]] %>% mutate(group = "0")

region_real_simu_noG <-
  map(1:20, ~ valr::bed_shuffle(
    x = real %>% select(chrom, start, end),
    genome = hg19.size,
    seed = .x
  ))

region_real_simu_noG_with_metric <-
  region_real_simu_noG %>%
  bind_rows(.id = "group") %>%
  CalculateG() %>%
  arrange(desc(simu_ng)) %>%
  mutate(simu_pg = simu_ng / length)

saveRDS(
  region_real_simu_noG_with_metric,
  sprintf(
    "data/first_case/region_noG%%_noChr_simu.%sx%s.rds",
    nrow(region_real_simu_noG[[1]]), length(region_real_simu_noG)
  )
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. (Main plot) Check the G percent of oxoG -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
col_grey = rep("grey60", 5) %>% set_names(1:5)
oxidative_all
real = oxidative_all[[1]] %>% mutate(group = "0")
simu = oxidative_all[2:6] %>% bind_rows(.id = "group")
data = bind_rows(
  real %>% select(group, gcontent = pG, length),
  simu %>% select(group, gcontent = simu_pg, length),
)
dist_len =
  ggplot(data, aes(length)) +
  xlim(0, 1000) +
  geom_density(aes(color = group)) +
  guides(color = guide_legend(title = "Obs.(0)\nExp.(>1)")) +
  ggtitle("Length distribution") +
  labs(x = "Length(bases)") +
  #scale_color_manual(values = c("0"="red",col_grey)) +
  scale_size_manual(values = c(2,1,1,1,1,1)) +
  theme(axis.title = element_text(face = "bold"))
dist_gcontent =
  ggplot(data, aes(gcontent)) +
  ylim(0, 8) +
  geom_density(aes(color = group)) +
  guides(color = guide_legend(title = "Obs.(0)\nExp.(>1)")) +
  ggtitle("G% distribution") +
  labs(x = "G%")+
  scale_color_manual(values = c("0"="red",col_grey)) +
  scale_x_continuous(labels = scales::percent, limits = c(0,0.6)) +
  theme(axis.title = element_text(face = "bold"))
(oxog_distribution_withG = dist_gcontent/dist_len)


col_grey = rep("grey60", 5) %>% set_names(1:5)
simu_length = read_rds("data/first_case/region_noG%_noChr_simu.50027x20.rds")
data = bind_rows(
  real %>% select(group, gcontent = pG, length),
  simu_length %>%
    filter(group %in% 1:5) %>%
    select(group, gcontent = simu_pg, length),
)
dist_len =
  ggplot(data, aes(length)) +
  xlim(0, 1000) +
  geom_density(aes(color = group)) +
  guides(color = guide_legend(title = "Obs.(0)\nExp.(>1)")) +
  ggtitle("Length distribution") +
  labs(x = "Length(bases)") +
  #scale_color_manual(values = c("0"="red",col_grey)) +
  theme(axis.title = element_text(face = "bold"))
dist_gcontent =
  ggplot(data, aes(gcontent)) +
  geom_density(aes(color = group)) +
  guides(color = guide_legend(title = "Obs.(0)\nExp.(>1)")) +
  ggtitle("G% distribution") +
  labs(x = "G%")+
  scale_color_manual(values = c("0"="red",col_grey)) +
  scale_x_continuous(labels = scales::percent, limits = c(0,0.6)) +
  theme(axis.title = element_text(face = "bold"))
oxog_distribution_withoutG = dist_gcontent/dist_len

oxog_distribution_withG |
  oxog_distribution_withoutG


