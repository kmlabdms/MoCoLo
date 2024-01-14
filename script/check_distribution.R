# This script is used check the length distribution of
# overlapped region
# for both case 1 and case2



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Case2: Non-B/8-oxoG ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1.1 Running the mapping ----
oxidative <- oxidative_all[1]
map_shuffle_tbl <- MapGenomewide(
  oxidative,
  nonbDB,
  n_core = 36L
)
host_feature = oxidative[[1]] %>% mutate(Type = "8-oxo-dG")
nest_feature = nonbDB %>% bind_rows() %>% mutate(length = end - start)

## 1.2 Distribution of Overlap ----
bt_intersect <- bed_intersect(host_feature, nest_feature) %>%
  filter(.overlap != 0) %>%
  mutate(nonb_length = end.y - start.y)
dplyr::count(bt_intersect, .overlap) %>% arrange(desc(n)) %>%
  dplyr::rename(overlap_legnth = .overlap)
p_inter <- ggplot(bt_intersect, aes(.overlap)) +
  geom_freqpoly(aes(color = Type.y), binwidth = 1) +
  guides(color = guide_legend(title = "Non-B types")) +
  ggtitle("Frequency of intersections between non-B and 8-oxo-dG  (full)") +
  labs(x = "overlap length")
p_inter_lim <- p_inter + xlim(0, 100) +
  guides(color = guide_legend(title = "Non-B types")) +
  ggtitle("Frequency of intersections between non-B and 8-oxo-dG (1,100)")+
  labs(x = "overlap length")
p_inter_fill <- ggplot(bt_intersect, aes(.overlap)) +
  geom_histogram(aes(fill = Type.y),
                 binwidth = 1, position = "fill") +
  xlim(0, 100) +
  guides(fill = guide_legend(title = "Non-B types")) +
  ggtitle("Conditional density plots (1,100)")+
  labs(x = "overlap length")

p_inter / p_inter_lim / p_inter_fill


## 1.3 non-B length distribution -----
p_host <- ggplot(host_feature, aes(length)) +
  xlim(0, 2000) +
  geom_freqpoly(color = "blue", binwidth = 1) +
  guides(color = guide_legend(title = "Feature")) +
  ggtitle("Length distribution of 8-oxo-dG") +
  theme(legend.position = "none") +
  labs(x = "length")
p_host
p_nest <- ggplot(nest_feature, aes(length)) +
  xlim(0, 100) +
  geom_freqpoly(color = "blue", binwidth = 1) +
  guides(color = guide_legend(title = "Feature")) +
  ggtitle("Length distribution of non-B forming motifs") +
  labs(x = "length")
p_inter <- ggplot(bt_intersect, aes(.overlap)) +
  geom_freqpoly(color = "red",binwidth = 1) +
  guides(color = guide_legend(title = "Non-B types"))  +
  xlim(0, 100) +
  ggtitle("Length distribution of intersections")+
  labs(x = "overlap length")

p_host / p_nest / p_inter





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Case2: Histone length distribution ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h3k9me3 <- read_csv("data/second_case/peaks/H3K9me3_R1_peaks.csv") %>% clean_names() %>%
  left_join(h3k9Anno) %>% dplyr::rename(chrom = chr) %>%
  mutate(chrom = paste0("chr", chrom))
h4k20me3 <- read_csv("data/second_case/peaks/H4K20me3_R1_peaks.csv") %>% clean_names() %>%
  left_join(h4k20Anno) %>% dplyr::rename(chrom = chr) %>%
  mutate(chrom = paste0("chr", chrom))

## 2.1 Distribution of Overlap ----
bt_intersect_histone <- bed_intersect(h3k9me3, h4k20me3) %>%
  filter(.overlap != 0) %>%
  mutate(Type.y = "histones" )
p_inter <- ggplot(bt_intersect, aes(.overlap)) +
  geom_freqpoly(aes(color = Type.y), binwidth = 1) +
  guides(color = guide_legend(title = "histones")) +
  ggtitle("Frequency of intersections between H3K9me3 and H4K20me3  (full)") +
  labs(x = "overlap length")
p_inter_lim <- p_inter + xlim(0, 2500) +
  guides(color = guide_legend(title = "histones")) +
  ggtitle("Frequency of intersections between H3K9me3 and H4K20me3  (1,2500)") +
  labs(x = "overlap length")
p_inter_xylim <- p_inter +   xlim(0, 2500) + ylim(0, 25) +
  guides(color = guide_legend(title = "histones")) +
  ggtitle("Frequency of intersections between H3K9me3 and H4K20me3  (1,2500)") +
  labs(x = "overlap length")

count(bt_intersect, .overlap) %>% arrange(desc(n)) %>%
  dplyr::rename(overlap_legnth = .overlap)

p_inter / p_inter_lim / p_inter_xylim

## 2.2 Distribution of histone length ----
p_host <- ggplot(h3k9me3, aes(length)) +
  xlim(0, 5000) + ylim(0, 50) +
  geom_freqpoly(color = "blue", binwidth = 1) +
  guides(color = guide_legend(title = "Feature")) +
  ggtitle("Length distribution of H3K9me3") +
  theme(legend.position = "none") +
  labs(x = "length")
p_host
p_nest <- ggplot(h4k20me3, aes(length)) +
  xlim(0, 5000) + ylim(0, 50) +
  geom_freqpoly(color = "blue", binwidth = 1) +
  guides(color = guide_legend(title = "Feature")) +
  ggtitle("Length distribution of H4K20me3") +
  labs(x = "length")
p_inter <- ggplot(bt_intersect_histone, aes(.overlap)) +
  geom_freqpoly(color = "red",binwidth = 1) +
  xlim(0, 5000) + ylim(0,30) +
  ggtitle("Length distribution of intersections")+
  labs(x = "overlap length")

p_host / p_nest / p_inter








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. check the Gpercent distribution of nonB -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
real = oxidative_all[[1]] %>% mutate(group = "0")

real_oxog_gdata = bind_rows(
  real %>% select(group, gcontent = pG, length)
)
nonbDB
all_nonb = nonbDB %>% bind_rows(.id = "group")

all_nonb_property = all_nonb %>%
  mutate(length = nchar(Sequence)) %>%
  mutate(gcontent = str_count(Sequence, "G") / length)
dist_len =
  ggplot(all_nonb_property, aes(length)) +
  geom_density(aes(color = group)) +
  guides(color = guide_legend(title = "group")) +
  ggtitle("Length distribution") +
  labs(x = "length") +
  xlim(0, 50)
dist_gcontent =
  ggplot(
    # filter G Content
    all_nonb_property %>% filter(gcontent>0),
    aes(gcontent)) +
  geom_density(aes(color = group)) +
  geom_density(data = real_oxog_gdata, color = "red", alpha =0.1, size = 1) +
  guides(color = guide_legend(title = "Feature types")) +
  ggtitle("gcontent distribution") +
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  theme(axis.title = element_text(face = "bold")) +
  labs(x = "G%", title = "G% distribution")
dist_gcontent/dist_len


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Others ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Number of non-B motifs with G% > 0
mean(all_nonb_property$gcontent>0)
all_nonb_property %>%
  summarise(mean(gcontent > 0), .by = "Type") %>%
  ggplot(aes(y = Type, x = `mean(gcontent > 0)`, fill = Type)) +
  xlim(0,1) +
  theme(legend.position = "none") +
  geom_col()
# check correlation between length and Gpercent (long time to numb)
if (F) {
  all_nonb_property %>%
    filter(group == "Z_DNA_Motif") %>%
    ggplot(aes(x = length, y = gcontent)) +
    geom_point(alpha = 1/20) +
    ggpubr::stat_cor(method = "spearman") +
    xlim(0, 100) +
    ylim(0, 1) +
    geom_smooth(method=lm, se=FALSE)
}








