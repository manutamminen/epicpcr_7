
library(tidyverse)


bact <-
  read_tsv(snakemake@input[[1]]) %>%
  # read_tsv("../../data/final/16S_bc_tax.txt") %>%
  select(-BC) %>%
  group_by(Sample, Taxonomy) %>%
  summarise(Count = sum(Count)) %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Type = "Bacteria",
         Sample = str_replace(Sample, "16S", ""))


euk <-
  read_tsv(snakemake@input[[2]]) %>%
  # read_tsv("../../data/final/18S_bc_tax.txt") %>%
  select(-BC) %>%
  group_by(Sample, Taxonomy) %>%
  summarise(Count = sum(Count)) %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Type = "Eukaryotes",
         Sample = str_replace(Sample, "18S", ""))


bind_rows(bact, euk) %>%
  ggplot(aes(x=Count)) +
  geom_density() +
  facet_grid(Sample~Type, scales="free") +
  scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
ggsave(snakemake@output[[1]])
