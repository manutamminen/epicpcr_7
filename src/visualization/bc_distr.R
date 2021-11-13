
library(tidyverse)


bact_bc <- read_tsv(snakemake@input[[1]])
euk_bc <- read_tsv(snakemake@input[[2]])

# bact_bc <- read_tsv("../../data/final/16S_bc_tax.txt")
# euk_bc <- read_tsv("../../data/final/18S_bc_tax.txt")


bact_bc_counts <-
  bact_bc %>%
  group_by(Sample, BC) %>%
  summarise(n = sum(Count)) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "16S",
         Sample = str_replace(Sample, "16S", ""))


euk_bc_counts <-
  euk_bc %>%
  group_by(Sample, BC) %>%
  summarise(n = sum(Count)) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "18S",
         Sample = str_replace(Sample, "18S", ""))


bc_counts <-
  bind_rows(bact_bc_counts, euk_bc_counts) %>%
  mutate(Sample = as.factor(Sample))


# levels(bc_counts$Sample) <-
#   c("Biolstd_Mock_Mag",
#   "Biolstd_Mock_Nomag",
#   "Biolstd_Nomock_Mag",
#   "Biolstd_Nomock_Nomag",
#   "BiolstdWW_Mock_Mag",
#   "BiolstdWW_Mock_Nomag",
#   "BiolstdWW_Nomock_Mag",
#   "BiolstdWW_Nomock_Nomag",
#   "WW_Mock_Mag",
#   "WW_Mock_Nomag",
#   "WW_Nomock_Mag",
#   "WW_Nomock_Nomag")


png(snakemake@output[[1]], units="in", width=5, height=5, res=300)
bc_counts %>%
  # separate(Sample, c("Smp", "Mock", "Mag"), sep="_") %>%
  ggplot(aes(x = id, y = n)) +
  geom_density(stat = "identity") +
  # facet_grid(Smp + Mock + Mag ~ Ribo, scales = "free") +
  facet_grid(Sample ~ Ribo, scales = "free") +
  # scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()
# ggsave(snakemake@output[[1]], last_plot())


bact_bc_taxa <-
  bact_bc %>%
  select(-Count) %>%
  count(Sample, BC) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "16S",
         Sample = str_replace(Sample, "16S", ""))


euk_bc_taxa <-
  euk_bc %>%
  select(-Count) %>%
  count(Sample, BC) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Junk"), sep="_") %>%
  select(-Junk) %>%
  mutate(Ribo = "18S",
         Sample = str_replace(Sample, "18S", ""))


bc_taxa <-
  bind_rows(bact_bc_taxa, euk_bc_taxa) %>%
  mutate(Sample = as.factor(Sample))


# levels(bc_taxa$Sample) <-
#   c("Biolstd_Mock_Mag",
#   "Biolstd_Mock_Nomag",
#   "Biolstd_Nomock_Mag",
#   "Biolstd_Nomock_Nomag",
#   "BiolstdWW_Mock_Mag",
#   "BiolstdWW_Mock_Nomag",
#   "BiolstdWW_Nomock_Mag",
#   "BiolstdWW_Nomock_Nomag",
#   "WW_Mock_Mag",
#   "WW_Mock_Nomag",
#   "WW_Nomock_Mag",
#   "WW_Nomock_Nomag")

png(snakemake@output[[2]], units="in", width=5, height=5, res=300)
bc_taxa %>%
  # separate(Sample, c("Smp", "Mock", "Mag"), sep="_") %>%
  ggplot(aes(x = id, y = n)) +
  geom_density(stat = "identity") +
  # facet_grid(Smp + Mock + Mag ~ Ribo, scales = "free") +
  facet_grid(Sample ~ Ribo, scales = "free") +
  # scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()
# ggsave(snakemake@output[[2]], last_plot())

