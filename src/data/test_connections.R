
library(tidyverse)


euks <-
  # read_tsv("../../tables/18S_all_abunds.txt") %>%
  read_tsv(snakemake@input[[1]]) %>%
  rename(Eukaryotic_taxonomy = Taxonomy,
         Eukaryotic_OTU_count = Count) %>%
  filter(str_detect(Eukaryotic_taxonomy, "Mock")) %>%
  separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
  group_by(Sample, Eukaryotic_taxonomy) %>%
  summarise(Eukaryotic_OTU_count = sum(Eukaryotic_OTU_count))


bacts <-
  # read_tsv("../../tables/16S_all_abunds.txt") %>%
  read_tsv(snakemake@input[[2]]) %>%
  rename(Bacterial_taxonomy = Taxonomy,
         Bacterial_OTU_count = Count) %>%
  filter(str_detect(Bacterial_taxonomy, "Mock")) %>%
  separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
  group_by(Bacterial_taxonomy) %>%
  summarise(Bacterial_OTU_count = sum(Bacterial_OTU_count))


conns <-
  # read_tsv("../../tables/all_euk_bact_connections.txt") %>%
  read_tsv(snakemake@input[[3]]) %>%
  rename(Connection_count = Count) %>%
  filter(str_detect(Bacterial_taxonomy, "Mock"),
         str_detect(Eukaryotic_taxonomy, "Mock")) %>%
  separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
  separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
  group_by(Sample, Eukaryotic_taxonomy, Bacterial_taxonomy) %>%
  summarise(Connection_count = sum(Connection_count))


png(snakemake@output[[1]], units="in", width=10, height=10, res=300)
conns %>%
  left_join(euks) %>%
  left_join(bacts) %>%
  mutate(Smp = Sample) %>%
  filter(!str_detect(Sample, "Nomock")) %>%
  unite(Eukaryotic_taxonomy,
        Sample,
        Eukaryotic_taxonomy,
        Eukaryotic_OTU_count) %>%
  unite(Bacterial_taxonomy,
        Smp,
        Bacterial_taxonomy,
        Bacterial_OTU_count) %>%
  ggplot(aes(Eukaryotic_taxonomy,
             Bacterial_taxonomy,
             fill=Connection_count)) +
  geom_tile() +
  geom_text(aes(label=Connection_count), size=0.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


