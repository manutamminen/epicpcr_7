options(tibble.print_max = 50, tibble.print_min = 50)

library(tidyverse)

euks <-
  read_tsv("../../data/final/18S_bc_tax.txt") %>%
  # read_tsv(snakemake@input[[1]]) %>%
  rename(Eukaryotic_taxonomy = Taxonomy) %>%
  separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
  filter(str_detect(Sample, "Nomock")) %>%
  select(-Count) %>%
  unique %>%
  group_by(Sample) %>%
  count(Eukaryotic_taxonomy, name = "Eukaryotic_OTU_count") %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Sample = str_replace(Sample, "18S", "")) %>%
  filter(!str_detect(Eukaryotic_taxonomy, "Mock")) %>%
  mutate(Eukaryotic_taxonomy = case_when(str_detect(Eukaryotic_taxonomy, "Rhodomonas") ~ "Rhodomonas",
                                        str_detect(Eukaryotic_taxonomy, "Cryptomonas") ~ "Chilomonas",
                                        str_detect(Eukaryotic_taxonomy, "Ochromonas") ~ "Ochromonas",
                                        TRUE ~ Eukaryotic_taxonomy)) %>%
  group_by(Sample, Eukaryotic_taxonomy) %>%
  summarise(Eukaryotic_OTU_count = sum(Eukaryotic_OTU_count)) %>%
  group_by(Eukaryotic_taxonomy) %>%
  mutate(OTU_sum = sum(Eukaryotic_OTU_count)) %>%
  filter(OTU_sum > 100) %>%
  select(-OTU_sum) %>%
  group_by(Sample) %>%
  arrange(desc(Eukaryotic_OTU_count), .by_group = TRUE)



bacts <-
  read_tsv("../../data/final/16S_bc_tax.txt") %>%
  # read_tsv(snakemake@input[[2]]) %>%
  rename(Bacterial_taxonomy = Taxonomy) %>%
  separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
  filter(str_detect(Sample, "Nomock")) %>%
  select(-Count) %>%
  unique %>%
  group_by(Sample) %>%
  count(Bacterial_taxonomy, name = "Bacterial_OTU_count") %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Sample = str_replace(Sample, "16S", "")) %>%
  filter(!str_detect(Bacterial_taxonomy, "Mock")) %>%
  mutate(Bacterial_taxonomy = case_when(str_detect(Bacterial_taxonomy, "hlorop") ~ "Chloroplast",
                                        TRUE ~ Bacterial_taxonomy)) %>%
  group_by(Sample, Bacterial_taxonomy) %>%
  summarise(Bacterial_OTU_count = sum(Bacterial_OTU_count)) %>%
  group_by(Bacterial_taxonomy) %>%
  mutate(OTU_sum = sum(Bacterial_OTU_count)) %>%
  arrange(desc(OTU_sum)) %>%
  filter(OTU_sum > 100) %>%
  select(-OTU_sum) %>%
  group_by(Sample) %>%
  arrange(desc(Bacterial_OTU_count), .by_group = TRUE)





connections <-
  read_tsv("../../tables/all_euk_bact_connections.txt") %>%
  # read_tsv(snakemake@input[[3]]) %>%
  filter(str_detect(Sample, "Nomock"),
         !str_detect(Bacterial_taxonomy, "Mock"),
         !str_detect(Eukaryotic_taxonomy, "Mock")) %>%
  rename(Connection_count = Count) %>%
  separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
  separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
  group_by(Sample, Eukaryotic_taxonomy, Bacterial_taxonomy) %>%
  summarise(Connection_count = sum(Connection_count)) %>%
  filter(!str_detect(Bacterial_taxonomy, "Mock")) %>%
  mutate(Eukaryotic_taxonomy = case_when(str_detect(Eukaryotic_taxonomy, "Rhodomonas") ~ "Rhodomonas",
                                         str_detect(Eukaryotic_taxonomy, "Cryptomonas") ~ "Chilomonas",
                                         str_detect(Eukaryotic_taxonomy, "Ochromonas") ~ "Ochromonas",
                                         TRUE ~ Eukaryotic_taxonomy)) %>%
  mutate(Bacterial_taxonomy = case_when(str_detect(Bacterial_taxonomy, "hlorop") ~ "Chloroplast",
                                        TRUE ~ Bacterial_taxonomy)) %>%
  inner_join(euks) %>%
  inner_join(bacts) %>%
  filter(str_detect(Sample, "Nomock")) %>%
  ungroup


connections %>%
  filter(Bacterial_taxonomy == "Chloroplast") %>%
  group_by(Sample, Eukaryotic_taxonomy) %>%
  summarise(CSum = sum(Connection_count)) %>%
  arrange(desc(CSum), .by_group=TRUE)




connections %>%
  arrange(desc(Connection_count)) %>%
  count(Sample, Bacterial_taxonomy) %>%
  arrange(desc(n))


connections %>%
  filter(!str_detect(Eukaryotic_taxonomy, "Mock"),
         !str_detect(Bacterial_taxonomy, "Mock")) %>%
  filter(str_detect(Bacterial_taxonomy, "lorop")) %>%
  group_by(Sample, Eukaryotic_taxonomy) %>%
  summarise(Sum = sum(Connection_count)) %>%
  mutate(Correct = case_when(str_detect(Eukaryotic_taxonomy, "Cryptomonas") ~ TRUE,
                             str_detect(Eukaryotic_taxonomy, "Rhodomonas") ~ TRUE,
                             str_detect(Eukaryotic_taxonomy, "Ochrophyta") ~ TRUE,
                             TRUE ~ FALSE)) %>%
  arrange(desc(Sum), .by_group = TRUE) %>%
  data.frame


  group_by(Sample, Correct) %>%
  summarise(Sum2 = sum(Sum))
