
library(tidyverse)
library(forcats)


biol_stds <-
  # read_tsv("../../data/final/18S_bc_tax.txt") %>%
read_tsv(snakemake@input[[1]]) %>%
  select(-BC, -Count) %>%
  filter(!str_detect(Taxonomy, "Mock")) %>%
  mutate(Taxonomy = case_when(str_detect(Taxonomy, "Rhodomonas") ~ "Rhodomonas",
                              str_detect(Taxonomy, "Cryptomonas") ~ "Chilomonas",
                              str_detect(Taxonomy, "Ochromonas") ~ "Ochromonas",
                              TRUE ~ "Other" ),
         Taxonomy = fct_relevel(Taxonomy,
                                c("Rhodomonas",
                                  "Chilomonas",
                                  "Ochromonas",
                                  "Other"))) %>%
  mutate(Sample = as.factor(Sample))


levels(biol_stds$Sample) <-
  c("Biolstd_Mock_Mag",
  "Biolstd_Mock_Nomag",
  "Biolstd_Nomock_Mag",
  "Biolstd_Nomock_Nomag",
  "BiolstdWW_Mock_Mag",
  "BiolstdWW_Mock_Nomag",
  "BiolstdWW_Nomock_Mag",
  "BiolstdWW_Nomock_Nomag",
  "WW_Mock_Mag",
  "WW_Mock_Nomag",
  "WW_Nomock_Mag",
  "WW_Nomock_Nomag")

biol_stds %>%
  count(Sample, Taxonomy) %>%
  separate(Sample, c("Smp", "Mock", "Mag"), sep="_") %>%
  ggplot(aes(x=Smp, y=n, fill=Taxonomy)) +
         geom_bar(stat="identity", position="dodge") +
         facet_grid(Mock + Mag ~ .) +
         coord_flip()

ggsave(snakemake@output[[1]])
