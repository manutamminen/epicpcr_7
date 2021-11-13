
library(tidyverse)

bact <-
  read_tsv(snakemake@input[[1]])
  # read_tsv("../../data/final/16S_bc_tax.txt")

euk <-
  read_tsv(snakemake@input[[2]])
  # read_tsv("../../data/final/18S_bc_tax.txt")


bact_bc_percs <-
  bact %>%
  separate(Taxonomy, "Mock", sep="_", remove=FALSE) %>%
  separate(Mock, c(NA, "Mock"), sep=":") %>%
  mutate(Mock = ifelse(str_detect(Mock, "Mock"), Mock, "Biol")) %>%
  select(-Taxonomy) %>%
  pivot_wider(id_cols=c("Sample", "BC"),
              names_from="Mock",
              values_from="Count",
              values_fn=sum,
              values_fill=0) %>%
  mutate(across(3:9, ~ .x > 0)) %>%
  pivot_longer(cols=-c("Sample", "BC"),
               names_to="Mock",
               values_to="Val") %>%
  group_by(Sample, Mock, Val) %>%
  summarise(Count = n()) %>%
  mutate(Tot = sum(Count), Prop = Count/Tot*100) %>%
  filter(Val) %>%
  select(-Tot, -Val) %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Sample = str_replace(Sample, "16S", ""),
         Type = "Bacteria")



euk_bc_percs <-
  euk %>%
  separate(Taxonomy, "Mock", sep="_", remove=FALSE) %>%
  separate(Mock, c(NA, "Mock"), sep=":") %>%
  mutate(Mock = ifelse(str_detect(Mock, "Mock"), Mock, "Biol")) %>%
  select(-Taxonomy) %>%
  pivot_wider(id_cols=c("Sample", "BC"),
              names_from="Mock",
              values_from="Count",
              values_fn=sum,
              values_fill=0) %>%
  mutate(across(3:9, ~ .x > 0)) %>%
  pivot_longer(cols=-c("Sample", "BC"),
               names_to="Mock",
               values_to="Val") %>%
  group_by(Sample, Mock, Val) %>%
  summarise(Count = n()) %>%
  mutate(Tot = sum(Count), Prop = Count/Tot*100) %>%
  filter(Val) %>%
  select(-Tot, -Val) %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Sample = str_replace(Sample, "18S", ""),
         Type = "Eukaryota")


percs <-
  bind_rows(bact_bc_percs, euk_bc_percs) %>%
  mutate(Sample = as.factor(Sample))

# levels(percs$Sample) <-
#   c("Biolstd_Mock_Mag", "Biolstd_Mock_Nomag", "Biolstd_Nomock_Mag",
#     "Biolstd_Nomock_Nomag", "BiolstdWW_Mock_Mag", "BiolstdWW_Mock_Nomag",
#     "BiolstdWW_Nomock_Mag", "BiolstdWW_Nomock_Nomag", "WW_Mock_Mag",
#     "WW_Mock_Nomag", "WW_Nomock_Mag", "WW_Nomock_Nomag")


percs %>%
  # separate(Sample, c("Smp_type", "Mock_type", "Mag_type"), sep="_") %>%
  # ggplot(aes(x=Smp_type, y=Count, fill=Mock)) +
  ggplot(aes(x=Sample, y=Count, fill=Mock)) +
  geom_bar(stat="identity", position="dodge") +
  facet_grid(. ~ Type) +
  # facet_grid(Mock_type + Mag_type ~ Type) +
  coord_flip()
ggsave(snakemake@output[[1]])

percs %>%
  # separate(Sample, c("Smp_type", "Mock_type", "Mag_type"), sep="_") %>%
  # ggplot(aes(x=Smp_type, y=Prop, fill=Mock)) +
  ggplot(aes(x=Sample, y=Prop, fill=Mock)) +
  geom_bar(stat="identity", position="dodge") +
  facet_grid(. ~ Type) +
  # facet_grid(Mock_type + Mag_type ~ Type) +
  coord_flip()
ggsave(snakemake@output[[2]])

