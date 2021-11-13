
library(tidyverse)

bact <-
  read_tsv("../../data/final/16S_bc_tax.txt")
  # read_tsv(snakemake@input[[1]])

euk <-
  read_tsv("../../data/final/18S_bc_tax.txt")
  # read_tsv(snakemake@input[[2]])

  bact %>%
  separate(Taxonomy, "Mock", sep="_", remove=FALSE) %>%
  separate(Mock, c(NA, "Mock"), sep=":") %>%
  mutate(Mock = ifelse(str_detect(Mock, "Mock"), Mock, "Biol")) %>%
  select(-Taxonomy)
