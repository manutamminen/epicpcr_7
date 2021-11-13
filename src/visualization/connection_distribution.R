
library(tidyverse)

conns <-
read_tsv(snakemake@input[[1]]) %>%
  # read_tsv("../../tables/nonmock_euk_bact_connections.txt") %>%
  arrange(Sample, desc(Count)) %>%
  group_by(Sample) %>%
  mutate(Ix = row_number(),
         Sample = as.factor(Sample))


levels(conns$Sample) <-
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

png(snakemake@output[[2]], units="in", width=5, height=5, res=300)
conns %>%
  separate(Sample, c("Smp", "Mock", "Mag"), sep="_") %>%
  ggplot(aes(x=Ix, y=Count)) +
  geom_density(stat = "identity") +
  facet_grid(Smp + Mock + Mag ~ ., scales="free") +
  scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))
dev.off()
# ggsave(snakemake@output[[1]])
