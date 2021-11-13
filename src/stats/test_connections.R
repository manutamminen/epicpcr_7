
library(tidyverse)
library(tidymodels)
library(broom.mixed)
library(dotwhisker)
library(ranger)

euks <-
  read_tsv("../../data/final/18S_bc_tax.txt") %>%
  # read_tsv(snakemake@input[[1]]) %>%
  rename(Eukaryotic_taxonomy = Taxonomy) %>%
  separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
  filter(!str_detect(Sample, "Nomock")) %>%
  select(-Count) %>%
  unique %>%
  group_by(Sample) %>%
  count(Eukaryotic_taxonomy, name = "Eukaryotic_OTU_count") %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Sample = str_replace(Sample, "18S", ""))


bacts <-
  read_tsv("../../data/final/16S_bc_tax.txt") %>%
  # read_tsv(snakemake@input[[2]]) %>%
  rename(Bacterial_taxonomy = Taxonomy) %>%
  separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
  filter(!str_detect(Sample, "Nomock")) %>%
  select(-Count) %>%
  unique %>%
  group_by(Sample) %>%
  count(Bacterial_taxonomy, name = "Bacterial_OTU_count") %>%
  separate(Sample, "Sample", sep="_") %>%
  mutate(Sample = str_replace(Sample, "16S", ""))


connections <-
  read_tsv("../../tables/all_euk_bact_connections.txt") %>%
  # read_tsv(snakemake@input[[3]]) %>%
  rename(Connection_count = Count) %>%
  separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
  separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
  group_by(Sample, Eukaryotic_taxonomy, Bacterial_taxonomy) %>%
  summarise(Connection_count = sum(Connection_count)) %>%
  left_join(euks) %>%
  left_join(bacts) %>%
  filter(!str_detect(Sample, "Nomock")) %>%
  ungroup


mock_connections <-
  connections %>%
  filter(Bacterial_taxonomy %in% c("d:Mock1", "d:Mock2", "d:Mock3"),
         Eukaryotic_taxonomy %in% c("d:Mock1", "d:Mock2", "d:Mock3"))


# biol_connections <-
#   connections %>%
#   filter(!str_detect(Bacterial_taxonomy, "Mock"),
#          !str_detect(Eukaryotic_taxonomy, "Mock"))


mock_lm_model <- function(mock_connections)
  linear_reg() %>%
    set_engine("lm") %>%
    fit(Connection_count ~ Bacterial_OTU_count * Eukaryotic_OTU_count, data = mock_connections)


lm_model_preds <-
  mock_connections %>%
  nest(data = c(Eukaryotic_taxonomy, Bacterial_taxonomy, Connection_count,
                Eukaryotic_OTU_count, Bacterial_OTU_count)) %>%
  mutate(Null_model = map(data, possibly(mock_lm_model, otherwise = NA))) %>%
  select(-data) %>%
  left_join(nest(connections,
                 Biol_data = c(Eukaryotic_taxonomy, Bacterial_taxonomy, Connection_count,
                               Eukaryotic_OTU_count, Bacterial_OTU_count)), by="Sample") %>%
  mutate(Predicted_means = map2(Null_model, Biol_data, ~predict(.x, .y)),
         Predicted_cis = map2(Null_model, Biol_data, ~predict(.x, .y, type = "conf_int"))) %>%
  select(-Null_model) %>%
  unnest(-Sample) %>%
  ungroup %>%
  mutate(Status = case_when(Connection_count > .pred_upper ~ "Above",
                            Connection_count < .pred_lower ~ "Below",
                            TRUE ~ "Non-significant"),
         Dist = abs(Connection_count - .pred))



> lm_model_preds %>% mutate(Dist = abs(Connection_count - .pred)) %>% arrange(desc(Dist)) %>% filter(!str_detect(Eukaryotic_taxonomy, "Mock"), !str_detec
t(Bacterial_taxonomy, "Mock")) %>% ggplot(aes(x=Dist)) + geom_density() + facet_grid(Sample~.)

lm_model_preds %>%
  mutate(Bacterial_taxonomy = case_when(str_detect(Bacterial_taxonomy, "Rhodomonas salina") ~ "Rhodo_chloroplast",
                                        str_detect(Bacterial_taxonomy, "Cryptomonas paramecium") ~ "Chilo_chloroplast",
                                        TRUE ~ Bacterial_taxonomy)) %>%
filter(Bacterial_taxonomy == "Chilo_chloroplast" | Bacterial_taxonomy == "Rhodo_chloroplast" | str_detect(Bacterial_taxonomy, "Mock")) %>%
group_by(Sample, Bacterial_taxonomy) %>%
summarise(sum(Bacterial_OTU_count)) %>%
data.frame


count(Eukaryotic_taxonomy) %>%
arrange(desc(n))



%>%
filter(str_detect(Bacterial_taxonomy, "lorop")) %>%
group_by(Bacterial_taxonomy) %>% summarise(Sum = sum(Bacterial_OTU_count)) %>% arrange(desc(Sum)) %>% data.frame


mock_rf_model <- function(mock_connections)
  rand_forest(mode = "regression") %>%
    set_engine("ranger") %>%
    fit(Connection_count ~ Bacterial_OTU_count * Eukaryotic_OTU_count, data = mock_connections)


rf_model_preds <-
  mock_connections %>%
  nest(data = c(Eukaryotic_taxonomy, Bacterial_taxonomy, Connection_count,
                Eukaryotic_OTU_count, Bacterial_OTU_count)) %>%
  mutate(Null_model = map(data, possibly(mock_rf_model, otherwise = NA))) %>%
  select(-data) %>%
  left_join(nest(connections,
                 Biol_data = c(Eukaryotic_taxonomy, Bacterial_taxonomy, Connection_count,
                               Eukaryotic_OTU_count, Bacterial_OTU_count)), by="Sample")

rf_model_preds %>%
  mutate(Predicted_means = map2(Null_model, Biol_data, ~predict(.x, .y))) %>%
  select(-Null_model) %>%
  unnest(-Sample) %>%
  ungroup

# mock_recipe <-
#   recipe(Connection_count ~ .,
#          data = select(mock_connections,
#                        Connection_count,
#                        Eukaryotic_OTU_count,
#                        Bacterial_OTU_count)) %>%
#   step_log(all_numeric()) %>%
#   step_normalize(all_numeric())




mock_glm_model <- function(mock_connections)
  linear_reg(penalty = 0.001, mixture = 0.5) %>%
    set_engine("glmnet") %>%
    fit(Connection_count ~ Bacterial_OTU_count * Eukaryotic_OTU_count, data = mock_connections)


glm_model_preds <-
  mock_connections %>%
  nest(data = c(Eukaryotic_taxonomy, Bacterial_taxonomy, Connection_count,
                Eukaryotic_OTU_count, Bacterial_OTU_count)) %>%
  mutate(Null_model = map(data, possibly(mock_glm_model, otherwise = NA))) %>%
  select(-data) %>%
  left_join(nest(connections,
                 Biol_data = c(Eukaryotic_taxonomy, Bacterial_taxonomy, Connection_count,
                               Eukaryotic_OTU_count, Bacterial_OTU_count)), by="Sample")

glm_model_preds %>%
  mutate(Predicted_means = map2(Null_model, Biol_data, ~predict(.x, .y))) %>%
  select(-Null_model) %>%
  unnest(-Sample) %>%
  ungroup


# tidy(mock_model) %>%
#   dwplot(dot_args = list(size = 2, color = "black"),
#          whisker_args = list(color = "black"),
#          vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2))


# biol_connections <-
#   connections %>%
#   filter(!str_detect(Bacterial_taxonomy, "Mock"),
#          !str_detect(Eukaryotic_taxonomy, "Mock"))


mean_pred <-
  predict(mock_model,
          new_data = connections)


conf_int_pred <-
  predict(mock_model,
          new_data = connections,
          type = "conf_int")


connections_predictions <-
  connections %>%
  ungroup %>%
  bind_cols(mean_pred)

%>%
  mutate(Dist = abs(Connection_count - .pred)) %>%
  arrange(desc(Dist)) %>%
  filter(!str_detect(Eukaryotic_taxonomy, "Mock"),
         !str_detect(Bacterial_taxonomy, "Mock"))


  %>%
  filter(str_detect(Bacterial_taxonomy, "loropl")) %>%
  arrange(desc(Connection_count))


  head %>%
  data.frame



# connections_predictions <-
#   connections %>%
#   bind_cols(mean_pred) %>%
#   bind_cols(conf_int_pred)


connections_predictions %>%
  ungroup %>%
  mutate(Status = case_when(Connection_count > .pred_upper ~ "Above",
                            Connection_count < .pred_upper ~ "Below",
                            TRUE ~ "Non-significant")) %>%
  filter(Status != "Below")

  count(Sample, Status)


  mutate(Above = Connection_count > .pred_upper,
         Below = Connection_count < .pred_lower)

  xyzk

# euk_mocks <-
#   read_tsv("../../data/final/18S_bc_tax.txt") %>%
#   # read_tsv(snakemake@input[[1]]) %>%
#   rename(Eukaryotic_taxonomy = Taxonomy) %>%
#   separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
#   filter(str_detect(Eukaryotic_taxonomy, "Mock"),
#     !str_detect(Sample, "Nomock")) %>%
#   select(-Count) %>%
#   unique %>%
#   group_by(Sample) %>%
#   count(Eukaryotic_taxonomy, name = "Eukaryotic_OTU_count") %>%
#   separate(Sample, "Sample", sep="_") %>%
#   mutate(Sample = str_replace(Sample, "18S", ""))
#
#
# bact_mocks <-
#   read_tsv("../../data/final/16S_bc_tax.txt") %>%
#   # read_tsv(snakemake@input[[2]]) %>%
#   rename(Bacterial_taxonomy = Taxonomy) %>%
#   separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
#   filter(str_detect(Bacterial_taxonomy, "Mock"),
#     !str_detect(Sample, "Nomock")) %>%
#   select(-Count) %>%
#   unique %>%
#   group_by(Sample) %>%
#   count(Bacterial_taxonomy, name = "Bacterial_OTU_count") %>%
#   separate(Sample, "Sample", sep="_") %>%
#   mutate(Sample = str_replace(Sample, "16S", ""))
#
#
# mock_connections <-
#   read_tsv("../../tables/all_euk_bact_connections.txt") %>%
#   # read_tsv(snakemake@input[[3]]) %>%
#   rename(Connection_count = Count) %>%
#   filter(str_detect(Bacterial_taxonomy, "Mock"),
#          str_detect(Eukaryotic_taxonomy, "Mock")) %>%
#   separate(Eukaryotic_taxonomy, "Eukaryotic_taxonomy", sep="_") %>%
#   separate(Bacterial_taxonomy, "Bacterial_taxonomy", sep="_") %>%
#   group_by(Sample, Eukaryotic_taxonomy, Bacterial_taxonomy) %>%
#   summarise(Connection_count = sum(Connection_count)) %>%
#   left_join(euk_mocks) %>%
#   left_join(bact_mocks) %>%
#   filter(!str_detect(Sample, "Nomock"))
#

mock_model <-
  linear_reg() %>%
  set_engine("lm") %>%
  fit(Connection_count ~ Bacterial_OTU_count * Eukaryotic_OTU_count, data = mock_connections)


tidy(mock_model) %>%
  dwplot(dot_args = list(size = 2, color = "black"),
         whisker_args = list(color = "black"),
         vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2))










smp_bool <-
  function(length, true_proport)
    sample(c(rep(TRUE, true_proport),
             rep(FALSE, (length - true_proport))),
           length, replace=FALSE)

emulsion_model <-
  function(n_droplets,
           bcs,
           mock1,
           mock2,
           mock3)
    tibble(bcs = smp_bool(n_droplets, bcs),
           mock1 = smp_bool(n_droplets, mock1),
           mock2 = smp_bool(n_droplets, mock2),
           mock3 = smp_bool(n_droplets, mock3))


emulsion_model(1e8, 1e4, 3e5, 3e4, 3e3) %>%
  filter(bcs, mock1)


conns %>%
  ungroup %>%
  mutate(Lambda = Eukaryotic_OTU_count * Bacterial_OTU_count)






library(tidymodels)
library(readr)
library(broom.mixed)
library(dotwhisker)


urchins <-
  # Data were assembled for a tutorial 
  # at https://www.flutterbys.com.au/stats/tut/tut7.5a.html
  read_csv("https://tidymodels.org/start/models/urchins.csv") %>% 
  # Change the names to be a little more verbose
  setNames(c("food_regime", "initial_volume", "width")) %>% 
  # Factors are very helpful for modeling, so we convert one column
  mutate(food_regime = factor(food_regime, levels = c("Initial", "Low", "High")))


urchins

lm_fit <-
  linear_reg() %>%
  set_engine("lm") %>%
  fit(width ~ initial_volume * food_regime, data = urchins)


