#####################################################################################
# Imputing missing trait data
################################################################################

## Project: Identifying Climate change vulnerable biodiversity hotspots
##
## Purpose: Impute missing functional trait data
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Load in trait data and put it in workable form
##           3) Set up parameters for the imputation function
##           4) Run missForest algorithm to impute missing traits
##           5) Perform diagnostics on imputed data
##           6) Illustrate diagnostics to ensure imputation was appropriate
##           7) Output imputed data frame
##
##
## Inputs:  data/clean/Trait_data_TRY_Diaz_2022_PNW.xlsx
##          - trait data from the TRY database used in Diaz 2022
##
##          code -
##
##
## Outputs: data/clean/imputed_traits.csv - full list of imputed traits
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023
## _____________________________________________________________________________

# 1) Set-up, load in necessary packages and data-sets ----
library(tidyverse)
library(cowplot)
library(readxl)
library(reshape2)
library(caret)
library(missForest)

# 2) Load in trait data and put it in workable form ----
data_in <- read_csv("data/clean/final_dataset.csv") %>% 
  select("Species_full","Family","LDMC (g/g)", "Plant height (m)", 
         "Nmass (mg/g)", "Leaf area (mm2)", "Woodiness") %>% 
  mutate(Genus = stringr::word(Species_full, 1, 1, sep = " ")) %>%
  mutate(Species = stringr::word(Species_full, 2, 2, sep = " ")) %>%
  mutate(Ssp_var = stringr::word(Species_full, 3, 4, sep = " ")) %>%
  relocate(Genus, Species, Ssp_var)

# Modify original data
traits_df <- data_in %>%
  # expand binomial
  mutate(LeafArea_log = log(`Leaf area (mm2)`)) %>%
  mutate(PlantHeight_log = log(`Plant height (m)`)) %>%
  # mutate(DiasporeMass_log = log(`Diaspore mass (mg)`)) %>%
  # Remove original variables replaced by "logged" versions
  select(-c("Leaf area (mm2)", "Plant height (m)")) #%>% #, "Diaspore mass (mg)")) %>%
  # Remove traits that we won't be using in our study
  # select(-c("LDMC (g/g)", "SSD observed (mg/mm3)", "SSD imputed (mg/mm3)"))

# * Get some basic information about the data set ----
# compute the coverage of each trait
coverage_vals <- traits_df %>% 
  select(-c("Genus":"Ssp_var")) %>%
  melt(id = c("Species_full")) %>% 
  group_by(variable) %>% 
  mutate(count = n()) %>%
  mutate(missing_num = sum(is.na(value))) %>%
  mutate(cover_num = sum(!is.na(value))) %>%
  mutate(coverage = cover_num / count) %>% 
  select(variable, missing_num, count, coverage) %>% 
  unique()

# 3) Set up parameters for the imputation function (missForest) ----

# # missForest parameters
# # maximum number of iterations to be performed given the stopping criterion isn't met
# maxiter = 100
# # number of trees to grow in each forest
# ntree = 1000
# # if 'TRUE', gives additional output between iterations
# verbose = TRUE
# # if 'TRUE', the OOB error is returned for each variable separately
# variablewise = TRUE

# 4) Run missForest algorithm to impute missing traits ----
# get taxonomic columns
taxa <- select(traits_df, Genus)#, Species, Ssp_var) # Get order and family level data
test <- rep(1, nrow(taxa))
taxa <- cbind(taxa, test)

# turn taxonomic groupings into dummy binary variables
dummies <- dummyVars(test ~ ., data = taxa)
taxa <- predict(dummies, newdata = taxa)

# select traits to be imputed
traits_to_impute <- as.data.frame(traits_df) %>%
  select(-c("Genus":"Species_full")) %>%
  mutate(Family = as.factor(Family)) %>%  #%>%
  select(-Family) %>% 
  mutate(Woodiness = as.factor(Woodiness)) #%>%
  # mutate(`Growth Form` = as.factor(`Growth Form`))

# combine dummy vars with traits (comment this line out to not use taxa in imputation)
traits_to_impute <- cbind(traits_to_impute, taxa)

#
# set.seed(82) # for debugging

# run missForest imputation
PNW_imp <- missForest(traits_to_impute,
                      maxiter = 1000, # maximum number of iterations to be performed given the stopping criterion isn't met
                      ntree = 1000, # number of trees to grow in each forest
                      verbose = TRUE, # if 'TRUE', gives additional output between iterations
                      variablewise = TRUE # if 'TRUE', the OOB error is returned for each variable separately
)

# add actual taxonomic columns back in and remove dummy variables
imputed_traits <- cbind(traits_df[, 1:5], PNW_imp$ximp[, 1:5]) %>% 
  mutate("Plant height (m)" = exp(PlantHeight_log), .keep = "unused") %>% 
  mutate("Leaf area (mm2)" = exp(LeafArea_log), .keep = "unused")

imputed_data <- right_join(
  select(final_data, -c("LDMC (g/g)":Woodiness)),
  select(imputed_traits, -c(Genus:Ssp_var))
)

write_csv(imputed_data, "data/clean/final_dataset_IMPUTED.csv")

# Make a tidy dataframe of all traits
full_df <- traits_df %>% 
  # keep track of which traits were original values
  mutate(type = "original") %>%
  rbind(imputed_traits %>% mutate(type = "imputed")) %>%
  select(-c("Genus", "Species", "Ssp_var", "Family")) %>%
  melt(id = c("Species_full", "type")) %>% 
  mutate(value = as.numeric(value))

# Collect the error OOB error for each variable

var_df <- full_df %>% 
  filter(type == "imputed") %>% 
  select(-c("Species_full", "type")) %>% 
  # filter(!variable %in% c("Woodiness", "Growth Form")) %>% 
  group_by(variable) %>% 
  summarise(variance = case_when(
    (variable %in% c("Woodiness", "Growth Form")) ~ 0,
    (!variable %in% c("Woodiness", "Growth Form")) ~ var(value)
    )) %>% 
  unique()

error_df <- as.list(PNW_imp$OOBerror[1:5]) %>% 
  as_tibble(.name_repair = "universal") %>% 
  melt() %>% 
  rename(error_value = value) %>% 
  mutate(error_type = substr(as.character(variable), start = 1, stop = 3),
         .keep = "unused") %>% 
  cbind(variable = colnames(PNW_imp$ximp[,1:5])) %>% 
  # compute NRMSE from MSE
  inner_join(var_df, by = "variable") %>% 
  mutate(RMSE = sqrt(error_value)) %>% 
  mutate(NRMSE = RMSE / variance)

# 5) Perform diagnostics on imputed data ----

# * Calculate distributions of trait values for the original and imputed data sets ----
# Continuous variables
traitdist_compare_df <- traits_df %>%
  select(-c("Woodiness")) %>%
  mutate(type = "original") %>%
  rbind(imputed_traits %>%
          mutate(type = "imputed") %>%
          select(-c("Woodiness"))
  ) %>%
  select(-c("Genus", "Species", "Ssp_var", "Family")) %>%
  melt(id = c("Species_full", "type"))

# Discrete variables
traitdist_compare_df_discrete<- traits_df %>%
  select(c("Species_full", "Woodiness")) %>%
  mutate(type = "original") %>%
  melt(id = c("Species_full", "type")) %>% 
  group_by(type, variable, value) %>% 
  filter(!is.na(value)) %>% 
  summarise(count = n()) %>% 
  mutate(density = count / sum(count)) %>% 
  rbind(imputed_traits %>%
          select(c("Species_full", "Woodiness")) %>%
          mutate(type = "imputed") %>% 
          melt(id = c("Species_full", "type"))%>% 
          filter(!is.na(value)) %>% 
          group_by(type, variable, value) %>% 
          summarise(count = sum(!is.na(value))) %>% 
          mutate(density = count / sum(count)) 
  )

# 6) Illustrate diagnostics to ensure imputation was appropriate ----

# * Figure: histograms of original trait data vs. imputed trait data ----
# Put together annotation data
annotate_df <- inner_join(error_df, coverage_vals)

# Continuous variables
traitdist_compare_plot <- traitdist_compare_df %>%
  # select(contains(select_trait), type) %>%
  ggplot(aes()) +
  # plot frequency polynomial
  geom_freqpoly(aes(x = value, y = after_stat(density), color = type), linewidth = 2) +
  # # uncomment this to show density curves instead
  # geom_density(linewidth = 2) +
  # annotate plot with coverage percentages of traits
  geom_text(data = filter(annotate_df, 
                          !(variable %in% c("Woodiness", "Growth Form"))),
            aes(x = Inf, y = Inf, 
                label = paste0("coverage  = ", 100*round(coverage, 2), "%\n",
                               "error (NRMSE) = ", signif(NRMSE, 3) )),
            hjust="right", vjust="top") +
  # iterate over all traits
  facet_wrap( ~ variable,
              ncol = 3,
              scales = "free") +
  ggtitle("Distributions of continuous traits") +
  theme_cowplot(16) +
  theme(axis.title.x = element_blank())


# Discrete  variables
traitdist_compare_plot_discrete <- traitdist_compare_df_discrete %>%
  group_by(type, variable) %>% 
  ggplot() +
  geom_col(aes(x = value, y = density, color = type, group = type),
           fill = NA, position = "dodge", lwd = 1) +
  # annotate plot with coverage percentages of traits
  geom_text(data = filter(annotate_df, 
                          variable %in% c("Woodiness", "Growth Form")),
            aes(x = Inf, y = Inf, 
                label = paste0("coverage  = ", 100*round(coverage, 2), "%\n",
                               "error (", error_type, ") = ", signif(error_value,3) )),
            hjust="right", vjust="top") +
  # iterate over all traits
  facet_wrap( ~ variable, ncol = 3, scales = "free") +
  ggtitle("Distributions of discrete traits") +
  theme_cowplot(16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank()
  )

# Save figures
ggsave2(
  filename = "figures/ImputedDistributions_Continuous.png",
  plot = traitdist_compare_plot,
  width = 16,
  height = 9
)

ggsave2(
  filename = "figures/ImputedDistributions_Discrete.png",
  plot = traitdist_compare_plot_discrete,
  width = 16,
  height = 9
)

# 7) Output imputed data frame ----
data_out <-imputed_traits %>%
  select(-c("Genus", "Species", "Ssp_var")) %>%
  mutate(LeafArea = exp(LeafArea_log), .keep = "unused") %>% 
  mutate(PlantHeight = exp(PlantHeight_log), .keep = "unused") %>% 
  mutate(DiasporeMass = exp(DiasporeMass_log), .keep = "unused") %>% 
  melt(id = "Species_full")
# !!! KD: this isn't in the exact same form as the other data sets. I can change this later

write_csv(data_out,
          file = "data/clean/imputed_traits.csv")