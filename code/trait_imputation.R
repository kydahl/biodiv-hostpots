################################################################################
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
## Inputs:  data -
##
##          code -
##
##
## Outputs: data -
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
data_in <- read_excel("data/clean/Trait_data_TRY_Diaz_2022_PNW.xlsx") %>%
  mutate(Genus = stringr::word(Species_name, 1, 1, sep = " ")) %>%
  mutate(Species = stringr::word(Species_name, 2, 2, sep = " ")) %>%
  mutate(Ssp_var = stringr::word(Species_name, 3, 4, sep = " ")) %>%
  relocate(Genus, Species, Ssp_var)

# Modify original data
traits_df <- data_in %>%
  # # put in "tidy" form with traits as row entries
  # melt(
  #   id = "Species_name", variable.name = "trait_name",
  #   value.name = "trait_value"
  # ) %>%
  # expand binomial
  mutate(LeafArea_log = log(`Leaf area (mm2)`)) %>%
  mutate(PlantHeight_log = log(`Plant height (m)`)) %>%
  mutate(DiasporeMass_log = log(`Diaspore mass (mg)`)) %>%
  # Remove original variables replaced by "logged" versions
  select(-c("Leaf area (mm2)", "Plant height (m)", "Diaspore mass (mg)")) %>%
  # Remove traits that we won't be using in our study
  select(-c("LDMC (g/g)", "SSD observed (mg/mm3)", "SSD imputed (mg/mm3)"))

# 3) Set up parameters for the imputation function (missForest) ----


# 4) Run missForest algorithm to impute missing traits ----
# get taxonomic columns
taxa <- select(traits_df, Genus, Species, Ssp_var) # Get order and family level data
test <- rep(1, nrow(taxa))
taxa <- cbind(taxa, test)

# turn taxonomic groupings into dummy binary variables
dummies <- dummyVars(test ~ ., data = taxa)
taxa <- predict(dummies, newdata = taxa)

# select traits to be imputed
traits_to_impute <- as.data.frame(traits_df) %>%
  select(-c("Genus":"Species_name")) %>%
  mutate(Woodiness = as.factor(Woodiness)) %>%
  mutate(`Growth Form` = as.factor(`Growth Form`))

# combine dummy vars with traits (comment this line out to not use taxa in imputation)
traits_to_impute <- cbind(traits_to_impute, taxa)

#
# set.seed(82) # for debugging

# run missForest imputation
PNW_imp <- missForest(traits_to_impute,
                      maxiter = 100, # maximum number of iterations to be performed given the stopping criterion isn't met
                      ntree = 1000, # number of trees to grow in each forest
                      verbose = TRUE, # if 'TRUE', gives additional output between iterations
                      variablewise = TRUE # if 'TRUE', the OOB error is returned for each variable separately
)

# add actual taxonomic columns back in and remove dummy variables
imputed_traits <- cbind(traits_df[, 1:4], PNW_imp$ximp[, 1:7])
# imputed_traits <- cbind(traits_df[,1:4],PNW_imp$ximp)


# 5) Perform diagnostics on imputed data ----

# Figure: histograms of original trait data vs. imputed trait data

traitdist_compare_df <- traits_df %>%
  select(-c("Woodiness", "Growth Form")) %>%
  mutate(type = "original") %>%
  rbind(imputed_traits %>%
          mutate(type = "imputed") %>%
          select(-c("Woodiness", "Growth Form"))
  ) %>%
  select(-c("Genus", "Species", "Ssp_var")) %>%
  melt(id = c("Species_name", "type"))

traitdist_compare_df_discrete<- traits_df %>%
  select(c("Genus":"Growth Form")) %>%
  mutate(type = "original") %>%
  select(-c("Genus", "Species", "Ssp_var")) %>%
  melt(id = c("Species_name", "type")) %>% 
  group_by(type, variable, value) %>% 
  filter(!is.na(value)) %>% 
  summarise(count = n()) %>% 
  mutate(density = count / sum(count)) %>% 
  rbind(imputed_traits %>%
          select(c("Genus":"Growth Form")) %>% 
          mutate(type = "imputed") %>% 
          select(-c("Genus", "Species", "Ssp_var")) %>%
          melt(id = c("Species_name", "type"))%>% 
          filter(!is.na(value)) %>% 
          group_by(type, variable, value) %>% 
          summarise(count = sum(!is.na(value))) %>% 
          mutate(density = count / sum(count)) 
  )

traitdist_compare_plot <- traitdist_compare_df %>%
  # select(contains(select_trait), type) %>%
  ggplot(aes(value, color = type)) +
  geom_freqpoly(linewidth = 2) +
  # geom_density(linewidth = 2) +
  # iterate over all traits
  facet_wrap( ~ variable, ncol = 3, scales = "free") +
  ggtitle("Distributions of continuous traits") +
  theme_cowplot(16)

traitdist_compare_plot_discrete <- traitdist_compare_df_discrete %>%
  group_by(type, variable) %>% 
  # select(contains(select_trait), type) %>%
  ggplot(aes(x = value, y = density, color = type, group = type)) +
  # geom_point() +
  # geom_freqpoly(linewidth = 2) +
  geom_col(fill = NA, position = "dodge") +
  # geom_density(linewidth = 2) +
  # iterate over all traits
  facet_wrap( ~ variable, ncol = 3, scales = "free") +
  ggtitle("Distributions of discrete traits") +
  theme_cowplot(16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

# ggsave2(filename = "figures/WTaxa_ImputeCompare.png",
#         width = 16,
#         height = 9)
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

# 6) Illustrate diagnostics to ensure imputation was appropriate ----


# 7) Output imputed data frame ----
