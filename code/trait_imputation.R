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
data_in <- read_excel("data/clean/Trait_data_TRY_Diaz_2022_PNW.xlsx")

traits_df <- data_in %>%
  # # put in "tidy" form with traits as row entries
  # melt(
  #   id = "Species_name", variable.name = "trait_name",
  #   value.name = "trait_value"
  # ) %>%
  # expand binomial
  mutate(Genus = stringr::word(Species_name, 1, 1, sep = " ")) %>%
  mutate(Species = stringr::word(Species_name, 2, 2, sep = " ")) %>%
  mutate(Ssp_var = stringr::word(Species_name, 3, 4, sep = " ")) %>%
  relocate(Genus, Species, Ssp_var)

# 3) Set up parameters for the imputation function ----


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
  # select("Leaf area (mm2)":"SSD observed (mg/mm3)") 
  select("Woodiness":"SSD observed (mg/mm3)") %>%
  mutate(Woodiness = as.factor(Woodiness)) %>%
  mutate(`Growth Form` = as.factor(`Growth Form`))
  # select("Leaf area (mm2)":"SSD observed (mg/mm3)")

# combine dummy vars with traits 
traits_to_impute<-cbind(traits_to_impute,taxa)

#
set.seed(82) # for debugging

# run missForest imputation
PNW_imp <- missForest(traits_to_impute,
  maxiter = 100, # maximum number of iterations to be performed given the stopping criterion isn't met
  ntree = 1000, # number of trees to grow in each forest
  verbose = TRUE, # if 'TRUE', gives additional output between iterations
  variablewise = TRUE # if 'TRUE', the OOB error is returned for each variable separately
)

# add back taxonomic columns back in
imputed_traits <- cbind(traits_df[,1:4],PNW_imp$ximp[,1:9])
# imputed_traits <- cbind(traits_df[,1:4],PNW_imp$ximp)


# 5) Perform diagnostics on imputed data ----

# Figure: histograms of original trait data vs. imputed trait data

traitdist_compare_df <- traits_df %>% 
  select(-c("Woodiness", "Growth Form","SSD imputed (mg/mm3)")) %>% 
  mutate(type = "original") %>% 
  rbind(imputed_traits %>% 
           mutate(type = "imputed") %>% 
           select(-c("Woodiness", "Growth Form"))) %>% 
  select(-c("Genus", "Species", "Ssp_var")) %>% 
  melt(id = c("Species_name", "type"))

traitdist_compare_plot <- traitdist_compare_df %>% 
  # select(contains(select_trait), type) %>% 
  ggplot(aes(value, color = type)) +
  # geom_histogram() +
  geom_density(linewidth = 2) +
  # iterate over all traits
  facet_wrap( ~ variable, ncol = 3, scales = "free") +
  ggtitle("Distributions of original and imputed traits (no taxonomic variables)") +
  theme_cowplot(16)

ggsave2(filename = "figures/WTaxa_ImputeCompare.png",
        width = 16,
        height = 9)

# 6) Illustrate diagnostics to ensure imputation was appropriate ----


# 7) Output imputed data frame ----
