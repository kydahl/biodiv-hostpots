################################################################################
# Data intake for plugging into model
################################################################################

# 0) Load in necessary packages, functions, settings ###########################
# Packages
require(tidyverse)
require(reshape2)

## Helper functions
# Function to determine if a given data column is empty
dataCheck_func <- function(x) {
  sum(!is.na(x)) > 0
}

# 1) Pull in raw data ##########################################################
clade_labels <- c(
  "Gymnosperms", "Ferns and fern-allies",
  "Flowering plants (angiosperms)",
  "Species extracted from", # !!! just being lazy with these two. could remove these rows by hand instead
  "There are also algae, lichens, fungi and bryophytes listed in the document. We could decide to omit these"
)

# Read in data set
base_data <- read_csv("data/PNW_SpeciesData.csv") %>%
  # Remove empty data columns
  select_if(dataCheck_func) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove NA rows
  filter(!is.na(Species))

TEK_data <- read_csv("data/TEKdata.csv") %>%
  # Remove empty data columns
  select_if(dataCheck_func) %>% # !!! removes `Number of names (common)` currently
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove NA rows
  filter(!is.na(Species))

full_data <- right_join(base_data, TEK_data) 

# 2) Clean the data for analysis ###############################################

data_in <- full_data %>%
  # Only keep columns used in analysis
  select(
    `Species`,
    # `Family`,
    `IUCN Status`,
    `Number of unique names (Indigenous)`,
    `Number of unique  languages (from Appendix 2B)`,
    # `Number of names (common)`,
    `Number of uses`
  ) %>%
  rename("Status" = `IUCN Status`) %>% 
  melt(id = c("Species", "Status"),
       variable.name = "Trait") %>%
  # entityID gives a unique number for each SPECIES
  mutate(entityID = as.numeric(factor(Species))) %>% 
  # Add genus and species labels for phylogenetic analysis
  mutate(Species = stringr::word(Species, 1,2, sep=" ")) %>% # get rid of spp. and var.
  mutate(Genus = stringr::word(Species, 1,1, sep=" "))

