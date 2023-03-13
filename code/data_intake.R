################################################################################
# Processing data from databases
################################################################################

## Project: Identifying Climate change vulnerable biodiversity hotspots
##
## Purpose: Process data to be used to inform simulations
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Load in main species list
##           2) Load in and process Diaz data set
##           3) Load in and process TRY data sets
##           4) Combine datasets
##           5) Put dataset into workable form for imputation and phylogeny
##              steps
##
##
## Inputs:  data - found in the data/raw folder
##          Trait_data_TRY_Diaz_2022\Dataset\Species_mean_traits.xlsx
##          - trait data from Diaz 2022 (originally from the TRY database)
##          TRY_data_Feb2023 folder
##          - trait data from a call to the TRY database
##          
##
## Outputs: data/clean/main_dataset.csv
##
## Written and by: Karen Castillioni, Elisa van Cleemput, Kyle Dahlin
## Maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023
## _____________________________________________________________________________


# 0) Set-up, load in necessary packages and data-sets ---------------------
# Packages
library(tidyverse)
library(reshape2)
library(readxl)

## Helper functions
# Function to determine if a given data column is empty
dataCheck_func <- function(x) {
  sum(!is.na(x)) > 0
}

# 1) Load in main species list --------------------------------------------

clade_labels <- c(
  "Gymnosperms", "Ferns and fern-allies",
  "Flowering plants (angiosperms)",
  "Species extracted from", # !!! just being lazy with these two. could remove these rows by hand instead
  "There are also algae, lichens, fungi and bryophytes listed in the document. We could decide to omit these"
)

# Read in data set
base_data <- read_csv("data/raw/PNW_Species_w_Metadata.csv") %>%
  # Remove anu empty data columns
  select_if(dataCheck_func) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove NA rows
  filter(!is.na(Species))

TEK_data <- read_csv("data/raw/TEKdata_20230222.csv") %>%
  # Remove any empty data columns
  select_if(dataCheck_func) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove NA rows
  filter(!is.na(Species))

full_data <- right_join(base_data, TEK_data) %>%
  # Create a column with the full species scientific species name
  # i.e. indicating whether a record is a ssp. or a var. or not
  dplyr::mutate(Var = str_replace_na(Var, replacement="")) %>%
  dplyr::mutate(Ssp = str_replace_na(Ssp, replacement="")) %>%
  dplyr::mutate(Ssp_author = str_replace_na(Ssp_author, replacement="")) %>%
  dplyr::mutate(Species_full = paste(Species, Var, Ssp)) %>%
  dplyr::mutate(Species_full_author = paste(Species, Var, author,  Ssp, Ssp_author)) %>%
  # replace double white spaces with one space
  dplyr::mutate(Species_full = str_replace_all(Species_full,"  "," ")) %>% 
  # remove white spaces at the end of a string
  dplyr::mutate(Species_full = trimws(Species_full, which=c("right"))) %>% 
  # replace double white spaces with one space
  dplyr::mutate(Species_full_author = str_replace_all(Species_full_author,"  "," ")) %>% 
  # remove white spaces at the end of a string
  dplyr::mutate(Species_full_author = trimws(Species_full_author, which=c("right"))) 
  
# Create synonym data frame
synonym_list <- select(full_data, Species_full) %>% 
  mutate(temp_name = case_when(
    Species_full == "Chamaecyparis nootkatensis" ~ "Cupressus nootkatensis",
    Species_full == "Equisetum hyemale ssp. affine" ~ "Equisetum hyemale",
    Species_full == "Adiantum aleuticum" ~ "Adiantum pedatum",
    Species_full == "Athyrium filix-femina ssp. cyclosorum" ~ "Athyrium filix-femina",
    Species_full == "Polypodium glycyrrhiza" ~ "Polypodium vulgare",
    Species_full == "Allium schoenoprasum var. sibiricum" ~ "Allium schoenoprasum",
    Species_full == "Alnus viridis ssp. crispa" ~ "Alnus alnobetula ssp. crispa",
    Species_full == "Alnus viridis ssp. sinuata" ~ "Alnus alnobetula ssp. sinuata",
    Species_full == "Arctous ruber" ~ "Arctous alpina",
    Species_full == "Argentina egedii" ~ "Potentilla anserina ssp. groenlandica",
    Species_full == "Argentina anserina" ~ "Potentilla anserina",
    Species_full == "Betula pumila var. glandulifera" ~ "Betula pumila",
    Species_full == "Cornus unalaschkensis" ~ "Cornus canadensis",
    Species_full == "Dodecatheon pauciflorum" ~ "Dodecatheon meadia", # Dodecatheon pauciflorum not in full data set
    Species_full == "Eurybia conspicua" ~ "Aster conspicuus", 
    Species_full == "Glaux maritima" ~ "Lysimachia maritima", 
    Species_full == "Lappula occidentalis" ~ "Lappula redowskii", 
    Species_full == "Hierochloe hirta" ~ "Hierochloe odorata", # Hierochloe hirta not in full data set
    Species_full == "Mahonia nervosa" ~ "Berberis nervosa", 
    Species_full == "Maianthemum racemosum ssp. amplexicaule" ~ "Maianthemum racemosum", 
    Species_full == "Phyllospadix scouleri" ~ "Phyllospadix iwatensis", 
    Species_full == "Platanthera stricta" ~ "Platanthera dilatata", 
    Species_full == "Populus balsamifera ssp. balsamifera" ~ "Populus balsamifera", 
    Species_full == "Populus balsamifera ssp. trichocarpa" ~ "Populus trichocarpa", 
    Species_full == "Pseudoroegneria spicata" ~ "Elymus spicatus", 
    Species_full == "Rhodiola rosea" ~ "Rhodiola rosea var. rosea", 
    Species_full == "Rhododendron groenlandicum" ~ "Ledum palustre ssp. groenlandicum", # Rhododendron groenlandicum not in full data set
    Species_full == "Rhododendron neoglandulosum" ~ "Ledum glandulosum", # Rhododendron neoglandulosum not in full data set
    Species_full == "Toxicodendron rydbergii" ~ "Toxicodendron radicans", 
    Species_full == "Vicia nigricans ssp. gigantea" ~ "Vicia nigricans", 
    Species_full == "Zigadenus venenosus" ~ "Toxicoscordion venenosum", 
    Species_full == "Zigadenus elegans" ~ "Anticlea elegans",
    TRUE ~ Species_full
    )) %>% 
  mutate(Synonym = ifelse(temp_name != Species_full, Species_full, NA)) %>% 
  mutate(Species_full = if_else(temp_name == Species_full, Species_full, temp_name)) %>% 
  select(-temp_name)

# Save synonym list for reference later  
write_csv(synonym_list,"data/clean/synonym_list.csv")


# 2) Load in and process Diaz data set ------------------------------------

# Load in trait data compiled by Diáz et al. 2020 (based on TRY)
# https://www.nature.com/articles/s41597-022-01774-9 
Diaz_data <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx", 
                         sheet = 1)
# Diaz_data_meta <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx", 
#                               sheet = 2)
# Diaz_data_ref <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/References.xlsx", 
#                              sheet = "References")
# Diaz_data_ref_meta <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/References.xlsx", 
#                                   sheet = 2)

# 3) Load in and process TRY data sets ------------------------------------



# 4) Combine data sets -----------------------------------------------------



# 5) Put data set into workable form for imputation and phylogeny steps --------

# Reduce dataset to only include data used in simulations
data_in <- full_data %>%
  # Only keep columns used in analysis
  select(
    # `Species`,
    `Species_full`, # EVC: use this instead of Species
    `Family`,
    `IUCN Status`,
    `Number of unique names (Indigenous)`,
    `Number of unique  languages (from Appendix 2B)`,
    # `Number of names (common)`,
    `Number of uses`
  ) %>%
  rename("Status" = `IUCN Status`) %>% 
  melt(id = c("Species_full", "Family", "Status"),
       variable.name = "Trait") %>%
  # entityID gives a unique number for each SPECIES
  mutate(entityID = as.numeric(factor(Species_full))) %>% 
  # Add genus and species labels for phylogenetic analysis
  mutate(Species = stringr::word(Species_full, 1,2, sep=" ")) %>% # get rid of spp. and var.
  mutate(Genus = stringr::word(Species_full, 1,1, sep=" "))



