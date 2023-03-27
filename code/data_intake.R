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


# 0) Set-up, load in necessary packages and data sets ---------------------
# Packages
library(tidyverse)
library(reshape2)
library(readxl)
library(taxize)

## Helper functions
# Function to determine if a given data column is empty
dataCheck_func <- function(x) {
  sum(!is.na(x)) > 0
}

# # Create synonym data frame
# synonym_list <- select(full_data, Species_full) %>%
#   mutate(Synonym = case_when(
#     Species_full == "Chamaecyparis nootkatensis" ~ "Cupressus nootkatensis",
#     Species_full == "Equisetum hyemale ssp. affine" ~ "Equisetum hyemale",
#     Species_full == "Adiantum aleuticum" ~ "Adiantum pedatum",
#     Species_full == "Athyrium filix-femina ssp. cyclosorum" ~ "Athyrium filix-femina",
#     Species_full == "Polypodium glycyrrhiza" ~ "Polypodium vulgare",
#     Species_full == "Allium schoenoprasum var. sibiricum" ~ "Allium schoenoprasum",
#     Species_full == "Alnus viridis ssp. crispa" ~ "Alnus alnobetula ssp. crispa",
#     Species_full == "Alnus viridis ssp. sinuata" ~ "Alnus alnobetula ssp. sinuata",
#     Species_full == "Arctous ruber" ~ "Arctous alpina",
#     Species_full == "Argentina egedii" ~ "Potentilla anserina ssp. groenlandica",
#     Species_full == "Argentina anserina" ~ "Potentilla anserina",
#     Species_full == "Betula pumila var. glandulifera" ~ "Betula pumila",
#     Species_full == "Cornus unalaschkensis" ~ "Cornus canadensis",
#     Species_full == "Dodecatheon pauciflorum" ~ "Dodecatheon meadia", # Dodecatheon pauciflorum not in full data set
#     Species_full == "Eurybia conspicua" ~ "Aster conspicuus",
#     Species_full == "Glaux maritima" ~ "Lysimachia maritima",
#     Species_full == "Lappula occidentalis" ~ "Lappula redowskii",
#     Species_full == "Hierochloe hirta" ~ "Hierochloe odorata", # Hierochloe hirta not in full data set
#     Species_full == "Mahonia nervosa" ~ "Berberis nervosa",
#     Species_full == "Maianthemum racemosum ssp. amplexicaule" ~ "Maianthemum racemosum",
#     Species_full == "Phyllospadix scouleri" ~ "Phyllospadix iwatensis",
#     Species_full == "Platanthera stricta" ~ "Platanthera dilatata",
#     Species_full == "Populus balsamifera ssp. balsamifera" ~ "Populus balsamifera",
#     Species_full == "Populus balsamifera ssp. trichocarpa" ~ "Populus trichocarpa",
#     Species_full == "Pseudoroegneria spicata" ~ "Elymus spicatus",
#     Species_full == "Rhodiola rosea" ~ "Rhodiola rosea var. rosea",
#     Species_full == "Rhododendron groenlandicum" ~ "Ledum palustre ssp. groenlandicum", # Rhododendron groenlandicum not in full data set
#     Species_full == "Rhododendron neoglandulosum" ~ "Ledum glandulosum", # Rhododendron neoglandulosum not in full data set
#     Species_full == "Toxicodendron rydbergii" ~ "Toxicodendron radicans",
#     Species_full == "Vicia nigricans ssp. gigantea" ~ "Vicia nigricans",
#     Species_full == "Zigadenus venenosus" ~ "Toxicoscordion venenosum",
#     Species_full == "Zigadenus elegans" ~ "Anticlea elegans",
#     TRUE ~ Species_full
#   ))

# Save synonym list for reference later
# write_csv(synonym_list,"data/clean/synonym_list.csv")
synonym_list <- read_csv("data/clean/synonym_list.csv")

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

TEK_data <- read_csv("data/raw/TEKdata.csv") %>%
  # Remove any empty data columns
  select_if(dataCheck_func) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove NA rows
  filter(!is.na(Species))
# species in base_data but not in TEK_dataL Hierochloe, Ledum groenlandicum, and Ledum glandulosum

full_data <- right_join(base_data, TEK_data) %>%
  # Create a column with the full species scientific species name
  # i.e. indicating whether a record is a ssp. or a var. or not
  dplyr::mutate(Var = str_replace_na(Var, replacement = "")) %>%
  dplyr::mutate(Ssp = str_replace_na(Ssp, replacement = "")) %>%
  dplyr::mutate(Ssp_author = str_replace_na(Ssp_author, replacement = "")) %>%
  dplyr::mutate(Species_full = paste(Species, Var, Ssp)) %>%
  dplyr::mutate(Species_full_author = paste(Species, Var, author, Ssp, Ssp_author)) %>%
  # replace double white spaces with one space
  dplyr::mutate(Species_full = str_replace_all(Species_full, "  ", " ")) %>%
  # remove white spaces at the end of a string
  dplyr::mutate(Species_full = trimws(Species_full, which = c("right"))) %>%
  # replace double white spaces with one space
  dplyr::mutate(Species_full_author = str_replace_all(Species_full_author, "  ", " ")) %>%
  # remove white spaces at the end of a string
  dplyr::mutate(Species_full_author = trimws(Species_full_author, which = c("right"))) %>%
  # deal with synonyms
  right_join(synonym_list) %>%
  # keep track of old name
  mutate(old_name = Species_full) %>%
  # rename to main synonym
  mutate(Species_full = Synonym) %>%
  # set binomial to synonym binomial
  mutate(Species = stringr::word(Synonym, 1, 2)) %>%
  group_by(Species_full) %>%
  # combine number of uses for synonymous species
  mutate(`Number of uses` = case_when(
    n() > 1 ~ sum(`Number of uses`)
  )) %>%
  # remove any repeated entries
  unique()

# 2) Load in and process Diaz data set ------------------------------------

# Load in trait data compiled by Diáz et al. 2020 (based on TRY)
# https://www.nature.com/articles/s41597-022-01774-9
Diaz_data <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx",
                        sheet = 1
)

###* Resolve Diaz species names using Global Names Resolver ----
# # Partition the list of species (otherwise it's too large to pass to gnr_resolve())
# Diaz_species_list <- Diaz_data$`Species name standardized against TPL`
# species_num <- length(Diaz_species_list)
# num_seq <- 1+0:(species_num %% 1000) * 1000
#
# # Select name sources
# src <- c("EOL", "The International Plant Names Index", #these are examples, so choose the database of interest
#          "Index Fungorum", "ITIS", "Catalogue of Life",
#          "Tropicos - Missouri Botanical Garden")
# # Initialize list of synonyms
# Diaz_synonyms <- tibble(name_in = as.character(),
#                         name_out = as.character())
# # Collect synonyms
# for (ii in 1:length(num_seq)) {
#   temp_seq <- num_seq[ii]:min(num_seq[ii+1], length(Diaz_species_list))
#   temp_species_list <- Diaz_species_list[temp_seq]
#   temp_species_out <- gnr_resolve(temp_species_list,
#                                   data_source_ids = c(1,2,5,150, 165,167),
#                                   with_canonical_ranks=T,
#                                   best_match_only = TRUE) %>%
#     select(name_in = user_supplied_name, name_out = matched_name2)
#   Diaz_synonyms <- rbind(Diaz_synonyms,temp_species_out)
# }
# # Save list of synonyms
# write_csv(Diaz_synonyms, 'data/clean/Diaz_synonyms.csv')

# Load list of synonyms
Diaz_synonyms <- read_csv("data/clean/Diaz_synonyms.csv", show_col_types = FALSE)

Diaz_data_renamed <- Diaz_synonyms %>%
  mutate(`Species name standardized against TPL` = name_in, .keep = "unused") %>%
  right_join(Diaz_data) %>%
  mutate(Species_full = name_out, .keep = "unused") %>%
  mutate(Species = stringr::word(Species_full, 1,2)) %>%
  select(-`Species name standardized against TPL`)

# Get list of species missing from Diaz data
Diaz_species_list <- Diaz_data_renamed$Species_full
Diaz_species_short <- Diaz_data_renamed$Species
full_species_list <- full_data$Species_full
full_species_short <- full_data$Species

# 23 missing if using full species epithet
missing_from_Diaz <- full_species_list[which(!(full_species_list %in% c(Diaz_species_list,Diaz_species_short)))]

# 39 missing is using shortened species epithet
missing_from_Diaz_short <- full_species_short[which(!(full_species_short %in% c(Diaz_species_list,Diaz_species_short)))]

# Reduce Diaz data to species present in the main dataset
Diaz_data_sel <- Diaz_data_renamed %>%
  # Remove non-focal species
  filter(Species %in% c(full_data$Species, full_data$Species_full)) # %>% # including regular species names only leads to the additional inclusion of Alnus incana
# # reduce to relevant set of traits
# select("TRY 30 AccSpecies ID":"Family",
#        "Plant height (m)", "Nmass (mg/g)", "LDMC (g/g)", "Leaf area (mm2)", "Woodiness", )
# missing: rooting depth, stem density, vegetative reproduction

# Join Diaz data and full data
new_data <- full_join(full_data, Diaz_data_sel, by = c("Species_full"))

# 3) Load in and process TRY data sets ------------------------------------

# We include these to obtain trait data for the 23 (or 39) missing from Diaz
# load TRY data
# TRY data is separated by trait (these data is in our Google Drive folder in databases)
# Note that csv files starting with "TRY" followed by a trait name are a subset of the data provided for all TRY species present in their database. Only species present in "Trait_data_TRY_Diaz_2022_PNW.csv" were kept, in addition to possible synonmys found in funct_div.R
# the full TRY data with all species is called TRY.zip

TRY_plant_height <- read_csv("data/raw/TRY_data_Feb2023/TRY_plant_height_growth.csv",
                             show_col_types = FALSE) %>% dplyr::select(-...1)
plant_height_trait_names <- unique(TRY_plant_height$TraitName)

TRY_leafN <- read_csv("data/raw/TRY_data_Feb2023/TRY_leafN.csv",
                      show_col_types = FALSE) %>% dplyr::select(-...1)
leafN_trait_names <- unique(TRY_leafN$TraitName)

TRY_LDMC <- read_csv("data/raw/TRY_data_Feb2023/TRY_LDMC.csv",
                     show_col_types = FALSE) %>% dplyr::select(-...1)
LDMC_trait_names <- unique(TRY_LDMC$TraitName)

# wait to handle this data because we need to know what Diaz used for leaf area measurements
TRY_leafarea <- read_csv("data/raw/TRY_data_Feb2023/TRY_leafarea.csv",
                         show_col_types = FALSE) %>% dplyr::select(-...1)
leafarea_trait_names <- unique(TRY_leafarea$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_rooting_depth <- read_csv("data/raw/TRY_data_Feb2023/TRY_rooting_depth.csv",
                              show_col_types = FALSE) %>% dplyr::select(-...1)
rooting_depth_trait_names <- unique(TRY_rooting_depth$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_stem_density <- read_csv("data/raw/TRY_data_Feb2023/TRY_stem_density.csv",
                             show_col_types = FALSE) %>% dplyr::select(-...1)
stem_density_trait_names <- unique(TRY_stem_density$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_plant_woodiness <- read_csv("data/raw/TRY_data_Feb2023/TRY_plant_woodiness.csv",
                                show_col_types = FALSE) %>% dplyr::select(-...1)
plant_woodiness_trait_names <- unique(TRY_plant_woodiness$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_plant_veg_reproduction <- read_csv("data/raw/TRY_data_Feb2023/TRY_vegetative_reproduction.csv",
                                       show_col_types = FALSE) %>% dplyr::select(-...1)
plant_veg_reproduction_trait_names <- unique(TRY_plant_veg_reproduction$TraitName)

TRY_data <- rbind(TRY_plant_height, TRY_leafN, TRY_LDMC, TRY_leafarea,
                  TRY_rooting_depth, TRY_stem_density, TRY_plant_woodiness,
                  TRY_plant_veg_reproduction) %>%
  # subset to just the species we're interested in
  filter(AccSpeciesName %in% c(full_data$Species, full_data$Species_full))

# Get list of species which are in base data but not TRY data set
TRY_species <- unique(TRY_data$AccSpeciesName)
full_species <- unique(c(full_data$Species_full,full_data$Species))
missing_TRY <- full_species[which(!(full_species %in% TRY_species))]

# Subset TRY_data in two ways:
# 1) traits not included in Diaz
TRY_trait_list <- c()
Diaz_trait_list <- c()
new_trait_list <- TRY_trait_list[which(!(TRY_trait_list %in% Diaz_trait_list))]
new_traits_data <- TRY_data %>%
  filter(TraitName %in% new_trait_list)

# 2) species not included in Diaz
new_species <- TRY_species[which(!(TRY_species %in% Diaz_species))]

TRY_newtraits <- TRY_data %>%
  dplyr::select()

test_df <- full_join(full_data, TRY_data, by = "Species")




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
  melt(
    id = c("Species_full", "Family", "Status"),
    variable.name = "Trait"
  ) %>%
  # entityID gives a unique number for each SPECIES
  mutate(entityID = as.numeric(factor(Species_full))) %>%
  # Add genus and species labels for phylogenetic analysis
  mutate(Species = stringr::word(Species_full, 1, 2, sep = " ")) %>% # get rid of spp. and var.
  mutate(Genus = stringr::word(Species_full, 1, 1, sep = " "))
