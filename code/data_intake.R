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
 
# Function: replace species names with those from the GNR resolver
synonym_func <- function(in_df) {
  out_df <- in_df %>% 
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
    dplyr::mutate(Species_full_author = trimws(Species_full_author, which = c("right"))) 
  
  # Partition the list of species (otherwise it's too large to pass to gnr_resolve())
  species_list <- out_df$Species_full
  species_num <- length(species_list)
  num_seq <- 1+0:(species_num %/% 1000) * 1000

  # Select name sources
  src <- c("EOL", "The International Plant Names Index", #these are examples, so choose the database of interest
           "Index Fungorum", "ITIS", "Catalogue of Life",
           "Tropicos - Missouri Botanical Garden")
  # Initialize list of synonyms
  synonyms_list <- tibble(name_in = as.character(),
                          name_out = as.character())
  # Collect synonyms
  for (ii in num_seq) {
    temp_seq <- ii+0:999
    temp_species_list <- species_list[temp_seq]
    temp_species_out <- gnr_resolve(temp_species_list,
                                    data_source_ids = c(1,2,5,150, 165,167),
                                    with_canonical_ranks=T,
                                    best_match_only = TRUE) %>%
      select(name_in = user_supplied_name, name_out = matched_name2)
    synonyms_list <- rbind(synonyms_list,temp_species_out)
  }
  
  # Add our custom curated listed of synonyms
  # Load synonym list created by Elisa Van Cleemput
  custom_synonyms_list <- read_csv("data/clean/synonym_list.csv", show_col_types = FALSE)
  
  # Check if there are disagreements between our custom list and GNR
  GNR_synonyms <- rename(synonyms_list, original_name = name_in, Synonym_GNR = name_out)
  compare_synonyms_df <- left_join(custom_synonyms_list, GNR_synonyms)
  syn_disagree_df <- filter(compare_synonyms_df, Synonym != Synonym_GNR)
  cat(paste0("There are ", dim(syn_disagree_df)[1], " synonym disagreements. \n
               Check the file data/clean/synonym_disagreements.csv for details."))
  write_csv(syn_disagree_df, "data/clean/synonym_disagreements.csv")
  
  out_df2 <- out_df %>% 
    # add column of synonym names
    left_join(rename(synonyms_list, Species_full = name_in, Synonym = name_out),
              by = "Species_full") %>%
    # keep track of original name
    mutate(original_name = Species_full) %>%
    # rename species to their main synonym
    mutate(Species_full = Synonym) %>%
    # set binomial to synonym binomial
    mutate(Species = stringr::word(Synonym, 1, 2))
  
  # add column of synonym names
  out_df2 <- out_df %>% 
    left_join(rename(synonym_list, Species_full = original_name),
            by = "Species_full") %>%
    # keep track of original name
    mutate(original_name = Species_full) %>%
    # rename species to their main synonym
    mutate(Species_full = Synonym) %>%
    # set binomial to synonym binomial
    mutate(Species = stringr::word(Synonym, 1, 2)) %>%
    # remove any repeated entries
    group_by(Species_full) %>% # Cornus canadensis and Platanthera dilatata (4 entries -> 2)
    distinct(Species_full, N_Names, N_Langs, N_Uses, .keep_all = TRUE)
  
}



# # # Create synonym data frame
# # # Put together by Elisa Van Cleemput
# synonym_list <- select(base_data, Species_full) %>%
#   rename(original_name = Species_full) %>% 
#   mutate(Synonym = case_when(
#     original_name == "Chamaecyparis nootkatensis" ~ "Cupressus nootkatensis",
#     original_name == "Equisetum hyemale ssp. affine" ~ "Equisetum hyemale",
#     original_name == "Adiantum aleuticum" ~ "Adiantum pedatum",
#     original_name == "Athyrium filix-femina ssp. cyclosorum" ~ "Athyrium filix-femina",
#     original_name == "Polypodium glycyrrhiza" ~ "Polypodium vulgare",
#     original_name == "Allium schoenoprasum var. sibiricum" ~ "Allium schoenoprasum",
#     original_name == "Alnus viridis ssp. crispa" ~ "Alnus alnobetula ssp. crispa",
#     original_name == "Alnus viridis ssp. sinuata" ~ "Alnus alnobetula ssp. sinuata",
#     original_name == "Arctous ruber" ~ "Arctous alpina",
#     original_name == "Argentina egedii" ~ "Potentilla anserina ssp. groenlandica",
#     original_name == "Argentina anserina" ~ "Potentilla anserina",
#     original_name == "Betula pumila var. glandulifera" ~ "Betula pumila",
#     original_name == "Cornus unalaschkensis" ~ "Cornus canadensis",
#     original_name == "Dodecatheon pauciflorum" ~ "Dodecatheon meadia", # Dodecatheon pauciflorum not in full data set
#     original_name == "Eurybia conspicua" ~ "Aster conspicuus",
#     original_name == "Glaux maritima" ~ "Lysimachia maritima",
#     original_name == "Lappula occidentalis" ~ "Lappula redowskii",
#     original_name == "Hierochloe hirta" ~ "Hierochloe odorata", # Hierochloe hirta not in full data set
#     original_name == "Mahonia nervosa" ~ "Berberis nervosa",
#     original_name == "Maianthemum racemosum ssp. amplexicaule" ~ "Maianthemum racemosum",
#     original_name == "Phyllospadix scouleri" ~ "Phyllospadix iwatensis",
#     original_name == "Platanthera stricta" ~ "Platanthera dilatata",
#     original_name == "Populus balsamifera ssp. balsamifera" ~ "Populus balsamifera",
#     original_name == "Populus balsamifera ssp. trichocarpa" ~ "Populus trichocarpa",
#     original_name == "Pseudoroegneria spicata" ~ "Elymus spicatus",
#     original_name == "Rhodiola rosea" ~ "Rhodiola rosea var. rosea",
#     original_name == "Rhododendron groenlandicum" ~ "Ledum palustre ssp. groenlandicum", # Rhododendron groenlandicum not in full data set
#     original_name == "Rhododendron neoglandulosum" ~ "Ledum glandulosum", # Rhododendron neoglandulosum not in full data set
#     original_name == "Toxicodendron rydbergii" ~ "Toxicodendron radicans",
#     original_name == "Vicia nigricans ssp. gigantea" ~ "Vicia nigricans",
#     original_name == "Zigadenus venenosus" ~ "Toxicoscordion venenosum",
#     original_name == "Zigadenus elegans" ~ "Anticlea elegans",
#     TRUE ~ original_name
#   ))
# # 
# # # Save synonym list for reference later
# write_csv(synonym_list,"data/clean/synonym_list.csv")

# Load synonym list created by Elisa Van Cleemput
synonym_list <- read_csv("data/clean/synonym_list.csv", show_col_types = FALSE)

# 1) Load in main species list --------------------------------------------
clade_labels <- c(
  "Gymnosperms", "Ferns and fern-allies",
  "Flowering plants (angiosperms)",
  "Species extracted from", # !!! just being lazy with these two. could remove these rows by hand instead
  "There are also algae, lichens, fungi and bryophytes listed in the document. We could decide to omit these"
)

# Read in data set
base_data <- read_csv("data/raw/PNW_Species_w_Metadata.csv", show_col_types = FALSE) %>%
  # Remove anu empty data columns
  select_if(dataCheck_func) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove rows with NA as species name
  filter(!is.na(Species))%>%
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
  dplyr::mutate(Species_full_author = trimws(Species_full_author, which = c("right"))) 

TEK_data <- read_csv("data/raw/TEKdata.csv", show_col_types = FALSE) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove rows with NA as species name
  filter(!is.na(Species))%>%
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
  # Rename TEK columns for easier reference
  rename(N_Names = `Number of unique names (Indigenous)`,
         N_Langs = `Number of unique  languages (from Appendix 2B)`,
         N_Uses = `Number of uses`)
  
# species in base_data but not in TEK_data: L Hierochloe, Ledum groenlandicum, and Ledum glandulosum

full_data <- full_join(base_data, TEK_data) %>%
  # select columns relevant for analysis
  dplyr::select(Species:Family, Species_full:N_Uses) %>% 
  # remove species for which we lack TEK data
  filter(!is.na(N_Names)) %>%
  # add column of synonym names
  left_join(rename(synonym_list, Species_full = original_name),
            by = "Species_full") %>%
  # keep track of original name
  mutate(original_name = Species_full) %>%
  # rename species to their main synonym
  mutate(Species_full = Synonym) %>%
  # set binomial to synonym binomial
  mutate(Species = stringr::word(Synonym, 1, 2)) %>%
  # remove any repeated entries
  group_by(Species_full) %>% # Cornus canadensis and Platanthera dilatata (4 entries -> 2)
  distinct(Species_full, N_Names, N_Langs, N_Uses, .keep_all = TRUE)

# Check for synonymous species with multiple entries
duplicate_df <- full_data %>% 
  group_by(Species_full) %>% 
  filter(n()>1)

# List which species were in our original TEK dataset that are missing from the full dataset
full_species <- full_data$Species_full
TEK_species <- TEK_data$Species_full

lost_species <- TEK_species[which(!(TEK_species %in% full_species))]

# Check that these were all renamed to synonyms
all(lost_species %in% full_data$original_name) # No
true_lost_species <- lost_species[which(!(lost_species %in% full_data$original_name))]

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
