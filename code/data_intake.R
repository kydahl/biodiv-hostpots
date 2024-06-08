################################################################################
# Processing data from databases
################################################################################

## Project: Comparing biodiversity hotspot identification
##
## Purpose: Process data to be used in simulations
##
## Contents: 0) Set-up, load in necessary packages, functions, and data-sets
##           1) Load in main species list
##           2) Load in and process Diaz data set
##           3) Load in and process TRY data sets
##           4) Combine datasets
##           5) Put dataset into workable form for imputation and phylogeny
##              steps
##
## Inputs:  .csv and .xlsx files found in the data/raw folder
##          
## Outputs: .csv and .rds files found in the data/clean folder
##
## Written and by: Karen Castillioni, Elisa Van Cleemput, Kyle Dahlin
## Maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023
## _____________________________________________________________________________

# 0) Set-up, load in necessary packages and data sets ---------------------
# Packages
library(tidyverse)
library(reshape2)
library(readxl)
library(taxize)
library(caret)
library(missForest)

## Helper functions
# Function: determine if a given data column is empty
dataCheck_func <- function(x) {
  sum(!is.na(x)) > 0
}

# Function: replace species names with those from the GNR resolver
synonym_func <- function(in_df) {
  # If data frame contains information at the subspecies and variety level
  # make sure the full species name includes this information
  if ("Var" %in% colnames(in_df)) {
    out_df <- in_df %>%
      # Create a column with the full species scientific species name
      # i.e. indicating whether a record is a ssp. or a var. or not
      dplyr::mutate(Var = str_replace_na(Var, replacement = "")) %>%
      dplyr::mutate(Ssp = str_replace_na(Ssp, replacement = "")) %>%
      dplyr::mutate(Ssp_author = str_replace_na(Ssp_author, replacement = "")) %>%
      # Write the full species epithet for name resolution
      dplyr::mutate(Species_full = paste(Species, Var, Ssp)) %>%
      # replace double white spaces with one space
      dplyr::mutate(Species_full = str_replace_all(Species_full, "  ", " ")) %>%
      # remove white spaces at the end of a string
      dplyr::mutate(Species_full = trimws(Species_full, which = c("right"))) %>%
      # Write the full species epithet with author name
      dplyr::mutate(Species_full_author = paste(Species, Var, author, Ssp, Ssp_author)) %>%
      # replace double white spaces with one space
      dplyr::mutate(Species_full_author = str_replace_all(Species_full_author, "  ", " ")) %>%
      # remove white spaces at the end of a string
      dplyr::mutate(Species_full_author = trimws(Species_full_author, which = c("right")))
    # otherwise, leave the data frame alone
  } else {
    out_df <- in_df
  }
  # Initialize the synonym list data frame
  synonyms_list <- tibble(
    name_in = as.character(),
    name_out = as.character()
  )
  # Partition the list of species into chunks of 1000 names
  # otherwise it may be too large to pass to gnr_resolve()
  species_list <- unique(out_df$Species_full)
  species_num <- length(species_list)
  num_seq <- 1 + 0:(species_num %/% 1000) * 1000
  
  # Collect synonyms
  for (ii in num_seq) {
    temp_seq <- ii + 0:999
    temp_species_list <- species_list[temp_seq]
    temp_species_out <- gnr_resolve(temp_species_list,
                                    data_source_ids = 198, # The Leipzig Catalogue of Vascular Plants
                                    with_canonical_ranks = T,
                                    best_match_only = TRUE
    ) %>%
      select(name_in = user_supplied_name, name_out = matched_name2)
    synonyms_list <- rbind(synonyms_list, temp_species_out)
  }
  
  num_name_changes <- dim(filter(synonyms_list, name_in != name_out))[1]
  print(paste0("There are ", num_name_changes, " necessary name changes."))
  
  # Add back in species that didn't have a synonym in the database
  synonyms_list <- rbind(synonyms_list,
                         tibble(name_in = species_list[!(species_list %in% synonyms_list$name_in)],
                                name_out = NA)
  )
  
  # Load synonym list created by Elisa Van Cleemput
  # Prefer name from GNR_resolver over custom list
  custom_synonyms_list <- read_csv("data/clean/synonym_list.csv", show_col_types = FALSE)
  
  # Check if there are disagreements between EVC's list and GNR
  GNR_synonyms <- rename(synonyms_list, original_name = name_in, Synonym_GNR = name_out)
  compare_synonyms_df <- left_join(custom_synonyms_list, GNR_synonyms)
  syn_disagree_df <- filter(compare_synonyms_df, Synonym != Synonym_GNR)
  cat(paste0("There are ", dim(syn_disagree_df)[1], " synonym disagreements between GNR_resolve and our custom synonym list. \n
               Check the file data/clean/synonym_disagreements.csv for details."))
  write_csv(syn_disagree_df, "data/clean/synonym_disagreements.csv")
  
  # If there are species not in GNR, use Synonym from custom list
  synonyms_list <- left_join(synonyms_list, 
                             rename(custom_synonyms_list, 
                                    name_in = original_name)) %>% 
    mutate(name_out = ifelse(is.na(name_out), Synonym, name_out)) %>% 
    select(-Synonym)
  
  final_df <- out_df %>%
    # add column of synonym names
    left_join(rename(synonyms_list, Species_full = name_in, Synonym = name_out),
              by = "Species_full"
    ) %>%
    # If synonym isn't a binomial, revert to original name
    mutate(Synonym = ifelse(is.na(stringr::word(Synonym, 2, 2)),
                            Species_full, Synonym)) %>% 
    # keep track of original name
    mutate(original_name = Species_full) %>%
    # rename species to their main synonym
    mutate(Species_full = Synonym) %>%
    # update the genus name with the synonym genus
    mutate(Genus = stringr::word(Species_full, 1, 1)) %>%
    # set binomial to synonym binomial
    mutate(Species = stringr::word(Species_full, 1, 2))
}


# 1) Load in main species list --------------------------------------------
clade_labels <- c(
  "Gymnosperms", "Ferns and fern-allies", "Flowering plants (angiosperms)",
  "Species extracted from",
  "There are also algae, lichens, fungi and bryophytes listed in the document. We could decide to omit these"
)

# Read in our base data set
base_data <- read_csv("data/raw/PNW_Species_w_Metadata.csv", show_col_types = FALSE) %>%
  # Remove any fully empty data columns
  select_if(dataCheck_func) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove rows with NA as species name
  filter(!is.na(Species)) %>%
  # Remove species for which we lack phylogenetic data: Pterospora andromedea
  # filter(Species != "Pterospora andromedea") %>% 
  synonym_func()

full_species_synonyms_list <- data.frame(
  original_name = base_data$original_name,
  Synonym = base_data$Synonym
) %>%
  unique()

# Read in TEK data set
TEK_data <- read_csv("data/raw/TEKdata.csv", show_col_types = FALSE) %>%
  # Remove rows corresponding to family labels
  filter(!Species %in% clade_labels) %>%
  # Remove rows with NA as species name
  filter(!is.na(Species)) %>%
  # Rename TEK columns for easier reference
  rename(
    N_Names = `Number of unique names (Indigenous)`,
    N_Langs = `Number of unique  languages (from Appendix 2B)`,
    N_Uses = `Number of uses`
  ) %>% 
  synonym_func()

full_species_synonyms_list <- rbind(full_species_synonyms_list,
                                    data.frame(
                                      original_name = TEK_data$original_name,
                                      Synonym = TEK_data$Synonym
                                    )
) %>%
  unique()

# species in base_data but not in TEK_data: L Hierochloe, Ledum groenlandicum, and Ledum glandulosum

full_data <- full_join(select(base_data, 
                              c("Synonym", "Family", "IUCN Status")), 
                       select(TEK_data,
                              c("Synonym", "N_Names":"N_Uses")), 
                       by = "Synonym") %>%
  # remove species for which we lack TEK data
  filter(!is.na(N_Names))

# # Uncomment and run the code below to take a closer look at which species were renamed
# # Check for synonymous species with multiple entries
# duplicate_df <- full_data %>%
#   group_by(Species_full) %>%
#   filter(n() > 1)
#
# # Find which species were in our original TEK dataset that are missing from the full dataset
# full_species <- full_data$Species
# TEK_species <- TEK_data$Species
#
# lost_species <- TEK_species[which(!(TEK_species %in% full_species))]
#
# # Check that these were all renamed to synonyms
# all(lost_species %in% full_data$original_name) # TRUE
# true_lost_species <- lost_species[which(!(lost_species %in% full_data$original_name))]

# 2) Load in and process Diaz data set ------------------------------------

# Load in trait data compiled by Diáz et al. 2020 (based on TRY)
# https://www.nature.com/articles/s41597-022-01774-9
Diaz_data <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx",
                        sheet = 1
)

Diaz_data_renamed <- Diaz_data %>%
  rename(Species_full = `Species name standardized against TPL`)

# Assign GNR synonyms to species in Diaz data
run_GNR_diaz <- FALSE # Change this to TRUE to re-run the synonym resolver (which takes a very long time!)
if (run_GNR_diaz) {
  # rename species with synonyms
  Diaz_synonyms <- Diaz_data_renamed %>% 
    synonym_func() %>%
    select(original_name, Synonym)
  
  # Save the synonyms so we don't have to poll the name resolver every time
  write_csv(Diaz_synonyms, "data/clean/Diaz_synonyms.csv")
} else {
  Diaz_synonyms <- read_csv("data/clean/Diaz_synonyms.csv", show_col_types = FALSE) %>%
    distinct(original_name, Synonym)
}

Diaz_data_joined <- Diaz_data_renamed %>%
  rename(original_name = Species_full) %>% 
  left_join(Diaz_synonyms, by = "original_name") %>%
  # update the genus name with the synonym genus
  mutate(Genus = stringr::word(Synonym, 1, 1)) %>%
  # set binomial to synonym binomial
  mutate(Species = stringr::word(Synonym, 1, 2))


# Add new original names to our synonyms list
# Only keep ones that have been assigned Synonyms which are in our main dataset
full_species_synonyms_list <- rbind(
  full_species_synonyms_list,
  filter(Diaz_synonyms, Synonym %in% full_data$Synonym)
) %>%
  unique()

# # Uncomment the code below to see which of our species are missing from the Diaz data set
# # Get list of species missing from Diaz data
# Diaz_species_list <- Diaz_data_renamed$Species_full
# Diaz_species_short <- Diaz_data_renamed$Species
# full_species_list <- full_data$Species_full
# full_species_short <- full_data$Species
#
# # 47 missing if using full species epithet
# missing_from_Diaz <- full_species_list[which(!(full_species_list %in% c(Diaz_species_list)))]
#
# # 38 missing if using shortened species epithet
# missing_from_Diaz_short <- full_species_short[which(!(full_species_short %in% c(Diaz_species_list, Diaz_species_short)))]

# Reduce Diaz data to species present in the main data set
Diaz_data_sel <- Diaz_data_joined %>%
  # Remove non-focal species
  filter(Synonym %in% c(full_data$Synonym)) %>% # NB: including shortened species names only leads to the additional inclusion of Alnus incana
  # reduce to relevant set of traits
  select(
    "Synonym", "original_name", "Genus", "Family", "Woodiness", "Growth Form", "Leaf area (mm2)",
    "Nmass (mg/g)", "LMA (g/m2)", "Plant height (m)", "Diaspore mass (mg)",
    "LDMC (g/g)"
  ) %>%
  rename(original_name_Diaz = original_name) %>%
  mutate(origin = "Diaz") %>% 
  # There is one case of a species having multiple entries for a trait. In this case, we use the average among the entries
  group_by(Synonym) %>% 
  mutate(`Plant height (m)` = ifelse(Synonym %in% c("Potamogeton natans", "Ericameria nauseosa"),
                                     mean(`Plant height (m)`), 
                                     `Plant height (m)`)) %>% 
  mutate(`Diaspore mass (mg)` = ifelse(Synonym == "Ericameria nauseosa",
                                       mean(`Diaspore mass (mg)`), 
                                       `Diaspore mass (mg)`)) %>% 
  # Remove the extra entries for these two species
  filter(!is.na(`Leaf area (mm2)`) | sum(!is.na(`Leaf area (mm2)`)) == 0)

# 3) Load in and process TRY data sets ------------------------------------

# We include these to obtain trait data for the 23 (or 39) missing from Diaz
# TRY data is separated by trait
TRY_plant_height <- read_csv("data/raw/TRY_data_Feb2023/TRY_plant_height_growth.csv",
                             show_col_types = FALSE
) %>%
  dplyr::select(-...1)
plant_height_trait_names <- unique(TRY_plant_height$TraitName)

TRY_leafN <- read_csv("data/raw/TRY_data_Feb2023/TRY_leafN.csv",
                      show_col_types = FALSE
) %>%
  dplyr::select(-...1)
leafN_trait_names <- unique(TRY_leafN$TraitName)

TRY_LDMC <- read_csv("data/raw/TRY_data_Feb2023/TRY_LDMC.csv",
                     show_col_types = FALSE
) %>%
  dplyr::select(-...1)
LDMC_trait_names <- unique(TRY_LDMC$TraitName)

TRY_leafarea <- read_csv("data/raw/TRY_data_Feb2023/TRY_leafarea.csv",
                         show_col_types = FALSE
) %>%
  dplyr::select(-...1)
leafarea_trait_names <- unique(TRY_leafarea$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_rooting_depth <- read_csv("data/raw/TRY_data_Feb2023/TRY_rooting_depth.csv",
                              show_col_types = FALSE
) %>%
  dplyr::select(-...1)
rooting_depth_trait_names <- unique(TRY_rooting_depth$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_stem_density <- read_csv("data/raw/TRY_data_Feb2023/TRY_stem_density.csv",
                             show_col_types = FALSE
) %>%
  dplyr::select(-...1)
stem_density_trait_names <- unique(TRY_stem_density$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_plant_woodiness <- read_csv("data/raw/TRY_data_Feb2023/TRY_plant_woodiness.csv",
                                show_col_types = FALSE
) %>%
  dplyr::select(-...1)
plant_woodiness_trait_names <- unique(TRY_plant_woodiness$TraitName)

# this trait was recommended for inclusion in the recommended traits paper
TRY_plant_veg_reproduction <- read_csv("data/raw/TRY_data_Feb2023/TRY_vegetative_reproduction.csv",
                                       show_col_types = FALSE
) %>%
  dplyr::select(-...1)
plant_veg_reproduction_trait_names <- unique(TRY_plant_veg_reproduction$TraitName)

# combine TRY trait data sets
TRY_data <- rbind(
  TRY_plant_height, TRY_leafN, TRY_LDMC, TRY_leafarea,
  TRY_rooting_depth, TRY_stem_density, TRY_plant_woodiness,
  TRY_plant_veg_reproduction
) %>%
  rename(Species_full = AccSpeciesName) %>%
  # deal with synonyms
  synonym_func() %>%
  # Remove non-focal species
  filter(Synonym %in% full_data$Synonym) %>%
  # Add Family names in
  right_join(select(full_data, Synonym, Family), by = "Synonym") %>% 
  filter(!is.na(Species_full))

full_species_synonyms_list <- rbind(
  full_species_synonyms_list,
  select(TRY_data, original_name, Synonym)
) %>%
  unique()

# # Uncomment this to see which species are missing from TRY data
# # Get list of species which are in base data but not TRY data set
# TRY_species <- unique(TRY_data$Species_full)
# full_species <- unique(c(full_data$Species_full, full_data$Species))
# # There are 170 species in our full dataset that are missing from the TRY data
# missing_TRY <- full_species[which(!(full_species %in% TRY_species))]

# Combine TRY data entries (i.e. take means etc.)
TRY_data_processed <- TRY_data %>%
  filter(
    # Exclude unused traits (see manuscript for rationale)
    !(TraitName %in% c(
      "Plant growth rate relative (plant relative growth rate, RGR)",
      "Plant vegetative reproduction: clonality of ramets",
      "Leaf area (in case of compound leaves: leaflet, petiole included)",
      "Leaf area (in case of compound leaves: leaf, petiole included)",
      "Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)",
      "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included",
      "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded"
    )),
    # Only include single values or measures of central tendency (this removes e.g. maximum, minimum estimates)
    ValueKindName %in% c(
      "Mean", "Best estimate", "Site specific mean",
      "Species mean", "Median", "Single", NA
    ),
    # Disregard rows with no values provided
    !is.na(OrigValueStr)
  ) %>%
  dplyr::select(DatasetID, Synonym, TraitName:ValueKindName, Comment:Family) %>%
  group_by(DatasetID, Synonym, TraitName, ValueKindName) %>%
  # Additional processing needed for the "Plant woodiness" categorical trait
  mutate(OrigValueStr = case_when(
    TraitName == "Plant woodiness" ~ case_when(
      OrigValueStr %in% c("non-woody", "h", "H", "N", "0") ~ 0, # "non-woody",
      OrigValueStr == "1" ~ ifelse(Comment == "0=nonwoody; 1=woody", 2, 1), # "woody", "semi-woody"),
      OrigValueStr %in% c("woody", "w", "W", "Y") ~ 2, # "woody",
      OrigValueStr %in% c("semi-woody") ~ 1, # "semi-woody",
      OrigValueStr == "2" ~ 2, # "soft wood",
      OrigValueStr == "3" ~ 2, # "dense wood"
    ),
    TRUE ~ as.double(OrigValueStr)
  )) %>%
  mutate(processed_value = case_when(
    # Leave the value alone if its already a measure of central tendency
    ValueKindName %in% c("Mean", "Site specific mean", "Species mean", "Median") ~ OrigValueStr,
    # Take the mean of single values
    ValueKindName == "Single" ~ (mean(as.double(OrigValueStr))),
    TRUE ~ OrigValueStr
  )) %>%
  select(-OrigValueStr) %>%
  distinct() %>%
  ungroup() %>%
  select(Synonym, TraitName, processed_value, Synonym, original_name, Family, Species) %>%
  arrange(Synonym, TraitName)

TRY_data_processed_wide <- TRY_data_processed %>%
  # Rename and combine "synonymous" traits
  # These trait names are chosen to be similar to those used in Diaz
  mutate(TraitName = case_when(
    TraitName %in% c("Plant woodiness") ~ "Woodiness",
    # Combine all the measures of leaf area
    TraitName %in% c(
      "Leaf area (in case of compound leaves: leaflet, petiole excluded)",
      "Leaf area (in case of compound leaves: leaf, petiole excluded)"
    ) ~ "Leaf area (mm2)",
    # Combine all the measures of leaf area per leaf dry mass
    TraitName %in% c(
      "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded"
    ) ~ "SLA",
    TraitName %in% c("Leaf nitrogen (N) content per leaf dry mass") ~ "Nmass (mg/g)",
    TraitName %in% c("Plant height vegetative") ~ "Plant height (m)",
    TraitName %in% c() ~ "Diaspore mass (mg)",
    TraitName %in% c("Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)") ~ "LDMC (g/g)",
    # Extra ones I added just to make the trait names shorter
    TraitName %in% c("Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density") ~ "SSD",
    TraitName %in% c("Root rooting depth") ~ "RRD",
    TraitName %in% c("Leaf nitrogen (N) content per leaf area") ~ "LNLA",
    TRUE ~ TraitName # Just "Root rooting depth" and "Leaf nitrogen (N) content per leaf dry mass" right now
  )) %>%
  mutate(processed_value = as.double(processed_value)) %>%
  pivot_wider(
    id_cols = c(Synonym, Family, Species, original_name),
    names_from = TraitName, values_from = processed_value,
    values_fn = mean
  ) %>%
  # Translate Woodiness back to a character
  mutate(Woodiness = case_when(
    round(Woodiness) == 0 ~ "non-woody",
    round(Woodiness) == 1 ~ "semi-woody",
    round(Woodiness) == 2 ~ "woody"
  )) %>%
  rename(original_name = original_name) %>%
  mutate(origin = "TRY")


# 4) Combine data sets -----------------------------------------------------

# Deal with any species with multiple synonyms
synonym_doubles <- full_species_synonyms_list %>%
  # We don't need to keep track of "synonyms" if they're the same as the original name
  filter(original_name != Synonym) %>% # 15 left
  # Check if there are multiple synonyms assigned to one original name
  group_by(original_name) %>%
  # filter(n() > 1) %>% # 16 left, these all seem to differ in how subspecies or variant is spelled
  mutate(original_sspvar = stringr::word(original_name, 3, 3, sep = " ")) %>%
  mutate(synonym_sspvar = stringr::word(Synonym, 3, 3, sep = " "))

write_csv(synonym_doubles, "data/clean/synonym_doubles.csv")

# Save a list of all the species synonyms we encountered in the datasets
custom_synonyms_list <- read_csv("data/clean/synonym_list.csv", show_col_types = FALSE)
full_species_synonyms_list %>%
  # Add in Elisa's custom list
  rbind(rename(custom_synonyms_list, Synonym = Synonym)) %>%
  unique() %>%
  write_csv("data/clean/all_synonyms.csv")

# Join Diaz and TRY data sets, preferring data entries from Diaz

# Get coverage of traits
coverage_df <-  rbind(
  Diaz_data_sel %>% rename(original_name = original_name_Diaz),
  TRY_data_processed_wide,
  # Add back in species which lack trait data along with TEK data
  full_data %>% mutate(origin = "this_study")
) %>% 
  ungroup() %>%
  summarise(across(everything(), ~ sum(!is.na(.)) / dim(TRY_data_processed_wide)[1]))

final_data <- rbind(
  Diaz_data_sel %>% rename(original_name = original_name_Diaz),
  TRY_data_processed_wide,
  # Add back in species which lack trait data along with TEK data
  full_data %>% mutate(origin = "this_study")
) %>%
  arrange(Synonym, origin) %>%
  # if there are data for a trait from both the Diaz and TRY data sets, just use the value from Diaz
  pivot_longer(
    cols = c(-"Synonym", -"origin"),
    names_to = "trait", values_to = "value",
    values_transform = list(value = as.character)
  ) %>%
  # for each species, for each trait,
  pivot_wider(
    names_from = "origin",
    values_from = "value"
  ) %>%
  # if value is not NA for origin = "TRY" and origin = "Diaz",
  # set value to value from "Diaz"
  group_by(Synonym, trait) %>%
  mutate(
    value = case_when(
      !is.na(Diaz) & !is.na(TRY) ~ Diaz,
      is.na(Diaz) & !is.na(TRY) ~ TRY,
      !is.na(Diaz) & is.na(TRY) ~ Diaz,
      is.na(Diaz) & is.na(TRY) ~ this_study,
      TRUE ~ NA
    ),
    # delete origin column
    .keep = "unused"
  ) %>%
  pivot_wider(names_from = "trait") %>%
  # Change numeric variables back to type double
  mutate_at(vars(`LDMC (g/g)`:`Leaf area (mm2)`, "Diaspore mass (mg)":"N_Uses"), as.double) %>%
  ungroup() %>%
  mutate(Species = stringr::word(Synonym, 1, 2, sep = " ")) %>%
  mutate(Ssp_var = stringr::word(Synonym, 3, 4, sep = " "))

# Look for repeats of species
subspecies_list <- final_data %>%
  group_by(Species) %>%
  filter(n() > 1)

# List of species for which we have zero numeric trait data besides TEK
test_data <- final_data %>%
  rowwise() %>%
  mutate(n_notNA = sum(!is.na(c_across(where(is.double))))) %>%
  filter(n_notNA == 3) # All species should have at least 3 trait values (from TEK data)

# Get coverage of traits
coverage_df <- final_data %>%
  ungroup() %>%
  summarise(across(everything(), ~ sum(!is.na(.)) / dim(final_data)[1]))

write_csv(coverage_df, "data/clean/coverage_pcts.csv")


# 5) Impute missing trait data --------------------------------------------

# Write the final data from before imputation
write_csv(final_data, "data/clean/dataset_no_imputation.csv")

# Put trait data in workable form for imputation
data_in <- final_data %>%
  select(
    "Synonym", "Family", "LDMC (g/g)", "Plant height (m)",
    "Nmass (mg/g)", "Leaf area (mm2)", "Woodiness", "Diaspore mass (mg)"
  ) %>%
  mutate(Genus = stringr::word(Synonym, 1, 1, sep = " ")) %>%
  mutate(Species = stringr::word(Synonym, 2, 2, sep = " ")) %>%
  mutate(Ssp_var = stringr::word(Synonym, 3, 4, sep = " ")) %>%
  relocate(Genus, Species, Ssp_var)

# Modify original data
traits_df <- data_in %>%
  mutate(LeafArea_log = log(`Leaf area (mm2)`)) %>%
  mutate(PlantHeight_log = log(`Plant height (m)`)) %>%
  mutate(DiasporeMass_log = log(`Diaspore mass (mg)`)) %>%
  mutate(LDMC_log = log(`LDMC (g/g)`)) %>%
  # Remove original variables replaced by "logged" versions
  select(-c("LDMC (g/g)", "Leaf area (mm2)", "Plant height (m)", "Diaspore mass (mg)"))

# Run missForest algorithm to impute missing traits
# get taxonomic columns
taxa <- select(traits_df, c(Genus, Family)) # Get order and family level data
test <- rep(1, nrow(taxa))
taxa <- cbind(taxa, test)

# turn taxonomic groupings into dummy binary variables
dummies <- dummyVars(test ~ ., data = taxa)
taxa <- predict(dummies, newdata = taxa) #%>% 

# select traits to be imputed
traits_to_impute <- as.data.frame(traits_df) %>%
  select(-c("Genus":"Family")) %>%
  mutate(Woodiness = as.factor(Woodiness))

# # Uncomment this to examine if there are any Genera for which there are species with no trait data
# emptyGenera_data <- traits_to_impute  %>%
#   rowwise() %>%
#   mutate(n_notNA = sum(!is.na(c_across(where(is.double))))) %>% 
#   relocate(n_notNA) %>% 
#   filter(n_notNA == 0) %>% 
#   pivot_longer(cols = GenusAbies:GenusZostera, names_to = "Genus") %>% 
#   filter(value == 1) %>% select(Genus)

# combine dummy vars with traits (comment this line out to not use taxa in imputation)
traits_to_impute <- cbind(traits_to_impute, taxa) %>% 
  mutate_if(is.character, function(x) {as.double(x) * 100})
true_df = traits_to_impute %>% filter_at(vars(`LDMC_log`:`DiasporeMass_log`), all_vars(!is.na(.)))

# run missForest imputation
PNW_imp <- missForest(traits_to_impute,
                      maxiter = 10000, # maximum number of iterations to be performed given the stopping criterion isn't met
                      ntree = 10000, # number of trees to grow in each forest
                      verbose = TRUE, # if 'TRUE', gives additional output between iterations
                      variablewise = TRUE # if 'TRUE', the OOB error is returned for each variable separately
)

imp_df = PNW_imp$ximp
# traits_to_impute
trait_names = c("LDMC_log", "Nmass (mg/g)", "LeafArea_log", "PlantHeight_log", "DiasporeMass_log")

# add actual taxonomic columns back in and remove dummy variables
imputed_traits <- cbind(traits_df[, 1:5], PNW_imp$ximp[, 1:6]) 

imputed_data <- right_join(
  select(final_data, -c("LDMC (g/g)":"Diaspore mass (mg)")),
  select(imputed_traits, -c(Genus:Ssp_var))
)

full_df <- traits_df %>%
  # keep track of which traits were original values and not imputed values
  mutate(type = "original") %>%
  rbind(imputed_traits %>% mutate(type = "imputed")) %>%
  mutate("LDMC (g/g)" = exp(LDMC_log), .keep = "unused") %>%
  mutate("Diaspore mass (mg)" = exp(DiasporeMass_log), .keep = "unused") %>%
  mutate("Plant height (m)" = exp(PlantHeight_log), .keep = "unused") %>%
  mutate("Leaf area (mm2)" = exp(LeafArea_log), .keep = "unused") %>%
  select(-c("Genus", "Species", "Ssp_var", "Family")) %>%
  melt(id = c("Synonym", "type")) %>%
  mutate(value = as.numeric(value))

# Calculate variance in each continuous variable
var_df <- imputed_traits %>%
  select(c("Nmass (mg/g)":"LDMC_log", -"Woodiness")) %>%
  summarise(across("Nmass (mg/g)":"LDMC_log", var)) %>% 
  cbind("Woodiness" = NA) %>% 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "variance")

error_df <- as.list(PNW_imp$OOBerror[1:6]) %>%
  as_tibble(.name_repair = "universal") %>%
  melt() %>%
  rename(error_value = value) %>%
  mutate(error_type = substr(as.character(variable), start = 1, stop = 3),
         .keep = "unused") %>%
  cbind(variable = colnames(PNW_imp$ximp[,1:6])) %>%
  # compute NRMSE from MSE
  inner_join(var_df, by = "variable") %>%
  mutate(RMSE = sqrt(error_value)) %>%
  mutate(NRMSE = RMSE / variance)

# 6) Put data set into workable form for imputation and phylogeny steps --------

write_csv(imputed_data, "data/clean/final_dataset.csv")