################################################################################
# Species presence and ecoregions
################################################################################

## Project: Comparing biodiversity hotspot identification
##
## Purpose: Determine species occurrence in the ecoregions of the PNW
##
## Contents: 0) Load packages
##           1) Load data for BIEN
##           2) Define the PNW
##           3) Extract the ecoregion(s) of each species
##           4) Extract species for PNW ecoregions
##           5) Organize and save species occurrence data
##
## Inputs:  - data/clean/all_synonyms.csv
##          - data/clean/final_dataset.csv
##          - data/raw/NA_Terrestrial_Ecoregions_v2_Level_III_Shapefile/NA_TerrestrialEcoregions_LIII/data/NA_Terrestrial_Ecoregions_v2_level3.shp
##          
## Outputs: - data/clean/species_occurrences.rds
##          - data/clean/missing_species_occurrences.rds
##
## Written and by: Elisa Van Cleemput
## Maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized April 2023
## _____________________________________________________________________________

# 0) Load packages --------------------------------------------------------
library(sp) # for "over"
library(tidyverse) # for "read_csv"
if (!require('BIEN')) install.packages('BIEN'); library('BIEN')
library(progress) # for keeping track of progress of long processes
library(doParallel) # to run code in parallel
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('maps')) install.packages('maps'); library('maps')
if (!require('grafify')) install.packages('grafify'); library('grafify')

# 1) Load data for BIEN ---------------------------------------------------

# Load in list of all the species synonyms we found in the original datasets
synonym_list <- read_csv("data/clean/all_synonyms.csv") %>%
  mutate(dropped_sspvar = stringr::word(original_name, 1,2))

all_species <- tibble(species = c(synonym_list$original_name,
                                  synonym_list$Synonym,
                                  synonym_list$dropped_sspvar)) %>%
  unique()

species_list <- BIEN_list_all() %>%
  filter(species %in% all_species) %>%
  mutate(Species_full = species) %>%
  right_join(rename(synonym_list, species = original_name), by = "species") %>%
  select(species = Synonym) %>%
  unique()

# 1. Restrict species we're looking for occurrences for to big list of species names
# 2. Find occurrence data for each of these species at Levels 1, 2, and 3
# 3. Then we need to go from the big list of species to our "synonymized" list
# 4. Assign a presence if ANY synonymous species are present at a level

## ---- Load list of species --------------
# Check which of our species are in the BIEN database
final_data_in_BIEN <- read_csv("data/clean/final_dataset.csv") %>%
  filter(Synonym %in% c(species_list$species, species_list$Synonym)) %>%
  select(Synonym) %>%
  unique()

## ---- Load Ecoregions map --------------
# Ecoregions II and III maps donwloaded on April 14 2023
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-ii/
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-iii/

Ecoregions3 <- read_sf("data/raw/NA_Terrestrial_Ecoregions_v2_Level_III_Shapefile/NA_TerrestrialEcoregions_LIII/data/NA_Terrestrial_Ecoregions_v2_level3.shp")
# Ecoregions dataset 3 contains information about level 1, 2 and 3


# 2) Define the PNW -------------------------------------------------------
# The Pacific Northwest of North America contains only the following Ecoregions: 7.1, 6.1,6.2,10.1,3.2,3.1,3.3,2.2,2.3,5.4

Ecoregions3_PNW <- Ecoregions3 %>%
  filter(LEVEL2 %in% c("2.2", "2.3", "3.1", "3.2", "3.3", "5.4","6.1", "6.2", "7.1", "10.1"))


# 3) Extract the ecoregion(s) of each species -----------------------------

# This idea is based on Echeverría-Londoño et al. 2018:
# "We overlaid the BIEN 2.0 plant species maps on a 100 × 100 km grid map with a
# Lambert Azimuthal Equal Area projection to obtain a presence/absence matrix of species for each grid cell."

# function to extract species occurence and PNW regions for each species
extract_species_occ <- function(species, Ecoregions){
  start = Sys.time()
  print(paste0("Extracting information for species ", species))
  
  # Extract species occurrences from BIEN database
  # This can take a while
  species_occ <- BIEN::BIEN_occurrence_species(species,
                                               cultivated=F, natives.only=F,
                                               new.world = T, # we only care about new world occurrences
                                               all.taxonomy = F) %>% 
    rename(species_full = scrubbed_species_binomial)
  
  if (nrow(species_occ) == 0){
    print("The species is has no occurrence data")
    occ_time <- Sys.time() - start
    output <- data.frame(species = species, LEVEL3 = NA, occ_time = occ_time)
  } else {
    
    species_occ <- species_occ %>%
      drop_na(longitude) %>%
      drop_na(latitude)
    
    if (nrow(species_occ) == 0){
      print("The species is has no records for longitude and latitude")
      occ_time <- Sys.time() - start
      output <- data.frame(species = species, LEVEL3 = NA, occ_time = occ_time)
    } else {
      
      # Project occurrence map to the same projection as the PNW map
      # The BIEN dataset is in WGS84 (EPSG: 4326) :
      # https://github.com/NiDEvA/R-protocols/blob/main/R_Protocol_get%26curate_data.R#L330
      species_occ_Lb <- species_occ %>%
        st_as_sf(coords = c("longitude","latitude")) %>%
        st_set_crs(4326) %>%
        st_transform(st_crs(Ecoregions))
      
      # Overlap BIEN occurrence data with PNW map: Extract ecoregion(s) where species occurs
      species_occ_eco <- data.frame(species = character(), LEVEL3 = double())
      for (index_species in species) {
        filter_occ = filter(species_occ_Lb, species_full == index_species)
        
        temp_eco <- over(as_Spatial(filter_occ), as_Spatial(Ecoregions)) %>%
          select(LEVEL3) %>%
          distinct() %>%
          filter(!is.na(LEVEL3)) %>%
          mutate(species = index_species)
        
        species_occ_eco <- rbind(temp_eco, species_occ_eco)
        
      }
      
      output <- species_occ_eco %>% mutate(occ_time = occ_time)# Levels 1 and 2 can be recreated later by subtracting suffixes
    }
    
  }
  return(output)
}

# 4) Extract species for PNW ecoregions -----------------------------------

# Start new cluster for parallel processing
cluster_size <- parallel::detectCores()

my.cluster <- parallel::makeCluster(
  cluster_size,
  type = "PSOCK"
)
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

# Function for slicing the data into memory-manageable chunks
slice_func<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

species_slice <- slice_func(species_list$species, 1)

# Set up progress bar
pb <- progress_bar$new(
  format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
  total = length(species_slice),
  width = 120)
progress <- function(n){
  pb$tick(tokens = list(system = n))
}
opts <- list(progress = progress)

# Calculate species occurrences
species_occurrences <- foreach(index_species = species_slice,
                               .combine = 'rbind',
                               .packages = c('tidyverse', 'BIEN', 'sp', 'sf'),
                               .options.snow = opts) %dopar% {
                                 extract_species_occ(index_species, Ecoregions3_PNW)
                               }


# 5) Organize and save species occurrence data ----------------------------

# Assign proper synonyms to each species in list
species_occurrences_out <- species_occurrences %>%
  select(-occ_time) %>%
  left_join(rename(synonym_list, species = original_name),
            relationship = "many-to-many") %>%
  select(Synonym, Level = LEVEL3) %>%
  unique()

write_rds(species_occurrences_out, "data/clean/species_occurrences.rds", compress = 'gz')

missing_occurrence_data <- filter(species_occurrences_out, is.na(Level))
write_rds(missing_occurrence_data, "data/clean/missing_species_occurrences.rds", compress = 'gz')
