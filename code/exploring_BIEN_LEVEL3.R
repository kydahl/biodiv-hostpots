##############################################################################
###                                                                        ###
###                       Species presence and ecoregions                  ###
###                                                                        ###
##############################################################################

# Elisa Van Cleemput, April 2023
######################################################## Explore BIEN data

# Load packages
library(sp) # for "over"
library(tidyverse) # for "read_csv"
if (!require('BIEN')) install.packages('BIEN'); library('BIEN')
library(progress) # for keeping track of progress of long processes
library(doParallel) # to run code in parallel
# vignette("BIEN")
# vignette("BIEN_tutorial")
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('maps')) install.packages('maps'); library('maps')
if (!require('grafify')) install.packages('grafify'); library('grafify')

# Specify map where species occurrence maps will be saved
dir_fig <- "figures/Species_occurrence/"

#######################################################
### 1) Load information
#######################################################
## ---- Load BIEN list of species --------------

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

dim(species_list)
head(species_list)


# 1. Restrict species we're looking for occurrences for to big list of species names
# 2. Find occurrence data for each of these species at Levels 1, 2, and 3
# 3. Then we need to go from the big list of species to our "synonymized" list
# 4. Assign a presence if ANY synonymous species are present at a level

## ---- Load list of species --------------
# Check which of our species are the BIEN database

final_data <- read_csv("data/clean/final_dataset.csv")
dim(final_data) # 318
length(unique(final_data$Species_full)) # 2 species have the same "base name": Alnus viridis and Populus balsamifera

final_data_in_BIEN <- final_data %>%
  filter(Species_full %in% c(species_list$species, species_list$Synonym)) %>%
  select(Species_full) %>%
  unique()
dim(final_data_in_BIEN) # 318, so 0 missing species

## ---- Load Ecoregions map --------------
# Ecoregions II and III maps donwloaded on April 14 2023
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-ii/
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-iii/

Ecoregions3 <- read_sf("data/raw/NA_Terrestrial_Ecoregions_v2_Level_III_Shapefile/NA_TerrestrialEcoregions_LIII/data/NA_Terrestrial_Ecoregions_v2_level3.shp")
# Ecoregions dataset 3 contains information about level 1, 2 and 3

# Explore metadata
# str(Ecoregions2)
# class(Ecoregions2)
# Ecoregions2$geometry # Projected CRS: Sphere_ARC_INFO_Lambert_Azimuthal_Equal_Area
# st_crs(Ecoregions2_PNW) 

# Quick visualization 
# plot(Ecoregions2["LEVEL2"])
# plot(st_geometry(Ecoregions2), axes=T)


######################################################
### 2) Define the PNW
#######################################################
# According to Lynette, the following Ecoregions consitute the PNW: 7.1, 6.1,6.2,10.1,3.2,3.1,3.3,2.2,2.3,5.4 (meeting notes April 5 2023)

Ecoregions3_PNW <- Ecoregions3 %>%
  filter(LEVEL2 %in% c("2.2", "2.3", "3.1", "3.2", "3.3", "5.4","6.1", "6.2", "7.1", "10.1"))
# plot(Ecoregions3_PNW["LEVEL3"])

#######################################################
# 3) Loop through all species and extract the locations and ecoregion(s) they occur in
#######################################################
# This idea is based on Echeverría-Londoño et al. 2018: 
# "We overlaid the BIEN 2.0 plant species maps on a 100 × 100 km grid map with a 
# Lambert Azimuthal Equal Area projection to obtain a presence/absence matrix of species for each grid cell."

# function to extract species occurence and PNW regions for each species
extract_species_occ <- function(species, Ecoregions){
  start = Sys.time()
  print(paste0("Extracting information for species ", species))
  
  # Extract species occurrences from BIEN database
  # This can take a while
  system.time(
  species_occ <- BIEN::BIEN_occurrence_species(species,
                                               cultivated=F, natives.only=F,
                                               new.world = T, # we only care about new world occurrences
                                               all.taxonomy = F) %>%  # KD: do we need all the taxonomy info? This might not actually speed things up, however
    rename(species_full = scrubbed_species_binomial)
  )
  
  if (nrow(species_occ) == 0){
    print("The species is has no occurrence data")
    occ_time <- Sys.time() - start
    output <- data.frame(species = species, LEVEL3 = NA, occ_time = occ_time)
    # output <- list("species_occ_Lb" = NA,
    #                "species_occ_eco" = NA)
  } else {
    
    species_occ <- species_occ %>%
      drop_na(longitude) %>%
      drop_na(latitude)
    
    if (nrow(species_occ) == 0){
      print("The species is has no records for longitude and latitude")
      occ_time <- Sys.time() - start
      output <- data.frame(species = species, LEVEL3 = NA, occ_time = occ_time)
      # output <- list("species_occ_Lb" = NA,
      #                "species_occ_eco" = NA)
    } else {
      
      # Optional: Quickly visualize points on a map
      # library(maps)
      # map('world', fill = TRUE, col= "grey", bg = "light blue") 
      # points(cbind(species_occ$longitude,
      #              species_occ$latitude),
      #        col = "red",
      #        pch = 20,
      #        cex = 1) 
      
      # Project occurrence map to the same projection as the PNW map 
      # I assume that the BIEN dataset is in WGS84 (EPSG: 4326) but I can't find documentation on that
      # According to this code it is in WGS84 indeed:
      # https://github.com/NiDEvA/R-protocols/blob/main/R_Protocol_get%26curate_data.R#L330
      species_occ_Lb <- species_occ %>% 
        st_as_sf(coords = c("longitude","latitude")) %>%
        st_set_crs(4326) %>%
        st_transform(st_crs(Ecoregions))
      
      # Overlap BIEN occurrence data with PNW map: Extract ecoregion(s) where species occurs 
      # This takes a while!
      # Slow option, using sf objects:
      # system.time (species_occ_eco <- st_join(species_occ_Lb, Ecoregions, join = st_within))
      # Faster option, using sp objects
      # system.time (species_occ_eco <- over(as_Spatial(species_occ_Lb), as_Spatial(Ecoregions)))
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

      # KD: don't we just need the levels the species occur at at this point?
      occ_time <- Sys.time() - start
      
      output <- species_occ_eco %>% mutate(occ_time = occ_time)# Levels 1 and 2 can be recreated later by subtracting suffixes
      
      # output <- list("species_occ_Lb" = species_occ_Lb,
      #                "species_occ_eco" = species_occ_eco)
      
    }
    
  }
  return(output)
}
# Error issues that came up while running this function:
# - Error for Chamaecyparis nootkatensis, because the species has no ocurrence data --> error issue solved
# - Error for Pinus albicaulis, because of missing values in coordinates --> error issue solved

# Extract species occurence and PNW regions for each species
# Start new cluster for doParallel
cluster_size <- parallel::detectCores()

my.cluster <- parallel::makeCluster(
  cluster_size, 
  type = "PSOCK"
)
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

slice_func<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

species_slice <- slice_func(species_list$species, 1) # ceiling(length(species_list$Species)/cluster_size))

pb <- progress_bar$new(
  format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
  total = length(species_slice),
  width = 120)          

progress <- function(n){
  pb$tick(tokens = list(system = n))
}
opts <- list(progress = progress)

#     should take the same amount of time on each core.
# Function: Slice data to optimize usage of parallel processing and memory


species_occurences <- foreach(index_species = species_slice,
                              .combine = 'rbind',
                              .packages = c('tidyverse', 'BIEN', 'sp', 'sf'),
                              .options.snow = opts) %dopar% {
  extract_species_occ(index_species, Ecoregions3_PNW)
                              }

#######################################################
# 4) Concatenate information from all species in a nice data frame
#######################################################

# Assign proper synonyms to each species in list
species_occurences_out <- species_occurences %>% 
  select(-occ_time) %>% 
  left_join(rename(all_synonyms, species = original_name),
            relationship = "many-to-many") %>% 
  select(species = Synonym, level = LEVEL3) %>% 
  unique()

write_rds(species_occurences_out, "data/clean/species_occurrences_LEVEL3_only.rds", compress = 'gz')