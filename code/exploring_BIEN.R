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

species_list <- BIEN_list_all() %>% 
  mutate(Species_full = species) %>%
  # filter to species in BIEN that match synonyms "original name" of our original species list
  synonym_func() # !!! re-run later when it lets me

dim(species_list)
head(species_list)

## ---- Load list of species --------------
# Check which of our species are the BIEN database

final_data <- read_csv("data/clean/final_dataset.csv")
dim(final_data) # 318
length(unique(final_data$Species_full)) # 2 species have the same "base name": Alnus viridis and Populus balsamifera

final_data_in_BIEN <- final_data %>%
  filter(Species %in% unique(species_list$species)) %>%
  select(Species) %>%
  unique()
dim(final_data_in_BIEN) # 315, so 1 missing species

species_list <- final_data_in_BIEN # We will use this dataframe in part 3)

## ---- Load Ecoregions map --------------
# Ecoregions II and III maps donwloaded on April 14 2023
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-ii/
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-iii/

Ecoregions2 <- read_sf("data/raw/NA_Terrestrial_Ecoregions_v2_Level_II_Shapefile/NA_TerrestrialEcoregions_LII/data/NA_Terrestrial_Ecoregions_v2_level2.shp")
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
# This is based on XXX
Ecoregions2_PNW <- Ecoregions2 %>%
  filter(LEVEL2 %in% c("2.2", "2.3", "3.1", "3.2", "3.3", "5.4","6.1", "6.2", "7.1", "10.1"))
# plot(Ecoregions2_PNW["LEVEL2"])

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
  
  print(paste0("Extracting information for species ", species))
  
  # Extract species occurrences from BIEN database
  # This can take a while
  species_occ <- BIEN::BIEN_occurrence_species(species,
                                               cultivated=F, natives.only=T)
  
  if (nrow(species_occ) == 0){
    print("The species is has no occurrence data")
    output <- list("species_occ_Lb" = NA,
                   "species_occ_eco" = NA)
  } else {
    
    species_occ <- species_occ %>%
      drop_na(longitude) %>%
      drop_na(latitude)
    
    if (nrow(species_occ) == 0){
      print("The species is has no records for longitude and latitude")
      output <- list("species_occ_Lb" = NA,
                     "species_occ_eco" = NA)
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
      system.time (species_occ_eco <- over(as_Spatial(species_occ_Lb), as_Spatial(Ecoregions)))
      
      output <- list("species_occ_Lb" = species_occ_Lb,
                     "species_occ_eco" = species_occ_eco)
      
    }
    
  }
}
# Error issues that came up while running this function:
# - Error for Chamaecyparis nootkatensis, because the species has no ocurrence data --> error issue solved
# - Error for Pinus albicaulis, because of missing values in coordinates --> error issue solved

# Extract species occurence and PNW regions for each species
# This takes a while!
library(doParallel)
numCores <- detectCores() - 2
cl <- makeCluster(numCores)
clusterExport(cl, c("Ecoregions3_PNW","extract_species_occ"))
clusterEvalQ(cl, {
  library(sp)
  library(tidyverse)
  library(BIEN)
  library(sf)
  library(ggplot2)
  library(maps)
  library(grafify)
})

system.time(species_occurences <- parApply(cl, species_list, 1, function(x) extract_species_occ(x, Ecoregions3_PNW)))
# system.time(species_occurences <- apply(species_list, 1, function(x) extract_species_occ(x, Ecoregions3_PNW)))
names(species_occurences) <- species_list$Species

# save the output
write_rds(species_occurences, "data/clean/species_occurences.rds", compress = 'gz')

# FYI:
# - NA for Chamaecyparis nootkatensis,  Pinus albicaulis, Asarum caudatum, Calypso bulbosa, 
#          Goodyera oblongifolia, Opuntia fragilis, Opuntia polyacantha, Platanthera dilatata, Platanthera stricta


#######################################################
# 4) Concatenate information from all species in a nice data frame
#######################################################
species_occurences <- read_rds("data/clean/species_occurences.rds")

species_occurences_eco = lapply(species_occurences, function(x) x[["species_occ_eco"]])
# Remove species that have no occurrence data
species_occurences_eco = species_occurences_eco[-(which(sapply(species_occurences_eco,is.logical),arr.ind=TRUE))]

# Create species list for each Ecoregion 1
species_occurences_L1 = lapply(species_occurences_eco, function(x) unique(x["LEVEL1"])) %>%
  bind_rows(.id = "groups") %>%
  rename(Species = groups) %>% 
  `rownames<-`( NULL ) %>%
  drop_na()

# Create species list for each Ecoregion 2
species_occurences_L2 = lapply(species_occurences_eco, function(x) unique(x["LEVEL2"])) %>%
  bind_rows(.id = "groups") %>%
  rename(Species = groups) %>% 
  `rownames<-`( NULL ) %>%
  drop_na()

# Create species list for each Ecoregion 3
species_occurences_L3 = lapply(species_occurences_eco, function(x) unique(x["LEVEL3"])) %>%
  bind_rows(.id = "groups") %>%
  rename(Species = groups) %>% 
  `rownames<-`( NULL ) %>%
  drop_na()

dim(species_occurences_L1)
unique(species_occurences_L1$LEVEL1)
dim(species_occurences_L2)
unique(species_occurences_L2$LEVEL2)
dim(species_occurences_L3)
unique(species_occurences_L3$LEVEL3)

# Remind from 1) and 3) that some species are missing!
length(unique(species_occurences_L2$Species))

# save the output
write_csv(species_occurences_L1, "data/clean/species_occurences_L1.csv")
write_csv(species_occurences_L2, "data/clean/species_occurences_L2.csv")
write_csv(species_occurences_L3, "data/clean/species_occurences_L3.csv")

#######################################################
# 5) Visualize where species occur in the PNW
#######################################################
species_occurences <- read_rds("data/clean/species_occurences.rds")

species_occurences_Lb = lapply(species_occurences, function(x) x[["species_occ_Lb"]])
# Remove species that have no occurrence data
species_occurences_Lb = species_occurences_Lb[-(which(sapply(species_occurences_Lb,is.logical),arr.ind=TRUE))]

# Crop world map to the PNW
world_cropped <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE)) %>%
  st_transform(st_crs(Ecoregions2_PNW))%>%
  st_crop((Ecoregions2_PNW))

# Visualize occurrence points of each species
for (i in 15:length(species_occurences_Lb)){
  print(paste0("Mapping ", names(species_occurences_Lb)[i]))
  
  # Crop species occurrences to the PNW (this takes a while)
  species_occurences_Lb_PNW <- species_occurences_Lb[[i]] %>%
    st_intersection((Ecoregions2_PNW))
  
  if (nrow(species_occurences_Lb_PNW) == 0){
    print("This species has no occurrence data in the PNW.")
  } else {
    g1 <- ggplot() +
      geom_sf(data=world_cropped) +
      geom_sf(data=Ecoregions2_PNW, aes(fill=LEVEL1)) +
      scale_fill_grafify(palette="bright") + 
      geom_sf(data=species_occurences_Lb_PNW, color="black") +
      theme_bw()
    g2 <- ggplot() +
      geom_sf(data=world_cropped) +
      geom_sf(data=Ecoregions2_PNW, aes(fill=LEVEL2)) +
      scale_fill_grafify(palette="all_grafify") + # long palettes: all_grafify, kelly, safe
      geom_sf(data=species_occurences_Lb_PNW, color="black") +
      theme_bw()
    g3 <- ggplot() +
      geom_sf(data=world_cropped) +
      geom_sf(data=Ecoregions3_PNW, aes(fill=LEVEL3)) +
      scale_fill_grafify(palette="kelly") + 
      geom_sf(data=species_occurences_Lb_PNW, color="black") +
      theme_bw()
    
    ggsave(paste0(dir_fig, "/Ecoregion_L1/", names(species_occurences_Lb)[i],".png"), g1, dpi=300,
           width = 7, height = 7, units = "in")
    ggsave(paste0(dir_fig, "/Ecoregion_L2/", names(species_occurences_Lb)[i],".png"), g2, dpi=300,
           width = 7, height = 7, units = "in")
    ggsave(paste0(dir_fig, "/Ecoregion_L3/", names(species_occurences_Lb)[i],".png"), g3, dpi=300,
           width = 9, height = 7, units = "in")
    
    ggsave(paste0(dir_fig, "/Ecoregion_L1/", names(species_occurences_Lb)[i],".pdf"), g1, dpi=300,
           width = 7, height = 7, units = "in")
    ggsave(paste0(dir_fig, "/Ecoregion_L2/", names(species_occurences_Lb)[i],".pdf"), g2, dpi=300,
           width = 7, height = 7, units = "in")
    ggsave(paste0(dir_fig, "/Ecoregion_L3/", names(species_occurences_Lb)[i],".pdf"), g3, dpi=300,
           width = 9, height = 7, units = "in")
  }
}


# Notes:
# Anaphalis margaritacea does not occur in the PNW according to occurrence data

#######################################################

