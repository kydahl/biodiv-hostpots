##############################################################################
###                                                                        ###
###                       Species presence and ecoregions                  ###
###                                                                        ###
##############################################################################

# Elisa Van Cleemput, April 2023
######################################################## Explore BIEN data

# Load packages
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

species_list <- BIEN_list_all()
dim(species_list)
head(species_list)

## ---- Load list of species --------------
# Check which of our species are the BIEN database

### OLD CODE, based on si dataframe, COMMENTED: can be deleted once we are finalizing the code
# load si dataframe from funct_div script
# library(dplyr)
# species_list_PNW <- species_list %>%
#   filter(species %in% unique(si$Species_full))
# dim(species_list_PNW) # 315
# length(unique(si$Species_full)) # 329
# # -> Not all species are in there. May be due to synonyms or full names (subsp.)
# 
# species_list_PNW_modif <- species_list %>%
#   filter(species %in% unique(si_modif$Species_full))
# dim(species_list_PNW_modif) # 321
# # A few more overlaps
# 
# # What happens if we don't use full species names?
# si_modif <- si_modif %>%
#   mutate(Species_short = paste(stringr::word(Species_full, 1,1, sep=" "),
#                           stringr::word(Species_full, 2,2, sep=" "),
#                           sep=" "))
# species_list_PNW_modif_short <- species_list %>%
#   filter(species %in% unique(si_modif$Species_short))
# dim(species_list_PNW_modif_short) # 325
# 
# species_list_PNW_modif_long <- si_modif %>%
#   filter(Species_short %in% unique(species_list$species)) %>%
#   select(Species_full) %>%
#   unique()
# dim(species_list_PNW_modif_long) # 327, this means that 2 full species names have the same short name
# # 2 missing species left
# species_list <- species_list_PNW_modif_long

final_data <- read_csv("data/clean/final_dataset.csv")
dim(final_data) # 318
length(unique(final_data$Species)) # 316, so 2 species have the same "base name"

final_data_in_BIEN <- final_data %>%
  filter(Species %in% unique(species_list$species)) %>%
  select(Species) %>%
  unique()
dim(final_data_in_BIEN) # 315, so 3 missing species

species_list <- final_data_in_BIEN # We will use this dataframe in part 3)

## ---- Load Ecoregions map --------------
# Ecoregions II and III maps donwloaded on April 14 2023
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-ii/
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-iii/

Ecoregions2 <- read_sf("data/raw/NA_Terrestrial_Ecoregions_v2_Level_II_Shapefile/NA_TerrestrialEcoregions_LII/data/NA_Terrestrial_Ecoregions_v2_level2.shp")
Ecoregions3 <- read_sf("data/raw/NA_Terrestrial_Ecoregions_v2_Level_III_Shapefile/NA_TerrestrialEcoregions_LIII/data/NA_Terrestrial_Ecoregions_v2_level3.shp")

# Explore metadata
# str(Ecoregions2)
# class(Ecoregions2)
# Ecoregions2$geometry # Projected CRS: Sphere_ARC_INFO_Lambert_Azimuthal_Equal_Area
# st_crs(Ecoregions2_PNW) 

# Quick visualization 
# plot(Ecoregions2["LEVEL2"])
# plot(st_geometry(Ecoregions2), axes=T)


#######################################################
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
# 3) Loop through all species and extract the ecoregion(s) the occur in
#######################################################
# This idea is based on Echeverría-Londoño et al. 2018: 
# "We overlaid the BIEN 2.0 plant species maps on a 100 × 100 km grid map with a 
# Lambert Azimuthal Equal Area projection to obtain a presence/absence matrix of species for each grid cell."

# function ot extract species occurence and PNW regions for each species
extract_species_occ <- function(species, Ecoregions3_PNW){
  
  print(paste0("Extracting information for species ", species))
  
  # Extract species occurrences from BIEN database
  species_occ <- BIEN::BIEN_occurrence_species(species,
                                         cultivated=F, natives.only=T)
  
  if (nrow(species_occ) == 0){
    print("The species is has no occurrence data")
    output <- list("species_occ_Lb" = NA,
                   "species_occ_eco3" = NA)
  } else {
    
    species_occ <- species_occ %>%
      drop_na(longitude) %>%
      drop_na(latitude)
    
    if (nrow(species_occ) == 0){
      print("The species is has no records for longitude and latitude")
      output <- list("species_occ_Lb" = NA,
                     "species_occ_eco3" = NA)
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
        st_transform(st_crs(Ecoregions3_PNW))
      
      # Overlap BIEN occurrence data with PNW map: Extract ecoregion(s) where species occurs 
      # This takes a while!
      # Slow option, using sf objects:
      # system.time (species_occ_eco2 <- st_join(species_occ_Lb, Ecoregions2_PNW, join = st_within))
      # Faster option, using sp objects
      # system.time (species_occ_eco2 <- over(as_Spatial(species_occ_Lb), as_Spatial(Ecoregions2_PNW)))
      system.time (species_occ_eco3 <- over(as_Spatial(species_occ_Lb), as_Spatial(Ecoregions3_PNW)))
      
      output <- list("species_occ_Lb" = species_occ_Lb,
                     "species_occ_eco3" = species_occ_eco3)
      
    }
    
  }
}
# Error issues that came up while running this function:
# - Error for Chamaecyparis nootkatensis, because the species has no ocurrence data --> error issue solved
# - Error for Pinus albicaulis, because of missing values in coordinates --> error issue solved
# - Error for Pinus contorta, 

# Extract species occurence and PNW regions for each species
# This takes a while!
system.time(species_occurences <- apply(species_list, 1, function(x) extract_species_occ(x, Ecoregions3_PNW)))
names(species_occurences) <- species_list$Species


# save

# Concatenate information from all species
# Create a dataframe showing in which ecoregions each species occurs

unique(species_occ_Lb$scrubbed_species_binomial)
unique(species_occ_eco3$LEVEL1)
unique(species_occ_eco3$LEVEL2)
# unique(species_occ_eco3$LEVEL3)

#######################################################
# 4) Visualize where species occur in the PNW
#######################################################

world_cropped <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE)) %>%
  st_transform(st_crs(Ecoregions2_PNW))%>%
  st_crop((Ecoregions2_PNW))


for (i in 1:nrow(species_list)){
 
  # Extract species occurrence:
  species_occ <- test(species_list[i,])
  
  # Visualize species occurrences
  g1 <- ggplot() +
    geom_sf(data=world_cropped) +
    geom_sf(data=Ecoregions2_PNW, aes(fill=LEVEL2)) +
    scale_fill_grafify(palette="all_grafify") + # long palettes: all_grafify, kelly, safe
    geom_sf(data=test2$`Abies amabilis`$species_occ_Lb, color="black") +
    theme_bw()
  g2 <- ggplot() +
    geom_sf(data=world_cropped) +
    geom_sf(data=Ecoregions3_PNW, aes(fill=LEVEL3)) +
    geom_sf(data=species_occ_sf, color="black") +
    theme_bw()
  
  # Compile ecoregion coverage in a dataframe
  #   
}



# started around 2:30 PM - end 2:40???
names(test2) <- species_list[1:2,]$Species

system.time(test3 <- lapply(
  split(species_list[1:2,],1:nrow(species_list[1:2,])),
  function(x) do.call(test, x)))

# Save this result and then proceede to mapping so that this data can always be visualized differently afterwards

result <- lapply(
  split(df,1:nrow(df)),
  function(x) do.call(power.t.test,x, Ecoregions3_PNW)
)








dir_fig

species_occ <- BIEN_occurrence_species(species=species_list[i,1],
                                       cultivated=F, natives.only=T)
# species_occ2 <- BIEN_occurrence_species(species=species_list_PNW_modif_long[i,1], all.taxonomy=T, new.world=T, native.status=T)
# str(species_occ)

# Quickly visualize points on a map
library(maps)
map('world', fill = TRUE, col= "grey", bg = "light blue") 
points(cbind(species_occ$longitude,
             species_occ$latitude),
       col = "red",
       pch = 20,
       cex = 1) 
# points(cbind(species_occ2$longitude,
#              species_occ2$latitude),
#        col = "blue",
#        pch = 20,
#        cex = 1) 

## ---- Overlap BIEN occurence with PNW map --------------

# I assume that the BIEN dataset is in WGS84 (EPSG: 4326) but I can't find documentation on that
# According to this code it is in WGS84 indeed:
# https://github.com/NiDEvA/R-protocols/blob/main/R_Protocol_get%26curate_data.R#L330
species_occ_sf <- species_occ %>% 
  st_as_sf(coords = c("longitude","latitude")) %>%
  st_set_crs(4326) %>%
  st_transform(st_crs(Ecoregions2_PNW))

# plot(Ecoregions2_PNW["LEVEL2"], axes=T)
# plot(species_occ_sf_Lb[1], pch=20, col="red", add=T)
# For some reason this visualization doesn't work. I don't know why

library(ggplot2)
library(maps)
library(grafify)
# plot_grafify_palette(palette="all_grafify")
world_cropped <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE)) %>%
  st_transform(st_crs(Ecoregions2_PNW))%>%
  st_crop((Ecoregions2_PNW))
ggplot() +
  geom_sf(data=world_cropped) +
  geom_sf(data=Ecoregions2_PNW, aes(fill=LEVEL2)) +
  scale_fill_grafify(palette="all_grafify") + # long palettes: all_grafify, kelly, safe
  geom_sf(data=species_occ_sf, color="black") +
  theme_bw()
ggplot() +
  geom_sf(data=world_cropped) +
  geom_sf(data=Ecoregions3_PNW, aes(fill=LEVEL3)) +
  geom_sf(data=species_occ_sf, color="black") + 
  theme_bw()

# Extract ecoregion(s) where species occurs 
# Slow option, using sf objects:
# system.time (species_occ_eco2 <- st_join(species_occ_sf, Ecoregions2_PNW, join = st_within))
# Fast option, using sp objects
system.time (species_occ_eco2 <- over(as_Spatial(species_occ_sf), as_Spatial(Ecoregions2_PNW)))
# system.time (species_occ_eco3 <- over(as_Spatial(species_occ_sf), as_Spatial(Ecoregions3_PNW)))

unique(species_occ_sf$scrubbed_species_binomial)
unique(species_occ_eco2$LEVEL1)
unique(species_occ_eco2$LEVEL2)
# unique(species_occ_eco3$LEVEL3)


#######################################################
## ---- --------------

