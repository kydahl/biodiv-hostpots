##############################################################################
###                                                                        ###
###                       Species presence and ecoregions                  ###
###                                                                        ###
##############################################################################

# Elisa Van Cleemput, April 2023
######################################################## Explore BIEN data

if (!require('BIEN')) install.packages('BIEN'); library('BIEN')
# vignette("BIEN")
# vignette("BIEN_tutorial")


#######################################################
## ---- Exploring BIEN dataset --------------

species_list <- BIEN_list_all()
head(species_list)
dim(species_list)

# Check which of our species are the BIEN database
# load si dataframe from funct_div script
library(dplyr)
species_list_PNW <- species_list %>%
  filter(species %in% unique(si$Species_full))
dim(species_list_PNW) # 315
length(unique(si$Species_full)) # 329
# -> Not all species are in there. May be due to synonyms or full names (subsp.)

species_list_PNW_modif <- species_list %>%
  filter(species %in% unique(si_modif$Species_full))
dim(species_list_PNW_modif) # 321
# A few more overlaps

# What happens if we don't use full species names?
si_modif <- si_modif %>%
  mutate(Species_short = paste(stringr::word(Species_full, 1,1, sep=" "),
                          stringr::word(Species_full, 2,2, sep=" "),
                          sep=" "))
species_list_PNW_modif_short <- species_list %>%
  filter(species %in% unique(si_modif$Species_short))
dim(species_list_PNW_modif_short) # 325

species_list_PNW_modif_long <- si_modif %>%
  filter(Species_short %in% unique(species_list$species)) %>%
  select(Species_full) %>%
  unique()
dim(species_list_PNW_modif_long) # 327, this means that 2 full species names have the same short name

# 2 missing species left


#######################################################
## ---- Load Ecoregions map --------------
# Ecoregions II and III maps donwloaded on April 14 2023
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-ii/
# http://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-iii/
library(sf)

Ecoregions2 <- read_sf("data/raw/NA_Terrestrial_Ecoregions_v2_Level_II_Shapefile/NA_TerrestrialEcoregions_LII/data/NA_Terrestrial_Ecoregions_v2_level2.shp")
Ecoregions3 <- read_sf("data/raw/NA_Terrestrial_Ecoregions_v2_Level_III_Shapefile/NA_TerrestrialEcoregions_LIII/data/NA_Terrestrial_Ecoregions_v2_level3.shp")

str(Ecoregions2)
class(Ecoregions2)
Ecoregions2$geometry # Projected CRS: Sphere_ARC_INFO_Lambert_Azimuthal_Equal_Area
st_crs(Ecoregions2_PNW) 


plot(Ecoregions2["LEVEL2"])
plot(st_geometry(Ecoregions2), axes=T)

# According to Lynette, the following Ecoregions consitute the PNW: 7.1, 6.1,6.2,10.1,3.2,3.1,3.3,2.2,2.3,5.4 (meeting notes April 5 2023)
Ecoregions2_PNW <- Ecoregions2 %>%
  filter(LEVEL2 %in% c("2.2", "2.3", "3.1", "3.2", "3.3", "5.4","6.1", "6.2", "7.1", "10.1"))
# plot(Ecoregions2_PNW["LEVEL2"])


Ecoregions3_PNW <- Ecoregions3 %>%
  filter(LEVEL2 %in% c("2.2", "2.3", "3.1", "3.2", "3.3", "5.4","6.1", "6.2", "7.1", "10.1"))
# plot(Ecoregions3_PNW["LEVEL3"])


## ---- Exploring BIEN occurence --------------
i <- 1

species_occ <- BIEN_occurrence_species(species=species_list_PNW_modif_long[i,1],
                                       cultivated=F, natives.only=T)
# species_occ2 <- BIEN_occurrence_species(species=species_list_PNW_modif_long[i,1], all.taxonomy=T, new.world=T, native.status=T)
# str(species_occ)

# Suickly visualize points on a map
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
# Echeverría-Londoño et al. 2018: "We overlaid the BIEN 2.0 plant species maps on a 100 × 100 km grid map with a 
# Lambert Azimuthal Equal Area projection to obtain a presence/absence matrix of species for each grid cell."

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

