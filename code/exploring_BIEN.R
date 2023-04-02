# Explore BIEN data

if (!require('BIEN')) install.packages('BIEN'); library('BIEN')
vignette("BIEN")
vignette("BIEN_tutorial")

species_list <- BIEN_list_all()
head(species_list)
dim(species_list)

# Check which of our species are the BIEN database
# si dataframe from funct_div script
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