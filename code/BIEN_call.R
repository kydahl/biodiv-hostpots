# BIEN database call

library("BIEN")
library("RPostgreSQL")
library("DBI")
library("ape")
library("maps")
library("sp")

#Set your directory

# setwd("~/Documents/Postdoc_Oak_Morton/Research_2022/")



#upload your Quercus species list
Species <- read.csv("Data/SpeciesList_BioHots.csv", sep = ",", # !!! change to the right fill
                    header = T)


head(Species)
dim(Species)

# To use the tutorial 

vignette("BIEN_tutorial")

Species_traits <- BIEN_trait_species (Species$Species, all.taxonomy = T, political.boundaries = T, source.citation = T)


unique(Species_traits$name_matched)
write.csv (Species_traits, "Species_traits.csv", sep = ",", col.names = TRUE)
