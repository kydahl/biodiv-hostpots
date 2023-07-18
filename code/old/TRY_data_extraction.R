################################################################################
# Processing data from the TRY database
################################################################################

## Project: Identifying Climate change vulnerable biodiversity hotspots
##
## Purpose: Process data to be used to inform simulations
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Load in trait data and put it in workable form
##           3) Set up parameters for the imputation function
##           4) Run missForest algorithm to impute missing traits
##           5) Perform diagnostics on imputed data
##           6) Illustrate diagnostics to ensure imputation was appropriate
##           7) Output imputed data frame
##
##
## Inputs:  data/clean/Trait_data_TRY_Diaz_2022_PNW.xlsx
##          - trait data from the TRY database used in Diaz 2022
##          
##
##          code -
##
##
## Outputs: data/clean/imputed_traits.csv - full list of imputed traits
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023
## _____________________________________________________________________________
#see next steps in the bottom of this script

#load the libraries
library(dplyr)

#load Diaz data 
Trait_Diaz <- read.csv("Trait_data_TRY_Diaz_2022_PNW.csv")

#load TRY data 
   #TRY data is separated by trait (these data is in our Google Drive folder in databases)
#Note that csv files starting with "TRY" followed by a trait name are a subset of the data provided for all TRY species present in their database. Only species present in "Trait_data_TRY_Diaz_2022_PNW.csv" were kept, in addition to possible synonmys found in funct_div.R 
#the full TRY data with all species is called TRY.zip

TRY_plant_height <- read.csv("TRY_plant_height_growth.csv")
TRY_leafN <- read.csv("TRY_leafN.csv") 
TRY_LDMC <- read.csv("TRY_LDMC.csv") 
TRY_leafarea <- read.csv("TRY_leafarea.csv") #wait to handle this data because we need to know what Diaz used for leaf area mesuarements
TRY_rooting_depth <- read.csv("TRY_rooting_depth.csv") #this trait was recommended for inclusion in the recommended traits paper
TRY_stem_density <- read.csv("TRY_stem_density.csv")  #this trait was recommended for inclusion in the recommended traits paper
TRY_plant_woodiness <- read.csv("TRY_plant_woodiness.csv")  #this trait was recommended for inclusion in the recommended traits paper
TRY_plant_veg_reproduction <- read_csv("TRY_vegetative_reproduction.csv")  #this trait was recommended for inclusion in the recommended traits paper

#---------------=== P   A   R   T    1 ===---------------------------------------
# C L E A N   U P   TRY  D A T A 
######---plant height ---######
TRY_plant_height <- read.csv("TRY_plant_height_growth.csv")
head(TRY_plant_height, 10)
unique(TRY_plant_height$TraitName) #look at the trait names that are under TraitName to make sure we are extracting data from the correct trait 
#here we have two types of traits "Plant height vegetative" and "Plant growth rate relative (plant relative growth rate, RGR)", but extract data only for "Plant height vegetative"

TRY_plant_height_2 <-TRY_plant_height %>% 
  select(AccSpeciesName, TraitName, ValueKindName, StdValue, UnitName) %>%
  filter(ValueKindName %in% c("Single", "Mean") & TraitName %in% "Plant height vegetative") %>% #Look at column "ValueKindName" keep either mean or single types
  group_by(AccSpeciesName, TraitName, UnitName) %>% 
  summarise(StdValue_mean = mean (StdValue)) #take the average of values for each species

######---leaf N ---######
TRY_leafN <- read.csv("TRY_leafN.csv")
head(TRY_leafN, 10)
unique(TRY_leafN$TraitName) #look at the trait names that are under TraitName to make sure we are extracting data from the correct trait 
#[1] "Leaf nitrogen (N) content per leaf dry mass" "Leaf nitrogen (N) content per leaf area"   #not sure which to keep, maybe the latter? 

TRY_leafN_2 <-TRY_leafN %>% 
  select(AccSpeciesName, TraitName, ValueKindName, StdValue, UnitName) %>%
  filter(ValueKindName %in% c("Single", "Mean") & TraitName %in% "Leaf nitrogen (N) content per leaf area") %>% #Look at column "ValueKindName" keep either mean or single types
  group_by(AccSpeciesName, TraitName, UnitName) %>% 
  summarise(StdValue_mean = mean (StdValue)) #take the average of values for each species

#---------------=== P   A   R   T    2 ===---------------------------------------

#---THESE TRAITS WERE NOT IN DIAZ DATA AND ARE ADDED HERE---######
######---rooting depth ---######
TRY_rooting_depth <- read.csv("TRY_rooting_depth.csv")
head(TRY_rooting_depth, 10)
unique(TRY_rooting_depth$TraitName) #look at the trait names that are under TraitName to make sure we are extracting data from the correct trait 

#columns to keep: "AccSpeciesName", "TraitName", "ValueKindName", "StdValue", "UnitName"
TRY_rooting_depth_2 <-TRY_rooting_depth %>% 
  select(AccSpeciesName, TraitName, ValueKindName, StdValue, UnitName) %>%
  filter(ValueKindName %in% c("Single", "Mean") & TraitName %in% "Root rooting depth") %>% #Look at column "ValueKindName" keep either mean or single types
  group_by(AccSpeciesName, TraitName, UnitName) %>% 
  summarise(StdValue_mean = mean (StdValue)) #take the average of values for each species

######---stem density ---######
TRY_stem_density <- read.csv("TRY_stem_density.csv")
head(TRY_stem_density, 10)
unique(TRY_stem_density$TraitName) #look at the trait names that are under TraitName to make sure we are extracting data from the correct trait 

#columns to keep: "AccSpeciesName", "TraitName", "ValueKindName", "StdValue", "UnitName"
TRY_stem_density_2 <-TRY_stem_density %>% 
  select(AccSpeciesName, TraitName, ValueKindName, StdValue, UnitName) %>%
  filter(ValueKindName %in% c("Single", "Mean") & TraitName %in% "Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density") %>% #Look at column "ValueKindName" keep either mean or single types
  group_by(AccSpeciesName, TraitName, UnitName) %>% 
  summarise(StdValue_mean = mean (StdValue)) #take the average of values for each species

######---TRY_vegetative_reproduction ---###### : this one has so few data that it is not even worth it adding to the master file

######---plant woodiness ---###### : we probably don't need this one either
TRY_plant_woodiness <- read.csv("TRY_plant_woodiness.csv")
head(TRY_plant_woodiness, 10)
unique(TRY_plant_woodiness$TraitName) #look at the trait names that are under TraitName to make sure we are extracting data from the correct trait 
unique(TRY_plant_woodiness$OriglName)

#columns to keep: "AccSpeciesName", "TraitName", "ValueKindName", "StdValue", "UnitName"
TRY_plant_woodiness_2 <-TRY_plant_woodiness %>% 
  select(AccSpeciesName, TraitName, OriglName ) %>% #we don't need "UnitName" and "ValueKindName" here because traits are categorical
  filter(OriglName %in% c("non-woody", "woody") & TraitName %in% "Plant woodiness") %>% #Look at column "ValueKindName" keep either mean or single types
  group_by(AccSpeciesName, TraitName) %>% 
  summarise(OriglName = OriglName) #take the average of values for each species
###########################################################################################################

##### Next action steps are:
#1) add the extra traits that are recommended (according to the paper) and keep their units to the Trait_Diaz dataframe
#we will need to decide if we need all these extra traits because for some we have very few data for species

#2) fill in the NAs in Trait_Diaz dataframe using the summarized TRY data

#3) we are left to deal with LMA and SLA traits once we know which trait data type to use (with petiole or without or both)





