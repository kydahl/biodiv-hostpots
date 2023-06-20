#### Create simulations of patches and calculate biodiversity ##################

## Load libraries---------------------------------------------------------------
library(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(picante)

# Register CPU cores for parallel processing
numCores <- parallel::detectCores()
registerDoParallel(cores = numCores - 2)

# Load in analysis functions
source("code/functions.R")

# Load in trait data
final_data <- read_csv("data/clean/final_dataset.csv") #%>% 
  # rename(Species = Species_full)

# Load in phylogenetic tree data
# tree <- read.tree(file = "data/clean/phylogenetic_tree.csv")
# Elisa: I commented this out, because we have been playing around with species names and hence, the tree might be different now


## Helper functions ------------------------------------------------------------
# Uniformly sample from interval [1, MaxNumEntities]
get.NumEntities <- function(MaxNumEntities, Level) {
  if (!(Level %in% MaxNumEntities$Level)) {
    warning("Incorrect level choice")
    break
  }
  
  MaxNum <- filter(MaxNumEntities, Level == Level)$count
  round(runif(1, 1, MaxNum))
}

## Parameters ------------------------------------------------------------------
## Set random seed (for debugging)
# set.seed(8797)

## Number of iterations to compute biodiversity over
numIterations <- 1000

# Load species occurrences ------------------------------------------------

# Output .csvs from "exploring_BIEN.R"
species_occurences_L1 <- read_csv("data/clean/species_occurences_L1.csv")
species_occurences_L2 <- read_csv("data/clean/species_occurences_L2.csv")
species_occurences_L3 <- read_csv("data/clean/species_occurences_L3.csv")

# level 2 seems like the best bet, but try exploring level 3. see what feels best when doing simulations
SpeciesOccs <- species_occurences_L3 %>% 
  rename(Level = LEVEL3)

# Create patches ----------------------------------------------------------

# Define total number of patches
NumPatches <- 400 # Just for testing purposes at this point. Need to settle on scale later


# Assign levels to patches ------------------------------------------------

# Each patch has an equal probability of being any of the levels
Levels <- sample(unique(SpeciesOccs$Level), NumPatches, replace = T)

Init_df <- tibble(Patch = 1:NumPatches, Level = Levels)

# Change to repeat each level N times

# Assign number of species to each patch ------------------------------

# Maximum needs to be set to the total number of species at that level
MaxNumEntities <- SpeciesOccs %>% 
  group_by(Level) %>% 
  summarise(count = n())


# Assign species to patches -----------------------------------------------

Patch_df <- tibble(Patch = as.integer(), Level = as.double(), Species = as.character())

# For each patch, get its level
for (index_patch in Init_df$Patch) {
  # Get the level of the patch
  level <- filter(Init_df, Patch == index_patch)$Level
  
  # Get species list from the right level
  species_list <- filter(SpeciesOccs, Level == level)$Species
  
  # Get number of entities in patch
  entity_count <- min(get.NumEntities(MaxNumEntities, level), length(species_list))
  
  # Sample entities without replacement from species list
  entities <- sample(species_list, entity_count, replace = FALSE)
  
  temp_df <- tibble(Patch = index_patch, Level = level, Species = entities)
  
  Patch_df <- rbind(Patch_df, temp_df)
  
}

# Add species traits to the dataframe
final_data <- rename(final_data)

full_df <- right_join(Patch_df, final_data, by = "Species") %>%
  filter(!is.na(Patch), !is.na(Level))
  # Elisa: I added this to remove the patch with NA as a number, but it turns out there are more incomplete cases. 
# I didn't check the reason for this, but shouldn't this dataset be complete (because it has imputed trait data)
# Kyle: I think na.omit removes rows where any of the columns are NA. So this is removing all the species with an NA in the subspecies or variant column. There may be species with NA in the Patch or Level column. These are the ones that weren't assigned to any of the patches. I made it so these get removed.

# I added trait_names to the input of get.biodiv_df
trait_names = c("LDMC (g/g)","Nmass (mg/g)", "Woodiness", "Plant height (m)", "Leaf area (mm2)" )
biodiv_df <- get.biodiv_df(full_df, trait_names)

#################################
# Elisa: I used full_df to make phylo_div and funct_div work. what is the difference between full_df and test_df?
#################################

get.full_df <- function(NumPatches, LevelOrder) {
  SpeciesOccs <- if (LevelOrder == 1) {
    read_csv("data/clean/species_occurences_L1.csv", show_col_types = FALSE) %>% 
      rename(Level = LEVEL1)
  } else if (LevelOrder == 2) {
    read_csv("data/clean/species_occurences_L2.csv", show_col_types = FALSE) %>% 
      rename(Level = LEVEL2)
  } else if (LevelOrder == 3) {
    read_csv("data/clean/species_occurences_L3.csv", show_col_types = FALSE) %>% 
      rename(Level = LEVEL3)
  } 
  
  # Assign levels to patches 
  
  # Equal numbers across each level (for now)
  Levels <- sample(SpeciesOccs$Level, NumPatches, replace = TRUE)
  
  Init_df <- tibble(Patch = 1:NumPatches, Level = Levels)
  
  # Assign number of species to each patch 
  
  # Maximum needs to be set to the total number of species at that level
  MaxNumEntities <- SpeciesOccs %>% 
    group_by(Level) %>% 
    summarise(count = n())
  
  
  # Assign species to patches -----------------------------------------------
  
  Patch_df <- tibble(Patch = as.integer(), Level = as.double(), Species = as.character())
  
  # For each patch, get its level
  for (index_patch in Init_df$Patch) {
    # Get the level of the patch
    level <- filter(Init_df, Patch == index_patch)$Level
    
    # Get species list from the right level
    species_list <- filter(SpeciesOccs, Level == level)$Species
    
    # Get number of entities in patch
    entity_count <- min(get.NumEntities(MaxNumEntities, level), length(species_list))
    
    # Sample entities without replacement from species list
    entities <- sample(species_list, entity_count, replace = FALSE)
    
    temp_df <- tibble(Patch = index_patch, Level = level, Species = entities)
    
    Patch_df <- rbind(Patch_df, temp_df)
    
  }
  
  # Add species traits to the dataframe
  full_df <- right_join(Patch_df, final_data, by = "Species", relationship = "many-to-many") %>% 
    filter(!is.na(Patch), !is.na(Level))
}

# Compare hotspots identified by different metrics
test_df <- get.full_df(NumPatches = 2, LevelOrder = 3)

biodiv_df <- get.biodiv_df(test_df, trait_names)

compare_biodiv_hotspots_df <- get.biodiv.compare_df(test_df)


compare_df <- tibble(iteration = as.integer(),
                      endemic = as.double(),
                      indig.name = as.double(),
                      indig.lang = as.double(),
                      use = as.double(),
                      CV = as.double(),
                      # SR = biodiv.compare_df$SR,
                      # PD_unrooted = biodiv.compare_df$PD_unrooted
)

# Build data frame
compare_df <- foreach(
  j = 1:numIterations,
  int = icount(),
  .combine = "rbind",
  .packages = c("tidyverse", "reshape2", "picante")
) %dopar% {
  # Run a simulation
  full_df <- get.full_df(400, 2)

  # # Get hotspot comparison values
  biodiv.compare_df <- get.biodiv.compare_df(full_df) %>%
    pivot_wider(names_from = variable) %>%
    unique()

  # Add to the list
  compare_df <- add_row(compare_df,
                        iteration = int,
                        endemic = biodiv.compare_df$NumEndemic,
                        indig.name = biodiv.compare_df$NumIndigName,
                        indig.lang = biodiv.compare_df$NumIndigLang,
                        use = biodiv.compare_df$NumUse,
                        CV = biodiv.compare_df$CoeffVar
                        # SR = biodiv.compare_df$SR,
                        # PD_unrooted = biodiv.compare_df$PD_unrooted
  )
  
}


