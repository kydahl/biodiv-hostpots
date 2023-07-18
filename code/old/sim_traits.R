################################################################################
# Initial lists used for simulation study
################################################################################

# Load libraries----------------------------------------------------------------
library(tidyverse)

# Parameters--------------------------------------------------------------------

## Random seed (to keep things consistent while coding at first)
# set.seed(8797)

## Number of entities
numEntities <- 1000 # got to 10 million before code took a while to run

## Number of traits
# TODO automate number of traits (low priority as we'll probably be taking traits from a database)
numTraits <- 4

## Number of patches
numPatches <- 10 # TODO: automate number of patches

## Maximum number of entities per patch
maxEntitiesPatch <- numEntities / 10 # very arbitrary choice for now.

## IUCN categories
IUCN_cats <- c(
  "Extinct",
  "Critically endangered",
  "Endangered",
  "Vulnerable",
  "Near threatened",
  "Least concern"
)
state_list <- tibble(nums = 1:6, cats = IUCN_cats)

# Helper functions--------------------------------------------------------------

# Sample from a discrete uniform distribution
# Source: https://stats.stackexchange.com/questions/3930/are-there-default-functions-for-discrete-uniform-distributions-in-r
dunif_sampleone <- function(n) sample(1:n, 1, replace = T)

# Assign entries from df to groups of size 1 to n (chosen uniformly randomly)
group_assign <- function(df, n) sample(df, size = dunif_sampleone(n), replace = FALSE)

# Functions for creating tables-----------------------------------------------------------------

## Create entities and assign trait values
get.entities <- function(numEntities, numTraits) {
  tibble(expand.grid(
    entityID = 1:numEntities,
    trait.num = 1:numTraits
  ),
  trait.val = rnorm(numEntities * numTraits)
  )
}


## Assign entities to patches
get.patches <- function(entity_df, numPatches, maxEntitiesPatch) {
  # initialize data.frame
  regional_df <- tibble(
    entityID = NA,
    patch = NA
  )

  # sample from entities list to fill patches
  for (i in 1:numPatches) {
    regional_df <- add_row(regional_df, patch = i, entityID = group_assign(entity_df$entityID, maxEntitiesPatch))
  }

  # assign uniqueIDs to entities in patches (and remove NAs we added in initialization)
  regional_df <- regional_df %>%
    mutate(uniqueID = paste("E", entityID, "_P", patch, sep = "")) %>%
    filter(!is.na(entityID))
}

## Assign states to each of the entities in each patch
assign.states <- function(regional_df) {
  regional_df <- regional_df %>%
    # Assign states (except for "Extinct")
    # these are numbers with higher = better
    mutate(state = sample(2:6, dim(regional_df)[1], replace = TRUE)) %>%
    mutate(year = 0) # initial year
}

## Create full data.frame of entities and their traits in patches with assigned states
get.full_df <- function(numEntities, numTraits, numPatches, maxEntitiesPatch) {
  # create entities data.frame and assign trait values
  entity_df <- get.entities(numEntities, numTraits)
  # assign patches and states to entities
  states_df <- entity_df %>% 
    # assign entities to patches
    get.patches(numPatches, maxEntitiesPatch) %>% 
    # assign states to entities within patches
    assign.states(.) 
  
  # add trait values back in
  full_df <- inner_join(entity_df, states_df) %>%
    # re-order the data.frame for easier viewing
    relocate(entityID, patch, uniqueID)
}

