################################################################################
# Initial lists used for simulation study
################################################################################

# Load libraries----------------------------------------------------------------
library(tidyverse)

# Parameters--------------------------------------------------------------------

## Random seed (to keep things consistent while coding at first)
#set.seed(8797)

## Number of entities
numEntities <- 1000 # got to 10 million before code took a while to run

## Number of traits
# TODO automate number of traits (low priority as we'll probably be taking traits from a database)
numTraits <- 4

## Number of patches
numPatches <- 10 # TODO: automate number of patches

## Maximum number of entities per patch
maxEntitiesPatch <- numEntities/2 # very arbitrary choice for now.

## IUCN categories
IUCN_cats <- c("Extinct",
               "Critically endangered",
               "Endangered", 
               "Vulnerable", 
               "Near threatened", 
               "Least concern")
state_list <- tibble(nums = 1:6, cats = IUCN_cats)

# Helper functions--------------------------------------------------------------

# Sample from a discrete uniform distribution
# Source: https://stats.stackexchange.com/questions/3930/are-there-default-functions-for-discrete-uniform-distributions-in-r
dunif_sampleone <- function(n) sample(1:n, 1, replace = T)

# Assign entries from df to groups of size 1 to n (chosen uniformly randomly)
group_assign <- function(df, n) sample(df, size = dunif_sampleone(n), replace = FALSE)

# Create tables-----------------------------------------------------------------

## Define all the entities to be assigned to patches
# Assign entity IDs and trait values
entity_df <- tibble(expand.grid(entityID = 1:numEntities, 
                         trait.num = 1:numTraits), 
                    trait.val = rnorm(numEntities * numTraits))

# # Assign traits ( these are currently IID normal(0, 1) )
# for (i in 1:numTraits) {
#   varname <- paste("trait", i , sep=".")
#   entity_df[[varname]] <- rnorm(numEntities)
# }


## Assign entities to patches (5 for now with a max of 25 entities in each)
# will be automated in future update
# uses a discrete uniform distribution to assign entities to patches for now
regional_df <- tibble(
  entityID = NA,
  patch = NA) 

for (i in 1:numPatches) {
  regional_df <- add_row(regional_df, patch = i, entityID = group_assign(entity_df$entityID, maxEntitiesPatch))
}

regional_df <- regional_df %>%
  mutate(uniqueID = paste("E", entityID, "_P", patch, sep = "")) %>% 
  filter(!is.na(entityID))

## Assign states to entities in patches
# done randomly for now

regional_df <- regional_df %>%
  # Assign states (except for "Extinct")
  # these are numbers with higher = better
  mutate(state = sample(2:6, dim(regional_df)[1], replace = TRUE)) %>% 
  mutate(year = 0) # initial year


## Expand the global entity data frame to include patch assignments
full_df <- inner_join(entity_df, regional_df) %>%
  relocate(entityID, patch, uniqueID)




