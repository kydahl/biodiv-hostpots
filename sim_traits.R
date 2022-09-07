# Load libraries----------------------------------------------------------------
library(tidyverse)

# Helper functions--------------------------------------------------------------

# Sample from a discrete uniform distribution
# Source: https://stats.stackexchange.com/questions/3930/are-there-default-functions-for-discrete-uniform-distributions-in-r
dunif_sampleone <- function(n) sample(1:n, n, replace = T)

# Assign entries from df to groups of size 1 to n (chosen uniformly randomly)
group_assign <- function(df, n) sort(sample(df, dunif_sampleone(n), replace = FALSE))

# Create tables-----------------------------------------------------------------
numEntities <- 100

## Define all the entities to be assigned to patches
entity_df <- tibble(
  # Assign traits (for now this is random) will be automated in future update
  trait1 = rnorm(numEntities),
  trait2 = rnorm(numEntities),
  trait3 = rnorm(numEntities),
  trait4 = rnorm(numEntities)
) %>%
  # Assign entityID (for now this is row number)
  mutate(entityID = row_number())

## Assign entities to patches (5 for now with a max of 25 entities in each)
# will be automated in future update
# uses a discrete uniform distribution to assign entities to patches for now
maxEntitiesPatch <- 25

regional_df <- tibble(
  entityID = c(),
  patch = c()
) %>%
  add_row(patch = 1, entityID = group_assign(entity_df$entityID, maxEntitiesPatch)) %>%
  add_row(patch = 2, entityID = group_assign(entity_df$entityID, maxEntitiesPatch)) %>%
  add_row(patch = 3, entityID = group_assign(entity_df$entityID, maxEntitiesPatch)) %>%
  add_row(patch = 4, entityID = group_assign(entity_df$entityID, maxEntitiesPatch)) %>%
  add_row(patch = 5, entityID = group_assign(entity_df$entityID, maxEntitiesPatch)) %>%
  mutate(uniqueID = paste("E", entityID, "_P", patch, sep = ""))

## Expand the global entity data frame to include patch assignments
full_df <- inner_join(entity_df, regional_df) %>%
  relocate(entityID, patch, uniqueID)

## Assign states to entities in patches
# done randomly for now
IUCN_cats <- c("Extinct", "Critically endangered", "Endangered", "Vulnerable", 
               "Near threatened", "Least concern")

full_df <- full_df %>%
  # Assign states (except for "Extinct")
  mutate(state = sample(IUCN_cats[2:6], dim(full_df)[1], replace = TRUE))