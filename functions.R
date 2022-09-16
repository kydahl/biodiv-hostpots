################################################################################
# Functions used in the simulation study
################################################################################

## Load libraries---------------------------------------------------------------
library(tidyverse)

## Helper functions-------------------------------------------------------------
# Sample from a discrete uniform distribution
# Source: https://stats.stackexchange.com/questions/3930/are-there-default-functions-for-discrete-uniform-distributions-in-r
dunif_sampleone <- function(n) sample(1:n, 1, replace = T)

# Assign entries from df to groups of size 1 to n (chosen uniformly randomly)
group_assign <- function(df, n) sample(df, size = dunif_sampleone(n), replace = FALSE)

## Functions for creating data frames-------------------------------------------
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


## Vulnerability function-------------------------------------------------------
# take traits, current state, and patch number and outputs future state

# simple function for now:
# if average of traits > 0.5, increase IUCN status
# if average of traits < -0.5, decrease IUCN status
# otherwise maintain status
get.vulnerability <- function(in_df, in_year) {
  out_df <- in_df %>%
    # select the correct year
    filter(year == in_year) %>%
    # group by all columns except traits
    group_by(across(c(-trait.num, -trait.val)), .drop = FALSE) %>%
    # compute average of traits
    mutate(trait_avg = mean(trait.val)) %>%
    # assign new states based on average of traits (in future, this will be a real function based on biology)
    mutate(state = state + case_when(
      state == 1 ~ ifelse(rnorm(1) > 3, 1, 0), # 2.326, 1, 0), # with small prob., extirpated species re-emerge
      trait_avg > 0.5 ~ 1,
      trait_avg < -0.5 ~ -1,
      TRUE ~ 0
    )) %>%
    # 6 is the highest IUCN level (Least concern), so don't go higher than 6
    mutate(state = ifelse(state > 6, 6, state)) %>%
    mutate(year = year + 1) %>%
    mutate(trait.val = case_when(
      trait.num == 1 ~ rnorm(1), # one trait changes randomly over time due to environmental changes, say
      TRUE ~ trait.val
    )) %>%
    select(-trait_avg) %>%
    ungroup()
}

## Biodiversity metric functions------------------------------------------------
# for a patch, look at entities and their states & traits and output
# biodiversity measure of the patch

# simple function for now:
# number of unique entityIDs in a patch
calc.biodiv_1 <- function(in_df) {
  out_df <- in_df %>%
    filter(state > 1) %>% # remove extinct entities
    # remove duplicates due to rows from trait.num
    select(patch, year, entityID) %>%
    distinct() %>%
    # count number of distinct entities
    count(patch, year, name = "biodiv1")
}

# Number of endemic entities
calc.biodiv_2 <- function(in_df) {
  out_df <- in_df %>%
    filter(state > 1) %>% # remove extinct entities
    select(patch, year, entityID) %>%
    # only keep rows where entityID and year are unique ( == endemics)
    group_by(entityID, year) %>%
    # remove duplicates due to rows from trait.num
    distinct() %>%
    filter(n() == 1) %>%
    ungroup() %>%
    # count number of endemics in each patch
    count(patch, year, name = "biodiv2")
}


## Hotspot identifier function--------------------------------------------------
# look at the biodiversity of all patches and determine which patches are
# "hotspots" relative to the other patches

# simple function for now:
# just choose the top patch in terms of biodiversity
find.hotspots <- function(in_df, biod_metric) {
  out_df <- in_df %>%
    select(patch, year, !!as.name(biod_metric)) %>%
    group_by(year) %>%
    slice_max(order_by = !!as.name(biod_metric), n = 1, with_ties = FALSE) %>%  # !!! this breaks ties arbitrarily
    ungroup()
}


## Hotspot comparison function--------------------------------------------------
# measure the difference in the hotspot maps produced under different
# biodiversity metrics for the same dataset

# simple function for now:
# count the number of overlaps in which patches were considered hotspots
calc.hotspot_compare <- function(df1, df2) {
  vec1 <- select(df1, patch)
  vec2 <- select(df2, patch)
  # get the Euclidean distance between the patch lists (!!! placeholder)
  norm(vec1-vec2, "2")
  
}

