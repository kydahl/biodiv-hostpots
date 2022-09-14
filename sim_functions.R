################################################################################
# Functions used in the simulation study
################################################################################

## Load libraries---------------------------------------------------------------
library(tidyverse)

## Helper functions-------------------------------------------------------------

## Vulnerability function-------------------------------------------------------
# take traits, current state, and patch number and outputs future state

# simple function for now: 
# if average of traits > 0.5, increase IUCN status
# if average of traits < -0.5, decrease IUCN status
# otherwise maintain status
vulner_func <- function(in_df, in_year) {
  out_df <- in_df %>%
    filter(year == in_year) %>% 
    rowwise() %>% 
    # compute average of traits
    # TODO replace to use starts_with("trait") somehow
    mutate(trait_avg = mean(c(trait1, trait2, trait3, trait4))) %>%
    mutate(trait1 = rnorm(1)) %>% # one trait changes over time due to environmental changes, say
    # assign new states
    mutate(state = state + case_when(
      state == 1       ~ ifelse(rnorm(1) > 3, 1, 0), #2.326, 1, 0), # with small prob., extirpated species re-emerge
      trait_avg > 0.5  ~ 1,
      trait_avg < -0.5 ~ -1,
      TRUE             ~ 0)
      ) %>% 
    mutate(state = ifelse(state > 6, 6, state)) %>% # 6 is the highest IUCN level (Least concern)
    mutate(year = year + 1) %>% 
    select(-trait_avg) %>% 
    ungroup()
}

## Biodiversity metric functions------------------------------------------------
# for a patch, look at entities and their states & traits and output 
# biodiversity measure of the patch

# simple function for now:
# number of unique entityIDs in a patch
biodiv_1 <- function(in_df) {
  out_df <- in_df %>%
    filter(state > 1) %>% # remove extinct entities
    select(patch, year, entityID) %>% 
    count(patch, year)
}

# Number of endemic entities
biodiv_2 <- function(in_df) {
  out_df <- in_df %>% 
    filter(state > 1) %>% # remove extinct entities
    select(patch, year, entityID) %>%  
    # only keep rows where entityID and year are unique ( == endemics)
    group_by(entityID, year) %>% 
    filter(n() == 1) %>% 
    ungroup() %>% 
    # count number of endemics in each patch
    count(patch, year)
}


## Hotspot identifier function--------------------------------------------------
# look at the biodiversity of all patches and determine which patches are 
# "hotspots" relative to the other patches

# simple function for now:
# just choose the top 5 patches in terms of biodiversity


## Hotspot comparison function--------------------------------------------------
# measure the difference in the hotspot maps produced under different 
# biodiversity metrics for the same dataset

# simple function for now:
# measure the overlap in which patches were considered hotspots



