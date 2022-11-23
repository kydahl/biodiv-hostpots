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
group_assign <- function(in_df, 
                         mean.NumEntities,
                         sd.NumEntities) {
  Num.Entities <- min(max(
    floor(rnorm(1, mean = mean.NumEntities, sd.NumEntities)),
    1),
    length(in_df))
  
  sample(in_df, size = Num.Entities, replace = FALSE)
}

## Functions for creating data frames-------------------------------------------
## Create entities and assign trait values
get.entities <- function(numEntities, numTraits) {
  tibble(
    expand.grid(
      entityID = 1:numEntities,
      trait.num = 1:numTraits
    ),
    trait.val = rnorm(numEntities * numTraits)
  )
}


## Assign entities to patches
get.patches <- function(entity_df,
                        numPatches, 
                        mean.NumEntities,
                        sd.NumEntities) {
  # initialize data.frame
  regional_df <- tibble(
    entityID = NA,
    patch = NA
  )
  
  # sample from entities list to fill patches
  for (i in 1:numPatches) {
    regional_df <- add_row(regional_df,
      patch = i,
      entityID = group_assign(
        unique(entity_df$entityID),
        mean.NumEntities,
        sd.NumEntities
      )
    )
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
    mutate(state = case_when(
      Status == "Increasing (Least Concern)" ~ 7,
      Status == "Stable (Least Concern)" ~ 6,
      Status == "Unknown (Least Concern)" ~ 5,
      Status == "Decreasing (Least Concern)" ~ 4,
      Status == "Decreasing (Near threatened)" ~ 3,
      Status == "Decreasing (Endangered)" ~ 2,
      Status == "/" ~ sample(2:6, 1) # !!! placeholder. just assign something random if it's blank...
    ))

  # mutate(state = sample(2:6, dim(regional_df)[1], replace = TRUE)) %>%
}

## Create full data.frame of entities and their traits in patches with assigned states
get.full_df <- function(entity_df, numTraits, numPatches, mean.NumEntities, sd.NumEntities) {
  # create entities data.frame and assign trait values
  # entity_df <- get.entities(numEntities, numTraits)
  # assign patches and states to entities
  regional_df <- get.patches(entity_df,
                             numPatches, 
                             mean.NumEntities,
                             sd.NumEntities) 

  states_df <- regional_df %>% 
    
  states_df <- full_join(entity_df, regional_df, by = "entityID") %>%
    filter(!is.na(patch))

  full_df <- states_df %>%
    # assign states to entities within patches%>%
    # Assign states (except for "Extinct")
    # these are numbers with higher = better
    mutate(state = ifelse(test = (Status == "/"),
      # !!! place holder. if we don't know status, just assign randomly
      yes = sample(
        x = 2:7,
        size = dim(filter(states_df, Status == "/"))[1],
        replace = TRUE
      ),
      no = case_when(
        Status == "Increasing (Least Concern)" ~ 7,
        Status == "Stable (Least Concern)" ~ 6,
        Status == "Unknown (Least Concern)" ~ 5,
        Status == "Decreasing (Least Concern)" ~ 4,
        Status == "Decreasing (Near threatened)" ~ 3,
        Status == "Decreasing (Endangered)" ~ 2,
        # Status == "/" ~ sample.int(6, 1, replace = TRUE) # !!! placeholder. just assign something random if it's blank...
      )
    )) %>%
    # re-order the data.frame for easier viewing
    relocate(entityID, patch, uniqueID)

  # states_df <- entity_df %>%
  #   # assign entities to patches
  #   get.patches(numPatches, maxEntitiesPatch) %>%
  #   # assign states to entities within patches
  #   assign.states(.)

  # add trait values back in
  # full_df <- inner_join(entity_df, states_df) %>%
  #   # re-order the data.frame for easier viewing
  #   relocate(entityID, patch, uniqueID)
}


## Vulnerability function-------------------------------------------------------
# take traits, current state, and patch number and outputs future state

# simple function for now:
# if average of traits > 0.5, increase IUCN status
# if average of traits < -0.5, decrease IUCN status
# otherwise maintain status
get.vulnerability <- function(in_df) {
  out_df <- in_df %>%
    # group by all columns except traits
    group_by(across(c(-Trait, -value)), .drop = FALSE) %>%
    # compute average of traits
    # mutate(trait_avg = mean(value)) %>%
    # assign new states based on average of traits (in future, this will be a real function based on biology)
    mutate(state = state + case_when(
      state == 1 ~ ifelse(rnorm(1) > 3, 1, 0), # 2.326, 1, 0), # with small prob., extirpated species re-emerge
      trait_avg > 0.5 ~ 1,
      trait_avg < -0.5 ~ -1,
      TRUE ~ 0
    )) %>%
    # 6 is the highest IUCN level (Least concern), so don't go higher than 6
    mutate(state = ifelse(state > 7, 7, state)) %>%
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

# number of unique entityIDs in a patch
num.unique.metric <- function(in_df) {
  out_df <- in_df %>%
    filter(state > 1) %>% # remove extinct entities
    # remove duplicates due to rows from trait.num
    select(patch, entityID) %>%
    distinct() %>%
    # count number of distinct entities
    count(patch, name = "biodiv")
}

# Number of endemic entities
num.endemic.metric <- function(in_df) {
  out_df <- in_df %>%
    filter(state > 1) %>% # remove extinct entities
    select(patch, entityID) %>%
    # only keep rows where entityID is unique ( == endemics)
    group_by(entityID) %>%
    # remove duplicates due to rows from trait.num
    distinct() %>%
    mutate(endemic_flag = ifelse(n() == 1, 1, 0)) %>% 
    ungroup() %>% 
    group_by(patch) %>% 
    summarise("biodiv" = sum(endemic_flag, na.rm = TRUE),
              .groups = "keep")
    # filter(n() == 1) %>%
    # ungroup() %>%
    # # count number of endemics in each patch
    # count(patch, name = "biodiv")
}

trait.count.metric <- function(in_df, trait_name) {
  out_df <- in_df %>%
    # remove extinct entities
    filter(state > 1) %>%
    filter(Trait == trait_name) %>%
    select(patch, entityID, value) %>%
    group_by(patch) %>%
    summarise("biodiv" = sum(value, na.rm = TRUE),
              .groups = "keep")
}


## Hotspot identifier function--------------------------------------------------
# look at the biodiversity of all patches and determine which patches are
# "hotspots" relative to the other patches

# simple function for now:
# give all the patches in the 95% quantile
find.hotspots <- function(in_df) {
  
  # Select top 5th percentile
  in_df[in_df$biodiv > quantile(in_df$biodiv,0.95), ]
}


## Hotspot comparison function--------------------------------------------------
# measure the difference in the hotspot maps produced under different
# biodiversity metrics for the same dataset

# simple function for now:
# count the number of overlaps in which patches were considered hotspots
calc.hotspot_compare <- function(hotspots.unique, hotspots.compare) {
  unique.hotspots <- hotspots.unique$patch
  compare.hotspots <- hotspots.compare$patch
  
  TP_count <- sum(compare.hotspots %in% unique.hotspots)
  
  # false_count <- 
  
  # Use precision as our quantifier:
  # of the identified hotspots, what proportion match with species diversity?
  comparison.quantifier <- TP_count / length(compare.hotspots)
  return(comparison.quantifier)
}
