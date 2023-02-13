#### Parameters for simulation study ####

# Kyle Dahlin, Jan 2023

## Parameters-------------------------------------------------------------------
## Set random seed (for debugging)
# set.seed(8797)

## Number of patches
numPatches <- 400 # !!! Decide on a realistic value for the scale we care about

## Number of entities per patch
mean.NumEntities <- 20
CV.NumEntities <- 1 / 2
sd.NumEntities <- mean.NumEntities * CV.NumEntities


# patchEntities_denominator <- 50
# maxEntitiesPatch <- ceiling(numEntities / patchEntities_denominator)
# Notes / Observations:
# * Species richness and num. endemics become more similar as
#   patchEntities_denominator is increased.
# * This is because if there are fewer entities/patch and the entities come from
#   the same (relatively large) pool of entities, then it is very likely that 
#   any given entity in a patch is endemic there.
# * Be careful when making patchEntities_denominator small because that can
#   substantially increase the size of 'full_df' and thus the overall
#   computation time for the simulations below

## IUCN categories
IUCN_cats <- c(
  "Extinct",
  "Critically endangered",
  "Endangered",
  "Vulnerable",
  "Near threatened",
  "Least concern"
)
# in case we want to convert numbers to categories
state_keys <- tibble(nums = 1:6, cats = IUCN_cats) 

## Number of iterations to compute biodiversity over
numIterations <- 1000