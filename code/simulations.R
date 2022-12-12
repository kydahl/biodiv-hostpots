#### Create simulations of patches and calculate biodiversity ##################

## Load libraries---------------------------------------------------------------
library(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(picante)

# Register CPU cores for parallel processing
registerDoParallel(cores = 12 - 2)

# Load in analysis functions
source("code/functions.R")
# Load in data
source("code/data_intake.R")

## Helper functions-------------------------------------------------------------


## Parameters-------------------------------------------------------------------
## Set random seed (for debuggin)
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

# Run simulations ---------------------------------------------------------

# Load in initial data
entity_df <- data_in

# Create a single simulation for visualizations
full_df <- full_df <- get.full_df(entity_df, 
                                  numPatches, 
                                  mean.NumEntities, 
                                  sd.NumEntities)


# Get tree data
tree <- read.tree(file = "data/clean/phylogenetic_tree.csv")

# Comparing the hotspots of each biodiversity metric against species richness

# Preallocate comparison data frame
# compare_df <- tibble(
#   iteration = numeric(),
#   endemic = numeric(),
#   indig.name = numeric(),
#   indig.lang = numeric(),
#   use = numeric(),
#   PD = numeric(),
#   SR = numeric(),
#   PD_unrooted = numeric()
# )

# Build data frame
compare_df <- foreach(
  j = 1:numIterations,
  int = icount(),
  .combine = "rbind",
  .packages = c("tidyverse", "reshape2", "picante")
) %dopar% {
  # Run a simulation
  full_df <- get.full_df(entity_df, 
                         numPatches, 
                         mean.NumEntities, 
                         sd.NumEntities)

  # Get hotspot comparison values
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
    PD = biodiv.compare_df$PD,
    SR = biodiv.compare_df$SR,
    PD_unrooted = biodiv.compare_df$PD_unrooted
  )
}

# !!! BUG: for some reason, foreach is adding repeated entries to the data frame
#          code below removes the repeats, but this could cause issues for
#          higher numbers of iterations
compare_df <- unique(compare_df)

