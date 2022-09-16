################################################################################
# Simulation trajectories over the time horizon
################################################################################

## Load libraries---------------------------------------------------------------
library(tidyverse)
library(doParallel)
registerDoParallel(cores = 7)
source("functions.R") # !!! combine these two later

## Helper functions-------------------------------------------------------------


## Parameters-------------------------------------------------------------------
## Random seed (to keep things consistent while coding at first)
# set.seed(8797)

## Number of entities
numEntities <- 1000

## Number of traits
numTraits <- 3

## Number of patches
numPatches <- 10

## Maximum number of entities per patch
patchEntities_denominator <- 2
maxEntitiesPatch <- numEntities / patchEntities_denominator # very arbitrary choice for now.
# Notes / Observations:
# * the difference between species richness and num. endemics decreases as 
#   patchEntities_denominator is increased. 
# * This is because if there are fewer entities/patch and the entities come from
#   the same pool of entities, then it is very likely that any given entity in a
#   patch is endemic there.
#   

## IUCN categories
IUCN_cats <- c(
  "Extinct",
  "Critically endangered",
  "Endangered",
  "Vulnerable",
  "Near threatened",
  "Least concern"
)
state_keys <- tibble(nums = 1:6, cats = IUCN_cats) # in case we want to convert numbers to categorgies

## Number of years to simulate
simYears <- 20

## Number of iterations to compute biodiversity over
numIterations <- 100

## Illustrative figures---------------------------------------------------------

### FIGURE 1 ###
### Comparing two biodiversity measures for a simulated set of patches over some
### number of years
full_df <- get.full_df(numEntities, numTraits, numPatches, maxEntitiesPatch)

# iterate over years to obtain states over time
for (i in 1:simYears) {
  full_df <- add_row(full_df, get.vulnerability(full_df, i - 1))
}

# compute biodiversity at each time point for each biodiversity metric
biod_df <- calc.biodiv_1(full_df) %>%
  full_join(calc.biodiv_2(full_df)) %>%
  # !!! placeholder for now, need to change calc.biodiv functions to just output this natively
  pivot_longer(
    cols = starts_with("biodiv"), names_to = "biodiv_metric",
    names_prefix = "biodiv", values_to = "biodiv_val"
  ) %>%
  mutate(biodiv_metric = case_when(
    biodiv_metric == 1 ~ "species_richness",
    biodiv_metric == 2 ~ "num_endemics"
  ))

# quick plot of biodiversity over time by patch
biod_plot <- biod_df %>%
  group_by(patch) %>%
  ggplot(aes(x = year, y = biodiv_val, color = as.factor(patch)), linetype = biodiv_metric) +
  geom_line(lwd = 1) +
  geom_point() +
  scale_y_continuous(name = "Biodiversity") +
  scale_color_discrete(name = "Patch") +
  # faceting:
  facet_wrap(~biodiv_metric, nrow = 1, scales = "free", labeller = labeller(variable = label_parsed)) +
  theme_minimal()

biod_plot

## Compare biodiversity metrics by iterating------------------------------------

biod_compare <- foreach(
  i = 1:numIterations,
  .combine = 'cbind',
  .packages = "tidyverse"
) %dopar% {
  full_df <- get.full_df(numEntities, numTraits, numPatches, maxEntitiesPatch)

  # iterate over years to obtain states over time
  for (i in 1:simYears) {
    full_df <- add_row(full_df, get.vulnerability(full_df, i - 1))
  }

  # at all times, identify hotspots under each biodiversity metric
  biod1_hotspots <- find.hotspots(calc.biodiv_1(full_df), "biodiv1")
  biod2_hotspots <- find.hotspots(calc.biodiv_2(full_df), "biodiv2")

  # compare identified hotspots between metrics
  calc.hotspot_compare(biod1_hotspots, biod2_hotspots)

}

# get quick distribution of comparison values, see how close they are to zero
hist(biod_compare)
