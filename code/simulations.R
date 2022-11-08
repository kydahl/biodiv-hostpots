################################################################################
# Simulation trajectories over the time horizon
################################################################################

## Load libraries---------------------------------------------------------------
library(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
registerDoParallel(cores = 6 - 1)
source("code/functions.R")
source("code/data_intake.R")

## Helper functions-------------------------------------------------------------


## Parameters-------------------------------------------------------------------
## Random seed (to keep things consistent while coding at first)
# set.seed(8797)

## Number of traits
numTraits <- 3 # !!! Is this still used?

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
#   the same (relatively large) pool of entities, then it is very likely that any given entity in a
#   patch is endemic there.
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
state_keys <- tibble(nums = 1:6, cats = IUCN_cats) # in case we want to convert numbers to categories


## Number of iterations to compute biodiversity over
numIterations <- 1000

## Create data frame
entity_df <- data_in
full_df <- get.full_df(entity_df, numTraits, numPatches, mean.NumEntities, sd.NumEntities)

## Basic statistics-------------------------------------------------------------

## Experiment 1:
## Generate 1000 "environments", in each, calculate biodiversity metrics and
## determine the hotspot patches, then measure the difference in how each
## metric assigns hotspots. Across those 1000 environments, so how each metric
## matches up. Which are most similar to each other? Which are most different?
## We can do this by creating a correlation matrix.

# for (i in 1 : numIterations) {
full_df <- get.full_df(entity_df, numTraits, numPatches, mean.NumEntities, sd.NumEntities)

# Number of unique entities
num.unique_df <- num.unique.metric(full_df)
hotspots.unique <- find.hotspots(num.unique_df)

# Number of endemic entities
num.endemic_df <- num.endemic.metric(full_df)
hotspots.endemic <- find.hotspots(num.endemic_df)

# Number of Indigenous names
num.indig.name_df <- trait.count.metric(
  full_df,
  "Number of unique names (Indigenous)"
)
hotspots.indig.name <- find.hotspots(num.indig.name_df)

# Number of Indigenous languages
num.indig.lang_df <- trait.count.metric(
  full_df,
  "Number of unique  languages (from Appendix 2B)"
)
hotspots.indig.lang <- find.hotspots(num.indig.lang_df)

# Number of uses
num.use_df <- trait.count.metric(
  full_df,
  "Number of uses"
)
hotspots.use <- find.hotspots(num.use_df)


# Quick analysis plots --------------------------------------------------------

## Comparing the distributions of each biodiversity metric


# Put together one big dataframe of biodiversity metrics of each patch
biodiv_df <- rename(num.unique_df, NumUnique = biodiv) %>% # Number of unique entities
  # Number of endemic entities
  right_join(rename(num.endemic_df, NumEndemic = biodiv)) %>%
  # Number of Indigenous names
  right_join(rename(num.indig.name_df, NumIndigName = biodiv)) %>%
  # Number of Indigenous languages
  right_join(rename(num.indig.lang_df, NumIndigLang = biodiv)) %>%
  # Number of uses
  right_join(rename(num.use_df, NumUse = biodiv))

hist(biodiv_df$NumUnique)

# Plot histograms in a single column
histogram_plot <- biodiv_df %>%
  melt(id = "patch") %>%
  ggplot() +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
    colour = "black",
    fill = "white"
  ) +
  # add a spline approximating the probability denstiy function
  geom_density(aes(x = value),
    alpha = .2,
    fill = "#FF6666"
  ) +
  # one histogram for each biodiversity metric
  facet_wrap(~variable, scales = "free") +
  theme_cowplot()

histogram_plot


# Comparing the hotspots of each biodiversity metric to number of species ------

# Preallocate comparison data frame
compare_df <- tibble(
  iteration = numeric(),
  endemic = numeric(),
  indig.name = numeric(),
  indig.lang = numeric(),
  use = numeric()
)

# Generate many datasets

compare_df <- foreach(
  j = 1:numIterations,
  int = icount(),
  .combine = "rbind",
  .packages = c("tidyverse")
) %dopar% {
  # Get the data frame
  full_df <- get.full_df(entity_df, numTraits, numPatches, mean.NumEntities, sd.NumEntities)

  # Calculate number of unique species in each patch
  num.unique_df <- num.unique.metric(full_df)
  hotspots.unique <- find.hotspots(num.unique_df)

  # Number of endemic entities
  num.endemic_df <- num.endemic.metric(full_df)
  hotspots.endemic <- find.hotspots(num.endemic_df)
  endemic.compare <- calc.hotspot_compare(hotspots.unique, hotspots.endemic)

  # Number of Indigenous names
  num.indig.name_df <- trait.count.metric(
    full_df,
    "Number of unique names (Indigenous)"
  )
  hotspots.indig.name <- find.hotspots(num.indig.name_df)
  indig.name.compare <- calc.hotspot_compare(hotspots.unique, hotspots.indig.name)

  # Number of Indigenous languages
  num.indig.lang_df <- trait.count.metric(
    full_df,
    "Number of unique  languages (from Appendix 2B)"
  )
  hotspots.indig.lang <- find.hotspots(num.indig.lang_df)
  indig.lang.compare <- calc.hotspot_compare(hotspots.unique, hotspots.indig.lang)

  # Number of uses
  num.use_df <- trait.count.metric(
    full_df,
    "Number of uses"
  )
  hotspots.use <- find.hotspots(num.use_df)
  use.compare <- calc.hotspot_compare(hotspots.unique, hotspots.use)

  compare_df <- add_row(compare_df,
    iteration = int,
    endemic = endemic.compare,
    indig.name = indig.name.compare,
    indig.lang = indig.lang.compare,
    use = use.compare
  )
}

compare_plot <- compare_df %>%
  select(-iteration) %>%
  melt() %>%
  ggplot() +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
    bins = 10
  ) +
  # add a spline approximating the probability denstiy function
  # geom_density(aes(x = value),
  #   alpha = .2,
  #   fill = "#FF6666"
  # ) +
  # one histogram for each biodiversity metric
  facet_wrap(~variable, scales = "free") +
  # x label
  xlab("Relative precision") +
  # title
  ggtitle("How similar is each biodiversity metric to species diversity?") +
  theme_cowplot()

compare_plot
