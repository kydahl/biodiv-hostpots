################################################################################
# Simulation trajectories over the time horizon
################################################################################

## Load libraries---------------------------------------------------------------
library(tidyverse)
library(gridExtra)
source("sim_traits.R")
source("sim_functions.R")

## Helper functions-------------------------------------------------------------


## Parameters-------------------------------------------------------------------

# Number of years to simulate
simYears <- 20

## Create trajectories----------------------------------------------------------

# iterate over years to obtain states over time
for (i in 1:simYears) {full_df <- add_row(full_df, vulner_func(full_df, i-1))}

# compute biodiversity at each time point for each biodiversity metric
biod1_df <- biodiv_1(full_df)
biod2_df <- biodiv_2(full_df)

# quick plot of biodiversity1 (species richness) over time by patch
biod1_plot <- biod1_df %>% 
  group_by(patch) %>% 
  ggplot(aes(x = year, y = n, color = as.factor(patch))) +
  geom_line(lwd = 1) +
  geom_point() +
  scale_y_continuous(name = "Species richness") +
  scale_color_discrete(name = "Patch") +
  theme_minimal()

# quick plot of biodiversity2 (# endemics) over time by patch
biod2_plot <- biod2_df %>% 
  group_by(patch) %>% 
  ggplot(aes(x = year, y = n, color = as.factor(patch))) +
  geom_line(lwd = 1, linetype =2) +
  geom_point() +
  scale_y_continuous(name = "Num. endemics") +
  scale_color_discrete(name = "Patch") +
  theme_minimal()

grid.arrange(biod1_plot + guides(color = "none"),biod2_plot, ncol = 2)

# at end time, identify hotspots under each biodiversity metric

# compare identified hotspots between metrics

# repeat this X times

