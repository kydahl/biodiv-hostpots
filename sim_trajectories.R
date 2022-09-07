################################################################################
# Simulation trajectories over the time horizon
################################################################################

## Load libraries---------------------------------------------------------------
library(tidyverse)
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

# quick plot of biodiversity1 over time by patch
biod1_plot <- biod1_df %>% 
  group_by(patch) %>% 
  ggplot(aes(x = year, y = n, color = as.factor(patch))) +
  geom_line()




# at end time, identify hotspots under each biodiversity metric



