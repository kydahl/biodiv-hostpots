################################################################################
# Simulation trajectories over the time horizon
################################################################################

## Load libraries---------------------------------------------------------------
library(tidyverse)
library(gridExtra)
source("sim_functions.R")

## Helper functions-------------------------------------------------------------


## Parameters-------------------------------------------------------------------

# Number of years to simulate
simYears <- 10

## Create trajectories----------------------------------------------------------

biod_compare <- c()
for (i in 1:100) {
source("sim_traits.R")
# iterate over years to obtain states over time
for (i in 1:simYears) {full_df <- add_row(full_df, get.vulnerability(full_df, i-1))}

# compute biodiversity at each time point for each biodiversity metric
biod1_df <- calc.biodiv_1(full_df)
biod2_df <- calc.biodiv_2(full_df)

# uncomment below if we want one big dataframe of biodiversity measures
biodiv_df <- full_join(biod1_df, biod2_df)

# # quick plot of biodiversity1 (species richness) over time by patch
# biod1_plot <- biod1_df %>% 
#   group_by(patch) %>% 
#   ggplot(aes(x = year, y = biodiv1, color = as.factor(patch))) +
#   geom_line(lwd = 1) +
#   geom_point() +
#   scale_y_continuous(name = "Species richness") +
#   scale_color_discrete(name = "Patch") +
#   theme_minimal()
# 
# # quick plot of biodiversity2 (# endemics) over time by patch
# biod2_plot <- biod2_df %>% 
#   group_by(patch) %>% 
#   ggplot(aes(x = year, y = biodiv2, color = as.factor(patch))) +
#   geom_line(lwd = 1, linetype =2) +
#   geom_point() +
#   scale_y_continuous(name = "Num. endemics") +
#   scale_color_discrete(name = "Patch") +
#   theme_minimal()
# 
# grid.arrange(biod1_plot + guides(color = "none"),biod2_plot, ncol = 2)

# at all times, identify hotspots under each biodiversity metric
biod1_hotspots <- find.hotspots(biod1_df, "biodiv1")
biod2_hotspots <- find.hotspots(biod2_df, "biodiv2")


# compare identified hotspots between metrics
biod_compare[i] <- calc.hotspot_compare(biod1_hotspots, biod2_hotspots)


# repeat this X times
rm(full_df)
i
}

# get distribution of comparison values, see how close they are to zero
hist(biod_compare)