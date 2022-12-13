#### Create visualizations from simulations ####################################

## Load libraries --------------------------------------------------------------
require(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(GGally) # used to make paired scatter plots

# Load in simulations
source("code/simulations.R")

# Figures -----------------------------------------------------------------

# Figure 0: Distribution of TEK traits (across all species)
data_in

# Plot histograms in a single column
trait_plot <- data_in %>%
  select(Trait, value) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  # one histogram for each biodiversity metric
  facet_wrap( ~ Trait, scales = "free") +
  theme_cowplot(24)

trait_plot

# Figure 1: Compare the distribution of biodiversity metric values across patches

# Get a single simulation
full_df_sample <- get.full_df(entity_df, numPatches, mean.NumEntities, sd.NumEntities)

# Calculate biodiversity metrics
biodiv_df <- get.biodiv_df(full_df_sample)

# Plot histograms in a single column
histogram_plot <- biodiv_df %>%
  melt(id = "patch") %>%
  ggplot() +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  # one histogram for each biodiversity metric
  facet_wrap(~variable, scales = "free") +
  theme_cowplot(24)

histogram_plot

# Figure 2: Scatterplots of biodiversity metrics against each other
biodiv_scatter <- biodiv_df %>% 
  select(-patch) %>% 
  ggpairs(aes()) +
  theme_cowplot(16)
  
biodiv_scatter

# Figure 3: Comparing biodiversity hotspots relative to species richness -------
compare_plot <- compare_df %>%
  select(-iteration) %>%
  melt() %>%
  group_by(variable) %>% 
  mutate(mean = mean(value)) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 bins = 30
  ) +
  # add a vertical line showing the mean
  geom_vline(aes(xintercept = mean, group = variable), 
             colour = "red") +  
  
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +

  facet_wrap( ~ variable, scales = "free") +
  # x axis
  scale_x_continuous(name = "Relative precision", 
                   breaks = seq(0, 1, by = 0.1),
                   limits = c(0,1)) +
  # title
  ggtitle("How similar is each biodiversity metric to species richness?") +
  theme_cowplot(font_size = 24)

compare_plot


