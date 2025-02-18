################################################################################
# Patch simulations and biodiversity hotspot calculations
################################################################################

## Project: Comparing biodiversity hotspot identification
##
## Purpose: Create simulations of patches, calculate biodiversity metrics, identify hotspots
##
## Contents: 0) Set-up, load in necessary packages, functions, and data-sets
##           1) Exploratory simulations
##           2) Calculate pairwise comparisons of biodiversity hotspot identifications
##
## Inputs:  - Source: code/functions.R
##          - data/clean/final_dataset.csv
##          - data/clean/full_tree.rds
##          
## Outputs: - results/full_comparisons.rds
##
## Written and by: Kyle Dahlin and Elisa Van Cleemput
## Maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized April 2023
## _____________________________________________________________________________


# 0) Set-up, load in necessary packages, functions, and data-sets ---------

## Load libraries ----
library(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(picante)
library(progress)
library(retry)
library(doFuture)
library(future.callr)
library(progressr)
library(GGally) # To nicely plot correlations among variables

# Load in analysis functions
source("code/functions.R")

# Load in dataset and phylogenetic tree
final_data <- read_csv("data/clean/final_dataset_no_taxa.csv") #%>% # read_csv("data/clean/final_dataset.csv") %>% 
# # log transform traits with large outliers
# mutate(LeafArea_log = log(`Leaf area (mm2)`), .keep = 'unused') %>%
# mutate(PlantHeight_log = log(`Plant height (m)`), .keep = 'unused') %>%
# mutate(DiasporeMass_log = log(`Diaspore mass (mg)`), .keep = 'unused') %>%
# mutate(LDMC_log = log(`LDMC (g/g)`), .keep = 'unused') %>%
# # Turn categorical traits into quantitative ones
# mutate(Woodiness = ifelse(Woodiness == "woody", 1, 0)) %>% 
# # Put traits at the end
# relocate(c(`Nmass (mg/g)`, Woodiness, LeafArea_log:LDMC_log), .after = last_col())

# tree <- readRDS("data/clean/full_tree.rds") # old tree
tree <- readRDS("data/clean/final_tree.rds")
# tree <- get.phylo_tree(final_data)

# Set the random seed used to generate figures in the manuscript 
set.seed(9523)

## Simulation parameters ----

# 1) Exploratory simulations ----------------------------------------------

# # Compare hotspots identified by different metrics
# full_df <- get.full_df(NumPatches = 1000)

# # I added trait_names to the input of get.biodiv_df
# trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", 
#                  "Leaf area (mm2)", "Diaspore mass (mg)")

#### Explore an example simulation with 40 patches ----

# # Plot correlations among biodiversity metrics
# pair_plot <- explore_df %>% 
#   select(-c(Patch)) %>% 
#   ggpairs()

# 2) Calculate pairwise comparisons of biodiversity hotspot identi --------

# List of trait names used in final analysis
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", 
                 "Leaf area (mm2)", "Diaspore mass (mg)")

# List of all metric names
metric_names <- c(
  # Taxonomic
  "NumUnique",
  # TEK
  "NumIndigName", "NumUse",
  # Phylogenetic
  "richness", "PD", "PSVs", "PSR",
  # Functional
  "FD", "FRic", "FDis"
)

# Set comparison parameters
numIterations <- 100
NumPatches <- 1000 # For scaling up: increase to 10000, for scaling down go to 100

# Set up progress bar
handlers(global = TRUE)
handlers("cli")

full_compare_df <- tibble(
  baseline = as.character(), # baseline biodiversity metric
  comparison = as.character(), # comparison biodiversity metric
  type = as.character(), # type (precision or list length)
  mean = as.double(),
  var = as.double()
)
metric_index = 0

for (metric_name in metric_names) {
  gc()
  metric_index = metric_index + 1
  print(paste0("Calculating ", metric_name,
               " ( Metric # ", metric_index, " of ", length(metric_names), " ):"))
  baseline_metric <- metric_name
  
  # Initialize comparison data frame
  compare_df <- tibble(
    # Taxonomic
    NumUnique = as.double(),
    # TEK
    NumIndigName = as.double(), NumUse = as.double(),
    # Phylogenetic
    PD = as.double(), richness = as.double(), PSVs = as.double(), PSR = as.double(),
    # Functional
    FD = as.double(), FRic = as.integer(), FDis = as.integer(),
    iteration = as.integer()
  ) %>%
    # Remove the focal metric
    select(-one_of(baseline_metric))
  
  # Calculate comparisons among diversity metrics
  get.temp_compare_df = function(slice_range) {
    p = progressor(along = slice_range)
    temp_compare_df = foreach(
      j = slice_range,
      .inorder = FALSE,
      .combine = "rbind",
      .options.future = list(seed = TRUE) # ensures true RNG among parallel processes
    ) %dofuture% {
      
      # Get biodiversity hotspot comparisons
      biodiv.compare_df = retry(
        biodiv_comp_helper_func(NumPatches, tree, baseline_metric),
        until = function(val, cnd) {
          !is.null(val)
        },
        silent = FALSE,
        interval = 0
      )
      
      # Set up output dataframe
      out_df = rbind(
        mutate(biodiv.compare_df %>%
                 select(-c(list_length, recall, jaccard)) %>% 
                 pivot_wider(names_from = variable) %>%
                 unique(),
               type = "precision"),
        mutate(biodiv.compare_df %>%
                 select(-c(value, recall, jaccard)) %>% 
                 pivot_wider(names_from = variable, values_from = list_length) %>%
                 unique(),
               type = "list_length"),
        mutate(biodiv.compare_df %>%
                 select(-c(value, list_length, jaccard)) %>% 
                 pivot_wider(names_from = variable, values_from = recall) %>%
                 unique(),
               type = "recall"),
        mutate(biodiv.compare_df %>%
                 select(-c(value, list_length, recall)) %>% 
                 pivot_wider(names_from = variable, values_from = jaccard) %>%
                 unique(),
               type = "jaccard")
      )
      
      # Iterate progress bar
      p(sprintf("j=%g", j))
      biodiv.compare_df
    }
  }
  
  sliceSize = numIterations/1 # Change this value to match memory/CPU availability
  sliceRange = 1:sliceSize
  for (i in 1:(numIterations/sliceSize)) {
    print(paste0("Collect simulation chunk # ", i, " of ", (numIterations/sliceSize), ":"))
    plan(multisession, workers = availableCores()-2, gc = TRUE)
    compare_df = rbind(compare_df, get.temp_compare_df(sliceRange))
    plan(sequential)
  }
  
  out_df <- compare_df %>% 
    rename(
      comparison = variable,
      precision = value
    ) %>% 
    pivot_longer(cols = precision:list_length, names_to = "type") %>% 
    group_by(comparison, type) %>%
    summarise(
      mean = mean(value),
      var = var(value),
      .groups = "keep"
    ) %>%
    mutate(baseline = baseline_metric)
  
  full_compare_df <- add_row(full_compare_df, out_df)
  rm(out_df)
  rm(compare_df)
}

saveRDS(full_compare_df, file = "results/final_comparisons.rds")