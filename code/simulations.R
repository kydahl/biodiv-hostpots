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
final_data <- read_csv("data/clean/final_dataset.csv")
tree <- readRDS("data/clean/full_tree.rds")

# Set the random seed used to generate figures in the manuscript 
set.seed(9523)

## Simulation parameters ----

## Number of iterations to compute biodiversity over
numIterations <- 100

# 1) Exploratory simulations ----------------------------------------------

# Compare hotspots identified by different metrics
full_df <- get.full_df(NumPatches = 1000)

# I added trait_names to the input of get.biodiv_df
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", 
                 "Leaf area (mm2)", "Diaspore mass (mg)")

#### Explore an example simulation with 40 patches ----


cl <- new_cluster(14)
cluster_library(cl, .packages())
cluster_copy(cl, lsf.str())
cluster_copy(cl, c("PD_dist_mat", "FD_dist_mat"))

explore_df <- get.full_df(1000) %>% 
  get.biodiv_df(., trait_names, tree)

# Plot correlations among biodiversity metrics
pair_plot <- explore_df %>% 
  select(-c(Patch)) %>% 
  ggpairs()

# 2) Calculate pairwise comparisons of biodiversity hotspot identi --------

# List of trait names used in final analysis
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", 
                 "Leaf area (mm2)", "Diaspore mass (mg)")

# List of all metric names
metric_names <- c("NumUnique", "NumIndigName", "NumUse",
                  "FD", "PD"
                  # "FRic", "FDiv", "FDis", "FEve", "Q", "richness", "GiniSimpson",
                  # "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh",
                  # "PSVs", "PSR"
)

# Set comparison parameters
numIterations <- 100
NumPatches <- 1000

# Set up parallel processing
plan(multisession, workers = 12, gc = TRUE)

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
    NumUnique = as.double(),
    NumIndigName = as.double(), NumUse = as.double(),
    FD = as.double(), PD = as.double(),
    # richness = as.double(), GiniSimpson = as.double(),
    # Simpson = as.double(), Shannon = as.double(),
    # Margalef = as.double(), Menhinick = as.double(),
    # McIntosh = as.double(), PSVs = as.double(),
    # PSR = as.double(), FRic = as.integer(),  FDiv = as.integer(),  
    # FDis = as.integer(),  FEve = as.integer(),  Q = as.integer(),
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
  # ) %dofuture% {
    ) %do% {
      # Iterate progress bar
      p(sprintf("j=%g", j))
      
      biodiv.compare_df = retry(
        biodiv_comp_helper_func(NumPatches, trait_names, tree, baseline_metric),
        until = function(val, cnd) {
          !is.null(val)
        },
        silent = FALSE,
        interval = 0
      )
      
      out_df = rbind(
        mutate(biodiv.compare_df %>%
                 select(-c(list_length, recall)) %>% 
                 pivot_wider(names_from = variable) %>%
                 unique(),
               type = "precision"),
        mutate(biodiv.compare_df %>%
                 select(-c(value, recall)) %>% 
                 pivot_wider(names_from = variable, values_from = list_length) %>%
                 unique(),
               type = "list_length"),
        mutate(biodiv.compare_df %>%
                 select(-c(value, list_length)) %>% 
                 pivot_wider(names_from = variable, values_from = recall) %>%
                 unique(),
               type = "recall")
      )
    }
  }
  
  sliceSize = 25
  sliceRange = 1:sliceSize
  for (i in 1:(numIterations/sliceSize)) {
    print(paste0("Collect simulation chunk # ", i, " of ", (numIterations/sliceSize), ":"))
    plan(multisession, workers = 12, gc = TRUE)
    compare_df = rbind(compare_df, get.temp_compare_df(sliceRange))
    plan(sequential)
  }
  
  out_df <- compare_df %>% 
    pivot_longer(cols = all_of(metric_names), names_to = "comparison") %>% 
    group_by(comparison, type) %>% 
    summarise(mean = mean(value), 
              var = var(value),
              .groups = "keep") %>% 
    mutate(baseline = baseline_metric)
  
  full_compare_df <- add_row(full_compare_df, out_df)
  rm(out_df)
  rm(compare_df)
}

saveRDS(full_compare_df, file = "results/full_comparisons_new.rds")
