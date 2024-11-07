# New biodiversity measuring functions following Doxa et al., 2020

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

library(gtools)
library(multidplyr)

# Calculate a version of Rao's quadratic entropy for taxonomic, functional, and phylogenetic diversity.

# Diversity is estimated as:
# Q = sum_i sum_j d_ij p_i p_j

# where p_i and p_j are the frequencies of species i and j
# d_ij is the distance between them.

# In our case, because we are only doing presence-absence, p_i and p_j will either be one or zero.

# Taxonomic diversity: just set d_ij = 1 to get the Simpson diversity index

# Functional diversity:
# Measured by the Euclidean distance between the species in PCoA space, after normalization to [0,1] by dividing by the maximum value.
#   Question: Is this the FULL PCoA space? Or the reduced space of just that community? Going with full space.

# Phylogenetic diversity:
# Calculate the ultrametric phylogenetic distances between species then divide by the maximum.
# Going with full space here as well


# They measured whether or not there were significant differences among the metrics by calculating the standardised effect size:
# SES = (Div_obs - mean(Div_sim)) / sd_sim

# Where in their case, simulations amounted to permuting the species names in a plot. 

# Step 1: Calculate distance matrices for FD and PD ----
# We will create "look up tables" for these FD and PD to quickly assess the value of d_ij for each pair of species.

# Load in analysis functions
source("code/functions.R")

# Load in dataset and phylogenetic tree
final_data <- read_csv("data/clean/final_dataset.csv")
tree <- readRDS("data/clean/full_tree.rds")

## Calculate functional distances ----

# Get principal components
PCoA = final_data %>% 
  mutate(
    Woodiness = case_when(
      Woodiness == "non-woody" ~ 0,
      Woodiness == "woody" ~ 1
    )) %>% 
  select("Nmass (mg/g)":"Leaf area (mm2)") %>%
  prcomp(scale = TRUE)
# Calculate distance matrix
FD_dist_mat = dist(PCoA$x) %>% as.matrix()
# Normalize values to 0, 1
FD_dist_mat = FD_dist_mat / max(FD_dist_mat)

saveRDS(FD_dist_mat, file = "data/clean/FD_dist_mat.rds")

# FD_distance_function <- function(species_id_i, species_id_j) {
#   return(FD_dist_mat[species_id_i, species_id_j])
# }


## Calculate phylogenetic distances ----
phylo_tree = tree$scenario.3

PD_dist_mat = cophenetic.phylo(phylo_tree)
PD_dist_mat = PD_dist_mat/max(PD_dist_mat)

saveRDS(PD_dist_mat, file = "data/clean/PD_dist_mat.rds")

# PD_distance_function <- function(species_name_i, species_name_j) {
#   return(PD_dist_mat[species_name_i, species_name_j])
# }

# Calculate all biodiversity types ----


# Given a community J

# TD = species richness
TD_function = function(community_df) {
  return(dim(community_df)[1])
}

# Efficient pairwise PD distance summation using matrix indexing
PD_dataframe_function <- function(species_names) {
  if (length(species_names) <= 1) return(0)
  
  # Match species to indices
  species_names <- gsub(" ", "_", species_names)
  species_indices <- match(species_names, colnames(PD_dist_mat))
  
  # Direct summation of pairwise distances from PD_dist_mat
  dist_values <- PD_dist_mat[species_indices, species_indices, drop = FALSE]
  sum(dist_values[lower.tri(dist_values)])
  
}

# Efficient pairwise FD distance summation using matrix indexing
FD_dataframe_function <- function(species_ids) {
  if (length(species_ids) <= 1) return(0)
  
  # Match species to indices
  species_indices <- match(species_ids, colnames(FD_dist_mat))
  
  # Direct summation of pairwise distances from FD_dist_mat
  dist_values <- FD_dist_mat[species_indices, species_indices, drop = FALSE]
  sum(dist_values[lower.tri(dist_values)])
  
}

# Refactored FDPD_function with optimized distance calculation
FDPD_function <- function(in_df) {
  in_df %>%
    select(Patch, species_id, Synonym) %>%
    group_by(Patch) %>%
    summarise(
      FD = FD_dataframe_function(species_id),
      PD = PD_dataframe_function(Synonym),
      .groups = 'drop'
    ) 
}










