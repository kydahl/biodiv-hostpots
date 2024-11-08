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

## Calculate phylogenetic distances ----
phylo_tree = tree$scenario.3

PD_dist_mat = cophenetic.phylo(phylo_tree)
PD_dist_mat = PD_dist_mat/max(PD_dist_mat)

saveRDS(PD_dist_mat, file = "data/clean/PD_dist_mat.rds")



