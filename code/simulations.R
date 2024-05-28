#### Create simulations of patches and calculate biodiversity ##################

## Load libraries---------------------------------------------------------------
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

# Load in analysis functions
source("code/functions.R")

# Load in trait data
final_data <- read_csv("data/clean/final_dataset.csv") # %>%
# rename(Species = Species_full)
tree <- readRDS("data/clean/full_tree.rds")

# Load in phylogenetic tree data

# tree <- read.tree(file = "data/clean/phylogenetic_tree.csv")
# Elisa: I commented this out, because we have been playing around with species names and hence, the tree might be different now

set.seed(9523)


## Parameters ------------------------------------------------------------------
## Set random seed (for debugging)
# set.seed(8797)

## Number of iterations to compute biodiversity over
numIterations <- 100

# Elisa: I added this to remove the patch with NA as a number, but it turns out there are more incomplete cases.
# I didn't check the reason for this, but shouldn't this dataset be complete (because it has imputed trait data)
# Kyle: I think na.omit removes rows where any of the columns are NA. So this is removing all the species with an NA in the subspecies or variant column. There may be species with NA in the Patch or Level column. These are the ones that weren't assigned to any of the patches. I made it so these get removed.

#################################
# Elisa: I used full_df to make phylo_div and funct_div work. what is the difference between full_df and test_df?
#################################

# Compare hotspots identified by different metrics
full_df <- get.full_df(NumPatches = 1000)

# I added trait_names to the input of get.biodiv_df
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", "Leaf area (mm2)")
# I set it up to just calculate the tree once
tree <- get.phylo_tree(final_data)

# KD: takes approximately 2 minutes, mostly because of phylogenetic diversity functions taking a long time


# Explore an example simulation -------------------------------------------
explore_df <- get.full_df(40) %>% 
  get.biodiv_df(., trait_names, tree)


# Plot correlations among biodiversity metrics
library(GGally)
pair_plot <- explore_df %>% 
  select(-c(Patch, NumEndemic, richness, GiniSimpson, Margalef, Menhinick, McIntosh)) %>% 
  ggpairs()


# Set up iterated simulations ---------------------------------------------

# # Set comparison parameters
# numIterations <- 100
# NumPatches <- 1000
# 
# baseline_metric <- "NumUnique"
# 
# # Initialize comparison data frame
# compare_df <- tibble(
#   NumUnique = as.double(),
#   NumEndemic = as.double(), NumIndigName = as.double(), NumUse = as.double(),
#   # NumIndigLang = as.double()
#   richness = as.double(), GiniSimpson = as.double(),
#   Simpson = as.double(), Shannon = as.double(),
#   Margalef = as.double(), Menhinick = as.double(),
#   McIntosh = as.double(), PSVs = as.double(),
#   PSR = as.double(), FRic = as.integer(),  FDiv = as.integer(),  
#   FDis = as.integer(),  FEve = as.integer(),  Q = as.integer(),
#   iteration = as.integer()
# ) %>%
#   # Remove the focal metric
#   select(-one_of(baseline_metric))
# 
# # Start new cluster for doParallel
# cluster_size <- parallel::detectCores() - 2
# my.cluster <- parallel::makeCluster(cluster_size, type = "PSOCK")
# # Register cluster for doParallel
# doSNOW::registerDoSNOW(cl = my.cluster)
# 
# # Set up progress bar
# pb <- progress_bar$new(
#   format = ":spin progress = :percent [:bar] elapsed: :elapsed | eta: :eta",
#   total = numIterations,
#   width = 100
# )
# progress <- function(n) {
#   pb$tick()
# }
# opts <- list(progress = progress)
# 
# # Calculate comparisons among diversity metrics
# compare_df <- foreach(
#   j = 1:numIterations,
#   # int = icount(),
#   .combine = "rbind",
#   .packages = c("tidyverse", "reshape2", "picante", "fundiversity", "adiv", "retry"),
#   .options.snow = opts
# ) %dopar%
#   {
#     j=j+1
#     gc()
#     biodiv.compare_df <- retry(
#       biodiv_comp_helper_func(NumPatches, trait_names, tree, baseline_metric),
#       until = function(val, cnd) {
#         !is.null(val)
#       }
#     ) %>% 
#       mutate(iteration = j)
#     
#     precision_df <- biodiv.compare_df %>%
#       select(-c(list_length, recall)) %>% 
#       pivot_wider(names_from = variable) %>%
#       unique() 
#     
#     list_length_df <- biodiv.compare_df %>%
#       select(-c(value, recall)) %>% 
#       pivot_wider(names_from = variable, values_from = list_length) %>%
#       unique() 
#     
#     recall_df <- biodiv.compare_df %>%
#       select(-c(value, list_length)) %>% 
#       pivot_wider(names_from = variable, values_from = recall) %>%
#       unique() 
#     
#     out_df <- rbind(
#       mutate(precision_df, type = "precision"),
#       mutate(list_length_df, type = "list_length"),
#       mutate(recall_df, type = "recall")
#     )
#     
#     # # Add to the list
#     # compare_df <- add_row(
#     #   compare_df,
#     #   biodiv.compare_df
#     # )
#   } %>% unique()
# 
# stopCluster(my.cluster)


# Pairwise precision comparisons -------------------------------------
# I added trait_names to the input of get.biodiv_df
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", 
                 "Leaf area (mm2)", "Diaspore mass (mg)")
# Calculate the full phylogenetic tree once
# tree <- get.phylo_tree(final_data)

# List of all metric names
metric_names <- c("NumUnique", "NumIndigName", "NumUse",
                  "FRic", "FDiv", "FDis", "FEve", "Q", "richness", "GiniSimpson",
                  "Simpson", "Shannon", "Margalef", "Menhinick", "McIntosh",
                  "PSVs", "PSR")
# Set comparison parameters
numIterations <- 100
NumPatches <- 1000

# plan(multisession, workers = 12, gc = TRUE)
# plan(callr)
# Set up progress bar
handlers(global = TRUE)
handlers("cli") #handlers("progress")

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
               " ( # ", metric_index, " of ", length(metric_names), " ):"))
  baseline_metric <- metric_name
  
  # Initialize comparison data frame
  compare_df <- tibble(
    NumUnique = as.double(),
    NumIndigName = as.double(), NumUse = as.double(),
    # NumIndigLang = as.double()
    richness = as.double(), GiniSimpson = as.double(),
    Simpson = as.double(), Shannon = as.double(),
    Margalef = as.double(), Menhinick = as.double(),
    McIntosh = as.double(), PSVs = as.double(),
    PSR = as.double(), FRic = as.integer(),  FDiv = as.integer(),  
    FDis = as.integer(),  FEve = as.integer(),  Q = as.integer(),
    iteration = as.integer()
  ) %>%
    # Remove the focal metric
    select(-one_of(baseline_metric))
  
  
  # for (slice_index in 1:slice_indices) {
  #   # Set up chunk of simulations to run and progress bar
  #   slice_range = (slice_index - 1)*sliceLength + 1:sliceLength
  
  
  # Calculate comparisons among diversity metrics
  get.temp_compare_df = function(slice_range) {
    p = progressor(along = slice_range)
    temp_compare_df = foreach(
      j = slice_range,
      .inorder = FALSE,
      .combine = "rbind",
      .options.future = list(seed = TRUE) # ensures true RNG among parallel processes
    ) %dofuture% {
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
  # slice_range = 1:(numIterations/2) #slice_range[slice_range <= numIterations]
  # 
  # temp_compare_df1 = get.temp_compare_df(slice_range)
  # temp_compare_df2 = get.temp_compare_df(slice_range)
  # 
  # compare_df <- rbind(temp_compare_df1, temp_compare_df2)
  
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

saveRDS(full_compare_df, file = "full_comparisons.rds")
