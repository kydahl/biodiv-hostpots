#### Create simulations of patches and calculate biodiversity ##################

## Load libraries---------------------------------------------------------------
library(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(picante)
library(progress)
library(retry)

# Load in analysis functions
source("code/functions.R")

# Load in trait data
final_data <- read_csv("data/clean/final_dataset.csv") # %>%
# rename(Species = Species_full)

# Load in phylogenetic tree data
# tree <- read.tree(file = "data/clean/phylogenetic_tree.csv")
# Elisa: I commented this out, because we have been playing around with species names and hence, the tree might be different now


## Helper functions ------------------------------------------------------------
# Uniformly sample from interval [1, MaxNumEntities]
get.NumEntities <- function(MaxNumEntities, Level) {
  if (!(Level %in% MaxNumEntities$Level)) {
    warning("Incorrect level choice")
    break
  }
  
  MaxNum <- filter(MaxNumEntities, Level == Level)$count
  round(runif(1, 1, MaxNum))
}

## Parameters ------------------------------------------------------------------
## Set random seed (for debugging)
# set.seed(8797)

## Number of iterations to compute biodiversity over
numIterations <- 1000

# Elisa: I added this to remove the patch with NA as a number, but it turns out there are more incomplete cases.
# I didn't check the reason for this, but shouldn't this dataset be complete (because it has imputed trait data)
# Kyle: I think na.omit removes rows where any of the columns are NA. So this is removing all the species with an NA in the subspecies or variant column. There may be species with NA in the Patch or Level column. These are the ones that weren't assigned to any of the patches. I made it so these get removed.

#################################
# Elisa: I used full_df to make phylo_div and funct_div work. what is the difference between full_df and test_df?
#################################

get.full_df <- function(NumPatches) {
  SpeciesOccs <- read_rds("data/clean/species_occurrences.rds")
  
  ### Assign levels to patches ###
  
  # Equal numbers across each level (for now)
  Levels <- sample(SpeciesOccs$Level, NumPatches, replace = TRUE)
  
  Init_df <- tibble(Patch = 1:NumPatches, Level = Levels)
  
  ### Assign number of species to each patch ###
  
  # Maximum needs to be set to the total number of species at that level
  MaxNumEntities <- SpeciesOccs %>%
    group_by(Level) %>%
    summarise(count = n())
  
  ### Assign species to patches ###
  
  Patch_df <- tibble(Patch = as.integer(), Level = as.double(), Synonym = as.character())
  
  # For each patch, get its level
  for (index_patch in Init_df$Patch) {
    # Get the level of the patch
    level <- filter(Init_df, Patch == index_patch)$Level
    
    # Get species list from the right level
    species_list <- filter(SpeciesOccs, Level == level)$Synonym
    
    # Get number of entities in patch
    entity_count <- min(get.NumEntities(MaxNumEntities, level), length(species_list))
    
    # Sample entities without replacement from species list
    entities <- sample(species_list, entity_count, replace = FALSE)
    
    temp_df <- tibble(Patch = index_patch, Level = level, Synonym = entities)
    
    Patch_df <- rbind(Patch_df, temp_df)
  }
  
  ### Add species traits to the dataframe ###
  full_df <- right_join(Patch_df, final_data, by = "Synonym", relationship = "many-to-many") %>%
    filter(!is.na(Patch), !is.na(Level))
}

# Compare hotspots identified by different metrics
full_df <- get.full_df(NumPatches = 400)

# I added trait_names to the input of get.biodiv_df
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", "Leaf area (mm2)")
# I set it up to just calculate the tree once
tree <- get.phylo_tree(full_df)

# KD: takes approximately 2 minutes, mostly because of phylogenetic diversity functions taking a long time


# Explore an example simulation -------------------------------------------
explore_df <- get.full_df(400) %>% 
  get.biodiv_df(., trait_names, tree)


# Plot correlations among biodiversity metrics
library(GGally)
pair_plot <- explore_df %>% 
  select(-c(Patch, NumEndemic, richness, GiniSimpson, Margalef, Menhinick, McIntosh)) %>% 
  ggpairs()


# Set up iterated simulations ---------------------------------------------

# Set comparison parameters
numIterations <- 1000
NumPatches <- 40

baseline_metric <- "NumUnique"

# Initialize comparison data frame
compare_df <- tibble(
  NumUnique = as.double(),
  NumEndemic = as.double(), NumIndigName = as.double(),
  NumIndigLang = as.double(), NumUse = as.double(),
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

# Start new cluster for doParallel
cluster_size <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(cluster_size, type = "PSOCK")
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

# Set up progress bar
pb <- progress_bar$new(
  format = ":spin progress = :percent [:bar] elapsed: :elapsed | eta: :eta",
  total = numIterations,
  width = 100
)
progress <- function(n) {
  pb$tick()
}
opts <- list(progress = progress)

# Calculate comparisons among diversity metrics
compare_df <- foreach(
  j = 1:numIterations,
  # int = icount(),
  .combine = "rbind",
  .packages = c("tidyverse", "reshape2", "picante", "fundiversity", "adiv", "retry"),
  .options.snow = opts
) %dopar%
  {
    gc()
    biodiv.compare_df <- retry(
      biodiv_comp_helper_func(NumPatches, trait_names, tree),
      until = function(val, cnd) {
        !is.null(val)
      }
    ) %>% 
      mutate(iteration = j)
    
    # # Add to the list
    # compare_df <- add_row(
    #   compare_df,
    #   biodiv.compare_df
    # )
  } %>% unique()

stopCluster(my.cluster)

# for (j in 1:numIterations) {
#   print(paste0("########### ITERATION # ", j, " ########### " ))
#   
#   biodiv.compare_df <- retry(
#     get.comparisons(NumPatches, trait_names, tree),
#     until = function(val, cnd) {
#       !is.null(val)
#     }
#   ) %>% 
#     mutate(iteration = j)
#   
#   # Add to the list
#   compare_df <- add_row(
#     compare_df,
#     biodiv.compare_df
#   )
# } 
