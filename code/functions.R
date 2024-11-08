################################################################################
# Functions used in the simulation study
################################################################################

## Project: Comparing biodiversity hotspot identification
##
## Purpose: Define functions to be used in simulations
##
## Contents: 0) Set-up, load in necessary packages, functions, and data-sets
##           1) Helper functions
##           2) Functions to create patches and communities
##           3) Biodiversity metric functions
##           4) Hotspot functions
##           5) Put dataset into workable form for imputation and phylogeny
##              steps
##
## Inputs:  "data/clean/final_dataset.csv" and "data/clean/species_occurrences.rds"
##          
## Outputs: none
##
## Written and by: Kyle Dahlin and Elisa Van Cleemput
## Maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized September 2022
## _____________________________________________________________________________


# 0) Set-up, load in necessary packages, functions, and data-sets ---------

##* Load libraries ----
library(tidyverse)
library(progress) # for progress bars for long-running tasks
library(fundiversity) # For functional diversity metrics
library(adiv) # For phylogenetic diversity metrics
library(stringr)
library(V.PhyloMaker2)
library(picante)
library(maditr)

##* Load necessary data ----
final_data <- read_csv("data/clean/final_dataset.csv", show_col_types = FALSE)

SpeciesOccs <- read_rds("data/clean/species_occurrences.rds")
PD_dist_mat <- readRDS(file = "data/clean/PD_dist_mat.rds")
FD_dist_mat <- readRDS(file = "data/clean/FD_dist_mat.rds")


# 1) Helper functions -----------------------------------------------------
# Sample from a discrete uniform distribution
dunif_sampleone <- function(n) sample(1:n, 1, replace = T)

# Assign entries from df to groups of size 1 to n (chosen uniformly randomly)
group_assign <- function(in_df,
                         mean.NumEntities,
                         sd.NumEntities) {
  Num.Entities <- min(
    max(
      floor(rnorm(1, mean = mean.NumEntities, sd.NumEntities)),
      2 # to prevent the creation of single species communities
    ),
    length(in_df)
  )
  
  sample(in_df, size = Num.Entities, replace = FALSE)
}


# 2) Functions to create patches and communities --------------------------

# Get the number of entities for a patch
get.NumEntities <- function(MaxNumEntities, Level) {
  if (!(Level %in% MaxNumEntities$Level)) {
    warning("Incorrect level choice")
    break
  }
  MaxNum <- filter(MaxNumEntities, Level == Level)$count
  round(runif(1, 1, MaxNum))
}

## Create full data.frame of entities and their traits in patches with assigned states
get.full_df <- function(NumPatches) {
  
  ### Assign levels to patches ###
  
  # Equal numbers across each level (for now)
  Levels = sample(SpeciesOccs$Level, NumPatches, replace = TRUE)
  
  Init_df = tibble(Patch = 1:NumPatches, Level = Levels)
  
  ### Assign number of species to each patch ###
  
  # Maximum needs to be set to the total number of species at that level
  MaxNumEntities = SpeciesOccs %>% 
    group_by(Level) %>% 
    summarise(count = n())
  
  ### Assign species to patches ###
  
  Patch_df = data.table(Patch = as.integer(), Level = as.character(), Synonym = as.character())
  
  # For each patch, get its level
  Patch_df = foreach(
    j = 1:length(Init_df$Patch),
    .combine = "rbind"
  ) %do% {
    level = Init_df[j,]$Level
    
    # Get species list from the right level
    species_list = filter(SpeciesOccs, Level == level)$Synonym
    
    # Get number of entities in patch
    entity_count = min(get.NumEntities(MaxNumEntities, level), length(species_list))
    
    # Sample entities without replacement from species list
    entities = sample(species_list, entity_count, replace = FALSE)
    
    out_df = data.table(Patch = j, Level = level, Synonym = entities)
  }
  
  ### Add species traits to the dataframe ###
  right_join(Patch_df, final_data, by = "Synonym", relationship = "many-to-many") %>% 
    filter(!is.na(Patch), !is.na(Level))
}

# 3) Biodiversity metric functions ----------------------------------------
# for a patch, look at entities and their states & traits and output
# biodiversity measure of the patch

# If we add abundance, we can compute: Shannon diversity, Evenness, Simpson's 
# index, rare species indices, ... 

#### Calculate taxonomic diversity metrics ----
# number of unique entityIDs in a patch
num.unique.metric <- function(in_df) {
  out_df <- in_df %>%
    # remove duplicates due to rows from trait.num
    select(Patch, Species) %>%
    distinct() %>%
    # count number of distinct entities
    count(Patch, name = "biodiv")
}

# Number of endemic entities
num.endemic.metric <- function(in_df) {
  out_df <- in_df %>%
    select(Patch, Species) %>%
    # only keep rows where Species is unique ( == endemics)
    group_by(Species) %>%
    # remove duplicates due to rows from trait.num
    distinct() %>%
    mutate(endemic_flag = ifelse(n() == 1, 1, 0)) %>%
    ungroup() %>%
    group_by(Patch) %>%
    summarise(
      "biodiv" = sum(endemic_flag, na.rm = TRUE),
      .groups = "keep"
    )
}

#### Calculate TEK diversity metrics ----
# General trait counting metric function
trait.count.metric <- function(in_df, trait_name) {
  out_df <- in_df %>%
    select(Patch, one_of(trait_name)) %>% 
    group_by(Patch) %>% 
    # Biodiversity measured as sum in a given trait
    summarise(biodiv = sum(!!as.name(trait_name)))
}

#### Calculate functional diversity metrics ----

# Trait coefficient of variation function (only works for numerical traits)
trait.coeffvariance.metric <- function(in_df, trait_names) {
  out_df <- in_df %>%
    select(Patch, all_of(trait_names)) %>% 
    melt(id = "Patch") %>% 
    group_by(Patch, variable) %>% 
    # Biodiversity measured as sum in a given trait
    summarise(coeff_var = sd(value)/mean(value)) %>% 
    group_by(Patch) %>% 
    summarise(biodiv = sum(coeff_var))
}

# Functional diversity metrics (work for numerical and categorical traits)
trait.fdiv.metrics <- function(in_df){
  # Create dataframe of species x traits
  SpXTraits <- in_df %>%
    select(Species, `Nmass (mg/g)`:LDMC_log) %>% 
    # Scale quantitative traits
    mutate(across(where(is.numeric), function(.x){scale(.x, center = TRUE, scale = TRUE)})) %>% 
    unique() %>%
    arrange(Species) %>% # order alphabetically to match with PatchXSp
    column_to_rownames(var = "Species")
  
  # Create dataframe of patch x species with presence/absence
  PatchXSp <- in_df %>% 
    mutate(Presence = 1) %>%
    select(Patch, Species, Presence) %>%
    arrange(Species) %>% 
    # keys become column names and values are the entries
    pivot_wider(names_from = Species, values_from = Presence, values_fn = unique) %>% 
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "Patch")
  
  # # The number of species (columns) in PatchXSp must match the number of species (rows) in SpXTraits. 
  # if(nrow(SpXTraits) == ncol(PatchXSp)) {
  #   print("Equal amount of species in the SpXTraits and PatchXSp dataframes, so good to proceed")
  # } else {
  #   print("Warning: Unequal amount of species in the SpXTraits and PatchXSp dataframes, so check fdiv script")
  # }
  # 
  # # the species labels in PatchXSp and SpXTraits must be identical and in the same order.
  # if(identical(rownames(SpXTraits), colnames(PatchXSp)) == TRUE) {
  #   print("The species in the SpXTraits and PatchXSp dataframes are in the same order, so good to proceed")
  # } else {
  #   print("Warning: The species in the SpXTraits and PatchXSp dataframes are NOT in the same order, so check fdiv script")
  # }
  
  # Calculate functional diversity
  fdiv_df <- fd_fric(SpXTraits, as.matrix(PatchXSp), stand = TRUE) %>% # functional richness
    # functional divergence
    # right_join(fd_fdiv(SpXTraits, as.matrix(PatchXSp)), by = "site") %>%
    # functional dispersion
    right_join(fd_fdis(SpXTraits, as.matrix(PatchXSp)), by = "site") %>%
    # functional evenness
    # right_join(fd_feve(SpXTraits, as.matrix(PatchXSp)), by = "site") %>%
    # Rao's entropy (Q)
    # right_join(fd_raoq(SpXTraits, as.matrix(PatchXSp)), by = "site") %>%
    mutate(Patch = as.integer(site), .keep = "unused")
  
  return(fdiv_df)
  # Remarks: 
  # - Feve, Fric and FDiv cannot be calculated for communities with <3 functionally singular species. 
}

#### Calculate phylogenetic diversity metrics ----

get.phylo_tree <- function(in_df) {
  print("Create and prune phylogenetic tree")
  # Create list of all species in the data set
  list_species_unique <- in_df %>%
    mutate(Genus = stringr::word(Synonym, 1, 1, sep = " ")) %>%
    dplyr::select(Synonym, Genus, Family) %>%
    mutate(Synonym = str_replace_all(Synonym, "ssp.", "subsp.")) %>%
    unique() %>%
    as.data.frame()
  # Generate a phylogeny for the species list
  tree <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
                                     tree = GBOTB.extended.LCVP, # botanical nomenclature of the Leipzig catalogue of vascular plants
                                     nodes = nodes.info.1.LCVP, 
                                     output.tree = TRUE,
                                     scenarios = "S3")
  
  return(tree)
}

phylodiv.metrics <- function(in_df, tree) {
  
  # Create a community matrix for pruning:
  # convert in_df to wide format
  # Attention: Species names should have a dash between genus and species name
  full_df_for_comm = in_df %>%
    mutate(Synonym = str_replace_all(Synonym, "ssp.", "subsp.")) %>%
    mutate(Synonym = str_replace_all(Synonym," ","_")) %>%
    dplyr::select(Synonym, Patch) %>% 
    unique()
  # Community matrix
  comm = dcast(full_df_for_comm, formula = Patch ~ Synonym, fun.aggregate = length) %>%
    column_to_rownames(var="Patch")
  # Pruned phylogenetic tree
  tree_pruned = prune.sample(comm,tree$tree.scenario.3)
  
  # output:
  # - tree_pruned
  # - comm
  
  # Calculate phylogenetic diversity
  # 1) sum of the total phylogenetic branch length (Faith, 1992)
  pdiv_length = as.data.table(evodiv(tree_pruned, comm, method = "richness"))
  
  # 2) phylogenetic species variability, richness and evenness (Helmus et al., 2007)
  pdiv_psv = picante::psv(comm, tree_pruned, compute.var = FALSE) # 2.6 seconds
  pdiv_psr = picante::psr(comm, tree_pruned, compute.var = FALSE) # 2.4 seconds
  
  # Combine results
  pdiv_all = pdiv_length %>%
    tibble::rownames_to_column("Patch") %>%
    mutate(Patch = as.double(Patch)) %>% 
    full_join(pdiv_psv %>% dplyr::select(PSVs)  %>%
                tibble::rownames_to_column("Patch")%>%
                mutate(Patch = as.double(Patch)),
              by = c("Patch")) %>%
    full_join(pdiv_psr %>% dplyr::select(PSR)  %>%
                tibble::rownames_to_column("Patch")%>%
                mutate(Patch = as.double(Patch)),
              by = c("Patch"))
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

#### Build biodiversity data frame ----
get.biodiv_df <- function(in_df, tree) {
  # Metric 1: Species richness
  num.unique_df = num.unique.metric(in_df)
  
  # Metric 2: Number of endemic entities (NB: this was not considered in the final analysis)
  # num.endemic_df = num.endemic.metric(in_df)
  
  # Metric 3: Number of Indigenous names
  num.indig.name_df = trait.count.metric(in_df, "N_Names")
  
  # Metric 4: Number of Indigenous languages (NB: this was not considered in the final analysis)
  # num.indig.lang_df <- trait.count.metric(in_df, "N_Langs")
  
  # Metric 5: Number of uses
  num.use_df = trait.count.metric(in_df, "N_Uses")
  
  # Metric 6: Combined variation of quantitative traits
  # FD_df = FD_function(in_df)
  FDPD_df = FDPD_function(in_df)
  fdiv_df = trait.fdiv.metrics(in_df)
  
  # Phylogenetic diversity metrics
  # PD_df = PD_function(in_df)
  
  pdiv_df = phylodiv.metrics(in_df, tree)
  
  # Put together one big dataframe of biodiversity metrics of each Patch
  biodiv_df = rename(num.unique_df, NumUnique = biodiv) %>% # Number of unique entities
    # Number of Indigenous names
    right_join(rename(num.indig.name_df, NumIndigName = biodiv), by = "Patch") %>%
    # Number of uses
    right_join(rename(num.use_df, NumUse = biodiv), by = "Patch") %>% 
    # Functional diversity
    right_join(fdiv_df, by = "Patch") %>%
    right_join(FDPD_df, by = "Patch") %>%
    # Phylogenetic diversity
    # right_join(PD_df, by = "Patch") %>% 
    right_join(pdiv_df, by = "Patch")
}


# 4) Hotspot functions ----------------------------------------------------

#### Hotspot identifier function ----
# look at the biodiversity of all patches and determine which patches are
# "hotspots" relative to the other patches

# Get all the patches in the 95% quantile of a given biodiversity metric
find.hotspots <- function(in_df) {
  # Select top 5th percentile
  hotspots <- in_df[in_df$biodiv > quantile(in_df$biodiv, 0.95, na.rm = TRUE), ]
  
  if (dim(hotspots)[1] == 0) {hotspots <- in_df[which.max(in_df$biodiv),]}
  
  return(hotspots)
}

#### Hotspot comparison function ----
# measure the difference in the hotspot maps produced under different
# biodiversity metrics for the same dataset

# count the number of overlaps in which patches were considered hotspots
calc.hotspot_compare <- function(hotspots.baseline, hotspots.compare) {
  hotspots.baseline <- hotspots.baseline$Patch
  hotspots.compare <- hotspots.compare$Patch
  
  # Number of true positives
  TP_count <- sum(hotspots.compare %in% hotspots.baseline)
  
  # Use precision as our quantifier:
  # of the identified hotspots, what proportion match with the baseline list?
  # also calculate recall 
  out <- tibble(precision = TP_count / length(hotspots.compare),
                recall = TP_count / length(hotspots.baseline))
  return(out)
}


#### Build hotspot comparison data frame ----
get.compare_df <- function(in_df, baseline_metric) {
  # Get biodiversity metrics
  temp_df <- in_df %>%
    melt(id = "Patch")
  
  # Pre-allocate table
  hotspot.compare_df <- tibble(
    variable = character(),
    value = numeric(), 
    recall = numeric(), 
    list_length = numeric()
  )
  
  # Get baseline hotspots
  hotspots.base <- temp_df %>%
    filter(variable == baseline_metric) %>%
    select(-variable) %>%
    rename(biodiv = value) %>%
    find.hotspots()
  
  num.hotspots.base <- dim(hotspots.base)[1]
  
  other_metrics <- unique(filter(temp_df, variable != baseline_metric)$variable)
  
  hotspot.compare_df <- add_row(hotspot.compare_df,
                                variable = baseline_metric,
                                value = 1,
                                list_length = num.hotspots.base
                                
  )
  
  # Build hotspot comparison data frame
  for (j in other_metrics) {
    # Get hotspot list
    hotspot_list <- temp_df %>%
      filter(variable == j) %>%
      select(-variable) %>%
      rename(biodiv = value) %>%
      mutate(biodiv = ifelse(is.na(biodiv), 0, biodiv)) %>% 
      find.hotspots()
    
    num.hotspots.compare <- dim(hotspot_list)[1]
    
    # Compare lists
    compare_values <- calc.hotspot_compare(hotspot_list, hotspots.base)
    
    # Add to dataframe
    hotspot.compare_df <- add_row(hotspot.compare_df,
                                  variable = j,
                                  value = compare_values$precision,
                                  recall = compare_values$recall,
                                  list_length = num.hotspots.compare
                                  
    )
  }
  
  return(hotspot.compare_df)
}

# Helper function that gets wrapped in the for-loop to collect biodiversity hotspot comparisons
biodiv_comp_helper_func <- function(NumPatches, tree, baseline_metric) {
  
  # Get hotspot comparison values
  # Run a simulation
  out_df = get.full_df(NumPatches) %>% 
    # Calculate biodiversity metrics
    # multidplyr::partition(., cl) %>%
    get.biodiv_df(., tree) %>%
    # collect() %>% 
    # Make comparisons
    get.compare_df(., baseline_metric) 
  
}