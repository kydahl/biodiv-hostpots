################################################################################
# Functions used in the simulation study
################################################################################

# Kyle Dahlin, September 2022
##* Load libraries---------------------------------------------------------------
library(tidyverse)

##* Load necessary data---------------------------------------------------------
final_data <- read_csv("data/clean/final_dataset_IMPUTED.csv", 
                       show_col_types = FALSE) %>% 
  rename(Species = Species_full)

##* Helper functions-------------------------------------------------------------
# Sample from a discrete uniform distribution
# Source: https://stats.stackexchange.com/questions/3930/are-there-default-functions-for-discrete-uniform-distributions-in-r
dunif_sampleone <- function(n) sample(1:n, 1, replace = T)

# Assign entries from df to groups of size 1 to n (chosen uniformly randomly)
group_assign <- function(in_df,
                         mean.NumEntities,
                         sd.NumEntities) {
  Num.Entities <- min(
    max(
      floor(rnorm(1, mean = mean.NumEntities, sd.NumEntities)),
      2 # !!! to prevent single species communities
    ),
    length(in_df)
  )

  sample(in_df, size = Num.Entities, replace = FALSE)
}

##* Functions for creating data frames-------------------------------------------

# Get the number of entities for a patch
get.NumEntities <- function(MaxNumEntities, Level) {
  if (!(Level %in% MaxNumEntities$Level)) {
    warning("Incorrect level choice")
    break
  }
  
  MaxNum <- filter(MaxNumEntities, Level == Level)$count
  round(runif(1, 1, MaxNum))
}

## Assign entities to patches
get.patches <- function(entity_df,
                        numPatches,
                        mean.NumEntities,
                        sd.NumEntities) {
  # initialize data.frame
  regional_df <- tibble(
    Species = NA,
    patch = NA
  )

  # sample from entities list to fill patches
  for (i in 1:numPatches) {
    regional_df <- add_row(regional_df,
      patch = i,
      Species = group_assign(
        unique(entity_df$Species),
        mean.NumEntities,
        sd.NumEntities
      )
    )
  }

  # assign uniqueIDs to entities in patches (and remove NAs we added in initialization)
  regional_df <- regional_df %>%
    mutate(uniqueID = paste("E", Species, "_P", patch, sep = "")) %>%
    filter(!is.na(Species))
}

## Assign states to each of the entities in each patch
assign.states <- function(regional_df) {
  regional_df <- regional_df %>%
    # Assign states (except for "Extinct")
    # these are numbers with higher = better
    mutate(state = case_when(
      Status == "Increasing (Least Concern)" ~ 7,
      Status == "Stable (Least Concern)" ~ 6,
      Status == "Unknown (Least Concern)" ~ 5,
      Status == "Decreasing (Least Concern)" ~ 4,
      Status == "Decreasing (Near threatened)" ~ 3,
      Status == "Decreasing (Endangered)" ~ 2,
      Status == "/" ~ 5 # !!! placeholder. The status being available is kind of like it being unknown?
    ))
}

## Create full data.frame of entities and their traits in patches with assigned states
get.full_df <- function(NumPatches, LevelOrder) {
  SpeciesOccs <- if (LevelOrder == 1) {
    read_csv("data/clean/species_occurences_L1.csv", show_col_types = FALSE) %>% 
      rename(Level = LEVEL1)
  } else if (LevelOrder == 2) {
    read_csv("data/clean/species_occurences_L2.csv", show_col_types = FALSE) %>% 
      rename(Level = LEVEL2)
  } else if (LevelOrder == 3) {
    read_csv("data/clean/species_occurences_L3.csv", show_col_types = FALSE) %>% 
      rename(Level = LEVEL3)
  } 
  
  # Assign levels to patches 
  
  # Equal numbers across each level (for now)
  Levels <- sample(SpeciesOccs$Level, NumPatches, replace = TRUE)
  
  Init_df <- tibble(Patch = 1:NumPatches, Level = Levels)
  
  # Assign number of species to each patch 
  
  # Maximum needs to be set to the total number of species at that level
  MaxNumEntities <- SpeciesOccs %>% 
    group_by(Level) %>% 
    summarise(count = n())
  
  
  # Assign species to patches -----------------------------------------------
  
  Patch_df <- tibble(Patch = as.integer(), Level = as.double(), Species = as.character())
  
  # For each patch, get its level
  for (index_patch in Init_df$Patch) {
    # Get the level of the patch
    level <- filter(Init_df, Patch == index_patch)$Level
    
    # Get species list from the right level
    species_list <- filter(SpeciesOccs, Level == level)$Species
    
    # Get number of entities in patch
    entity_count <- min(get.NumEntities(MaxNumEntities, level), length(species_list))
    
    # Sample entities without replacement from species list
    entities <- sample(species_list, entity_count, replace = FALSE)
    
    temp_df <- tibble(Patch = index_patch, Level = level, Species = entities)
    
    Patch_df <- rbind(Patch_df, temp_df)
    
  }
  
  # Add species traits to the dataframe
  full_df <- right_join(Patch_df, final_data, by = "Species") %>% 
    filter(!is.na(Patch), !is.na(Level))
}


##* Vulnerability function (NYI) -----------------------------------------------
# take traits, current state, and patch number and outputs future state

# simple function for now:
# if average of traits > 0.5, increase IUCN status
# if average of traits < -0.5, decrease IUCN status
# otherwise maintain status
get.vulnerability <- function(in_df) {
  out_df <- in_df %>%
    # group by all columns except traits
    group_by(across(c(-Trait, -value)), .drop = FALSE) %>%
    # compute average of traits
    # mutate(trait_avg = mean(value)) %>%
    # assign new states based on average of traits (in future, this will be a real function based on biology)
    mutate(state = state + case_when(
      state == 1 ~ ifelse(rnorm(1) > 3, 1, 0), # 2.326, 1, 0), # with small prob., extirpated species re-emerge
      trait_avg > 0.5 ~ 1,
      trait_avg < -0.5 ~ -1,
      TRUE ~ 0
    )) %>%
    # 6 is the highest IUCN level (Least concern), so don't go higher than 6
    mutate(state = ifelse(state > 7, 7, state)) %>%
    mutate(trait.val = case_when(
      trait.num == 1 ~ rnorm(1), # one trait changes randomly over time due to environmental changes, say
      TRUE ~ trait.val
    )) %>%
    select(-trait_avg) %>%
    ungroup()
}

##* Biodiversity metric functions------------------------------------------------
# for a patch, look at entities and their states & traits and output
# biodiversity measure of the patch

# number of unique entityIDs in a patch
num.unique.metric <- function(in_df) {
  out_df <- in_df %>%
    # filter(state > 1) %>% # remove extinct entities
    # remove duplicates due to rows from trait.num
    select(Patch, Species) %>%
    distinct() %>%
    # count number of distinct entities
    count(Patch, name = "biodiv")
}

# Number of endemic entities
num.endemic.metric <- function(in_df) {
  out_df <- in_df %>%
    # filter(state > 1) %>% # remove extinct entities
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
  # filter(n() == 1) %>%
  # ungroup() %>%
  # # count number of endemics in each Patch
  # count(Patch, name = "biodiv")
}

# General trait counting metric function
trait.count.metric <- function(in_df, trait_name) {
  out_df <- in_df %>%
    select(Patch, one_of(trait_name)) %>% 
    group_by(Patch) %>% 
    # Biodiversity measured as sum in a given trait
    summarise(biodiv = sum(!!as.name(trait_name)))
}

# Trait coefficient of variation function
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

# If we add abundance, we can compute: Shannon diversity, Evenness, Simpson's 
# index, rare species indices, ... 


# Calculate phylogenetic diversity metrics --------------------------------
phylodiv.metrics <- function(in_df) {
  # Pruning happens based on a community matrix, so either
  # - create a fake community matrix from list_species_unique 
  # - convert full_df to wide format
  # Attention: Species names should have a dash between genus and species name
  # KD: I think the goal is to have each row be a patch and each column a
  #     species, with a 1 if the species is in that patch. Is that correct?
  #     If so, we can use the code below!
  
  si <- in_df %>%
    mutate(Species2 = paste(stringr::word(Species, 1,1, sep=" "),
                            stringr::word(Species, 2,2, sep=" "),
                            sep="_")) %>%
    dplyr::select(Species2, Patch) %>% 
    unique()
  
  comm <- dcast(si, formula = Patch ~ Species2, fun.aggregate = length, 
                value.var = "Patch")
  
  ## ---- Calculate phylogenetic diversity --------------
  # this takes a while!
  
  # We have to remove species missing in the phylogeny in order for this to work
  comm <- comm %>%
    dplyr::select(-Pterospora_andromedea)
  
  tree_pruned <- prune.sample(comm,tree)
  
  # 1) sum of the total phylogenetic branch length (Faith, 1992)
  pdiv_length <- pd(comm, tree_pruned, include.root=TRUE) 
  pdiv_length2 <- pd(comm, tree_pruned, include.root=FALSE)
  # warnings because this cannot be calculated for communities with 1 species only
  # !!! BUG: when there is a single species community, pd with include.root = TRUE
  #          throws an error. Solution for now: prevent single species communities
  
  # 2) phylogenetic species variability richness and evenness (Helmus et al., 2007)
  pdiv_psv <- psv(comm, tree_pruned)
  pdiv_psr <- psr(comm, tree_pruned)
  pdiv_pse <- pse(comm, tree_pruned) 
  # pdiv_pse only makes sense when working with relative species abundances, which is not the case here.
  # I set the abundance of each species to 1
  
  # Combine results
  pdiv_all <- pdiv_length %>%
    tibble::rownames_to_column("Patch") %>%
    full_join(pdiv_length2 %>% dplyr::select(PD) %>%
                rename(PD_unrooted = PD) %>%
                tibble::rownames_to_column("Patch"),
              by = c("Patch")) %>%
    full_join(pdiv_psv %>% dplyr::select(PSVs)  %>%
                tibble::rownames_to_column("Patch"),
              by = c("Patch")) %>%
    full_join(pdiv_psr %>% dplyr::select(PSR)  %>%
                tibble::rownames_to_column("Patch"),
              by = c("Patch")) %>%
    full_join(pdiv_pse %>% dplyr::select(PSEs)  %>%
                tibble::rownames_to_column("Patch"),
              by = c("Patch"))
}


# Build biodiversity data frame -------------------------------------------
get.biodiv_df <- function(in_df) {
  # Metric 1: Species richness
  num.unique_df <- num.unique.metric(in_df)
  # hotspots.unique <- find.hotspots(num.unique_df)

  # Metric 2: Number of endemic entities
  num.endemic_df <- num.endemic.metric(in_df)
  # hotspots.endemic <- find.hotspots(num.endemic_df)

  # Metric 3: Number of Indigenous names
  num.indig.name_df <- trait.count.metric(in_df, "N_Names")
  # hotspots.indig.name <- find.hotspots(num.indig.name_df)

  # Metric 4: Number of Indigenous languages
  num.indig.lang_df <- trait.count.metric(in_df, "N_Langs")
  # hotspots.indig.lang <- find.hotspots(num.indig.lang_df)

  # Metric 5: Number of uses
  num.use_df <- trait.count.metric(in_df, "N_Uses")
  # hotspots.use <- find.hotspots(num.use_df)
  
  # Metric 6: Combined variation of quantitative traits
  coeffvar_df <- trait.coeffvariance.metric(in_df, c("LDMC (g/g)","Nmass (mg/g)", "Plant height (m)", "Leaf area (mm2)" ))
  
  # Phylogenetic diversity metrics
  # PD_df <- phylodiv.metrics(in_df)
  # PD_df$Patch <- as.integer(PD_df$Patch)
  
  # !!! BUG: Getting NA for a couple metrics
  # PD_df <- select(PD_df, -PSVs, -PSR, -PSEs)

  # Put together one big dataframe of biodiversity metrics of each Patch
  biodiv_df <- rename(num.unique_df, NumUnique = biodiv) %>% # Number of unique entities
    # Number of endemic entities
    right_join(rename(num.endemic_df, NumEndemic = biodiv), by = "Patch") %>%
    # Number of Indigenous names
    right_join(rename(num.indig.name_df, NumIndigName = biodiv), by = "Patch") %>%
    # Number of Indigenous languages
    right_join(rename(num.indig.lang_df, NumIndigLang = biodiv), by = "Patch") %>%
    # Number of uses
    right_join(rename(num.use_df, NumUse = biodiv), by = "Patch") %>% 
    # Coefficient of variation
    right_join(rename(coeffvar_df, CoeffVar = biodiv), by = "Patch")
  
  biodiv_df <- biodiv_df #%>% 
    # # add Phylogenetic diversity metrics
    # right_join(PD_df, by = "Patch")
}



# Hotspot identifier function ---------------------------------------------
# look at the biodiversity of all patches and determine which patches are
# "hotspots" relative to the other patches

# simple function for now:
# give all the patches in the 95% quantile
find.hotspots <- function(in_df) {
  # Select top 5th percentile
  hotspots <- in_df[in_df$biodiv > quantile(in_df$biodiv, 0.95, na.rm = FALSE), ]
  
  if (dim(hotspots)[1] == 0) {hotspots <- in_df[which.max(in_df$biodiv),]}
  
  return(hotspots)
}


# Hotspot comparison function ---------------------------------------------
# measure the difference in the hotspot maps produced under different
# biodiversity metrics for the same dataset

# simple function for now:
# count the number of overlaps in which patches were considered hotspots
calc.hotspot_compare <- function(hotspots.unique, hotspots.compare) {
  unique.hotspots <- hotspots.unique$Patch
  compare.hotspots <- hotspots.compare$Patch

  TP_count <- sum(compare.hotspots %in% unique.hotspots)

  # false_count <-

  # Use precision as our quantifier:
  # of the identified hotspots, what proportion match with species diversity?
  comparison.quantifier <- TP_count / length(compare.hotspots)
  
  # Use Jaccard similarity coefficient instead:
  #  = number of hotspots shared in both lists / total number of hotspots identified
  # measures the amount of overlap, without considering one list the "true" list
  # This is not substantially different from above - you just also include 
  # "false negatives" in the denominator. But may be more familiar to ecology
  # and biology folks.
  total.hotspots <- length(compare.hotspots)+length(unique.hotspots)-TP_count
  Jaccard.sim.coef <- TP_count / total.hotspots
  
  return(comparison.quantifier)
}


# Build hotspot comparison data frame ------------------------------------------
get.biodiv.compare_df <- function(in_df) {
  # Get biodiversity metrics
  biodiv_df <- get.biodiv_df(in_df) %>%
    melt(id = "Patch")

  # Pre-allocate table
  hotspot.compare_df <- tibble(
    variable = character(),
    value = numeric()
  )

  # Build hotspot comparison data frame
  for (j in unique(biodiv_df$variable)) {
    # Comparisons are made relative to species richness ("NumUnique")
    if (j == "NumUnique") {
      hotspots.unique <- biodiv_df %>%
        filter(variable == j) %>%
        select(-variable) %>%
        rename(biodiv = value) %>%
        find.hotspots()
      
    } else {
      hotspot_list <- biodiv_df %>%
        filter(variable == j) %>%
        select(-variable) %>%
        rename(biodiv = value) %>%
        mutate(biodiv = ifelse(is.na(biodiv), 0, biodiv)) %>% 
        find.hotspots()

      compare_value <- calc.hotspot_compare(hotspot_list, hotspots.unique)

      hotspot.compare_df <- add_row(hotspot.compare_df,
        variable = j,
        value = compare_value
      )
    }
  }
  return(hotspot.compare_df)
}
