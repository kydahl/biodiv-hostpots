################################################################################
# Functions used in the simulation study
################################################################################

##* Load libraries---------------------------------------------------------------
library(tidyverse)

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

## Assign entities to patches
get.patches <- function(entity_df,
                        numPatches,
                        mean.NumEntities,
                        sd.NumEntities) {
  # initialize data.frame
  regional_df <- tibble(
    entityID = NA,
    patch = NA
  )

  # sample from entities list to fill patches
  for (i in 1:numPatches) {
    regional_df <- add_row(regional_df,
      patch = i,
      entityID = group_assign(
        unique(entity_df$entityID),
        mean.NumEntities,
        sd.NumEntities
      )
    )
  }

  # assign uniqueIDs to entities in patches (and remove NAs we added in initialization)
  regional_df <- regional_df %>%
    mutate(uniqueID = paste("E", entityID, "_P", patch, sep = "")) %>%
    filter(!is.na(entityID))
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
get.full_df <- function(entity_df, numPatches, mean.NumEntities, sd.NumEntities) {
  # assign patches and states to entities
  regional_df <- get.patches(
    entity_df,
    numPatches,
    mean.NumEntities,
    sd.NumEntities
  )

  states_df <- full_join(entity_df, regional_df, by = "entityID") %>%
    filter(!is.na(patch))

  full_df <- states_df %>%
    # assign states to entities within patches%>%
    # Assign states (except for "Extinct")
    # these are numbers with higher = better
    mutate(state = ifelse(test = (Status == "/"),
      # !!! place holder. if we don't know status, just assign randomly
      yes = sample(
        x = 2:7,
        size = dim(filter(states_df, Status == "/"))[1],
        replace = TRUE
      ),
      no = case_when(
        Status == "Increasing (Least Concern)" ~ 7,
        Status == "Stable (Least Concern)" ~ 6,
        Status == "Unknown (Least Concern)" ~ 5,
        Status == "Decreasing (Least Concern)" ~ 4,
        Status == "Decreasing (Near threatened)" ~ 3,
        Status == "Decreasing (Endangered)" ~ 2,
        # Status == "/" ~ sample.int(6, 1, replace = TRUE) # !!! placeholder. just assign something random if it's blank...
      )
    )) %>%
    # re-order the data.frame for easier viewing
    relocate(entityID, patch, uniqueID)

  # states_df <- entity_df %>%
  #   # assign entities to patches
  #   get.patches(numPatches, maxEntitiesPatch) %>%
  #   # assign states to entities within patches
  #   assign.states(.)

  # add trait values back in
  # full_df <- inner_join(entity_df, states_df) %>%
  #   # re-order the data.frame for easier viewing
  #   relocate(entityID, patch, uniqueID)
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
    filter(state > 1) %>% # remove extinct entities
    # remove duplicates due to rows from trait.num
    select(patch, entityID) %>%
    distinct() %>%
    # count number of distinct entities
    count(patch, name = "biodiv")
}

# Number of endemic entities
num.endemic.metric <- function(in_df) {
  out_df <- in_df %>%
    filter(state > 1) %>% # remove extinct entities
    select(patch, entityID) %>%
    # only keep rows where entityID is unique ( == endemics)
    group_by(entityID) %>%
    # remove duplicates due to rows from trait.num
    distinct() %>%
    mutate(endemic_flag = ifelse(n() == 1, 1, 0)) %>%
    ungroup() %>%
    group_by(patch) %>%
    summarise(
      "biodiv" = sum(endemic_flag, na.rm = TRUE),
      .groups = "keep"
    )
  # filter(n() == 1) %>%
  # ungroup() %>%
  # # count number of endemics in each patch
  # count(patch, name = "biodiv")
}

# General trait counting metric function
trait.count.metric <- function(in_df, trait_name) {
  out_df <- in_df %>%
    # remove extinct entities
    filter(state > 1) %>%
    filter(Trait == trait_name) %>%
    # KD: This calculates TOTAL number of traits (which perhaps makes sense for
    #     uses or names, but probably not for biological traits)
    select(patch, entityID, value) %>%
    group_by(patch) %>%
    summarise(
      "biodiv" = sum(value, na.rm = TRUE),
      .groups = "keep"
    )
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
    dplyr::select(Species2, patch) %>% 
    unique()
  
  comm <- dcast(si, formula = patch ~ Species2, fun.aggregate = length, 
                value.var = "patch")
  
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
    tibble::rownames_to_column("patch") %>%
    full_join(pdiv_length2 %>% dplyr::select(PD) %>%
                rename(PD_unrooted = PD) %>%
                tibble::rownames_to_column("patch"),
              by = c("patch")) %>%
    full_join(pdiv_psv %>% dplyr::select(PSVs)  %>%
                tibble::rownames_to_column("patch"),
              by = c("patch")) %>%
    full_join(pdiv_psr %>% dplyr::select(PSR)  %>%
                tibble::rownames_to_column("patch"),
              by = c("patch")) %>%
    full_join(pdiv_pse %>% dplyr::select(PSEs)  %>%
                tibble::rownames_to_column("patch"),
              by = c("patch"))
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
  num.indig.name_df <- trait.count.metric(
    in_df,
    "Number of unique names (Indigenous)"
  )
  # hotspots.indig.name <- find.hotspots(num.indig.name_df)

  # Metric 4: Number of Indigenous languages
  num.indig.lang_df <- trait.count.metric(
    in_df,
    "Number of unique  languages (from Appendix 2B)"
  )
  # hotspots.indig.lang <- find.hotspots(num.indig.lang_df)

  # Metric 5: Number of uses
  num.use_df <- trait.count.metric(
    in_df,
    "Number of uses"
  )
  # hotspots.use <- find.hotspots(num.use_df)
  
  # Phylogenetic diversity metrics
  PD_df <- phylodiv.metrics(in_df)
  PD_df$patch <- as.integer(PD_df$patch)
  
  # !!! BUG: Getting NA for a couple metrics
  PD_df <- select(PD_df, -PSVs, -PSR, -PSEs)

  # Put together one big dataframe of biodiversity metrics of each patch
  biodiv_df <- rename(num.unique_df, NumUnique = biodiv) %>% # Number of unique entities
    # Number of endemic entities
    right_join(rename(num.endemic_df, NumEndemic = biodiv), by = "patch") %>%
    # Number of Indigenous names
    right_join(rename(num.indig.name_df, NumIndigName = biodiv), by = "patch") %>%
    # Number of Indigenous languages
    right_join(rename(num.indig.lang_df, NumIndigLang = biodiv), by = "patch") %>%
    # Number of uses
    right_join(rename(num.use_df, NumUse = biodiv), by = "patch") 
  
  biodiv_df <- biodiv_df %>% 
    # add Phylogenetic diversity metrics
    right_join(PD_df, by = "patch")
}



# Hotspot identifier function ---------------------------------------------
# look at the biodiversity of all patches and determine which patches are
# "hotspots" relative to the other patches

# simple function for now:
# give all the patches in the 95% quantile
find.hotspots <- function(in_df) {
  # Select top 5th percentile
  in_df[in_df$biodiv > quantile(in_df$biodiv, 0.95, na.rm = FALSE), ]
}


# Hotspot comparison function ---------------------------------------------
# measure the difference in the hotspot maps produced under different
# biodiversity metrics for the same dataset

# simple function for now:
# count the number of overlaps in which patches were considered hotspots
calc.hotspot_compare <- function(hotspots.unique, hotspots.compare) {
  unique.hotspots <- hotspots.unique$patch
  compare.hotspots <- hotspots.compare$patch

  TP_count <- sum(compare.hotspots %in% unique.hotspots)

  # false_count <-

  # Use precision as our quantifier:
  # of the identified hotspots, what proportion match with species diversity?
  comparison.quantifier <- TP_count / length(compare.hotspots)
  return(comparison.quantifier)
}


# Build hotspot comparison data frame ------------------------------------------
get.biodiv.compare_df <- function(in_df) {
  # Get biodiversity metrics
  biodiv_df <- get.biodiv_df(in_df) %>%
    melt(id = "patch")

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
