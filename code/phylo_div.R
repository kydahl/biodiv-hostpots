##############################################################################
###                                                                        ###
###                 Calculating Phylogenetic diversity                     ###
###                                                                        ###
##############################################################################

# Elisa Van Cleemput, Nov 2022
#######################################################

# Important notes:

# 1) I removed spp. and var. from the analyses

# 1. Restrict species we're looking for occurrences for to big list of species names
# 2. Find phylogenetic data for these species
# 3. Build a tree
# 3. Then we need to go from the big list of species to our "synonymized" list
# 4. Decide how to deal with possible conflicts


#######################################################

## ---- Prepare species list to create a phylogenetic tree --------------
# --> phylo.maker needs a dataframe with three columns
# 1) species: without authors
# 2) genus
# 3) family

library(stringr)

list_species_unique <- full_df %>%
  mutate(Genus = stringr::word(Species, 1, 1, sep = " ")) %>%
  dplyr::select(Species, Genus, Family) %>%
  unique()
# dim(list_species_unique)
# head(list_species_unique)

# !!!! Pterospera andromedea is still in there? I thought it would have been removed.
# It is part of the species occurrrence list, so maybe that's the reason why it's in there?
# This list shoudl hence be updated!

# I did not trim the subsp. and var. from the species names and that doesn't seem to be a problem

## ---- Generate a phylogeny for the species list --------------
# library("devtools")
# devtools::install_github("jinyizju/V.PhyloMaker2")
library("V.PhyloMaker2")

# The following function takes +- 1 minute to run
system.time(
  tree.WP.S3 <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
                                        tree = GBOTB.extended.WP, # botanical nomenclature of the World Plants database
                                        nodes = nodes.info.1.WP, 
                                        output.tree = TRUE,
                                        scenarios = "S3")
) # 7 taxa fail to bind:
# "Zigadenus_venenosus"      "Zigadenus_elegans"        "Glaux_maritima"           "Streptopus_amplexifolius"  "Rumex_aquaticus"  "Pterospora_andromedea"    "Frasera_montana"
system.time(
  tree.LCVP.S3 <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
                                           tree = GBOTB.extended.LCVP, # botanical nomenclature of the Leipzig catalogue of vascular plants
                                           nodes = nodes.info.1.LCVP, 
                                           output.tree = TRUE,
                                           scenarios = "S3")
) # 6 taxa fail to bind
# Same species, except for Glaux_maritima
system.time(
  tree.TPL.S3 <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
                                             tree = GBOTB.extended.TPL, # botanical nomenclature of The Plant List (static since 2013)
                                             nodes = nodes.info.1.TPL, 
                                             output.tree = TRUE,
                                             scenarios = "S3")
) # 7 taxa fail to bind:
# Same species

# Before, I only had a problem for 1 species!!!!

# The functions tells us whether all species are present on the tree or not
# Warnings:
# - "Taxonomic classification not consistent between sp.list and tree."
#   e.g.,  genus family_in_sp.list family_in_tree
#          Abies            Abies       Pinaceae
#           Acer        Aceraceae    Sapindaceae
# --> which information does the algorithm use in these cases? I don't know...


tree <- tree.LCVP.S3

## ---- Add missing species to the tree --------------
if (!require('pez')) install.packages('pez'); library('pez')
# missing_species <- c("Pterospora_andromedea")
missing_species <- c("Zigadenus_venenosus", "Zigadenus_elegans", "Glaux_maritima", 
                     "Streptopus_amplexifolius", "Rumex_aquaticus", "Pterospora_andromedea", "Frasera_montana")

tree.WP.S3_complete <- congeneric.merge(tree.WP.S3$tree.scenario.3, missing_species)
tree.LCVP.S3_complete <- congeneric.merge(tree.LCVP.S3$tree.scenario.3, missing_species)
tree.TPL.S3_complete <- congeneric.merge(tree.TPL.S3$tree.scenario.3, missing_species)
# I don't think this worked. Is this because the genus is not present in the tree?
(tree.LCVP.S3$tree.scenario.3$Nnode)
(tree.LCVP.S3_complete$Nnode)
which(str_detect(tree.LCVP.S3_complete$tip.label, "Frasera_montana") == TRUE)
# --> Rumex_aquaticus, Pterospora_andromedea, Frasera_montana is not present in the tree


# save tree
# write.tree(tree.WP.S3_complete$tree.scenario.3, file = "data/clean/phylogenetic_tree.csv")

# Option 1: Consider changing species names before adding them to the tree with congeneric functions? 
custom_synonyms_list <- read_csv("data/clean/synonym_list.csv", show_col_types = FALSE)
Diaz_synonyms_list <- read_csv("data/clean/Diaz_synonyms.csv", show_col_types = FALSE)
missing_species_custom <- missing_species %>%
  as.data.frame() %>%
  rename("Species" = ".") %>%
  mutate(original_name = paste(stringr::word(Species, 1, 1, sep = "_"), 
                               stringr::word(Species, 2, 2, sep = "_"), 
                               sep = " ")) %>%
  left_join(custom_synonyms_list, by = "original_name") %>%
  left_join(Diaz_synonyms_list %>% rename(original_name = name_in), by = "original_name")
tree.WP.S3_complete <- congeneric.merge(tree.WP.S3$tree.scenario.3, missing_species_custom$Synonym)
# This didn't solve the issue

# Option 2: Consider changing species names before creating tree? 
list_species_unique2 <- list_species_unique %>%
  as.data.frame() %>%
  select(-Genus) %>%
  rename(original_name = Species) %>%
  # left_join(custom_synonyms_list %>% rename(Species = original_name), by = "Species")
  left_join(missing_species_custom %>% 
              select(original_name, Synonym), by = "original_name") %>%
  mutate(Species = coalesce(Synonym, original_name)) %>%
  mutate(Genus = stringr::word(Species, 1, 1, sep = " ")) %>%
  dplyr::select(Species, Genus, Family) %>%
  unique()
# !!!! Family equals Genus which should not be the case!!! Something wrong in the synonym_func function in the “data_intake.R”!


system.time(
  tree.LCVP.S3.test <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique2, 
                                             tree = GBOTB.extended.LCVP, # botanical nomenclature of the Leipzig catalogue of vascular plants
                                             nodes = nodes.info.1.LCVP, 
                                             output.tree = TRUE,
                                             scenarios = "S3")
) # 4 taxa fail to bind
#"Streptopus_amplexifolius" "Rumex_aquaticus"          "Pterospora_andromedea"    "Frasera_montana" 


## ---- Prune tree to species list --------------
if (!require('picante')) install.packages('picante'); library('picante')
# An introduction to the picante package:
# https://cran.r-project.org/web/packages/picante/vignettes/picante-intro.pdf
# Other package and function: https://search.r-project.org/CRAN/refmans/abdiv/html/faith_pd.html

# Create our simulated patch communities
full_df <- data_in %>%  
  filter(Species != "Pterospora andromedea") %>% 
  get.full_df(numPatches, mean.NumEntities, sd.NumEntities)

# Pruning happens based on a community matrix, so either
# - create a fake community matrix from list_species_unique 
# - convert full_df to wide format
# Attention: Species names should have a dash between genus and species name
# KD: I think the goal is to have each row be a patch and each column a
#     species, with a 1 if the species is in that patch. Is that correct?
#     If so, we can use the code below!
# EVC: yes, thanks for improving the code! I added a line turning the first column in rownames
si <- full_df %>%
  mutate(Species2 = paste(stringr::word(Species, 1,1, sep=" "),
                          stringr::word(Species, 2,2, sep=" "),
                          sep="_")) %>%
  mutate(Species_full2 = str_replace_all(Species_full," ","_")) %>%
  dplyr::select(Species_full2, patch) %>% 
  unique()

comm <- dcast(si, formula = patch ~ Species_full2, fun.aggregate = length) %>%
    column_to_rownames(var="patch")

# comm <- full_df %>%
#   mutate(Species2 = paste(stringr::word(Species, 1,1, sep=" "),
#                           stringr::word(Species, 2,2, sep=" "),
#                           sep="_")) %>%
#   pivot_wider(names_from = Species2, values_from = value) %>%
#   dplyr::select(-Species, -entityID, -patch, -Status, -Trait, -state) %>%
#   column_to_rownames(var="uniqueID")
# # --> the values do not make sense here, but that is okay
# comm[!is.na(comm)] <- 1 # replace values with 1
# comm[is.na(comm)] <- 0 


# I am not sure  what the communities' ID is: uniqueID or patch.
# KD: Community ID should be the patch number. UniqueID is just the way of 
#     referring to a specific species in a specific patch (eg to differentiate
#     between Abies amabilis in patches 3 and 300)
# When using patch, there are duplicate species present in a patch,
# and that causes problems when pivoting to a longer format
# comm <- full_df %>%
#   dplyr::mutate(Species2 = paste(stringr::word(Species, 1,1, sep=" "),
#                           stringr::word(Species, 2,2, sep=" "),
#                           sep="_")) %>%
#   pivot_wider(id_cols=patch, names_from = Species2, values_from = value) %>%
#   column_to_rownames(var="patch")

# KD: This looks for duplicates in the full dataframe. We only get 
#     Alnus_viridis and Populus_balsamifera now, as expected (see above)
duplicates_df <- full_df %>%
  dplyr::mutate(Species2 = paste(stringr::word(Species, 1,1, sep=" "),
                          stringr::word(Species, 2,2, sep=" "),
                          sep="_")) %>%
  mutate(Species_full2 = str_replace_all(Species_full," ","_")) %>%
  dplyr::group_by(patch, Species_full2) %>% # EVC: there are no duplicates when using Species_full2 instead of Species2
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) %>% 
  # There will be repeats for each unique entry for "Trait", divide by that
  # number to get the actual number of duplicates
  dplyr::mutate(actual_num = n / length(unique(full_df$Trait))) %>% 
  filter(actual_num > 1)


tree_complete_pruned <- prune.sample(comm,tree_complete)
tree_pruned <- prune.sample(comm,tree$tree.scenario.3)


# check if number of species corresponds with number of species in the species list
length(tree_complete_pruned$tip.label)
length(tree_pruned$tip.label)
dim(list_species_unique)
# As we know, one species is missing... 

# Check out the species list --> subspecies and variants are present, so this worked?
tree_pruned$tip.label

### plot the phylogenies with node ages displayed.
plot.phylo(tree_pruned, type="fan", cex = 0.5)
plot.phylo(tree_pruned, type="tidy", cex = 0.5)


## ---- Calculate phylogentic diversity --------------
# this takes a while!

# We have to remove species missing in the phylogeny in order for this to work
# In the updated version, this is already done in "simulations.R"
# comm <- comm %>%
#   dplyr::select(-Pterospora_andromedea)


# 1) sum of the total phylogenetic branch length (Faith, 1992)
pdiv_length <- pd(comm, tree_pruned, include.root=TRUE) 
pdiv_length2 <- pd(comm, tree_pruned, include.root=FALSE)
# warnings because this cannot be calculated for communities with 1 species only
# KD: Warnings resolved! I think because I dealt with the duplicate issues
# EVC: The second calculation still throws warnings, but that is normal because 
#     there are patches with 1 species only. In these cases PD = NA

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



plot(pdiv_all[,"PD"], pdiv_all[,"SR"])
plot(pdiv_all[,"PSVs"], pdiv_all[,"SR"])
plot(pdiv_all[,"PSR"], pdiv_all[,"SR"])
plot(pdiv_all[,"PSR"], pdiv_all[,"PD"])
# plot(pdiv_all[,"PSEs"], pdiv_all[,"SR"])
