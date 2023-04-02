##############################################################################
###                                                                        ###
###                 Calculating Phylogenetic diversity                     ###
###                                                                        ###
##############################################################################

# Elisa Van Cleemput, Nov 2022
#######################################################

# Important notes:
# 1)  Check 
# Populus balsamifera
#  --> the first record has spp. in the author column, is this an error?
# Some species have spp. or var. in the author column and others in the species column!
# Fix this!

# 2) I removed spp. and var. from the analyses

# 3) Pterospora_andromedea is not present in the tree, so delete from other analyses

# 4) I am not sure what the level is at which you calculate diversity. 
# Below I used uniqueID, but they have only 1 species so this doesn't seem right
# Within 1 patch there are sometimes duplicates of species, so I am not sure if that is correct
#  KD: Thanks for pointing this out! I fixed it

#######################################################

# Data used in this script:
# - full_data or data_in to create phylogenetic tree
# - full_df to calculate phylogenetic diversity for patches
source("code/data_intake.R") # full_data or data_in
source("code/functions.R") # to create community data frame
source("code/parameters.R") # to get parameters for creating community data frame
# KD: don't need to run "simulations.R" anymore. full_df is created in this script


## ---- Prepare species list to create a phylogenetic tree --------------
# --> phylo.maker needs a dataframe with three columns
# 1) species: without authors
# 2) genus
# 3) family

library(stringr)

# full_data <- base_data

# list_species_unique <- full_data %>%
#   as.data.frame()  %>%
#   dplyr::mutate(Species = stringr::word(Species, 1,2, sep=" ")) %>% # get rid of spp. and var.
#   dplyr::mutate(Genus = stringr::word(Species, 1,1, sep=" ")) %>%
#   dplyr::select(Species, Genus, Family) %>%
#   # group_by(Species) %>% # uncomment to check for how many species are "repeated"
#   # filter(n()>1)
#   unique()

list_species_unique <- data_in %>%
  dplyr::select(Species_full, Genus, Family) %>%
  # group_by(Species) %>% # uncomment to check for how many species are "repeated"
  # filter(n()>1)
  unique()

# KD: "Alnus viridis" and "Populus balsamifera" have duplicate entries because
#     they are listed by their subspecies in the original data set:
#     Alnus viridis ssp. crispa and Alnus viridis ssp. sinuata and
#     Populus balsamifera and Populus balsamifera ssp. trichocarpa
#     *** Should we trim off the "ssp." part of the species name before    ***
#     *** creating the full_data dataframe?                                ***
# EVC: These duplicates have the same nb of unique names and languages, supposedly they are the same?
#     They only differ in number of uses.
#     I am not sure if the tree has spp. and var. available.
#     I could check this, but first would have to clean up the dataset, i.e. add spp. to species names
#     in cases where this is not done already. 
#     E.g. Populus balsamifera should be Populus balsamifera ssp. balsamifera
#     It is probably also worthwile checking the trait database information to see if they provide
#     trait infomation for spp. or not


# dim(list_species_unique)
# head(list_species_unique)

## ---- Generate a phylogeny for the species list --------------
# library("devtools")
# devtools::install_github("jinyizju/V.PhyloMaker2")
library("V.PhyloMaker2")

# The following function takes +- 1 minute to run
system.time(
  tree.S3 <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
                                        tree = GBOTB.extended.WP,
                                        nodes = nodes.info.1.WP, 
                                        output.tree = TRUE,
                                        scenarios = "S3")
)
# The functions tells us whether all species are present on the tree or not
# Warnings:
# - "Taxonomic classification not consistent between sp.list and tree."
#   e.g.,  genus family_in_sp.list family_in_tree
#          Acer         Aceraceae    Sapindaceae
# --> which information does the algorithm use in these cases? I don't know...
# - Note: 1 taxa fail to be binded to the tree
#   "Pterospora_andromedea"

tree <- tree.S3

# Add missing species to the tree
if (!require('pez')) install.packages('pez'); library('pez')
missing_species <- c("Pterospora_andromedea")
tree_complete <- congeneric.merge(tree$tree.scenario.3, missing_species)
# I don't think this worked. Is this because the genus is not present in the tree?
(tree$tree.scenario.3$Nnode)
(tree_complete$Nnode)
which(str_detect(tree_complete$tip.label, "Pterospora_andromedea") == TRUE)
# --> it is not present in the tree


# save tree
write.tree(tree$tree.scenario.3, file = "data/clean/phylogenetic_tree.csv")

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
