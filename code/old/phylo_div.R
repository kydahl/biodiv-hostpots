##############################################################################
###                                                                        ###
###                 Calculating Phylogenetic diversity                     ###
###                                                                        ###
##############################################################################

# Elisa Van Cleemput, Nov 2022
#######################################################
## ---- Prepare species list to create a phylogenetic tree --------------
# --> phylo.maker needs a dataframe with three columns
# 1) species: without authors
# 2) genus
# 3) family

library(stringr)

list_species_unique <- full_df %>%
  mutate(Genus = stringr::word(Synonym, 1, 1, sep = " ")) %>%
  # dplyr::select(Species, Genus, Family) %>%
  dplyr::select(Synonym, Genus, Family) %>%
  mutate(Synonym = str_replace_all(Synonym, "ssp.", "subsp.")) %>%
  unique() %>%
  as.data.frame()
# dim(list_species_unique)
# View(list_species_unique)

# I did not trim the subsp. and var. from the species names and that doesn't seem to be a problem

## ---- Generate a phylogeny for the species list --------------
# library("devtools")
# devtools::install_github("jinyizju/V.PhyloMaker2")
library("V.PhyloMaker2")

# The following functions takes +- 1 minute to run
# Choose one nomenclature to work with
# system.time(
#   tree.WP.S3 <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
#                                         tree = GBOTB.extended.WP, # botanical nomenclature of the World Plants database
#                                         nodes = nodes.info.1.WP, 
#                                         output.tree = TRUE,
#                                         scenarios = "S3")
# ) 
system.time(
  tree.LCVP.S3 <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
                                           tree = GBOTB.extended.LCVP, # botanical nomenclature of the Leipzig catalogue of vascular plants
                                           nodes = nodes.info.1.LCVP, 
                                           output.tree = TRUE,
                                           scenarios = "S3")
) 
# system.time(
#   tree.TPL.S3 <- V.PhyloMaker2::phylo.maker(sp.list = list_species_unique, 
#                                              tree = GBOTB.extended.TPL, # botanical nomenclature of The Plant List (static since 2013)
#                                              nodes = nodes.info.1.TPL, 
#                                              output.tree = TRUE,
#                                              scenarios = "S3")
# ) 
# The functions tells us whether all species are present on the tree or not
# Warnings:
# - "Taxonomic classification not consistent between sp.list and tree."
#   e.g.,  genus family_in_sp.list family_in_tree
#           Acer        Aceraceae    Sapindaceae
# --> which information does the algorithm use in these cases? I don't know...


tree <- tree.LCVP.S3

## ---- Prune tree to species list --------------
if (!require('picante')) install.packages('picante'); library('picante')
# An introduction to the picante package:
# https://cran.r-project.org/web/packages/picante/vignettes/picante-intro.pdf
# Other package and function: https://search.r-project.org/CRAN/refmans/abdiv/html/faith_pd.html

# Pruning happens based on a community matrix, so convert full_df to wide format
# Attention: Species names should have a dash between genus and species name
full_df_for_comm <- full_df %>%
  mutate(Species_full = str_replace_all(Species_full, "ssp.", "subsp.")) %>%
  mutate(Species_full2 = str_replace_all(Species_full," ","_")) %>%
  dplyr::select(Species_full2, Patch) %>% 
  unique()

if (!require('maditr')) install.packages('maditr'); library('maditr')
comm <- dcast(full_df_for_comm, formula = Patch ~ Species_full2, fun.aggregate = length) %>%
    column_to_rownames(var="Patch")

tree_pruned <- prune.sample(comm,tree$tree.scenario.3)


# check if number of species corresponds with number of species in the species list
print(paste0("Number of tips in original tree is ", length(tree$tree.scenario.3$tip.label)))
print(paste0("Number of tips in pruned tree is ", length(tree_pruned$tip.label)))
print(paste0("Number of species in specieslist is ", nrow(list_species_unique)))


# Check out the species list --> subspecies and variants are present, so this seemsto have worked
# tree_pruned$tip.label

### Optional: plot the phylogenies with node ages displayed.
plot.phylo(tree_pruned, type="fan", cex = 0.5)
plot.phylo(tree_pruned, type="tidy", cex = 0.5)



