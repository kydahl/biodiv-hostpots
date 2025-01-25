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
final_df <- read_csv("data/clean/dataset_no_imputation.csv")

list_species_unique <- final_df %>%
  mutate(Genus = stringr::word(Synonym, 1, 1, sep = " ")) %>%
  # dplyr::select(Species, Genus, Family) %>%
  dplyr::select(Synonym, Genus, Family) %>%
  # mutate(Synonym = str_replace_all(Synonym, "ssp.", "subsp.")) %>%
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


write_rds(tree, "data/clean/final_tree.rds")

