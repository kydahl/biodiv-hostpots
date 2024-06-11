# biodiv-hostpots
Code for "Comparisons of traditional ecological knowledge, taxonomic, functional, and phylogenetic biodiversity metrics reveal dissimilarities in biodiversity hotspot identification"

The scripts below are listed in order of dependency: run the first scripts before the latter ones.

Run `data_intake.R` to load all the necessary data into an appropriate format for analysis. This script also includes the imputation of missing trait data.

`functions.R` contains all the special functions necessary for the study.

`species_occurrences.R` calculates the occurrence ecoregions in the PNW then the occurrence of PNW species within those ecoregions.

`simulations.R` runs the actual simulations. This script is likely to use a large amount of memory and time to complete.

`visualizations.R` produces the figures and tables seen in the manuscript.