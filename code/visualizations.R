################################################################################
# Visualization of simulation results
################################################################################

## Project: Comparing biodiversity hotspot identification
##
## Purpose: Produce figures illustrating features of simulation results
##
## Contents: Set-up, load in necessary packages, functions, and data-sets
##           Figure 0: Distribution of traits (across all species)
##           Figure 1: Compare the distribution of biodiversity metric values
##           Table 1: Numbers of hot spots identified
##           Figure ?: Biodiversity metric dendrogram clustered by precision
##           Figure S3: Pairwise comparison heatmaps
##           Figure S?: Joint density of TEK data
##            
## Inputs:  - Source: code/functions.R
##          - data/clean/final_dataset.csv
##          - data/clean/full_tree.rds
##          
## Outputs: - results/full_comparisons.rds
##
## Written and by: Kyle Dahlin and Elisa Van Cleemput
## Maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized April 2023
## _____________________________________________________________________________


# Set-up, load in necessary packages, functions, and data-sets ------------
#### Load libraries ----
require(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(GGally) # used to make paired scatter plots
library(retry)
library(cols4all)
library(gridExtra)
library(ggh4x)
library(scales)
library(ggpubr)
library(ggdendro) # to plot dendrograms

# Load in functions
source("code/functions.R")

# Load in dataset and phylogenetic tree
final_data <- read_csv("data/clean/final_dataset.csv") %>% #read_csv("data/clean/final_dataset.csv") %>% 
  # log transform traits with large outliers
  mutate(LeafArea_log = log(`Leaf area (mm2)`), .keep = 'unused') %>%
  mutate(PlantHeight_log = log(`Plant height (m)`), .keep = 'unused') %>%
  mutate(DiasporeMass_log = log(`Diaspore mass (mg)`), .keep = 'unused') %>%
  mutate(LDMC_log = log(`LDMC (g/g)`), .keep = 'unused') %>%
  # Turn categorical traits into quantitative ones
  mutate(Woodiness = ifelse(Woodiness == "woody", 1, 0)) %>% 
  # Put traits at the end
  relocate(c(`Nmass (mg/g)`, Woodiness, LeafArea_log:LDMC_log), .after = last_col())

# Assign colors to biodiversity metrics according to type
x_label_color_function <- function(string) {
  out <- case_match(string,
                    "Taxonomic" ~ "black",
                    "TEK" ~ c4a("brewer.set1", 3)[1],
                    "Phylogenetic" ~ c4a("brewer.set1", 3)[2],
                    "Functional" ~ c4a("brewer.set1", 3)[3],)
  out <- as.character(out)
}

# # Figure 1: Compare the distribution of biodiversity metric values --------
# # Define the relevant trait names
# trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", 
#                  "Leaf area (mm2)", "Diaspore mass (mg)")
# 
# # Initialize the full phylogenetic tree
# tree <- readRDS("data/clean/full_tree.rds")
# 
# 
# 
# 
# # Get a single simulation
# numPatches <- 1000
# full_df_sample <- get.full_df(numPatches)
# 
# # Calculate biodiversity metrics
# biodiv_df <- get.biodiv_df(full_df_sample, tree)
# 
# # Set up label colors to indicate what type of index it is
# biodiv_plot_df <- biodiv_df %>% 
#   pivot_longer(cols = -Patch) %>% 
#   filter(!is.na(value)) %>%
#   add_row(Patch = 1, name = "dummy")
# 
# biodiv_plot_df$group <- case_match(biodiv_plot_df$name,
#                                    c("NumUnique", "NumEndemic") ~ "Taxonomic",
#                                    c("NumIndigName", "NumUse") ~ "TEK",
#                                    c("FD", "FRic", "FDis") ~ "Functional",
#                                    c("PD","richness", "PSVs", "PSR") ~ "Phylogenetic"
#                                    # c("richness", "GiniSimpson", "Simpson",
#                                    # "Shannon", "Margalef", "Menhinick",
#                                    # "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
#                                    # c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"         
# )
# 
# 
# biodiv_plot_df$color <- x_label_color_function(biodiv_plot_df$group)
# 
# # Make labels more descriptive
# biodiv_plot_df$label <-  case_match(biodiv_plot_df$name,
#                                     "NumUnique" ~ "Species richness",
#                                     "NumEndemic" ~ "Number of 'endemic' species",
#                                     "NumIndigName" ~ "Indigenous names (TEK)",
#                                     "NumUse" ~ "Recorded Indigenous uses (TEK)",
#                                     "FD" ~ "Rao's entropy (FD)",
#                                     "PD" ~ "Rao's entropy (PD)",
#                                     "FRic" ~ "Species richness (FD)",
#                                     # "FDiv" ~ "Functional sp. divergence",
#                                     "FDis" ~ "Species dispersion (FD)",
#                                     # "FEve" ~ "Functional sp. evenness",
#                                     # "Q" ~ "Rao's entropy (Q)",
#                                     "richness" ~ "Faith's (PD)",
#                                     # "GiniSimpson" ~ "Gini and Simpson",
#                                     # "Simpson" ~ "Simpson",
#                                     # "Shannon" ~ "Shannon",
#                                     # "Margalef" ~ "Margalef",
#                                     # "Menhinick" ~ "Menhinick",
#                                     # "McIntosh" ~ "McIntosh",
#                                     "PSVs" ~ "Species variability (PD)",
#                                     "PSR" ~ "Species richness (PD)",
#                                     "dummy" ~ ""
# )
# 
# biodiv_plot_df$label <-  factor(biodiv_plot_df$label, levels = c(
#   "Species richness", "Indigenous names (TEK)",   "Recorded Indigenous uses (TEK)",
#   "Number of 'endemic' species", 
#   "Rao's entropy (FD)", "Species richness (FD)", "Species dispersion (FD)", 
#   "Rao's entropy (PD)", "Faith's (PD)", "Species variability (PD)", "Species richness (PD)"
#   
#   # "Functional sp. divergence", 
#   # 
#   # "Functional sp. evenness", "Rao's entropy (Q)", "Gini and Simpson", 
#   # "Margalef", "McIntosh", "Menhinick",  "Shannon", 
#   # "Simpson", ""
# )
# )
# 
# # Get mean and standard deviation values for each index
# meanSD_df <- biodiv_plot_df %>% 
#   filter(label != "Number of 'endemic' species") %>%
#   group_by(label) %>% 
#   summarise(mean = mean(value, na.rm = TRUE),
#             stdev = sd(value, na.rm = TRUE),
#             upper = mean + 1.96 * stdev,
#             hotspot_cutoff = quantile(value, c(0.95), na.rm = TRUE)
#   )
# 
# # Plot histograms in three columns
# # histogram_colors = c("black", rep("#E41A1C",2), rep("#4DAF4A", 5), rep("#377EB8", 9), NA)
# histogram_colors = c("black", rep("#E41A1C",2), rep("#4DAF4A", 3), rep("#377EB8", 4), NA)
# 
# strip_theme = strip_themed(background_x = elem_list_rect(fill = histogram_colors),
#                            text_x = element_text(color = "white",face = "bold", size = 8))
# 
# biodiv_plot_df2 <- right_join(biodiv_plot_df, meanSD_df, by = c("label")) %>% 
#   distinct()
# 
# histogram_plot <- biodiv_plot_df2 %>%
#   filter(label != "Number of 'endemic' species") %>%
#   ggplot(aes(x = value)) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value, ..scaled..),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   # add a vertical line showing where the top 5% are
#   stat_summary(aes(xintercept = after_stat(x), y = 0),
#                fun = quantile, fun.args = list(0.95),
#                geom = "vline", orientation = "y",
#                color = "blue", lwd = 1) +
#   geom_vline(aes(xintercept = mean), linetype = 2) +
#   ylab('Density') +
#   xlab('Value') +
#   # one histogram for each biodiversity metric
#   facet_wrap2( ~ label, ncol = 3, scales = "free",
#                strip = strip_theme) +
#   theme_cowplot(8)
# 
# histogram_plot
# 
# ggsave("figures/Figure1_BiodivDistributions.pdf", histogram_plot, width = 6.5, height = 6, units = "in")
# ggsave("figures/Figure1_BiodivDistributions.png", histogram_plot, width = 6.5, height = 6, units = "in")
# 
# 
# # Figure 2: Numbers of hot spots identified --------------------------------
# hotspot_nums <- readRDS('results/final_comparisons.rds') %>% 
#   filter(type == "list_length") %>% 
#   filter(baseline == "NumUnique") %>% 
#   rowwise() %>%  
#   mutate(stdev = sqrt(var/100)) %>% 
#   select(-c(baseline, type)) %>%   
#   mutate(low_95 = mean - 1.96 * stdev,
#          high_95 = min(50, mean + 1.96 * stdev))
# 
# hotspot_nums$metric_label <-  case_match(hotspot_nums$comparison,
#                                          "NumUnique" ~ "Species richness (TD)",
#                                          "NumEndemic" ~ "Number of 'endemic' species",
#                                          "NumIndigName" ~ "Indigenous names (TEK)",
#                                          "NumUse" ~ "Recorded Indigenous uses (TEK)",
#                                          "FD" ~ "Rao's entropy (FD)",
#                                          "PD" ~ "Rao's entropy (PD)",
#                                          "FRic" ~ "Species richness (FD)",
#                                          # "FDiv" ~ "Functional sp. divergence",
#                                          "FDis" ~ "Species dispersion (FD)",
#                                          # "FEve" ~ "Functional sp. evenness",
#                                          # "Q" ~ "Rao's entropy (Q)",
#                                          "richness" ~ "Faith's (PD)",
#                                          # "GiniSimpson" ~ "Gini and Simpson",
#                                          # "Simpson" ~ "Simpson",
#                                          # "Shannon" ~ "Shannon",
#                                          # "Margalef" ~ "Margalef",
#                                          # "Menhinick" ~ "Menhinick",
#                                          # "McIntosh" ~ "McIntosh",
#                                          "PSVs" ~ "Species variability (PD)",
#                                          "PSR" ~ "Species richness (PD)",
#                                          "dummy" ~ ""
# )
# 
# hotspot_nums$type <- case_match(hotspot_nums$comparison,
#                                 c("NumUnique", "NumEndemic") ~ "Taxonomic",
#                                 c("NumIndigName", "NumUse") ~ "TEK",
#                                 c("FD", "FRic", "FDis") ~ "Functional",
#                                 c("PD","richness", "PSVs", "PSR") ~ "Phylogenetic"
#                                 # c("richness", "GiniSimpson", "Simpson",
#                                 # "Shannon", "Margalef", "Menhinick",
#                                 # "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
#                                 # c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"         
# )
# 
# hotspot_nums$color <- x_label_color_function(hotspot_nums$type)
# 
# hotspot_nums$type <- factor(hotspot_nums$type, levels = c("Taxonomic", "TEK", "Phylogenetic", "Functional"))
# 
# # Plot bar graph showing mean and standard error of the numbers of hotspots
# hotspot_bars = hotspot_nums %>% 
#   arrange(metric_label) %>% 
#   ggplot(aes(x = metric_label, color = type, fill = type)) +
#   geom_col(aes(y = mean)) +
#   geom_errorbar(aes(ymin = high_95, ymax = low_95), color = "purple", lwd = 0.5) +
#   scale_color_manual(breaks = c("Taxonomic", "TEK", "Phylogenetic", "Functional"),
#                      values = c("black", c4a("brewer.set1", 3))) +
#   scale_fill_manual(breaks = c("Taxonomic", "TEK", "Phylogenetic", "Functional"),
#                     values = c("black", c4a("brewer.set1", 3))) +
#   scale_y_continuous("Mean number of identified hotspots", breaks = 40:50) +
#   facet_grid(. ~ type, scales = "free_x", space = "free_x") +
#   coord_cartesian(ylim = c(44.9, 50.1)) +
#   theme_minimal(8) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust=1),
#     legend.position = "none"
#   )
# 
# hotspot_bars
# 
# ggsave("figures/Figure2_HotspotCounts.pdf", hotspot_bars, width = 6.5, height = 3, units = "in")
# ggsave("figures/Figure2_HotspotCounts.png", hotspot_bars, width = 6.5, height = 3, units = "in")
# 
# 
# # Figure S1: Distribution of traits (across all species) -------------------
# # Plot histograms in a single column
# trait_no_impute <- read_csv("data/clean/dataset_no_imputation.csv") %>%
#   relocate("Synonym", "Family", "Species", "Ssp_var",
#            "N_Langs", "N_Names", "N_Uses") %>% 
#   mutate(Woodiness = Woodiness == "woody") %>% 
#   mutate(Leaf_area_log = log(`Leaf area (mm2)`), .keep = "unused") %>% 
#   mutate(Plant_height_log = log(`Plant height (m)`), .keep = "unused") %>% 
#   mutate(LDMC_log = log(`LDMC (g/g)`), .keep = "unused") %>% 
#   pivot_longer(cols = c("N_Names":"LDMC_log")) %>% 
#   select(name, value)
# 
# trait_no_impute$label <-  case_match(
#   trait_no_impute$name,
#   "N_Names" ~ "Number of Indigenous names (TEK)",
#   "N_Uses" ~ "Number of Recorded Indigenous uses (TEK)",
#   "Nmass (mg/g)" ~ "Leaf nitrogen content per dry mass (mg/g)",
#   "Woodiness" ~ "Woodiness",
#   "Leaf_area_log" ~ "log(Leaf area (mm2))",
#   "Plant_height_log" ~ "log(Plant height (m))",
#   "LDMC_log" ~ "log(Leaf dry matter content(g/g))"
# )
# 
# logged_plots <- trait_no_impute %>% 
#   filter(label %in% c("log(Leaf area (mm2))",
#                       "log(Plant height (m))",
#                       "log(Leaf dry matter content(g/g))")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value") +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12)
# 
# LNPDM_plot <- trait_no_impute %>% 
#   filter(label %in% c("Leaf nitrogen content per dry mass (mg/g)")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value") +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12)
# 
# Woodiness_plot <- trait_no_impute %>% 
#   filter(label %in% c("Woodiness")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value", limits = c(-1,2)) +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12)
# 
# remaining_plots <- grid.arrange(LNPDM_plot, Woodiness_plot, nrow = 1)
# 
# Indig_plots <- trait_no_impute %>% 
#   filter(label %in% c(
#     "Number of Indigenous names (TEK)",
#     "Number of Recorded Indigenous uses (TEK)")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value") +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12)
# 
# plot_trait_no_impute <- grid.arrange(logged_plots, remaining_plots, Indig_plots, nrow = 3)
# 
# # Plot histograms in a single column
# trait_impute <- read_csv("data/clean/final_dataset.csv") %>%
#   relocate("Synonym", "Family", "Species", "Ssp_var",
#            "N_Langs", "N_Names", "N_Uses") %>% 
#   mutate(Woodiness = Woodiness == "woody") %>% 
#   mutate(Leaf_area_log = log(`Leaf area (mm2)`), .keep = "unused") %>% 
#   mutate(Plant_height_log = log(`Plant height (m)`), .keep = "unused") %>% 
#   mutate(LDMC_log = log(`LDMC (g/g)`), .keep = "unused") %>% 
#   mutate(DiasporeMass_log = log(`Diaspore mass (mg)`), .keep = "unused") %>% 
#   pivot_longer(cols = c("N_Names":"DiasporeMass_log")) %>% 
#   select(name, value)
# 
# trait_impute$label <-  case_match(
#   trait_impute$name,
#   "N_Names" ~ "Number of Indigenous names (TEK)",
#   "N_Uses" ~ "Number of Recorded Indigenous uses (TEK)",
#   "Nmass (mg/g)" ~ "Leaf nitrogen content per dry mass (mg/g)",
#   "Woodiness" ~ "Woodiness",
#   "Leaf_area_log" ~ "log(Leaf area (mm2))",
#   "Plant_height_log" ~ "log(Plant height (m))",
#   "LDMC_log" ~ "log(Leaf dry matter content(g/g))",
#   "DiasporeMass_log" ~ "log(Diaspore mass (mg))"
# )
# 
# logged_plots <- trait_impute %>% 
#   filter(label %in% c("log(Leaf area (mm2))",
#                       "log(Plant height (m))",
#                       "log(Leaf dry matter content(g/g))")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value") +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12) +
#   theme(axis.title.x = element_blank())
# 
# LNPDM_plot <- trait_impute %>% 
#   filter(label %in% c("Leaf nitrogen content per dry mass (mg/g)")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value") +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12) +
#   theme(axis.title.y = element_blank()) +
#   theme(axis.title.x = element_blank())
# 
# Woodiness_plot <- trait_impute %>% 
#   filter(label %in% c("Woodiness")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value", limits = c(-1,2)) +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12) +
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_blank())
# 
# DiasporeMass_plot <- trait_impute %>% 
#   filter(label %in% c("log(Diaspore mass (mg))")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value", limits = c(-1,2)) +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12) +
#   theme(axis.title.x = element_blank())
# 
# remaining_plots <- grid.arrange(DiasporeMass_plot, LNPDM_plot, Woodiness_plot, nrow = 1)
# 
# Indig_plots <- trait_impute %>% 
#   filter(label %in% c(
#     "Number of Indigenous names (TEK)",
#     "Number of Recorded Indigenous uses (TEK)")) %>% 
#   ggplot(aes(x = value)) +
#   # plot histogram using density instead of count
#   geom_histogram(aes(x = value, y = ..density..),
#                  colour = "black",
#                  fill = "white"
#   ) +
#   # add a spline approximating the probability density function
#   geom_density(aes(x = value),
#                alpha = .2,
#                fill = "#FF6666"
#   ) +
#   scale_x_continuous(name = "Value") +
#   scale_y_continuous(name = "Density") +
#   # one histogram for each biodiversity metric
#   facet_wrap(~label, scales = "free") +
#   theme_cowplot(12) +
#   theme(axis.title.x = element_blank())
# 
# plot_trait_impute <- grid.arrange(logged_plots, remaining_plots, Indig_plots, nrow = 3)
# 
# ggsave("figures/FigS1_TraitDistributions.pdf", plot_trait_impute, width = 10.5, height = 6, units = "in")
# ggsave("figures/FigS1_TraitDistributions.png", plot_trait_impute, width = 10.5, height = 6, units = "in")


# Figure 3: Biodiversity metric dendrogram clustered by precision ---------

# Make more descriptive labels for biodiversity metrics
full_comp_df <- readRDS('results/final_comparisons.rds') #%>% 
  # filter(baseline != "FRic", comparison != "FRic")

full_comp_df$baseline_label <-  case_match(full_comp_df$baseline,
                                           "NumUnique" ~ "Species richness (TD)",
                                           "NumEndemic" ~ "Number of 'endemic' species",
                                           "NumIndigName" ~ "Indigenous names (TEK)",
                                           "NumUse" ~ "Recorded Indigenous uses (TEK)",
                                           "FD" ~ "Rao's entropy (FD)",
                                           "PD" ~ "Rao's entropy (PD)",
                                           "FRic" ~ "Species richness (FD)",
                                           # "FDiv" ~ "Functional sp. divergence",
                                           "FDis" ~ "Species dispersion (FD)",
                                           # "FEve" ~ "Functional sp. evenness",
                                           # "Q" ~ "Rao's entropy (Q)",
                                           "richness" ~ "Faith's (PD)",
                                           # "GiniSimpson" ~ "Gini and Simpson",
                                           # "Simpson" ~ "Simpson",
                                           # "Shannon" ~ "Shannon",
                                           # "Margalef" ~ "Margalef",
                                           # "Menhinick" ~ "Menhinick",
                                           # "McIntosh" ~ "McIntosh",
                                           "PSVs" ~ "Species variability (PD)",
                                           "PSR" ~ "Species richness (PD)",
                                           "dummy" ~ ""
)

full_comp_df$comparison_label <-  case_match(full_comp_df$comparison,
                                             "NumUnique" ~ "Species richness (TD)",
                                             "NumEndemic" ~ "Number of 'endemic' species",
                                             "NumIndigName" ~ "Indigenous names (TEK)",
                                             "NumUse" ~ "Recorded Indigenous uses (TEK)",
                                             "FD" ~ "Rao's entropy (FD)",
                                             "PD" ~ "Rao's entropy (PD)",
                                             "FRic" ~ "Species richness (FD)",
                                             # "FDiv" ~ "Functional sp. divergence",
                                             "FDis" ~ "Species dispersion (FD)",
                                             # "FEve" ~ "Functional sp. evenness",
                                             # "Q" ~ "Rao's entropy (Q)",
                                             "richness" ~ "Faith's (PD)",
                                             # "GiniSimpson" ~ "Gini and Simpson",
                                             # "Simpson" ~ "Simpson",
                                             # "Shannon" ~ "Shannon",
                                             # "Margalef" ~ "Margalef",
                                             # "Menhinick" ~ "Menhinick",
                                             # "McIntosh" ~ "McIntosh",
                                             "PSVs" ~ "Species variability (PD)",
                                             "PSR" ~ "Species richness (PD)",
                                             "dummy" ~ ""
)


full_comp_df$baseline_type <- case_match(full_comp_df$baseline,
                                         c("NumUnique", "NumEndemic") ~ "Taxonomic",
                                         c("NumIndigName", "NumUse") ~ "TEK",
                                         c("FD", "FRic", "FDis") ~ "Functional",
                                         # c("FD", "FDis") ~ "Functional",
                                         c("PD","richness", "PSVs", "PSR") ~ "Phylogenetic"
                                         # c("richness", "GiniSimpson", "Simpson",
                                         # "Shannon", "Margalef", "Menhinick",
                                         # "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
                                         # c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"            
)

full_comp_df$comparison_type <- case_match(full_comp_df$comparison,
                                           c("NumUnique", "NumEndemic") ~ "Taxonomic",
                                           c("NumIndigName", "NumUse") ~ "TEK",
                                           # c("FD", "FDis") ~ "Functional",
                                           c("FD", "FRic", "FDis") ~ "Functional",
                                           c("PD","PSVs", "PSR") ~ "Phylogenetic",
                                           # c("PD","richness", "PSVs", "PSR") ~ "Phylogenetic"
                                           # c("richness", "GiniSimpson", "Simpson",
                                           # "Shannon", "Margalef", "Menhinick",
                                           # "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
                                           # c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"           
)

full_comp_df$baseline_color <- x_label_color_function(full_comp_df$baseline_type)
full_comp_df$comparison_color <- x_label_color_function(full_comp_df$comparison_type)

# Cluster metrics based on the Euclidean distance of comparison values

# Jaccard
pair_jaccard_mat <- full_comp_df %>%
  filter(type == "jaccard") %>%
  select(-c(var, type, baseline, comparison, baseline_type, comparison_type,
            baseline_color, comparison_color)) %>%
  pivot_wider(values_from = c("mean"),
              names_sort = F,
              names_from = "comparison_label")

m_jaccard <- as.matrix((pair_jaccard_mat[, -1]), ncol = 10)
pair_jaccard_cluster <- hclust(dist(t(m_jaccard)), method = "ward.D2")

# Precision
pair_prec_mat <- full_comp_df %>%
  filter(type == "precision") %>% 
  select(-c(var, type, baseline, comparison, baseline_type, comparison_type,
            baseline_color, comparison_color)) %>% 
  pivot_wider(values_from = c("mean"),
              names_sort = F,
              names_from = "comparison_label")

m_prec <- as.matrix((pair_prec_mat[, -1]), ncol = 18)
pair_prec_cluster <- hclust(dist(t(m_prec)), method = "ward.D2")

# Recall
pair_recall_mat <- full_comp_df %>%
  filter(type == "recall") %>% 
  select(-c(var, type, baseline, comparison, baseline_type, comparison_type,
            baseline_color, comparison_color)) %>% 
  pivot_wider(values_from = c("mean"),
              names_sort = F,
              names_from = "comparison_label")

m_recall <- as.matrix((pair_recall_mat[, -1]), ncol = 18)
pair_recall_cluster <- hclust(dist(t(m_recall)), method = "ward.D2")

ddata <- dendro_data(pair_jaccard_cluster, type="rectangle")

ddata_labels = label(ddata)
# ddata_labels$color <- rev(c("#E41A1C", "#E41A1C", "#377EB8", "black", "#377EB8", "#377EB8", 
#                             "#4DAF4A", "#377EB8", "#377EB8", "#377EB8", "#4DAF4A", "#4DAF4A", 
#                             "#4DAF4A", "#4DAF4A", "#377EB8", "#377EB8", "#377EB8"))
ddata_labels$color <- rev(c("#E41A1C", "#E41A1C", "black", rep("#4DAF4A", 3), rep("#377EB8", 4)))
# ddata_labels$color <- rev(c("#E41A1C", "#E41A1C", "black", rep("#4DAF4A", 2), rep("#377EB8", 4)))


ddata_labels$group <- case_match(
  ddata_labels$label,
  "Species richness (TD)" ~ "Taxonomic",
  c("Indigenous names (TEK)", "Recorded Indigenous uses (TEK)") ~ "TEK",
  # c("Rao's entropy (FD)", "Species dispersion (FD)") ~ "Functional",
  c("Rao's entropy (FD)", "Species richness (FD)", "Species dispersion (FD)") ~ "Functional",
  c("Rao's entropy (PD)", "Faith's (PD)", "Species variability (PD)", "Species richness (PD)") ~ "Phylogenetic",
  # c("Faith's (PD)", "Gini and Simpson", "Simpson",
  #   "Shannon", "Margalef", "Menhinick",
  #   "McIntosh", "Species variability (PD)", "Species richness (PD)") ~ "Phylogenetic",
  # c("Species richness (FD)", "Functional sp. divergence", "Species dispersion (FD)", "Functional sp. evenness", "Rao's entropy (Q)") ~ "Functional"         
)

ddata_labels$color <- x_label_color_function(ddata_labels$group)

ddata_labels <- arrange(ddata_labels, label)

ggplot() +
  geom_segment(data = segment(ddata),
               aes_string(x = "x", y = "y", xend = "xend", yend = "yend"),
               linewidth = 0.1) +
  geom_text(data = ddata_labels, 
            aes(x = x, y = y, label = label, color = label),
            fontface = "bold", hjust = 0, angle = 0,
            size = .75) + 
  scale_color_manual(values = as.character(ddata_labels$color)) +
  scale_y_continuous(expand=c(0.5, 0),
                     trans = "reverse") + 
  coord_flip() +
  guides(color = "none") +
  theme_cowplot(4) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
  )

ggsave("figures/Figure3_Dendrogram_Jaccard.pdf", width = 3.5625, height = 1, units = "in",
       dpi = 600)
ggsave("figures/Figure3_Dendrogram_Jaccard.png", width = 3.5625, height = 1, units = "in",
       dpi = 600)

# Figure S2: Pairwise comparison heatmaps ---------------------------------

# Save and load the name order (rev(colnames(m_prec)[pair_prec_cluster$order])) from the "actual" analysis to use for the others to make comparisons 
# standard_name_order = (rev(colnames(m_prec)[pair_prec_cluster$order]))
# write.csv(standard_name_order, "results/standard_name_order.csv")

standard_name_order = (read_csv("results/standard_name_order.csv"))$x

# Get colors assigned to biodiversity metrics to use for text labels
plot_colours <- full_comp_df %>%
  select(baseline_color, baseline_label) %>%
  unique() %>%
  slice(match(colnames(m_jaccard)[pair_jaccard_cluster$order], baseline_label)) %>%
  select(baseline_color)

# Jaccard clustered heatmap
jaccard_heatmap <- full_comp_df %>%
  arrange(desc(baseline), desc(comparison)) %>%
  filter(type == "jaccard") %>%
  mutate(mean = ifelse(baseline_label == comparison_label, NA, mean)) %>%
  ggplot(aes(baseline_label, comparison_label, fill = mean, colour = baseline_type)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(mean, 2)), color = "black", size = 2.5) +
  scale_fill_gradient("Mean value", low = "white", high = "blue") +
  scale_y_discrete("Comparison index", limits = standard_name_order) +
  scale_x_discrete("Baseline index", limits = standard_name_order) +
  theme_cowplot(10) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1,
                               color = rev(plot_colours$baseline_color)),
    axis.text.y = element_text(color = rev(plot_colours$baseline_color))
  ) +
  guides(fill = guide_colorbar(barheight = 25,
                               barwidth = 0.5,
                               title.hjust = 0.5))

# Get colors assigned to biodiversity metrics to use for text labels
plot_colours <- full_comp_df %>% 
  select(baseline_color, baseline_label) %>% 
  unique() %>% 
  slice(match(colnames(m_prec)[pair_prec_cluster$order], baseline_label)) %>% 
  select(baseline_color)

# Precision clustered heatmap
precision_heatmap <- full_comp_df %>% 
  arrange(desc(baseline), desc(comparison)) %>% 
  filter(type == "precision") %>% 
  mutate(mean = ifelse(baseline_label == comparison_label, NA, mean)) %>% 
  ggplot(aes(baseline_label, comparison_label, fill = mean, colour = baseline_type)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(mean, 2)), color = "black", size = 2.5) +
  scale_fill_gradient("Mean value", low = "white", high = "blue") +
  scale_y_discrete("Comparison index", limits = standard_name_order) +
  scale_x_discrete("Baseline index", limits = standard_name_order) +
  theme_cowplot(10) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, 
                               color = rev(plot_colours$baseline_color)),
    axis.text.y = element_text(color = rev(plot_colours$baseline_color))
    ) +
  guides(fill = guide_colorbar(barheight = 25,
                               barwidth = 0.5,
                               title.hjust = 0.5))

# Recall clustered heatmap
recall_heatmap <- full_comp_df %>% 
  arrange(desc(baseline), desc(comparison)) %>% 
  filter(type == "recall") %>% 
  mutate(mean = ifelse(baseline_label == comparison_label, NA, mean)) %>% 
  ggplot(aes(baseline_label, comparison_label, fill = mean, colour = baseline_type)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(mean, 2)), color = "black", size = 2.5) +
  scale_fill_gradient("Mean value", low = "white", high = "blue") +
  scale_y_discrete("Comparison index", limits = standard_name_order) +
  scale_x_discrete("Baseline index", limits = standard_name_order) +
  theme_cowplot(10) +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, 
                               color = rev(plot_colours$baseline_color)),
    axis.text.y = element_text(color = rev(plot_colours$baseline_color))
    
  ) +
  guides(fill = guide_colorbar(barheight = 25,
                               barwidth = 0.5,
                               title.hjust = 0.5))

# Combine heatmaps into a single plot
heatmaps <- ggarrange(NULL, precision_heatmap, NULL, recall_heatmap, nrow = 4, 
                      heights = c(0.1, 1, 0.05, 1),
                      labels = c(NA, "A) Precision", NA, "B) Sensitivity"), 
                      label.y = c(1.1, 1.1),
                      legend = "right", common.legend = TRUE)

ggsave("figures/FigS2_PairwiseHeatmaps_Jaccard.png", jaccard_heatmap, width = 6.5, height = 3.5, units = "in")


ggsave("figures/FigS2_PairwiseHeatmaps.pdf", heatmaps, width = 6.5, height = 7, units = "in")
ggsave("figures/FigS2_PairwiseHeatmaps.png", heatmaps, width = 6.5, height = 7, units = "in")
