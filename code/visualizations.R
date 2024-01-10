#### Create visualizations from simulations ####################################

## Load libraries --------------------------------------------------------------
require(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(GGally) # used to make paired scatter plots
library(retry)
library(cols4all)
library(gridExtra)
library(ggh4x)

# Load in simulations
source("code/functions.R")

# Figures -----------------------------------------------------------------

# Figure 0: Distribution of TEK traits (across all species) ----

# Plot histograms in a single column
trait_no_impute <- read_csv("data/clean/dataset_no_imputation.csv") %>%
  relocate("Synonym", "Family", "Species", "Ssp_var",
           "N_Langs", "N_Names", "N_Uses") %>% 
  mutate(Woodiness = Woodiness == "woody") %>% 
  mutate(Leaf_area_log = log(`Leaf area (mm2)`), .keep = "unused") %>% 
  mutate(Plant_height_log = log(`Plant height (m)`), .keep = "unused") %>% 
  mutate(LDMC_log = log(`LDMC (g/g)`), .keep = "unused") %>% 
  pivot_longer(cols = c("N_Names":"LDMC_log")) %>% 
  select(name, value)

trait_no_impute$label <-  case_match(
  trait_no_impute$name,
  "N_Names" ~ "Number of Indigenous names",
  "N_Uses" ~ "Number of recorded Indigenous uses",
  "Nmass (mg/g)" ~ "Leaf nitrogen content per dry mass (mg/g)",
  "Woodiness" ~ "Woodiness",
  "Leaf_area_log" ~ "log(Leaf area (mm2))",
  "Plant_height_log" ~ "log(Plant height (m))",
  "LDMC_log" ~ "log(Leaf dry matter content(g/g))"
)

logged_plots <- trait_no_impute %>% 
  filter(label %in% c("log(Leaf area (mm2))",
                      "log(Plant height (m))",
                      "log(Leaf dry matter content(g/g))")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value") +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

LNPDM_plot <- trait_no_impute %>% 
  filter(label %in% c("Leaf nitrogen content per dry mass (mg/g)")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value") +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

Woodiness_plot <- trait_no_impute %>% 
  filter(label %in% c("Woodiness")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value", limits = c(-1,2)) +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

remaining_plots <- grid.arrange(LNPDM_plot, Woodiness_plot, nrow = 1)

Indig_plots <- trait_no_impute %>% 
  filter(label %in% c(
    "Number of Indigenous names",
    "Number of recorded Indigenous uses")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value") +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

plot_trait_no_impute <- grid.arrange(logged_plots, remaining_plots, Indig_plots, nrow = 3)


# Plot histograms in a single column
trait_impute <- read_csv("data/clean/final_dataset.csv") %>%
  relocate("Synonym", "Family", "Species", "Ssp_var",
           "N_Langs", "N_Names", "N_Uses") %>% 
  mutate(Woodiness = Woodiness == "woody") %>% 
  mutate(Leaf_area_log = log(`Leaf area (mm2)`), .keep = "unused") %>% 
  mutate(Plant_height_log = log(`Plant height (m)`), .keep = "unused") %>% 
  mutate(LDMC_log = log(`LDMC (g/g)`), .keep = "unused") %>% 
  pivot_longer(cols = c("N_Names":"LDMC_log")) %>% 
  select(name, value)

trait_impute$label <-  case_match(
  trait_impute$name,
  "N_Names" ~ "Number of Indigenous names",
  "N_Uses" ~ "Number of recorded Indigenous uses",
  "Nmass (mg/g)" ~ "Leaf nitrogen content per dry mass (mg/g)",
  "Woodiness" ~ "Woodiness",
  "Leaf_area_log" ~ "log(Leaf area (mm2))",
  "Plant_height_log" ~ "log(Plant height (m))",
  "LDMC_log" ~ "log(Leaf dry matter content(g/g))"
)

logged_plots <- trait_impute %>% 
  filter(label %in% c("log(Leaf area (mm2))",
                      "log(Plant height (m))",
                      "log(Leaf dry matter content(g/g))")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value") +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

LNPDM_plot <- trait_impute %>% 
  filter(label %in% c("Leaf nitrogen content per dry mass (mg/g)")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value") +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

Woodiness_plot <- trait_impute %>% 
  filter(label %in% c("Woodiness")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value", limits = c(-1,2)) +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

remaining_plots <- grid.arrange(LNPDM_plot, Woodiness_plot, nrow = 1)

Indig_plots <- trait_impute %>% 
  filter(label %in% c(
    "Number of Indigenous names",
    "Number of recorded Indigenous uses")) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, y = ..density..),
                 colour = "black",
                 fill = "white"
  ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  scale_x_continuous(name = "Value") +
  scale_y_continuous(name = "Density") +
  # one histogram for each biodiversity metric
  facet_wrap(~label, scales = "free") +
  theme_cowplot(12)

plot_trait_impute <- grid.arrange(logged_plots, remaining_plots, Indig_plots, nrow = 3)

# Compare distributions of traits before and after imputation (only functional traits)

compare_dists <- read_csv("data/clean/final_dataset.csv") %>%
  mutate(origin = "imputed") %>% 
  rbind(read_csv("data/clean/dataset_no_imputation.csv") %>% 
          mutate(origin = "original")) %>% 
  select(-c(Synonym:Ssp_var)) %>% 
  relocate(origin) %>% 
  mutate(Woodiness = Woodiness == "woody") %>% 
  pivot_longer(cols = c("Nmass (mg/g)":"LDMC (g/g)")) 

plot_compare_dists <- compare_dists %>% 
  ggplot(aes(x = value, fill = origin)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, after_stat(density))) +
  # # add a spline approximating the probability density function
  # geom_density(alpha = .2) +
  # one histogram for each biodiversity metric
  facet_wrap( ~ name, scales = "free") +
  theme_cowplot(12)

plot_compare_dists


# Figure 1: Compare the distribution of biodiversity metric values --------
# Define the relevant trait names
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", "Leaf area (mm2)")

# Initialize the full phylogenetic tree
tree <- get.phylo_tree(read_csv("data/clean/final_dataset.csv"))

# Get a single simulation
numPatches <- 1000
full_df_sample <- get.full_df(numPatches)

# Calculate biodiversity metrics
biodiv_df <- get.biodiv_df(full_df_sample, trait_names, tree)

# Set up label colors to indicate what type of index it is
biodiv_plot_df <- biodiv_df %>% 
  pivot_longer(cols = -Patch)

biodiv_plot_df$group <- case_match(biodiv_plot_df$name,
                                   c("NumUnique", "NumEndemic") ~ "Basic",
                                   c("NumIndigName", "NumUse") ~ "TEK",
                                   c("richness", "GiniSimpson", "Simpson",
                                     "Shannon", "Margalef", "Menhinick",
                                     "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
                                   c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"         
)

# Remove number of "endemic" from analysis
x_label_color_function <- function(string) {
  out <- #paste("<span style = 'color: ",
    case_match(string,
               "Basic" ~ "black",
               "TEK" ~ c4a("brewer.set1", 3)[1],
               "Phylogenetic" ~ c4a("brewer.set1", 3)[2],
               "Functional" ~ c4a("brewer.set1", 3)[3],)#,
  # ";'>",
  # string,
  # "</span>", sep = "")
  out <- as.character(out)
}

biodiv_plot_df$color <- x_label_color_function(biodiv_plot_df$group)

# Make labels more descriptive
biodiv_plot_df$label <-  case_match(biodiv_plot_df$name,
                                    "NumUnique" ~ "Species richness",
                                    "NumEndemic" ~ "Number of 'endemic' species",
                                    "NumIndigName" ~ "Number of Indigenous names",
                                    "NumUse" ~ "Number of recorded Indigenous uses",
                                    "FRic" ~ "Functional species richness",
                                    "FDiv" ~ "Functional species divergence",
                                    "FDis" ~ "Functional species dispersion",
                                    "FEve" ~ "Functional species evenness",
                                    "Q" ~ "Rao's entropy (Q)",
                                    "richness" ~ "Faith",
                                    "GiniSimpson" ~ "Gini and Simpson",
                                    "Simpson" ~ "Simpson",
                                    "Shannon" ~ "Shannon",
                                    "Margalef" ~ "Margalef",
                                    "Menhinick" ~ "Menhinick",
                                    "McIntosh" ~ "McIntosh",
                                    "PSVs" ~ "Phylogenetic species variability",
                                    "PSR" ~ "Phylogenetic species richness"
)

biodiv_plot_df$label <-  factor(biodiv_plot_df$label, levels = c(
  "Species richness","Number of Indigenous names",   "Number of recorded Indigenous uses", 
  "Number of 'endemic' species",
  "Functional species richness", "Functional species divergence",  "Functional species dispersion",  "Functional species evenness", "Rao's entropy (Q)",
  "Faith", "Gini and Simpson", "Margalef", "McIntosh", "Menhinick",   "Phylogenetic species variability",  "Phylogenetic species richness", "Shannon", "Simpson"
)
)

# biodiv_plot_df <- biodiv_plot_df %>%
#   select(Patch, value, label) %>%
#   pivot_wider(names_from = label, values_from = value)

# Get mean and standard deviation values for each index
meanSD_df <- biodiv_plot_df %>% 
  filter(label != "Number of 'endemic' species") %>%
  # select(-`Number of 'endemic' species`) %>% 
  # pivot_longer(cols = `Species richness`:`Phylogenetic species richness`,
  #              names_to = "variable") %>% 
  # filter(!is.na(value)) %>% 
  group_by(label) %>% 
  mutate(mean = mean(value, na.rm = TRUE),
         stdev = sd(value, na.rm = TRUE),
         upper = mean + 1.96 * stdev,
         hotspot_cutoff = quantile(value, c(0.95), na.rm = TRUE)
  )

hotspot_stats <- meanSD_df %>% 
  filter(value > hotspot_cutoff) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            stdev = sd(value, na.rm = TRUE),
            lower = mean - 1.96 * stdev
  )

nonhotspot_stats <- meanSD_df %>% 
  filter(value <= hotspot_cutoff) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            stdev = sd(value, na.rm = TRUE),
            upper = mean + 1.96 * stdev
  )

# Plot histograms in a single column
histogram_colors = c("black", rep("#E41A1C",2), rep("#4DAF4A", 5), rep("#377EB8", 9))

strip_theme = strip_themed(background_x = elem_list_rect(fill = histogram_colors),
                           text_x = element_text(color = "white",face = "bold", size = 14))

histogram_plot <- biodiv_plot_df %>%
  filter(label != "Number of 'endemic' species") %>%
  ggplot(aes(x = value)) +
  # # plot histogram using density instead of count
  # geom_histogram(aes(x = value, y = ..density..),
  #                colour = "black",
  #                fill = "white"
  # ) +
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  # add a vertical line showing where the top 5% are
  stat_summary(aes(xintercept = after_stat(x), y = 0),
               fun = quantile, fun.args = list(0.95),
               geom = "vline", orientation = "y",
               color = "blue", lwd = 1) + 
  geom_vline(data = meanSD_df, aes(xintercept = mean), linetype = 2) +
  # geom_vline(data = nonhotspot_stats, aes(xintercept = upper), 
  #            color = "green", linetype = 2) +
  # geom_vline(data = hotspot_stats, aes(xintercept = lower), 
  #            color = "blue", linetype = 2) +
  ylab('Density') +
  xlab('Value') +
  # one histogram for each biodiversity metric
  facet_wrap2( ~ label, ncol = 3, scales = "free",
               strip = strip_theme) +
  theme_cowplot(12)

histogram_plot



# Figure 2: Scatterplots of selected biodiversity metrics -----------------
biodiv_scatter <- biodiv_df %>% 
  select(-c(Patch, NumEndemic, 
            FDiv, FDis, FEve, Q,
            richness, GiniSimpson, Margalef, Menhinick, McIntosh)) %>% 
  ggpairs(aes()) +
  theme_cowplot(11)

biodiv_scatter



# Figure 3: Comparing biodiversity hotspots -------------------------------

# Number of patches to simulate
NumPatches <- 1000

# Number of sets of simulated patches to create
numIterations <- 100

# Names of traits to use for functional diversity
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", "Leaf area (mm2)")

# Tree to use for phylogenetic diversity
tree <- get.phylo_tree(final_data)

# Biodiversity metric to make comparisons with
baseline_metric <- "FDiv" # The full list is right below, the column names for compare_df

# Initialize comparison data frame
compare_df <- tibble(
  NumUnique = as.double(),
  NumEndemic = as.double(), NumIndigName = as.double(), NumUse = as.double(),
  # NumIndigLang = as.double()
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
      biodiv_comp_helper_func(NumPatches, trait_names, tree, baseline_metric),
      until = function(val, cnd) {
        !is.null(val)
      }
    ) %>% 
      mutate(iteration = j)
    
    precision_df <- biodiv.compare_df %>%
      select(-list_length) %>% 
      pivot_wider(names_from = variable) %>%
      unique() 
    
    list_length_df <- biodiv.compare_df %>%
      select(-"value") %>% 
      pivot_wider(names_from = variable, values_from = list_length) %>%
      unique() 
    
    gc()
    out_df <- rbind(
      mutate(precision_df, type = "precision"),
      mutate(list_length_df, type = "list_length")
    )
    
    # # Add to the list
    # compare_df <- add_row(
    #   compare_df,
    #   biodiv.compare_df
    # )
    
  } %>% unique()

stopCluster(my.cluster)

write_rds(compare_df, "results/compare_biodiv.rds")

compare_df <- read_rds("results/compare_biodiv.rds")

# Create figure comparing precision
compare_plot <- compare_df %>%
  filter(type == "precision") %>% 
  select(-iteration) %>%
  melt() %>%
  group_by(variable) %>% 
  mutate(mean = mean(value)) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, after_stat(density)),
                 bins = 30
  ) +
  # add a vertical line showing the mean
  geom_vline(aes(xintercept = mean, group = variable), 
             colour = "red") +  
  
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  
  facet_wrap( ~ variable, scales = "free") +
  # x axis
  scale_x_continuous(name = "Relative precision", 
                     breaks = seq(0, 1, by = 0.1),
                     limits = c(0,1)) +
  # title
  ggtitle("How similar is each biodiversity metric to functional divergence?") +
  theme_cowplot(font_size = 11)

compare_plot

# Create figure comparing list lengths
list_length_plot <- compare_df %>%
  filter(type == "list_length") %>% 
  select(-iteration) %>%
  melt() %>%
  group_by(variable) %>% 
  mutate(mean = mean(value)) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, after_stat(density)),
                 bins = 30
  ) +
  # add a vertical line showing the mean
  geom_vline(aes(xintercept = mean, group = variable), 
             colour = "red") +  
  facet_wrap( ~ variable, scales = "free") +
  # x axis
  scale_x_continuous(name = "Number of hotspots identified", 
                     limits = c(0,NA)) +
  # title
  ggtitle("How similar are the numbers of hotspots identified?") +
  theme_cowplot(font_size = 11)
list_length_plot

# Figure 4: Comparisons between functional/phylogenetic measures ----------

# Initialize the full phylogenetic tree
tree <- get.phylo_tree(read_csv("data/clean/final_dataset.csv"))
# Names of traits to use for functional diversity
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", "Leaf area (mm2)")

# Get a single simulation
numPatches <- 1000
full_df_sample <- get.full_df(numPatches)

# Calculate biodiversity metrics
biodiv_df <- get.biodiv_df(full_df_sample, trait_names, tree)

# Functional diversity measures
func_scatter <- biodiv_df %>% 
  select(FDiv, FRic, FDis, FEve, Q) %>% 
  ggpairs(aes()) +
  ggtitle('Comparison of functional trait biodiversity measures') +
  theme_cowplot(11)

# Phylogenetic diversity measures
phylo_scatter <- biodiv_df %>% 
  select(richness, GiniSimpson, Simpson, Shannon, Margalef, Menhinick, McIntosh, PSVs, PSR) %>% 
  ggpairs(aes()) +
  ggtitle('Comparison of phylogenetic trait biodiversity measures') +
  theme_cowplot(11)

# compare_df <- read_rds("results/compare_biodiv.rds")

# Functional diversity measures
func_compare_plot <- compare_df %>%
  filter(type == "precision") %>% 
  select(-iteration) %>%
  select(FDiv, FRic, FDis, FEve, Q) %>% 
  melt() %>%
  group_by(variable) %>% 
  mutate(mean = mean(value)) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, after_stat(density)),
                 bins = 30
  ) +
  # add a vertical line showing the mean
  geom_vline(aes(xintercept = mean, group = variable), 
             colour = "red") +  
  
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  
  facet_wrap( ~ variable, scales = "free") +
  # x axis
  scale_x_continuous(name = "Relative precision", 
                     breaks = seq(0, 1, by = 0.1),
                     limits = c(0,1)) +
  # title
  ggtitle("How similar are the functional diversity measures?") +
  theme_cowplot(font_size = 11)

func_compare_plot

# Phylogenetic diversity measures
phylo_compare_plot <- compare_df %>%
  filter(type == "precision") %>% 
  select(-iteration) %>%
  select(richness, GiniSimpson, Simpson, Shannon, Margalef, Menhinick, McIntosh, PSVs, PSR) %>% 
  melt() %>%
  group_by(variable) %>% 
  mutate(mean = mean(value)) %>% 
  ggplot(aes(x = value)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, after_stat(density)),
                 bins = 30
  ) +
  # add a vertical line showing the mean
  geom_vline(aes(xintercept = mean, group = variable), 
             colour = "red") +  
  
  # add a spline approximating the probability density function
  geom_density(aes(x = value),
               alpha = .2,
               fill = "#FF6666"
  ) +
  
  facet_wrap( ~ variable, scales = "free") +
  # x axis
  scale_x_continuous(name = "Relative precision", 
                     breaks = seq(0, 1, by = 0.1),
                     limits = c(0,1)) +
  # title
  ggtitle("How similar are the phylogenetic diversity measures?") +
  theme_cowplot(font_size = 11)

phylo_compare_plot

# Figure 5: Pairwise precision heatmap ------------------------------------

pair_prec_mean <- readRDS('full_comparisons.rds') %>% 
  filter(type == "precision") %>% 
  select(-c(var, type)) %>% 
  pivot_wider(values_from = c("mean"),
              names_sort = F,
              names_from = "comparison")

pair_prec_mean_mat <- matrix(unlist(pair_prec_mean[1:18, 2:19]), ncol = 18)

# Check how far off this matrix is from being symmetric
isSymmetric.matrix(pair_prec_mean_mat) # It's not symmetric
# How far off is it? This value would be near zero if matrix is near symmetriic
max(abs(pair_prec_mean_mat-t(pair_prec_mean_mat))) 

pair_prec_plot <- readRDS('full_comparisons.rds') %>% 
  filter(type == "precision") %>% 
  ggplot(mapping = aes(x = baseline, y = comparison, fill = mean
                       #, alpha = rev(log(1+prec_var)))
  )) +
  geom_tile() +
  geom_text(aes(label = round(mean, 2)))  +
  scale_fill_gradient(low = "white", high = "blue")

# Make more descriptive labels for biodiversity metrics
full_comp_df <- readRDS('full_comparisons.rds') 

full_comp_df$baseline_label <-  case_match(full_comp_df$baseline,
                                           "NumUnique" ~ "Species richness",
                                           "NumEndemic" ~ "Number of 'endemic' species",
                                           "NumIndigName" ~ "Number of Indigenous names",
                                           "NumUse" ~ "Number of recorded Indigenous uses",
                                           "FRic" ~ "Functional species richness",
                                           "FDiv" ~ "Functional species divergence",
                                           "FDis" ~ "Functional species dispersion",
                                           "FEve" ~ "Functional species evenness",
                                           "Q" ~ "Rao's entropy (Q)",
                                           "richness" ~ "Faith",
                                           "GiniSimpson" ~ "Gini and Simpson",
                                           "Simpson" ~ "Simpson",
                                           "Shannon" ~ "Shannon",
                                           "Margalef" ~ "Margalef",
                                           "Menhinick" ~ "Menhinick",
                                           "McIntosh" ~ "McIntosh",
                                           "PSVs" ~ "Phylogenetic species variability",
                                           "PSR" ~ "Phylogenetic species richness"
)

full_comp_df$comparison_label <-  case_match(full_comp_df$comparison,
                                             "NumUnique" ~ "Species richness",
                                             "NumEndemic" ~ "Number of 'endemic' species",
                                             "NumIndigName" ~ "Number of Indigenous names",
                                             "NumUse" ~ "Number of recorded Indigenous uses",
                                             "FRic" ~ "Functional species richness",
                                             "FDiv" ~ "Functional species divergence",
                                             "FDis" ~ "Functional species dispersion",
                                             "FEve" ~ "Functional species evenness",
                                             "Q" ~ "Rao's entropy (Q)",
                                             "richness" ~ "Faith",
                                             "GiniSimpson" ~ "Gini and Simpson",
                                             "Simpson" ~ "Simpson",
                                             "Shannon" ~ "Shannon",
                                             "Margalef" ~ "Margalef",
                                             "Menhinick" ~ "Menhinick",
                                             "McIntosh" ~ "McIntosh",
                                             "PSVs" ~ "Phylogenetic species variability",
                                             "PSR" ~ "Phylogenetic species richness"
)


full_comp_df$baseline_type <- case_match(full_comp_df$baseline,
                                         c("NumUnique", "NumEndemic") ~ "Basic",
                                         c("NumIndigName", "NumUse") ~ "TEK",
                                         c("richness", "GiniSimpson", "Simpson",
                                           "Shannon", "Margalef", "Menhinick",
                                           "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
                                         c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"         
)

full_comp_df$comparison_type <- case_match(full_comp_df$comparison,
                                           c("NumUnique", "NumEndemic") ~ "Basic",
                                           c("NumIndigName", "NumUse") ~ "TEK",
                                           c("richness", "GiniSimpson", "Simpson",
                                             "Shannon", "Margalef", "Menhinick",
                                             "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
                                           c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"         
)

# Remove number of "endemic" from analysis
full_comp_df <- filter(full_comp_df, baseline != "NumEndemic", comparison != "NumEndemic")

x_label_color_function <- function(string) {
  out <- #paste("<span style = 'color: ",
    case_match(string,
               "Basic" ~ "black",
               "TEK" ~ c4a("brewer.set1", 3)[1],
               "Phylogenetic" ~ c4a("brewer.set1", 3)[2],
               "Functional" ~ c4a("brewer.set1", 3)[3],)#,
  # ";'>",
  # string,
  # "</span>", sep = "")
  out <- as.character(out)
}

full_comp_df$baseline_color <- x_label_color_function(full_comp_df$baseline_type)
full_comp_df$comparison_color <- x_label_color_function(full_comp_df$comparison_type)

# Cluster based on Euclidean distance

pair_prec_mat <- full_comp_df %>%
  filter(type == "precision") %>% 
  select(-c(var, type, baseline, comparison, baseline_type, comparison_type,
            baseline_color, comparison_color)) %>% 
  pivot_wider(values_from = c("mean"),
              names_sort = F,
              names_from = "comparison_label")

m <- as.matrix((pair_prec_mat[, -1]), ncol = 18)

pair_prec_cluster <- hclust(dist(t(m)), method = "ward.D2")

# Cluster dendrogram
require(ggdendro)

dhc <- as.dendrogram(pair_prec_cluster)

ddata <- dendro_data(pair_prec_cluster, type="rectangle")

ddata_labels = label(ddata)
ddata_labels$color <- rev(c("#E41A1C", "#E41A1C", "#377EB8", "black", "#377EB8", "#377EB8", 
                            "#4DAF4A", "#377EB8", "#377EB8", "#377EB8", "#4DAF4A", "#4DAF4A", 
                            "#4DAF4A", "#4DAF4A", "#377EB8", "#377EB8", "#377EB8"))

ddata_labels$group <- case_match(
  ddata_labels$label,
  "Species richness" ~ "Basic",
  c("Number of Indigenous names", "Number of recorded Indigenous uses") ~ "TEK",
  c("Faith", "Gini and Simpson", "Simpson",
    "Shannon", "Margalef", "Menhinick",
    "McIntosh", "Phylogenetic species variability", "Phylogenetic species richness") ~ "Phylogenetic",
  c("Functional species richness", "Functional species divergence", "Functional species dispersion", "Functional species evenness", "Rao's entropy (Q)") ~ "Functional"         
)

x_label_color_function <- function(string) {
  out <- #paste("<span style = 'color: ",
    case_match(string,
               "Basic" ~ "black",
               "TEK" ~ c4a("brewer.set1", 3)[1],
               "Phylogenetic" ~ c4a("brewer.set1", 3)[2],
               "Functional" ~ c4a("brewer.set1", 3)[3],)#,
  # ";'>",
  # string,
  # "</span>", sep = "")
  out <- as.character(out)
}

ddata_labels$color <- x_label_color_function(ddata_labels$group)

ddata_labels <- arrange(ddata_labels, label)

ggplot() +
  geom_segment(data = segment(ddata),
               aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
  geom_text(data = ddata_labels, 
            aes(x = x, y = y, label = label, color = label),
            fontface = "bold",hjust = 0, angle = 0) + 
  scale_color_manual(values = as.character(ddata_labels$color)) +
  scale_y_continuous(expand=c(0.2, 0),
                     trans = "reverse") + 
  coord_flip() +
  guides(color = "none") +
  theme_cowplot(12) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    )


ggdendrogram(ddata, rotate = FALSE) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 12, color = dendrogram_colors,
                                   face = "bold", hjust = 1))

# Clustered heatmap

plot_colours <- full_comp_df %>% 
  select(baseline_color, baseline_label) %>% 
  unique() %>% 
  # left_join(tibble(baseline_label = colnames(m)[pair_prec_cluster$order]))
  slice(match(colnames(m)[pair_prec_cluster$order], baseline_label)) %>% 
  select(baseline_color)


full_comp_df %>% 
  filter(type == "precision") %>% 
  ggplot(aes(baseline_label, comparison_label, fill = mean, colour = baseline_type)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(mean, 2)), color = "black") +
  # geom_blank() +
  scale_fill_gradient("Mean precision", low = "white", high = "blue") +
  scale_y_discrete("Comparison index", limits = colnames(m)[pair_prec_cluster$order]) +
  scale_x_discrete("Baseline index", limits = colnames(m)[pair_prec_cluster$order]) +
  # scale_color_manual("test", values = c(
  #   "Basic" = "black", "TEK" = "blue", "Phylogenetic" = "orange", "Functional" = "red")) +
  # guides("test", colour = guide_legend(override.aes = 
  #                                        list(fill = "white"))) +
  # ggtitle("Precision of comparison biodiversity metrics") +
  theme_cowplot() +
  theme(
    
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, 
                               color = plot_colours$baseline_color),
    axis.text.y = element_text(color = plot_colours$baseline_color)
    
  ) +
  guides(fill = guide_colorbar(barheight = 20,
                               title.hjust = 0.5))

# Table 1: Numbers of hot spots identified ---------------------------------

hotspot_nums <- readRDS('full_comparisons.rds') %>% 
  filter(type == "list_length") %>% 
  filter(baseline == "NumUnique") %>% 
  select(-c(baseline, type))

hotspot_nums$metric_label <-  case_match(hotspot_nums$comparison,
                                         "NumUnique" ~ "Species richness",
                                         "NumEndemic" ~ "Number of 'endemic' species",
                                         "NumIndigName" ~ "Number of Indigenous names",
                                         "NumUse" ~ "Number of recorded Indigenous uses",
                                         "FRic" ~ "Functional species richness",
                                         "FDiv" ~ "Functional species divergence",
                                         "FDis" ~ "Functional species dispersion",
                                         "FEve" ~ "Functional species evenness",
                                         "Q" ~ "Rao's entropy (Q)",
                                         "richness" ~ "Faith",
                                         "GiniSimpson" ~ "Gini and Simpson",
                                         "Simpson" ~ "Simpson",
                                         "Shannon" ~ "Shannon",
                                         "Margalef" ~ "Margalef",
                                         "Menhinick" ~ "Menhinick",
                                         "McIntosh" ~ "McIntosh",
                                         "PSVs" ~ "Phylogenetic species variability",
                                         "PSR" ~ "Phylogenetic species richness"
)



# Figure S1: Pairwise *recall* --------------------------------------------
pair_prec_mean <- readRDS('full_comparisons.rds') %>% 
  filter(type == "recall") %>% 
  select(-c(var, type)) %>% 
  pivot_wider(values_from = c("mean"),
              names_sort = F,
              names_from = "comparison")

pair_prec_mean_mat <- matrix(unlist(pair_prec_mean[1:18, 2:19]), ncol = 18)

# Check how far off this matrix is from being symmetric
isSymmetric.matrix(pair_prec_mean_mat) # It's not symmetric
# How far off is it? This value would be near zero if matrix is near symmetriic
max(abs(pair_prec_mean_mat-t(pair_prec_mean_mat))) 

pair_prec_plot <- readRDS('full_comparisons.rds') %>% 
  filter(type == "recall") %>% 
  ggplot(mapping = aes(x = baseline, y = comparison, fill = mean
                       #, alpha = rev(log(1+prec_var)))
  )) +
  geom_tile() +
  geom_text(aes(label = round(mean, 2)))  +
  scale_fill_gradient(low = "white", high = "blue")

# Make more descriptive labels for biodiversity metrics
full_comp_df <- readRDS('full_comparisons.rds') 

full_comp_df$baseline_label <-  case_match(full_comp_df$baseline,
                                           "NumUnique" ~ "Species richness",
                                           "NumEndemic" ~ "Number of 'endemic' species",
                                           "NumIndigName" ~ "Number of Indigenous names",
                                           "NumUse" ~ "Number of recorded Indigenous uses",
                                           "FRic" ~ "Functional species richness",
                                           "FDiv" ~ "Functional species divergence",
                                           "FDis" ~ "Functional species dispersion",
                                           "FEve" ~ "Functional species evenness",
                                           "Q" ~ "Rao's entropy (Q)",
                                           "richness" ~ "Faith",
                                           "GiniSimpson" ~ "Gini and Simpson",
                                           "Simpson" ~ "Simpson",
                                           "Shannon" ~ "Shannon",
                                           "Margalef" ~ "Margalef",
                                           "Menhinick" ~ "Menhinick",
                                           "McIntosh" ~ "McIntosh",
                                           "PSVs" ~ "Phylogenetic species variability",
                                           "PSR" ~ "Phylogenetic species richness"
)

full_comp_df$comparison_label <-  case_match(full_comp_df$comparison,
                                             "NumUnique" ~ "Species richness",
                                             "NumEndemic" ~ "Number of 'endemic' species",
                                             "NumIndigName" ~ "Number of Indigenous names",
                                             "NumUse" ~ "Number of recorded Indigenous uses",
                                             "FRic" ~ "Functional species richness",
                                             "FDiv" ~ "Functional species divergence",
                                             "FDis" ~ "Functional species dispersion",
                                             "FEve" ~ "Functional species evenness",
                                             "Q" ~ "Rao's entropy (Q)",
                                             "richness" ~ "Faith",
                                             "GiniSimpson" ~ "Gini and Simpson",
                                             "Simpson" ~ "Simpson",
                                             "Shannon" ~ "Shannon",
                                             "Margalef" ~ "Margalef",
                                             "Menhinick" ~ "Menhinick",
                                             "McIntosh" ~ "McIntosh",
                                             "PSVs" ~ "Phylogenetic species variability",
                                             "PSR" ~ "Phylogenetic species richness"
)


full_comp_df$baseline_type <- case_match(full_comp_df$baseline,
                                         c("NumUnique", "NumEndemic") ~ "Basic",
                                         c("NumIndigName", "NumUse") ~ "TEK",
                                         c("richness", "GiniSimpson", "Simpson",
                                           "Shannon", "Margalef", "Menhinick",
                                           "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
                                         c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"         
)

full_comp_df$comparison_type <- case_match(full_comp_df$comparison,
                                           c("NumUnique", "NumEndemic") ~ "Basic",
                                           c("NumIndigName", "NumUse") ~ "TEK",
                                           c("richness", "GiniSimpson", "Simpson",
                                             "Shannon", "Margalef", "Menhinick",
                                             "McIntosh", "PSVs", "PSR") ~ "Phylogenetic",
                                           c("FRic", "FDiv", "FDis", "FEve", "Q") ~ "Functional"         
)

# Remove number of "endemic" from analysis
full_comp_df <- filter(full_comp_df, baseline != "NumEndemic", comparison != "NumEndemic")

x_label_color_function <- function(string) {
  out <- #paste("<span style = 'color: ",
    case_match(string,
               "Basic" ~ "black",
               "TEK" ~ c4a("brewer.set1", 3)[1],
               "Phylogenetic" ~ c4a("brewer.set1", 3)[2],
               "Functional" ~ c4a("brewer.set1", 3)[3],)#,
  # ";'>",
  # string,
  # "</span>", sep = "")
  out <- as.character(out)
}

full_comp_df$baseline_color <- x_label_color_function(full_comp_df$baseline_type)
full_comp_df$comparison_color <- x_label_color_function(full_comp_df$comparison_type)

# Cluster based on Euclidean distance

pair_prec_mat <- full_comp_df %>%
  filter(type == "recall") %>% 
  select(-c(var, type, baseline, comparison, baseline_type, comparison_type,
            baseline_color, comparison_color)) %>% 
  pivot_wider(values_from = c("mean"),
              names_sort = F,
              names_from = "comparison_label")

m <- as.matrix((pair_prec_mat[, -1]), ncol = 18)

pair_prec_cluster <- hclust(dist(t(m)), method = "ward.D2")

# Cluster dendrogram
require(ggdendro)

dhc <- as.dendrogram(pair_prec_cluster)

ddata <- dendro_data(pair_prec_cluster, type="rectangle")

ggdendrogram(ddata, rotate = TRUE) +
  xlab('Biodiversity index') +
  theme(axis.text = element_text(size = 12)) 

dendrogram <- plot(pair_prec_cluster)

# Clustered heatmap

plot_colours <- full_comp_df %>% 
  select(baseline_color, baseline_label) %>% 
  unique() %>% 
  # left_join(tibble(baseline_label = colnames(m)[pair_prec_cluster$order]))
  slice(match(colnames(m)[pair_prec_cluster$order], baseline_label)) %>% 
  select(baseline_color)


full_comp_df %>% 
  filter(type == "recall") %>% 
  ggplot(aes(baseline_label, comparison_label, fill = mean, colour = baseline_type)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(mean, 2)), color = "black") +
  # geom_blank() +
  scale_fill_gradient("Mean recall", low = "white", high = "blue") +
  scale_y_discrete("Comparison index", limits = colnames(m)[pair_prec_cluster$order]) +
  scale_x_discrete("Baseline index", limits = colnames(m)[pair_prec_cluster$order]) +
  # scale_color_manual("test", values = c(
  #   "Basic" = "black", "TEK" = "blue", "Phylogenetic" = "orange", "Functional" = "red")) +
  # guides("test", colour = guide_legend(override.aes = 
  #                                        list(fill = "white"))) +
  # ggtitle("Precision of comparison biodiversity metrics") +
  theme_cowplot() +
  theme(
    
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, 
                               color = plot_colours$baseline_color),
    axis.text.y = element_text(color = plot_colours$baseline_color)
    
  ) +
  guides(fill = guide_colorbar(barheight = 20,
                               title.hjust = 0.5))
