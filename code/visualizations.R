#### Create visualizations from simulations ####################################

## Load libraries --------------------------------------------------------------
require(tidyverse)
library(doParallel)
library(MetBrewer)
require(cowplot)
require(GGally) # used to make paired scatter plots
library(retry)

# Load in simulations
source("code/functions.R")

# Figures -----------------------------------------------------------------

# Figure 0: Distribution of TEK traits (across all species)

# Plot histograms in a single column
trait_plot_no_impute <- read_csv("data/clean/dataset_no_imputation.csv") %>%
  relocate("Synonym", "Family", "Species", "Ssp_var",
           "N_Langs", "N_Names", "N_Uses") %>% 
  mutate(Woodiness = Woodiness == "woody") %>% 
  mutate(Leaf_area_log = log(`Leaf area (mm2)`), .keep = "unused") %>% 
  mutate(Plant_height_log = log(`Plant height (m)`), .keep = "unused") %>% 
  pivot_longer(cols = c("N_Langs":"Plant_height_log")) %>% 
  select(name, value) %>% 
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
  # one histogram for each biodiversity metric
  facet_wrap( ~ name, scales = "free") +
  theme_cowplot(24)

trait_plot_no_impute


# Plot histograms in a single column
trait_plot <- read_csv("data/clean/final_dataset.csv") %>%
  relocate("Synonym", "Family", "Species", "Ssp_var",
           "N_Langs", "N_Names", "N_Uses") %>% 
  mutate(Woodiness = Woodiness == "woody") %>% 
  mutate(Leaf_area_log = log(`Leaf area (mm2)`), .keep = "unused") %>% 
  mutate(Plant_height_log = log(`Plant height (m)`), .keep = "unused") %>% 
  pivot_longer(cols = c("N_Langs":"Plant_height_log")) %>% 
  select(name, value) %>% 
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
  # one histogram for each biodiversity metric
  facet_wrap( ~ name, scales = "free") +
  theme_cowplot(24)

trait_plot

# Compare distributions of traits before and after imputation (only functional traits)

compare_dists_plot <- read_csv("data/clean/final_dataset.csv") %>%
  mutate(origin = "imputed") %>% 
  rbind(read_csv("data/clean/dataset_no_imputation.csv") %>% 
          mutate(origin = "original")) %>% 
  select(-c(Synonym:Ssp_var)) %>% 
  relocate(origin) %>% 
  mutate(Woodiness = Woodiness == "woody") %>% 
  mutate(Leaf_area_log = log(`Leaf area (mm2)`), .keep = "unused") %>% 
  mutate(Plant_height_log = log(`Plant height (m)`), .keep = "unused") %>% 
  pivot_longer(cols = c("LDMC (g/g)":"Plant_height_log")) %>% 
  ggplot(aes(x = value, fill = origin)) +
  # plot histogram using density instead of count
  geom_histogram(aes(x = value, after_stat(density))) +
  # # add a spline approximating the probability density function
  # geom_density(alpha = .2) +
  # one histogram for each biodiversity metric
  facet_wrap( ~ name, scales = "free") +
  theme_cowplot(24)
compare_dists_plot


# Figure 1: Compare the distribution of biodiversity metric values --------
# Define the relevant trait names
trait_names <- c("LDMC (g/g)", "Nmass (mg/g)", "Woodiness", "Plant height (m)", "Leaf area (mm2)")

# Initialize the full phylogenetic tree
tree <- get.phylo_tree(read_csv("data/clean/final_dataset.csv"))

# Get a single simulation
numPatches <- 400
full_df_sample <- get.full_df(numPatches)

# Calculate biodiversity metrics
biodiv_df <- get.biodiv_df(full_df_sample, trait_names, tree)

# Plot histograms in a single column
histogram_plot <- biodiv_df %>%
  melt(id = "Patch") %>%
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
  # add a vertical line showing where the top 5% are
  stat_summary(aes(xintercept = after_stat(x), y = 0),
               fun = quantile, fun.args = list(0.95),
               geom = "vline", orientation = "y",
               color = "blue", lwd = 1) +
  # one histogram for each biodiversity metric
  facet_wrap(~variable, scales = "free") +
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
numPatches <- 400
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
