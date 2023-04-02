##############################################################################
###                                                                        ###
###                       Calculating Functional diversity                 ###
###                                                                        ###
##############################################################################

# Elisa Van Cleemput, Jan 2023
#######################################################

# notes:


#######################################################

# Data used in this script:
# - full_df 
source("code/data_intake.R")
source("code/functions.R") 
source("code/parameters.R") 
full_df <- data_in %>%  
  filter(Species != "Pterospora andromedea") %>% 
  get.full_df(numPatches, mean.NumEntities, sd.NumEntities)

# - trait data 

# Trait data extracted from TRY
# (I added the date to not confuse different versions version)
# trait_data <- read_csv("data/Trait_data_20221130.csv")

# Trait data compiled by Diáz et al. 2020 (based on TRY)
# https://www.nature.com/articles/s41597-022-01774-9 
library(readxl)
trait_data <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx", 
                         sheet = 1)
trait_data_meta <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx", 
                              sheet = 2)
trait_data_ref <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/References.xlsx", 
                             sheet = "References")
trait_data_ref_meta <- read_excel("data/raw/Trait_data_TRY_Diaz_2022/Dataset/References.xlsx", 
                             sheet = 2)

#######################################################
## ---- Clean trait dataset (in the case of data extracted from TRY) --------------

#######################################################
## ---- Match the community dataset with the Diaz trait dataset --------------

# Create species list from the community dataset
si <- full_df %>%
  dplyr::select(Species_full, patch) %>% 
  dplyr::mutate(Species_full = str_replace_all(Species_full,"ssp.","subsp.")) %>% 
  unique()

# Check which species are present in the trait dataset
trait_data_sel <- trait_data %>%
  filter(`Species name standardized against TPL` %in% unique(si$Species_full))

length(unique(si$Species_full))
dim(trait_data_sel)
# --> there are quite a lot of species missing
# Are they actually missing or do they have a different name in the Diaz dataset? Let's check

missing_spp <- si %>%
  filter(!Species_full %in% trait_data$`Species name standardized against TPL`) %>%
  select(Species_full) %>%
  unique()
missing_spp # 54 missing species

# 1                   Chamaecyparis nootkatensis    synonym "Cupressus nootkatensis" is present
# 23             Equisetum hyemale subsp. affine    species without subsp. is present
# 44                          Adiantum aleuticum    synonym "Adiantum pedatum" (without var. aleuticum) is present
# 59     Athyrium filix-femina subsp. cyclosorum    species without subsp. is present  
# 82                      Polypodium glycyrrhiza    synonym "Polypodium vulgare" (without var. falcatum or var. occidentale) is present
# 101        Allium schoenoprasum var. sibiricum    species without subsp. is present  
# 115                Alnus viridis subsp. crispa    synonym "Alnus alnobetula subsp. crispa" is present
# 134               Alnus viridis subsp. sinuata    synonym "Alnus alnobetula subsp. sinuata" is present
# 159                         Angelica genuflexa    missing, many other Angelica spp. are present. No synonyms present
# 187                              Arctous ruber    synonym "Arctous alpina" (without var. rubra) is present
# 207                           Argentina egedii    synonym "Potentilla anserina subsp. groenlandica" is present
# 228                         Argentina anserina    synonym "Potentilla anserina" is present
# 259                            Asarum caudatum    missing, 4 other Asarum spp. are present. No synonyms present
# 290            Betula pumila var. glandulifera    species without subsp. is present  
# 313                        Boschniakia hookeri    missing. No synonyms present
# 340                        Boschniakia rossica    missing. No synonyms present
# 360                              Cirsium edule    missing, many other Cirsium spp. are present. No synonyms present
# 374                        Cirsium hookerianum    missing, many other Cirsium spp. are present. No synonyms present
# 402                      Conioselinum gmelinii    missing, 1 other Conioselinum sp. is present. No synonyms present
# 423                      Cornus unalaschkensis    synonym "Cornus canadensis" is present
# 443                    Dodecatheon pauciflorum    synonym "Dodecatheon meadia" is present
# 478                     Eriophorum chamissonis    missing, many other Eriophorum spp. are present. No synonyms present.
# 495                          Eurybia conspicua    synonym "Aster conspicuus" is present
# 530                            Frasera montana    missing, 4 other Frasera spp. are present. No synonyms present
# 558                 Fritillaria camschatcensis    missing, many other Fritillaria spp. are present. No synonyms present
# 578                             Glaux maritima    synonym "Lysimachia maritima" present
# 608                           Hackelia diffusa    missing, many other Hackelia spp. are present. No synonyms present
# 637                       Lappula occidentalis    synonym "Lappula redowskii" is present 
# 656                        Heuchera chlorantha    missing, many other Heuchera spp. are present. No synonyms present
# 683                           Hierochloe hirta    synonym "Hierochloe odorata" is present
# 706                          Ligusticum canbyi    missing, many other Lingusticum spp. are present. No synonyms present
# 729                            Lomatium geyeri    missing, many other Lomatium spp. are present. No synonyms present
# 754                       Lomatium macrocarpum    synonym "Peucedanum macrocarpum" is present. No synonyms present
# 782                         Mahonia aquifolium    missing, many other Mahonia spp. are present. No synonyms present
# 804                            Mahonia nervosa    synonym "Berberis nervosa" is present
# 823  Maianthemum racemosum subsp. amplexicaule    species without subsp. is present
# 849                      Phyllospadix scouleri    synonym "Phyllospadix iwatensis" is present
# 876                        Platanthera stricta    synoym "Platanthera dilatata" (without var. gracilis or var. viridiflora) is present
# 900                     Platanthera hyperborea    missing, many other Plantanthera spp. are present. No synonyms present
# 928     Populus balsamifera subsp. balsamifera    species without subsp. is present
# 949     Populus balsamifera subsp. trichocarpa    synonym "Populus trichocarpa" is present
# 974                    Pseudoroegneria spicata    synonym "Elymus spicatus" is present
# 1006                            Rhodiola rosea    variant "Rhodiola rosea var. rosea" is present
# 1030                Rhododendron groenlandicum    synonym "Ledum palustre subsp. groenlandicum" is present
# 1053               Rhododendron neoglandulosum    synonym "Ledum glandulosum" is present
# 1083                              Ribes lobbii    missing, many other Ribes spp. are present. No synonyms present
# 1115                             Rubus pedatus    missing, many other Rubus spp. are present. No synonyms present
# 1142                           Sedum divergens    missing, many other Sedum spp. are present. No synonyms present
# 1166                            Sedum oreganum    missing, many other Sedum spp. are present. No synonyms present
# 1192                   Toxicodendron rydbergii    synonym "Toxicodendron radicans" without (var. rydbergii) is present
# 1216                       Trillium petiolatum    missing, many other Trillium spp. are present  
# 1244           Vicia nigricans subsp. gigantea    species without subsp. is present
# 1267                       Zigadenus venenosus    synonym "Toxicoscordion venenosum" is present
# 1290                         Zigadenus elegans    synonym "Anticlea elegans" is present

# Replace species names with synonyms if possible
si_modif <- si %>%
  mutate(Species_full = ifelse(Species_full == "Chamaecyparis nootkatensis", "Cupressus nootkatensis", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Equisetum hyemale subsp. affine", "Equisetum hyemale", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Adiantum aleuticum", "Adiantum pedatum", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Athyrium filix-femina subsp. cyclosorum", "Athyrium filix-femina", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Polypodium glycyrrhiza", "Polypodium vulgare", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Allium schoenoprasum var. sibiricum", "Allium schoenoprasum", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Alnus viridis subsp. crispa", "Alnus alnobetula subsp. crispa", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Alnus viridis subsp. sinuata", "Alnus alnobetula subsp. sinuata", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Arctous ruber", "Arctous alpina", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Argentina egedii", "Potentilla anserina subsp. groenlandica", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Argentina anserina", "Potentilla anserina", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Betula pumila var. glandulifera", "Betula pumila", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Cornus unalaschkensis", "Cornus canadensis", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Dodecatheon pauciflorum", "Dodecatheon meadia", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Eurybia conspicua", "Aster conspicuus", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Glaux maritima", "Lysimachia maritima", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Lappula occidentalis", "Lappula redowskii", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Hierochloe hirta", "Hierochloe odorata", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Mahonia nervosa", "Berberis nervosa", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Maianthemum racemosum subsp. amplexicaule", "Maianthemum racemosum", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Phyllospadix scouleri", "Phyllospadix iwatensis", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Platanthera stricta", "Platanthera dilatata", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Populus balsamifera subsp. balsamifera", "Populus balsamifera", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Populus balsamifera subsp. trichocarpa", "Populus trichocarpa", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Pseudoroegneria spicata", "Elymus spicatus", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Rhodiola rosea", "Rhodiola rosea var. rosea", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Rhododendron groenlandicum", "Ledum palustre subsp. groenlandicum", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Rhododendron neoglandulosum", "Ledum glandulosum", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Toxicodendron rydbergii", "Toxicodendron radicans", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Vicia nigricans subsp. gigantea", "Vicia nigricans", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Zigadenus venenosus", "Toxicoscordion venenosum", Species_full)) %>%
  mutate(Species_full = ifelse(Species_full == "Zigadenus elegans", "Anticlea elegans", Species_full))
  

trait_data_sel <- trait_data %>%
  filter(`Species name standardized against TPL` %in% unique(si_modif$Species_full))

length(unique(si_modif$Species_full))
dim(trait_data_sel)

missing_spp <- si_modif %>%
  filter(!Species_full %in% trait_data$`Species name standardized against TPL`) %>%
  select(Species_full) %>%
  unique()
missing_spp # 22 missing species left

# Eriophorum chamissonis is present in self-collected TRY dataset

# Save the list of missing species
library(writexl)
trait_data <- write_xlsx(missing_spp,
                         "data/clean/Trait_data_TRY_Diaz_2022_PNW_missing.xlsx", col_names=TRUE)

#######################################################
# Create community matrix (patch x species)
comm <- dcast(si_modif, formula = patch ~ Species_full, fun.aggregate = length) %>%
  column_to_rownames(var="patch")


#######################################################
## ---- Calculate functional diversity --------------
library(FD)

# Select which traits to use for calculations
trait_data_sel_reduced <- trait_data_sel %>%
  rename(Species_name = `Species name standardized against TPL`) %>%
  select(Species_name, `Woodiness`, `Growth Form`, `Leaf area (mm2)`, `Nmass (mg/g)`, `LMA (g/m2)`, `Plant height (m)`, `Diaspore mass (mg)`, `LDMC (g/g)`, `SSD observed (mg/mm3)`, `SSD imputed (mg/mm3)`) %>%
  # select(Species_name, `Plant height (m)`, `Diaspore mass (mg)`) %>% # trying the dbFD algorithm with traits that are the most complete
  # select(Species_name, `Growth Form`, `Diaspore mass (mg)`) %>%
  column_to_rownames(var = "Species_name")

# Save a raw version of the cleaned dataset
library(writexl)
trait_data <- write_xlsx(trait_data_sel_reduced <- trait_data_sel %>%
                           rename(Species_name = `Species name standardized against TPL`) %>%
                           select(Species_name, `Woodiness`, `Growth Form`, `Leaf area (mm2)`, `Nmass (mg/g)`, `LMA (g/m2)`, `Plant height (m)`, `Diaspore mass (mg)`, `LDMC (g/g)`, `SSD observed (mg/mm3)`, `SSD imputed (mg/mm3)`),
                         "data/clean/Trait_data_TRY_Diaz_2022_PNW.xlsx", col_names=TRUE)

# Since we are experiencing problems with NAs below, we will delete incomplete rows for now
# This removes many species so we should try to solve this problem differently!
# trait_data_sel_reduced <- trait_data_sel_reduced[complete.cases(trait_data_sel_reduced), ] 
# KD: I think this removes all species!

# KD: I added this because I think the function below can only work with numeric
#     variables
trait_data_sel_reduced <- select(trait_data_sel_reduced, -c("Woodiness", "Growth Form"))

# Remove species from comm matrix that are not present in the trait dataset, 
# just for exploratory purposes at this time, but would be better to extend the dataset
comm_sel <- comm[, colnames(comm) %in% rownames(trait_data_sel_reduced)]
dim(comm_sel)
dim(trait_data_sel_reduced)

# order alphabetically to make sure both datasets have the same order of species
comm_sel <- comm_sel[, order(colnames(comm_sel))]
trait_data_sel_reduced <- trait_data_sel_reduced[order(rownames(trait_data_sel_reduced)),]

# Calculate functional diversity
# FRic = functional richness = convex hull volume (Villéger et al. 2008)
# Fdiv = functional divergence (Villéger et al. 2008)
# FDis = functional dispersion: weighted average distance to centroid (Laliberté and Legendre 2010). 
#        For communities composed of only one species, dbFD returns a FDis value of 0.
fdiv <- dbFD(trait_data_sel_reduced, as.matrix(comm_sel), w.abun=F, stand.x=T, 
           calc.FRic=T, m="max", 
           # calc.FGR=T,  clust.type="ward.D2", # this will ask for a manual decision on how to cut the tree
           calc.FDiv=T,
           calc.CWM= F)
# if not removing incomplete rows --> Error: NA's in the distance matrix.
# According to the help file, NAs are tolerated, but there are many NAs!
# The NA error message is related to NAs in the distance matrix which is calculated internally like gowdis(trait_data_sel_reduced)
# missing trait values can be imputed using the MICE (Multivariate Imputation by Chained Equations) algorithm with 
# predictive mean matching in the “mice” R package (van Buuren and Groothuis‐Oudshoorn, 2011).
# KD: I've also used the 'missForest' package to impute traits for primate species. 
#     With this package, it's one line of code to use random forests to fill in
#     missing traits.

str(fdiv)

#######################################################
