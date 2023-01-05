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
# Run first lines of "simulations.R" to create full_df

# - trait data 

# Trait data extracted from TRY
# (I added the date to not confuse different versions version)
# trait_data <- read_csv("data/Trait_data_20221130.csv")

# Trait data compiled by Diáz et al. 2020 (based on TRY)
# https://www.nature.com/articles/s41597-022-01774-9 
library(readxl)
trait_data <- read_excel("data/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx", 
                         sheet = 1)
trait_data_meta <- read_excel("data/Trait_data_TRY_Diaz_2022/Dataset/Species_mean_traits.xlsx", 
                              sheet = "ReadMe")
trait_data_ref <- read_excel("data/Trait_data_TRY_Diaz_2022/Dataset/References.xlsx", 
                             sheet = "References")
trait_data_ref_meta <- read_excel("data/Trait_data_TRY_Diaz_2022/Dataset/References.xlsx", 
                             sheet = "ReadMe")

#######################################################
# Clean trait dataset (in the case of data extracted from TRY)

#######################################################
# Create community matrix (patch x species)

si <- full_df %>%
  # mutate(Species2 = paste(stringr::word(Species, 1,1, sep=" "),
  #                         stringr::word(Species, 2,2, sep=" "),
  #                         sep="_")) %>%
  # mutate(Species_full2 = str_replace_all(Species_full," ","_")) %>%
  # dplyr::select(Species_full2, patch) %>% 
  dplyr::select(Species_full, patch) %>% 
  dplyr::mutate(Species_full = str_replace_all(Species_full,"ssp.","subsp.")) %>% 
  unique()

# comm <- dcast(si, formula = patch ~ Species_full2, fun.aggregate = length) %>%
#   column_to_rownames(var="patch")
comm <- dcast(si, formula = patch ~ Species_full, fun.aggregate = length) %>%
  column_to_rownames(var="patch")

#######################################################
# Match the community dataset with trait dataset

trait_data_sel <- trait_data %>%
  filter(`Species name standardized against TPL` %in% unique(si$Species_full))

length(unique(si$Species_full))
dim(trait_data_sel)
# --> there are quite a lot of species missing
# Are they actually missing or do they have a different name in the Diaz dataset?

missing_spp <- si %>%
  filter(!Species_full %in% trait_data$`Species name standardized against TPL`) %>%
  select(Species_full) %>%
  unique()
missing_spp

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
# 949     Populus balsamifera subsp. trichocarpa    species without subsp. is present
# 974                    Pseudoroegneria spicata    synonym "Elymus spicatus" is present
# 1006                            Rhodiola rosea    variant "Rhodiola rosea var. rosea" is present
# 1030                Rhododendron groenlandicum    synonym "Ledum palustre subsp. groenlandicum" is present
# 1053               Rhododendron neoglandulosum    synonym "Ledum glandulosum" is present
# 1083                              Ribes lobbii    missing, many other Ribes spp. are present. No synonyms present
# 1115                             Rubus pedatus    missing, many other Rubus spp. are present. No synonyms present
# 1142                           Sedum divergens    missing, many other Sedum spp. are present. No synonyms present
# 1166                            Sedum oreganum    missing, many other Sedum spp. are present. No synonyms present
# 1192                   Toxicodendron rydbergii    missing, many other Toxicodendron spp. are present. No synonyms present
# 1216                       Trillium petiolatum    missing, many other Trillium spp. are present  
# 1244           Vicia nigricans subsp. gigantea    species without subsp. is present
# 1267                       Zigadenus venenosus    synonym "Toxicoscordion venenosum" is present
# 1290                         Zigadenus elegans    synonym "Anticlea elegans" is present

#######################################################