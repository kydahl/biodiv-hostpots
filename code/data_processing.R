###   Find typos in species names   ###

#install.packages("BIOMASS")
library(BIOMASS)
#This function corrects typos for a given taxonomic name using the Taxonomic Name Resolution Service (TNRS).
typos <- correctTaxo(veg.df$Species_name)
#no typos found

###   Match species names with databases   ###
library(taxize)
veg.df <- read.csv("Trait_data.csv")

#inspect the full list of available sources: 
View(gnr_datasources())

#specify the reference databases you want to use. If you don’t indicate databases, GNR will match your names against all sources.
src <- c("World Checklist of Vascular Plant", "The International Plant Names Index",
          "USDA NRCS PLANTS Database", "The International Plant Names Index")

#Show specified sources:
subset(gnr_datasources(), title %in% src)

#Run Global Names Resolver:
result.long <- veg.df$Species %>%
  gnr_resolve(data_source_ids = c(1,2,5,150, 165,167), 
              with_canonical_ranks=T)
result.long 

write.table(result.long,
            "result.long.txt", 
            sep="\t", row.names = F, quote = F)
#write.csv(result.long, "result.long.csv")

#create a short version of your results that is free of duplicates and that omits the names of the reference databases:
result.short <- result.long %>%
  select(submitted_name, matched_name2, score)%>%
  distinct()
result.short

#Write results to file:
write.table(result.short,
            "result.short.txt", 
            sep="\t", row.names = F, quote = F)

