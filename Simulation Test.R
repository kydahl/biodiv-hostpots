rm(list = ls(all = TRUE))
graphics.off()

trait_number <- 4
global_species_number <- 10
regional_pool_num <- 5

regional_pool_prob <- 0.5

global_species_matrix <- matrix(nrow = global_species_number, ncol = trait_number)
regional_pool_matrix <-  matrix(nrow = regional_pool_num, ncol = global_species_number)

for(i in 1:global_species_number){
    
    for(j in 1:trait_number){
        
        global_species_matrix[i, j] <- rnorm(n = 1, mean = 0, sd = 1)
        
    }
    
}

for(i in 1:regional_pool_num){
    
    for(j in 1:global_species_number){
        
        regional_pool_chance <- runif(n = 1, min = 0, max = 1)
        
        regional_pool_matrix[i, j] <- regional_pool_chance < regional_pool_prob
        
    }
    
}