library(hypervolume)
library(foreach)
library(doParallel)
library(data.table)
library(tidyverse)

remove(list = ls(all.names = TRUE))

sampled_data <- fread("20210624_random_points_sampled_wPreds.csv") %>% 
  select(-Pixel_Lat, -Pixel_Long, -Resolve_Biome) %>% 
  rename_with(~str_remove(., '_mean_probability'))

taxonList <- sampled_data %>% select(starts_with("VTX"), starts_with(("SH"))) %>% names()

my.cluster <- parallel::makeCluster(
  72,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

output <- list() 

output <- foreach(i = taxonList, 
                  .combine="c",
                  .inorder=TRUE, 
                  .packages = c("hypervolume", "tidyverse"),
                  .errorhandling="pass",
                  .final = function(i) setNames(i, taxonList)
                  ) %dopar% {
                       taxon <- i
                       subsetData <- sampled_data %>% 
                         select(taxon, CHELSA_BIO_Annual_Mean_Temperature, SG_SOC_Content_005cm, CHELSA_BIO_Max_Temperature_of_Warmest_Month, CGIAR_PET, CHELSA_BIO_Precipitation_Seasonality, SG_Soil_pH_H2O_005cm, GlobBiomass_AboveGroundBiomass, SG_Sand_Content_005cm, EarthEnvTopoMed_Elevation) %>% 
                         unique() %>%
                         filter(!!sym(taxon) >= 0.9) %>% 
                         na.omit() %>% 
                         select(-taxon) %>% 
                         top_n(2500)
                       
                       if (nrow(subsetData) > 10){
                         hypervol <- hypervolume_gaussian(subsetData, name = taxon, verbose = FALSE)
                         
                         volume <- as.data.frame(hypervol@Volume) %>% rename(hypervolume = 1) 
                       } else {
                         volume <- NA
                       }
                       
                  }

output_df <- as.data.frame(output) %>% pivot_longer(everything()) %>% rename(volume = value, taxon = name)
write_csv(output_df, '20210705_hypervolume_results_predictions_mask90pct.csv')
parallel::stopCluster(cl = my.cluster)
