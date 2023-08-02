### Code to generate landscape attributes
### Uses the HUC codes associated with each stock
### Then grabs the landscape characteristics of each HUC
### Then gets a area-weighted total for each stock and each land use category
### This is repeated for 2016 and 1996
require(readxl)
require(tidyr)
require(dplyr)
source("R/PS-HERRING-FUNS.R")

## get the HUCS associated with each stocklet ####
stock.hucs <- get_stock_huc(filename = "data/Landscape characteristics by HUC12.xlsx")


## get the 2016 and 1996 HUC land cover ####
lc_1996 <- as.data.frame(read_excel(path = "data/Landscape characteristics by HUC12.xlsx", 
                                    sheet = "PIVOT C-CAP 1996", 
                                    range = "AB4:BA152",
                                    col_names = T))

lc_1996 <- dplyr::rename(lc_1996, HUC12 = `HUC12 ID`)


lc_2016 <- as.data.frame(read_excel(path = "data/Landscape characteristics by HUC12.xlsx", 
                                    sheet = "PIVOT C-CAP 2016", 
                                    range = "AA4:AZ152",
                                    col_names = T))
lc_2016 <- dplyr::rename(lc_2016, HUC12 = `HUC12 ID`)

## create a  loop to assign land use over to stockletsto generate area-weighted proportions ####
stocklets <- unique(stock.hucs$StockletID)

### create data frame to hold output ####
stocklet_lc_2016 <- stocklet_lc_1996 <- stocklet_change <- data.frame(
                              stockletID = stocklets,
                               BareLand = NA,
                               Cultivated = NA,
                               DeciduousForest= NA,
                               DevelopedOpenSpace = NA,
                               EstuarineAquaticBed = NA,
                               EstuarineEmergentWetland = NA,
                               EstuarineForestedWetland = NA,
                               EstuarineScrub_ShrubWetland = NA,
                               EvergreenForest = NA,
                               Grassland = NA,
                               HighIntensityDeveloped = NA,
                               LowIntensityDeveloped = NA,
                               MediumIntensityDeveloped = NA,
                               MixedForest = NA,
                               PalustrineAquaticBed = NA,
                               PalustrineEmergentWetland = NA,
                               PalustrineForestedWetland = NA,
                               PalustrineScrub_ShrubWetland = NA,
                               Pasture_Hay = NA,                   
                               Scrub_Shrub = NA,
                               Snow_Ice = NA,
                               UnconsolidatedShore = NA,
                               Water = NA,
                               NoDATA = NA)
### for each stocklet, find HUC codes and calculate proportion by land use category ####
for (i in 1:length(stocklets)) {
  
  # find the HUC codes associated with this stocklet
  huc_codes <- stock.hucs[stock.hucs$StockletID == stocklets[i],1]$HUC12
  
  if (length(huc_codes) == 1){
    stocklet_lc_2016[i,2:ncol(stocklet_lc_2016)] <- lc_2016[lc_2016$HUC12==huc_codes,3:ncol(lc_2016)]
    stocklet_lc_1996[i,2:ncol(stocklet_lc_1996)] <- lc_1996[lc_1996$HUC12==huc_codes,3:ncol(lc_1996)]
  } else {
    # extract areas
    areas <- lc_2016$`Area (km2)`[lc_2016$HUC12 %in% huc_codes]
    tot.area <- sum(areas)
    p.areas <- areas / tot.area
    for (j in 2:ncol(stocklet_lc_2016)) {
      stocklet_lc_2016[i,j] <- sum(p.areas * lc_2016[lc_2016$HUC12%in% huc_codes,j+1])
      stocklet_lc_1996[i,j] <- sum(p.areas * lc_1996[lc_1996$HUC12%in% huc_codes,j+1])
    }
  }
}


# 
# ## Aggregated land cover categories ####
lc_types <- colnames(stocklet_change)
forest <- grep("forest", lc_types, ignore.case= T)
developed <- grep("developed", lc_types, ignore.case = T)
agriculture <- which(lc_types %in% c("Cultivated", "Pasture_Hay"))
estuarine  <- grep("estuarine", lc_types, ignore.case= T)
palustrine <- grep("palustrine", lc_types, ignore.case = T)

#### Remove overlap between estuarine / pallustrine and forest, remove from forest ####
forest <- forest[-which(forest == intersect(forest, estuarine))]
forest <- forest[-which(forest == intersect(forest, palustrine))]


#### other land use cover types to keep ####
lc_2_keep <- c("BareLand" ,             
"Grassland",
 "Scrub_Shrub",
"Snow_Ice",
"UnconsolidatedShore",
"Water")

#### Condense 2016 land cover ####
stocklet_lc_2016$Forest <- rowSums(stocklet_lc_2016[,forest])
stocklet_lc_2016$Developed <- rowSums(stocklet_lc_2016[,developed])
stocklet_lc_2016$Agriculture <- rowSums(stocklet_lc_2016[,agriculture])
stocklet_lc_2016$Estuarine <- rowSums(stocklet_lc_2016[,estuarine])
stocklet_lc_2016$Palustrine <- rowSums(stocklet_lc_2016[,palustrine])

#### Condense 1996 land cover ####
stocklet_lc_1996$Forest <- rowSums(stocklet_lc_1996[,forest])
stocklet_lc_1996$Developed <- rowSums(stocklet_lc_1996[,developed])
stocklet_lc_1996$Agriculture <- rowSums(stocklet_lc_1996[,agriculture])
stocklet_lc_1996$Estuarine <- rowSums(stocklet_lc_1996[,estuarine])
stocklet_lc_1996$Palustrine <- rowSums(stocklet_lc_1996[,palustrine])

rm(forest, developed, agriculture, estuarine, palustrine)


red_stocklet_lc_2016 <-  stocklet_lc_2016 %>%
  select(stockletID, Forest, Developed, Agriculture, Estuarine, Palustrine, all_of(lc_2_keep)) %>%
  rename(`Scrub / Shrub` = Scrub_Shrub, 
         `Snow / Ice` = Snow_Ice,
         `Unconsolidated Shore` = UnconsolidatedShore)

red_stocklet_lc_1996 <-  stocklet_lc_1996 %>%
  select(stockletID, Forest, Developed, Agriculture, Estuarine, Palustrine, all_of(lc_2_keep))
red_stocklet_lc_1996 <- dplyr::rename(red_stocklet_lc_1996, `Scrub / Shrub` = Scrub_Shrub, 
         `Snow / Ice` = Snow_Ice,
         `Unconsolidated Shore` = UnconsolidatedShore)


