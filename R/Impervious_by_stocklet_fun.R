### Code to generate impervious levels for each HUC
### Then grabs the landscape characteristics of each HUC
### Then gets a area-weighted total for each stock and each land use category
### This is repeated for 2001 and 2016
impervious_by_stocklet <- function() {
source("R/PS-HERRING-FUNS.R")

## get the HUCS associated with each stocklet ####
stock.hucs <- get_stock_huc(filename = "data/Landscape characteristics by HUC12.xlsx")


## get the 2016 and 2001 HUC land cover ####
imp.2001 <- as.data.frame(read_excel(path = "data/Landscape characteristics by HUC12.xlsx", 
                                    sheet = "Impervious 2001 by HUC12 VAT", 
                                    range = "A1:D13623",
                                    col_names = T))


imp.2001 <- dplyr::rename(imp.2001, HUC12 = HUC12_ID)

imp.2001 <- imp.2001 %>%
  mutate(ImpervArea = Imperv2001 * Count) %>%
  group_by(HUC12) %>%
  summarise(ImpervPer = sum(ImpervArea)/ sum(Count), TotalArea = sum(Count))


imp.2016 <- as.data.frame(read_excel(path = "data/Landscape characteristics by HUC12.xlsx", 
                                     sheet = "Impervious 2016 by HUC12 VAT", 
                                     range = "A1:D13664",
                                     col_names = T))


imp.2016 <- dplyr::rename(imp.2016, HUC12 = HUC12_ID)

imp.2016 <- imp.2016 %>%
  mutate(ImpervArea = Imperv2016 * Count) %>%
  group_by(HUC12) %>%
  summarise(ImpervPer = sum(ImpervArea)/ sum(Count), TotalArea = sum(Count))


## create a  loop to assign land use over to stockletsto generate area-weighted proportions ####
stocklets <- unique(stock.hucs$StockletID)

### create data frame to hold output ####
imp_by_stocklet <- data.frame( stockletID = stocklets,
                               ImpervPer = NA,
                            changeImperv = NA
                            )
### for each stocklet, find HUC codes and calculate proportion by land use category ####
for (i in 1:length(stocklets)) {
  
  # find the HUC codes associated with this stocklet
  huc_codes <- stock.hucs[stock.hucs$StockletID == stocklets[i],1]$HUC12
  if (length(huc_codes) == 1){
    imp_by_stocklet[i,2] <- imp.2016[imp.2016$HUC12==huc_codes,2]
    imp_by_stocklet[i,3] <-  imp_by_stocklet[i,2] - imp.2001[imp.2001$HUC12==huc_codes,2]
  } else {
    # extract areas
    areas <- imp.2016$TotalArea[imp.2016$HUC12 %in% huc_codes]
    tot.area <- sum(areas)
    p.areas <- areas / tot.area
   imp_by_stocklet[i,2] <- sum(p.areas * imp.2016[imp.2016$HUC12%in% huc_codes,2])
   imp_by_stocklet[i,3] <- imp_by_stocklet[i,2] - sum(p.areas * imp.2001[imp.2001$HUC12%in% huc_codes,2])
}
}
return(imp_by_stocklet)
}
