# get human population density for each stocklet
human_by_stocklet <- function() {
# ge the HUCS associated with each stocklet
stock.hucs <- get_stock_huc("data/Landscape characteristics by HUC12.xlsx")

# get the human population size by size by HUC

human.huc.1998 <- as.data.frame(read_excel(path = "data/Landscape characteristics by HUC12.xlsx", 
                                      sheet = "LSPOP 1998 by HUC12 VAT", 
                                      range = "A1:D6907",
                                      col_names = T))

human.huc.2016 <- as.data.frame(read_excel(path = "data/Landscape characteristics by HUC12.xlsx", 
                                           sheet = "LSPOP 2016 by HUC12 VAT", 
                                           range = "A1:D10085",
                                           col_names = T))


# summarize by HUC12_ID
human.huc.1998 <- human.huc.1998 %>%
  mutate(Pop = Lspop1998 * Count) %>%
  group_by(HUC12_ID) %>%
  summarise(Population = sum(Pop))

human.huc.2016 <- human.huc.2016 %>%
  mutate(Pop = Lspop2016 * Count) %>%
  group_by(HUC12_ID) %>%
  summarise(Population = sum(Pop))


## ENDED HERE ON 4/19 #####
# need to get area of each HUC

huc_area <- as.data.frame(read_excel(path = "data/Landscape characteristics by HUC12.xlsx", 
                                    sheet = "PIVOT C-CAP 2016", 
                                    range = "AA4:AB152",
                                    col_names = T))


stocklets <- unique(stock.hucs$StockletID)
density_by_stocklet <- data.frame(stockletID = stocklets,
                                census = NA,
                                density = NA,
                                change = NA)
for (i in 1:length(stocklets)) {
  
  # find the HUC codes associated with this stocklet
  huc_codes <- stock.hucs[stock.hucs$StockletID == stocklets[i],1]
  
  human_count <- sum(human.huc.2016$Population[human.huc.2016$HUC12_ID %in% huc_codes$HUC12])
  human_density <- human_count / sum(huc_area$`Area (km2)`[huc_area$`HUC12 ID` %in% huc_codes$HUC12])
  human_change <- log(human_count / sum(human.huc.1998$Population[human.huc.1998$HUC12_ID %in% huc_codes$HUC12]))/(2016 - 1998)
  density_by_stocklet[i,2:4] <- c(human_count, human_density, human_change)
}

return(density_by_stocklet)
}

  

