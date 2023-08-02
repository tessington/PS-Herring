### Functions used in PS-Herring project
get_stock_huc <- function(filename) {
  ## get the HUCS associated with each stocklet ####
  stock.huc.matrix <- as.data.frame(read_excel(path = filename, 
                                               sheet = "Sites with HUC12 IDs Matrix", 
                                               range = "A4:W38", 
                                               col_names = TRUE))
  #### Pivot data ####
  stock.hucs <- pivot_longer(data =stock.huc.matrix,
                             cols = all_of(colnames(stock.huc.matrix)[-1]),
                             names_to = "stocklet")
  
  #### Remove NA and drop value column  ####
  stock.hucs <- stock.hucs %>%
    filter(value == 1) %>%
    select(HUC12, stocklet) %>%
    rename(StockletID = stocklet)
  return(stock.hucs)
}


### Function to match long stocklet to short stocklet names

convertNames <- function(longname) {
stockletNames <- tibble(long = c("Semiahmoo Bay", "Cherry Point" ,              
                                             "Samish/Portage Bay",  "Fidalgo Bay",                
                                             "Interior San Jan Islands",   "NW San Juan Islands",        
                                             "South Hood Canal" ,           "Quilcene Bay",               
                                             "Port Gamble"  ,               "Skagit",                     
                                             "Holmes Harbor" ,              "Port Susan",                 
                                             "Elliott Bay"    ,             "Quartermaster Harbor",       
                                             "Port Orchard – Port Madison", "Squaxin Pass",               
                                             "Purdy",                       "Wollochet",                  
                                             "Kilisut Harbor",              "Discovery Bay",              
                                             "Dungeness",                   "Sequim Bay"),
                           short = c(  "SEMIAHMOO" ,     "CHERRYPT",       "SAMISH_PORTAGE", "FIDALGO"  ,      "INT.SJI" ,      
                                               "NWSJI"      ,    "SOUTH_HC",       "QUILCENE"   ,    "PTGAMBLE" ,      "SKAGIT" ,       
                                               "HOLMES"      ,   "PTSUSAN",        "ELLIOTT"    ,    "QM"       ,      "PO.PM"  ,       
                                               "SQUAXIN"    ,    "PURDY",          "WOLLOCHET"  ,    "KILISUT"  ,      "DISCOVERY",     
                                               "DUNGENESS"  ,    "SEQUIM")
)

shortname <- stockletNames$short[stockletNames$long==longname]
return(shortname)
}

asinTransform <- function(p) asin(sqrt(p)) 
get_ratios <- function(x) {
  x <- x[which(!is.na(x))] # remove NAs
  nx <- length(x)
  first_to_last <- log(x[1] +1) - log(x[nx]+1) # adding 1 to deal with zeros
  ave_to_highest <- log(max(x) +1) -  mean(x[(nx-10):nx]+1) # mean over the last ten years of data, adding 1
  return(list(first_to_last = first_to_last, ave_to_highest = ave_to_highest))
}

calc_lb <- function(data_ts) {
  # first get the first year of data
  years <- data_ts$YEAR[which(!is.na(data_ts$biomass))]
  if (min(years) >= 2000) {
    lb.val = "no"
  } else {
    min.b <- min(data_ts$biomass, na.rm = T)
    min.b.recent <- min(data_ts$biomass[data_ts$YEAR >2010], na.rm = T)
    lb.val <- ifelse(min.b == min.b.recent, "yes", "no")
  }
  return(lb.val)
}

calc_zeros <- function(data_ts) {
  # get data past 1995
  data_ts <- dplyr::filter(data_ts, YEAR >= 1995)
  n.years <- length(which(!is.na(data_ts$biomass))) # how many years of data
  n.zeros <- length(which(data_ts$biomass ==0)) # how many years had 0 spawning biomass
  return(n.zeros / n.years)
  
  
}
run_regression <- function(wdfw.data) {
  
  ## stocknames
  stocklet.names <- c("SQUAXIN", 
                      "PURDY", 
                      "WOLLOCHET", 
                      "QM",
                      "ELLIOTT", 
                      "PO.PM",
                      "SOUTH_HC", 
                      "QUILCENE", 
                      "PTGAMBLE",
                      "KILISUT", 
                      "DISCOVERY", 
                      "SEQUIM", 
                      "DUNGENESS",
                      "PTSUSAN", 
                      "HOLMES", 
                      "SKAGIT",
                      "FIDALGO",
                      "SAMISH_PORTAGE",
                      "SEMIAHMOO",
                      "CHERRYPT",
                      "NWSJI",
                      "INT.SJI"
                      )
  
  
  # dataframe to hold output
  regression.output <- data.frame(
    stocklet = stocklet.names,
    slope = NA,
    se = NA,
    tval = NA
  )
  
  
  n.stocks <- length(stocklet.names)
  # loop through each stock, each time fitting a linear model, save the fitted slope and p value in dataframe
  for (i in 1:n.stocks) {
    stocklet.2.use <- stocklet.names[i]
    dat <-
      dplyr::filter(wdfw.data, stocklet == stocklet.2.use, YEAR >= earliest.year)
    # remove NA's from dat
    dat <- dat[complete.cases(dat), ]
    dat$YEAR <- dat$YEAR - 1996
    
    lm.output <- try(lm(biomass ~ YEAR, data = dat))
    se <- sqrt(diag(vcov(lm.output)))
    if (class(lm.output) != "try-error") {
      regression.output[i, 2] <- lm.output$coefficients[2]
      regression.output[i, 3] <- se[2]
      regression.output[i,4] <- summary(lm.output)$coef[2,3]
    }
  }
  return(regression.output)
}
calc_aicc <-function(fit) {
  k <- length(fit$coefficients)
  n <- length(fit$fitted.values)
  AIC <- AIC(fit)
  return(AIC + 2* k * (k+1) / (n - k - 1))
}
