library(reshape2)
library(ggplot2)
library(dplyr)
library(tweedie)
library(TMB)


earliest.year <- 1996

wdfw.data.raw <- read.csv(file = "data/wdfw-data.csv", header = T)
thedata <- reshape2::melt(wdfw.data.raw, id = "YEAR", variable.name = "stocklet", value.name = "biomass")
thedata$basin <- rep(NA, times = nrow(thedata))

## assign basins                     
south <- c("SQUAXIN", "PURDY", "WOLLOCHET", "QM")
central <- c("ELLIOTT", "PO.PM")
hoodcanal <- c("SOUTH_HC", "QUILCENE", "PTGAMBLE")
sanjuandefuca <- c("KILISUT", "DISCOVERY", "SEQUIM", "DUNGENESS")
whidbey <- c("PTSUSAN", "HOLMES", "SKAGIT")
north <- c("FIDALGO", "SAMISH_PORTAGE", "SEMIAHMOO", "CHERRYPT", "NWSJI", "INT.SJI")

south.index <- which(thedata$stocklet %in% south)
central.index <- which(thedata$stocklet %in% central) 
hoodcanal.index <- which(thedata$stocklet %in% hoodcanal)
sanjuandefuca.index <- which(thedata$stocklet %in% sanjuandefuca)
whidbey.index <- which(thedata$stocklet %in% whidbey)
north.index <- which(thedata$stocklet %in% north)
thedata$basin[south.index] <- "south"
thedata$basin[central.index] <- "central"
thedata$basin[hoodcanal.index] <- "hoodcanal"
thedata$basin[sanjuandefuca.index] <- "sanjaundefuca"
thedata$basin[whidbey.index] <- "whidbey"
thedata$basin[north.index] <- "north"

thedata$basin <- as.factor(thedata$basin)

# configure regression analysis
stocklet.names <- c(south, central, hoodcanal, sanjuandefuca, whidbey, north)
# dataframe to hold output


# configure TMB model

compile("R/exponential_state_space.cpp")
dyn.load(dynlib("R/exponential_state_space"))

n.stocks <- length(stocklet.names)
regression.output <- data.frame(stocklet = stocklet.names,
                                beta = NA,
                                beta_se = NA,
                                sigma_proc = NA)

manage_data <- function(dat) {
  # function to deal with NAs in the data from for each stock, returns a list containing array of observations, and year indices
  not.na.index <- which(!is.na(dat$biomass))
  years <- not.na.index[1] : max(not.na.index) 
  y <- dat$biomass[years]
  
  t <- years - years[1] # index for TMB to look up the year of each observation
  # now, remove any NAs that happen to be in the middle of the vector
  
  if(any(is.na(y))) {
    keep.index <- which(!is.na(y))
    y <- y[keep.index]
    t <- t[keep.index]
  }
  return(list(y = y, t = t))
}

# loop through each stock, each time fitting a linear model, save the fitted slope and p value in dataframe
for (i in 1:n.stocks) {
stocklet.2.use <- stocklet.names[i]
dat <- dplyr::filter(thedata, stocklet == stocklet.2.use, YEAR>=earliest.year)
data.2.use <- manage_data(dat)

data <- list(y = data.2.use$y,
             t = data.2.use$t)
parameters <- list(beta = 0,
           u = rep(log(mean(data.2.use$y+1)), max(data.2.use$t)+1),
           logsigma_proc = 0.5,
           logphi = log(20),
           thetaf = 0)

obj <- MakeADFun(data, parameters,random = "u", DLL = "exponential_state_space")
obj$hessian <- FALSE
opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
rep <- sdreport(obj)

u <- summary(rep, "random")[, "Estimate"]

u_se <- summary(rep, "random")[, "Std. Error"]
plot(data$t, data$y, col = "#00000080", pch = 21, bg = "#00000030", las = 1,
     ylab = "abundance", xlab = "time", main = stocklet.2.use)
all.years <- data$t[1]:max(data$t)
lines(all.years, exp(u), col = "red", lty = 1, lwd = 1.5)
polygon(c(all.years, rev(all.years)), c(exp(u - 2 * u_se), rev(exp(u + 2 * u_se))), 
        border = NA, col = "#FF000050")
# extract fixed effects:
fixed <- summary(rep, "fixed")
# get fixed effects of transformed variables
transformed <- summary(rep, "report")

print(fixed)
print(transformed)
regression.output$beta[i] <- fixed[1,1]
regression.output$beta_se[i] <- fixed[1,2]
regression.output$sigma_proc[i] <- transformed[1,1]
}

saveRDS(regression.output, file = "state_space_output.RDS")

