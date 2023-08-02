library(reshape2)
library(ggplot2)
library(dplyr)
library(tweedie)
library(TMB)
library(Matrix)

# Set up actual data
earliest.year <- 1995

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


basin.2.use <- "south"
basin.dat <- dplyr::filter(thedata, basin == basin.2.use, YEAR>=earliest.year)
stocks.in.basin <- unique(basin.dat$stocklet)
nstocks <- length(stocks.in.basin)
ymat <- umat <- tmat <- matrix(NA, nrow = 27, ncol = length(stocks.in.basin))
nobs <- nyear <- numeric(nstocks)
put.in.matrix <- function(x, X, colno) {
  nx <- length(x)
  X[1:nx, i] <- x
  return(X)
}
# create spares mat
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

for (i in 1:length(stocks.in.basin)) {
  stocklet.2.use <- stocks.in.basin[i]
  dat <- dplyr::filter(thedata, stocklet == stocklet.2.use, YEAR>=earliest.year)
  data.2.use <- manage_data(dat)
  nyear[i] <- max(data.2.use$t)+1
  nobs[i] <- length(data.2.use$y)
  ymat <- put.in.matrix(x = data.2.use$y,X = ymat, colno= i)
  tmat <- put.in.matrix(data.2.use$t, tmat, i)
  umat <- put.in.matrix(rep(log(mean(data.2.use$y+1)), nyear[i]), umat, i)
}

# replace real data with simulated data with specified betas and process error
set.seed(124)

real_sigma_proc <- exp(rnorm(n = 4, mean = log(0.1), sd = 0.2))
real_beta <- rnorm(n = 4, mean = -0.15, sd = 0.1)
p <- 1.2
phi <- 5

for (i in 1:length(stocks.in.basin)) {
  nyear <- length(ymat[,i][(!is.na(ymat[,i]))==TRUE])
  u <- numeric(nyear)
  u[1] <- log(ymat[1,i])
  for (j in 2:nyear) u[j] <- u[j-1] + real_beta[i] + rnorm(1, mean = 0, sd = real_sigma_proc[i])
  fake_ys <- rtweedie(n = nyear, mu = exp(u), power = p, phi = phi)
  ymat[1:nyear,i] <- fake_ys
  umat[1:nyear,i] <- log(fake_ys+1)
}

  
    

# create sparse matrix for t
#sp.tmat <- create.sparse(tmat)
#sp.yobs <- create.sparse(ymat)

# get a vector of u's

non.na <- which(!is.na(umat), arr.ind = T)
uvec <- umat[non.na]
u.stock.id <- non.na[,2] -1
for (i in 2:ncol(tmat)) {
  tmat[,i] <- tmat[,i] + max(tmat[,i-1], na.rm = T)+1
}

# and a vector of ys
non.na <- which(!is.na(ymat), arr.ind = T)
yvec <- ymat[which(!is.na(ymat))]
y.stock.ids <- non.na[,2] -1
# make a vector of ts
tvec <- tmat[non.na]


# configure TMB model
setwd("R")
compile("exponential_state_space_bybasin.cpp")
dyn.load(dynlib("exponential_state_space_bybasin"))
regression.output <- readRDS("state_space_output.RDS")
data <- list(y = yvec,
             t = tvec,
             nstock = 4,
             ustockid =u.stock.id,
             ystockid = y.stock.ids)
parameters <- list(beta = regression.output[1:4,2],
                   logsigma_beta = log(sd(regression.output[1:4,2])),
                   mu_beta = mean(regression.output[1:4,2]),
                   u = uvec,
                   logsigma_proc = rep(log(1), nstocks),
                   logphi = rep(log(20), nstocks),
                   thetaf = rep(0, nstocks)
)
obj <- MakeADFun(data, parameters,random = c("u","beta"), DLL = "exponential_state_space_bybasin")
obj$hessian <- FALSE
opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
rep <- sdreport(obj)

fixed <- summary(rep, "fixed")
random <- summary(rep, "random")[,"Estimate"]

u <-  random[names(random)=="u"]
beta <-  random[names(random) == "beta"]
 
random_se <- summary(rep, "random")[,"Std. Error"]
u_se <- random_se[names(random_se)=="u"]
beta_se <- random_se[names(random_se) == "beta"]

 
 for (i in 1:nstocks) {
   stock.index <- which(y.stock.ids == i-1)
   ys <- yvec[stock.index]
   ts <- tvec[stock.index]
   plot(ts, ys, col = "#00000080", pch = 21, bg = "#00000030", las = 1,
              ylab = "abundance", xlab = "time", main = stocks.in.basin[i])
   stock.index <- which(u.stock.id == i-1)
   us <- u[stock.index]
   uses <- u_se[stock.index]
   all.years <- ts[1]:max(ts)
    lines(all.years, exp(us), col = "red", lty = 1, lwd = 1.5)
    polygon(c(all.years, rev(all.years)), c(exp(us - 2 * uses), rev(exp(us + 2 * uses))), 
            border = NA, col = "#FF000050")
   
 } 
 


# loop through each stock, each time fitting a linear model, save the fitted slope and p value in dataframe
# for (i in 1:n.stocks) {
# stocklet.2.use <- stocklet.names[i]
# dat <- dplyr::filter(thedata, stocklet == stocklet.2.use, YEAR>=earliest.year)
# data.2.use <- manage_data(dat)
# 
# data <- list(y = data.2.use$y,
#              t = data.2.use$t)
# parameters <- list(beta = 0,
#            u = rep(log(mean(data.2.use$y+1)), max(data.2.use$t)+1),
#            logsigma_proc = 0.5,
#            logphi = log(20),
#            thetaf = 0)
# 
# obj <- MakeADFun(data, parameters,random = "u", DLL = "exponential_state_space")
# obj$hessian <- FALSE
# opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
# rep <- sdreport(obj)
# 
# u <- summary(rep, "random")[, "Estimate"]
# 
# u_se <- summary(rep, "random")[, "Std. Error"]
# plot(data$t, data$y, col = "#00000080", pch = 21, bg = "#00000030", las = 1,
#      ylab = "abundance", xlab = "time", main = stocklet.2.use)
# all.years <- data$t[1]:max(data$t)
# lines(all.years, exp(u), col = "red", lty = 1, lwd = 1.5)
# polygon(c(all.years, rev(all.years)), c(exp(u - 2 * u_se), rev(exp(u + 2 * u_se))), 
#         border = NA, col = "#FF000050")
# # extract fixed effects:
# fixed <- summary(rep, "fixed")
# # get fixed effects of transformed variables
# transformed <- summary(rep, "report")
# 
# print(fixed)
# print(transformed)
# regression.output$beta[i] <- fixed[1,1]
# regression.output$beta_se[i] <- fixed[1,2]
# }
# 
# saveRDS(regression.output, file = "state_space_output.RDS")
# 
