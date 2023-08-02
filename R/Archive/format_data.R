library(reshape2)
library(ggplot2)
library(dplyr)

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
# plot the data
p1 <- ggplot(data = thedata, aes(x = YEAR, y = biomass)) + 
  geom_point() + 
  geom_line()  + 
  facet_wrap(~stocklet, scales = "free")

p2 <- ggplot(data = thedata, aes(x = YEAR, y = biomass, col = stocklet)) + 
  geom_point() + 
  geom_line()  + 
  facet_wrap(~basin, scales = "free") +
  scale_color_discrete("viridis")
p1
p2

# filter out an individual stocklet, example using "ELLIOT" (elliot bay)
stocklet.2.use <- "ELLIOT" # you can make this anything
# create reduced data fram "dat" that only contains the stocklet that you want
dat <- dplyr::filter(thedata, stocklet == stocklet.2.use)

# you can do the same thing by basin
basin.2.use <- "south"
basin.dat <- dplyr::filter(thedata, basin == basin.2.use)

