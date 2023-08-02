library(reshape)
library(tidyr)
library(dplyr)


# Make plot of select state space estimates
# Load Herring data
wdfw.data.raw <- read.csv(file = "data/wdfw-data.csv", header = T)
thedata <- reshape2::melt(wdfw.data.raw, id = "YEAR", variable.name = "stocklet", value.name = "biomass")
thedata$basin <- rep(NA, times = nrow(thedata))



# Load state space results
stateSpaceResult <- readRDS(file = "analysis/state_space_estimates.RDS")

## Choose species for plotting

stock.2.plot <- c("PTGAMBLE", "QUILCENE", "SAMISH_PORTAGE", "CHERRYPT")
ylims <- matrix(c(0, 3000, 0,  7500, 0, 950, 0, 3800), nrow = 4, ncol = 2, byrow = T)
ytext <- 0.95*ylims[,2]
plottext <- c("Port Gamble", "Quilcene", "Samish / Portage", "Cherry Point")

png(filename = "graphics/example_trends.png",
    width = 6,
    height = 6,
    units = "in",
    res = 300)
earliest.year <- 1996
par(mfrow = c(2,2), las = 1, omi = c(.25, .5, .25, .25), mai = c(0.5, 0.5, 0.1, 0.1))
for (i in 1:length(stock.2.plot)) {
  stocklet.2.use <- stock.2.plot[i]
  stockletStateSpace <- as.data.frame(stateSpaceResult[[which(names(stateSpaceResult) == stock.2.plot[i])]])
  dat <- dplyr::filter(thedata, stocklet == stocklet.2.use, YEAR>=earliest.year)
  if(i ==1) ylab = "Spawning Biomass"
  if(i >1) ylab = ""
  
  plot(dat$YEAR, dat$biomass,
       type = "p",
       pch = 21,
       bg = "black",
       xlim = c(1995, 2021),
       ylim = ylims[i,],
       xlab = "",
       ylab = "",
       cex.axis = 1.2,
       cex.lab = 1.2
  )
  mtext(side = 2, text = ylab, line = 0, cex = 1, las = 0, outer = T)
  mtext(side = 1, text = "Year", line = 0, cex = 1, las = 0, outer = T)
  
  years <- stockletStateSpace$Year + earliest.year
  u <- stockletStateSpace$Est
  u_se <- stockletStateSpace$SE
  
  polygon(c(years, rev(years)), c(exp(u - 2 * u_se), rev(exp(u + 2 * u_se))), 
          border = NA, col = "gray")
  lines(years, exp(u), col = "black", lty = 1, lwd = 1.5)
  
  points(dat$YEAR, dat$biomass, pch = 21, bg = "black")
  text(x = 1996, y = ytext[i], labels = plottext[i],pos = 4)
}
dev.off()
