# Install Packages ###
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(vegan)
library(gridExtra)
library(viridis)


# set up ggplot theme ####
theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             axis.line = element_line(color = "black"),
             axis.text = element_text(size = 14),
             axis.title= element_text(size = 16),
             strip.text = element_text(size = 14),
             legend.text = element_text(size = 12),
             strip.background = element_blank()
             )
             
             

# Functions  ####
source("R/PS-HERRING-FUNS.R")

# Load herring sub stock data ####
earliest.year <- 1996
wdfw.data.raw <- read.csv(file = "data/wdfw-data.csv", header = T)
wdfw.data <-
  reshape2::melt(
    wdfw.data.raw,
    id = "YEAR",
    variable.name = "stocklet",
    value.name = "biomass"
  )


# Run regression ####
regression.output <- run_regression(wdfw.data)
stocklet.names <- unique(wdfw.data$stocklet)
n.stocks <- length(stocklet.names)

# load state space fits #####
state.space.output <- readRDS("analysis/state_space_output.RDS")

# put state space output into regression output dataframe
regression.output$beta <- state.space.output$beta
regression.output$beta_se <- state.space.output$beta_se

#  Create df called "thedata" that contains sub stock metrics ####
thedata <- data.frame(stocklet = stocklet.names,
                      ave_to_highest = NA,
                      first_to_last = NA,
                      freq_zeros = NA,
                      lb = NA)



for (i in 1:n.stocks) {
  stocklet.2.use <- stocklet.names[i]
  data.2.use <- dplyr::filter(wdfw.data, stocklet == stocklet.2.use)
  ratios <- get_ratios(data.2.use$biomass)
  thedata$ave_to_highest[i] <- ratios$ave_to_highest
  thedata$first_to_last[i] <- ratios$first_to_last
  thedata$lb[i] <- calc_lb(data.2.use)
  thedata$freq_zeros[i] <- calc_zeros(data.2.use)
  
}
#  add to "thedata" the regression output, including state space beta

thedata <- inner_join(thedata, regression.output)

## Analyse covariance in response variables
response_vars <- dplyr::select(thedata,
                               tval, freq_zeros, ave_to_highest, first_to_last, beta)


# calculate and print out the correlation in response variables
print(cor(response_vars, use = "complete.obs"))


## Plot Metrics ####
saveplot <- T
if (saveplot) png(filename = "metric_summary.png", width = 800, 
                  height = 500,
                  units = "px")
col <- "gray"
thedata$lb <- as.factor(thedata$lb)
tvalplot <- ggplot(thedata, aes(x = tval)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 10,
                 color = "black", fill = col) +
  labs(x = "Regression Slope", y = "Count") 
  

zerosplot <- ggplot(thedata, aes(x = freq_zeros)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 10,
                 color = "black", fill = col) + 
  labs(x = "Frequency of 0's", y = "Count") 
  
ave_highestplot <- ggplot(thedata, aes(x = ave_to_highest)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 10, color = "black", fill = col) + 
  labs(x = "log(mean / maximum )", y = "Count") 
  
first_lastplot <- ggplot(thedata, aes(x = first_to_last)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 10, color = "black", fill = col) + 
  labs(x = "log(initial / final) ", y = "Count") 
  
betaplot <- ggplot(thedata, aes(x = beta)) +
  geom_histogram(position = "identity", alpha = 1.0, bins = 10, color = "black", fill = col) + 
  labs(x = "Population Growth Rate", y = "Count") 
 
lowestplot <- ggplot(thedata, aes(lb)) +
  geom_bar(colour = "black", fill = col) + 
  labs(x = "Lowest Biomass?", y = "Count") 
  
grid.arrange(tvalplot, zerosplot, ave_highestplot, first_lastplot, betaplot, lowestplot, nrow = 2) 
if (saveplot) dev.off()


# Run Watershed Analysis ####
# Load and run source files to load functions  ####

source(file = "R/landuse_by_stocklet_fun.R")
source(file = "R/human_by_stocklet_fun.R")
source(file = "R/Impervious_by_stocklet_fun.R")
stock.hucs <- get_stock_huc(filename = "data/Landscape characteristics by HUC12.xlsx")

lcuc.list <- landuse_by_stocklet()
red_lc_2016_tf <- asinTransform(lcuc.list$red_stocklet_lc_2016[,-1])
red_lc_1996_tf <- asinTransform(lcuc.list$red_stocklet_lc_1996[,-1])


# remove snow/ ice from change because there no change there
red_lc_1996_tf <- red_lc_1996_tf %>%
  select(-`Snow / Ice`)
red_lc_2016_tf <- red_lc_2016_tf %>%
  select(-`Snow / Ice`)

## Ordination ####
#### PCoA on 2016 Land cover / Land  use
distance <- vegdist(red_lc_2016_tf, method = "bray")
ord_predictor_2016 <- stats::cmdscale(distance, k = 2, eig = T)
colnames(ord_predictor_2016$points) <- c("axs1", "axs2")
ax1 <- ord_predictor_2016$points[,1]
ax2 <- ord_predictor_2016$points[,2]
ord_2016 <- tibble(pcoa.1 = ax1,
                   pcoa.2 = ax2)
ord.labels <- colnames(red_lc_2016_tf)

# get correlation between each transformed land cover type and each axis
getcor <- function(x,y) cor(x,y)
corr.ax1 <- apply(X = red_lc_2016_tf, MAR = 2, FUN = getcor, y = ax1)
corr.ax2 <- apply(X = red_lc_2016_tf, MAR = 2, FUN = getcor, y = ax2)

mult <- 0.2 # adjusts overall size of arrows
arrow.dat <- tibble(lc = ord.labels,
                    x.0 = 0,
                    y.0 = 0,
                    x.end = corr.ax1* mult ,
                    y.end = corr.ax2 * mult )


landcover_2016 <- ggplot(data = ord_2016, aes(x = pcoa.1, y = pcoa.2)) +
  geom_point(size = 2.5) + 
  xlab("PCoA Axis 1") + 
  ylab("PCoA Axis 2") +
  geom_segment(
    data = arrow.dat,
    aes(
      x = x.0,
      y = y.0,
      xend = x.end,
      yend = y.end,
      colour = lc
    ),
    linewidth = 1.25,
    # col = cols,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  scale_color_viridis_d(option = "turbo") + 
  labs(col = "") +
  theme_classic()+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12, colour = "black"), 
        axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        legend.text= element_text(size = 12),
        legend.title= element_text(size = 18),
        legend.position = c(0.15, 0.75),
        legend.background = element_blank()
  )


##### Ordination on landcover change ####
red_change_tf <- red_lc_2016_tf - red_lc_1996_tf
ord_predictor_change <- princomp(x = red_change_tf, scores = T, cor = T)
ord.cors <- ord_predictor_change$loadings

ord.labels <- colnames(red_change_tf)

mult <- 3 # adjusts overall size of arrows
arrow.dat <- tibble(lc = rownames(ord_predictor_change$loadings),
                    x.0 = 0,
                    y.0 = 0,
                    x.end = ord_predictor_change$loadings[,1] * mult,
                    y.end = ord_predictor_change$loadings[,2] * mult)
pca.scores <- as.data.frame(ord_predictor_change$scores)

landcover_change <- ggplot(data = pca.scores, aes(x = Comp.1, y = Comp.2)) +
  xlab("PCA Axis 1") + 
  ylab("PCA Axis 2") +
  geom_segment(
    data = arrow.dat,
    aes(
      x = x.0,
      y = y.0,
      xend = x.end,
      yend = y.end,
      colour = lc
    ),
    size = 1.25,
    # col = cols,
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  scale_color_viridis_d(option = "turbo") + 
  labs(col = "Landcover") +
  theme_classic()+
  theme(plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(size=12, colour = "black"), 
        axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        legend.position = "none"
  ) +
  geom_point(size = 2.5) 

library(gridExtra)
ordinations <- gridExtra::grid.arrange(landcover_2016, landcover_change, nrow = 1, ncol = 2)
ggsave("graphics/ordinations.png", ordinations, width = 12, height = 6, units = "in", device = "png")


# Relate stock status to watersheds ####
imp_by_stocklet <- impervious_by_stocklet()
density_by_stocklet <- human_by_stocklet()
# Run analyses on ordinated predictor variable
# create new, revised dataset that has coordinate 1 and coordinate 2 values together
# with the two uncorrelated response variables
# complicated code but it is just putting data all together
## Change long stocklet names to short stocklet names ####

longstocklet <- unique(stock.hucs$StockletID)
shortstocklet <- sapply(X = longstocklet, FUN = convertNames)

# new.df contains watershed features and stock status metrics for each sub stock
new.df <- inner_join(data.frame(stocklet = shortstocklet,
                                coord1 = ax1,
                                coord2 = ax2,
                                changecoord1 =  ord_predictor_change$scores[,1],
                                changecoord2 =  ord_predictor_change$scores[,2],
                                changecoord3 = ord_predictor_change$scores[,3],
                                changecoord4 = ord_predictor_change$scores[,4]),
                     dplyr::select(thedata,
                                   stocklet,
                                   lb,
                                   beta,
                                   beta_se))

# now add human density
density_by_stocklet <- rename(density_by_stocklet, stocklet = stockletID)
density_by_stocklet$stocklet <- sapply(X = density_by_stocklet$stocklet,
                                     FUN = convertNames)
new.df <- inner_join(new.df,
                     density_by_stocklet)

# add imperviousness
imp_by_stocklet <- rename(imp_by_stocklet, stocklet = stockletID)
imp_by_stocklet$stocklet <- sapply(X = imp_by_stocklet$stocklet,
                                   FUN = convertNames)

new.df <- inner_join(new.df,
                     imp_by_stocklet)

# Fit linear models ####

## For binary response variable "lowest biomass" ####
# convert Yes / No in "lowest_biomass" to 1 and 0
new.df$lb_num <- 1
new.df$lb <- as.character(new.df$lb) # remove factors from low biomass
new.df$lb_num[new.df$lb=="no"] <- 0
new.df$wt <- 1/new.df$beta_se


# model selection on the probability that lowest biomass is in recent time period

fit.null <- glm(lb_num ~ 1 , family = binomial, data = new.df)
fit.lc <- glm(lb_num ~ ImpervPer, family = binomial, data = new.df)
fit.imp <- glm(lb_num ~ ImpervPer, family = binomial, data = new.df)
fit.human <- glm(lb_num ~ density, family = binomial, data = new.df)
fit.lc.change <- glm(lb_num ~  changecoord1 + changecoord2 + changecoord3, family = binomial, data = new.df)
fit.imp.change <- glm(lb_num ~ changeImperv, family = binomial, data = new.df)
fit.human.change <- glm(lb_num ~ change, family = binomial, data = new.df)

AIClb <- matrix(NA, nrow = 7, ncol = 1)
AIClb[1,1] <- calc_aicc(fit.null)
AIClb[2,1] <- calc_aicc(fit.lc)
AIClb[3,1] <- calc_aicc(fit.imp)
AIClb[4,1] <- calc_aicc(fit.human)
AIClb[5,1] <- calc_aicc(fit.lc.change)
AIClb[6,1] <- calc_aicc(fit.imp.change)
AIClb[7,1] <- calc_aicc(fit.human.change)

rownames(AIClb) <- c("null", "landcover", "impervious", "human density", "change landcover", "change impervious", "change human density")
DAIClb <- AIClb - min(AIClb)
colnames(DAIClb) <- "delta AIC"
print(DAIClb)

## For Population Growth Rate ####

fit.null <- glm(beta ~ 1, data = new.df, weights = wt)
fit.lc <- glm(beta ~ ImpervPer,  data = new.df, weights = wt)
fit.imp <- glm(beta~ ImpervPer, data = new.df, weights = wt)
fit.human <- glm(beta ~ density,  data = new.df, weights = wt)
fit.lc.change <- glm(beta ~  changecoord1 + changecoord2 + changecoord3,  data = new.df, weights = wt)
fit.imp.change <- glm(beta ~ changeImperv, data = new.df, weights = wt)
fit.human.change <- glm(beta ~ change,  data = new.df, weights = wt)

AICbeta <- matrix(NA, nrow = 7, ncol = 1)
rownames(AICbeta) <- c("null", "landcover", "impervious", "human density", "change landcover", "change impervious", "change human density")

AICbeta[1,1] <- calc_aicc(fit.null)
AICbeta[2,1] <- calc_aicc(fit.lc)
AICbeta[3,1] <- calc_aicc(fit.imp)
AICbeta[4,1] <- calc_aicc(fit.human)
AICbeta[5,1] <- calc_aicc(fit.lc.change)
AICbeta[6,1] <- calc_aicc(fit.imp.change)
AICbeta[7,1] <- calc_aicc(fit.human.change)


DAICbeta <- AICbeta - min(AICbeta)
colnames(DAICbeta) <- "delta AIC"
print(DAICbeta)

## Make plots showing supported relationships ####
r1 <- ggplot(data = new.df, aes(x = density, y = lb_num)) +
  geom_point() +
  labs(x = expression(paste('2016 Human Density (# / ',km^2, ')')), y = "Lowest Biomass post 2010?")
  
r2 <- ggplot(data = new.df, aes(x = change, y = beta)) +
  geom_point() +
  labs(x = expression(paste('Change in Human Density (# / ',km^2, ')')), y = expression(paste('Population Growth Rate (',yr^-1, ')')))

saveplot <- F
reg_plot <- grid.arrange(r1, r2, nrow = 1, ncol = 2)
if (saveplot) ggsave(reg_plot,
                     filename = "graphics/plot_regression.png",
                  height = 350,
                  width = 700)

