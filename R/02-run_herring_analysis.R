# Install Packages ###
library(reshape2)
library(ggplot2)
library(dplyr)
library(readxl)
library(vegan)
library(ape)
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

# Load data ####
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
state.space.output <- readRDS("R/state_space_output.RDS")

# put state space output into regression output dataframe
regression.output$beta <- state.space.output$beta
regression.output$beta_se <- state.space.output$beta_se

# loop through stocks, calculate metrics ####
thedata <- data.frame(stocklet = stocklet.names,
                      ave_to_highest = NA,
                      first_to_last = NA,
                      freq_zeros = NA,
                      lb = NA)


# now add to "thedata" the regression output, including state space beta
thedata <- inner_join(thedata, regression.output)
for (i in 1:n.stocks) {
  stocklet.2.use <- stocklet.names[i]
  data.2.use <- dplyr::filter(wdfw.data, stocklet == stocklet.2.use)
  ratios <- get_ratios(data.2.use$biomass)
  thedata$ave_to_highest[i] <- ratios$ave_to_highest
  thedata$first_to_last[i] <- ratios$first_to_last
  thedata$lb[i] <- calc_lb(data.2.use)
  thedata$freq_zeros[i] <- calc_zeros(data.2.use)
  
}
thedata <- inner_join(thedata, regression.output)

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
# Load and run source files to describe watersheds ####

source(file = "R/landuse_by_stocklet.R")
source(file = "R/human_by_stocklet.R")
source(file = "R/Impervious_by_stocklet.R")

red_lc_2016_tf <- asinTransform(red_stocklet_lc_2016[,-1])
red_lc_1996_tf <- asinTransform(red_stocklet_lc_1996[,-1])
red_change_tf <- red_lc_2016_tf - red_lc_1996_tf

# remove snow/ ice from change because there is little or no change there

red_change_tf <- red_change_tf %>%
  select(-`Snow / Ice`)
red_lc_2016_tf <- red_lc_2016_tf %>%
  select(-`Snow / Ice`)

distance <- vegdist(red_lc_2016_tf, method = "bray")
PCOA <- pcoa(distance)
ord_predictor_2016 <- stats::cmdscale(distance, k = 2, eig = T)
colnames(ord_predictor_2016$points) <- c("axs1", "axs2")
ax1 <- ord_predictor_2016$points[,1]
ax2 <- ord_predictor_2016$points[,2]

# get correlation between each transformed land cover type and each axis
getcor <- function(x,y) cor(x,y)
cor_ax1 <- apply(X = red_lc_2016_tf, MAR = 2, FUN = getcor, y = ax1)
cor_ax2 <- apply(X = red_lc_2016_tf, MAR = 2, FUN = getcor, y = ax2)
lcuc <- names(cor_ax1)
mult <- 0.3
plot_df2016 <- tibble(lcuc = lcuc,
                  ax1start = 0,
                  ax2start = 0,
                  ax1load = cor_ax1  * mult,
                  ax2load = cor_ax2 * mult)
pcoa_points <- tibble(stocklet = red_stocklet_lc_2016$stockletID,
                      axs1 = ax1,
                      axs2 = ax2
)
p <- ggplot(plot_df2016, aes(x = ax1load, y = ax2load, col = lcuc)) +
  geom_segment(aes(x = (ax1start), y = (ax2start), xend = (ax1load), yend = (ax2load)), arrow = arrow(),
               linewidth = 1.5,
               show.legend = FALSE) + 
  geom_line(aes(group = lcuc), linewidth = 1.5) + 
  scale_color_viridis_d(option = "turbo") +
  xlim(c(-0.5, 0.6)) + 
  ylim(c(-0.5,0.6)) +
  xlab("PCoA Axis 1") +
  ylab("PCoA Axis 2") +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8))
p1 <- p + geom_point(data = pcoa_points, aes(x = axs1, y = axs2),
               col = "black",
             show.legend = F)  

# repeat for change in land cover. Note here, distance measure is the euclidian distance of the correlation matrix of transformed proportional change, and we use PCA on this
# remove snow / ice because it doesn't change

ord_predictor_change <- princomp(x = red_change_tf, scores = T, cor = T)
ax1_loading<- ord_predictor_change$loadings[,1]
ax2_loading<- ord_predictor_change$loadings[,2]
x_change <- ord_predictor_change$scores[,1]
y_change <- ord_predictor_change$scores[,2]
mult = 5
plot_df <- tibble(lcuc = rownames(ord_predictor_change$loadings),
                  ax1start = 0,
                  ax2start = 0,
                  ax1load = ax1_loading * mult,
                  ax2load = ax2_loading * mult)
ord_points <- tibble(ord1 = ord_predictor_change$scores[,1],
                     ord2 = ord_predictor_change$scores[,2])

p2 <- ggplot(plot_df, aes(x = ax1load, y = ax2load, col = lcuc)) +
  geom_segment(aes(x = (ax1start), y = (ax2start), xend = (ax1load), yend = (ax2load)), arrow = arrow(),
               linewidth = 1.5,
               show.legend = FALSE) + 
  scale_color_viridis_d(option = "turbo") +
  xlab("PCA Axis 1") +
  ylab("PCA Axis 2") +
  geom_point(data = ord_points, aes(x = ord1, y = ord2),
             col = "black",
             show.legend = F)  
saveplot <- T
allplot <- grid.arrange(p1, p2, nrow = 1, ncol = 2)
if (saveplot) ggsave(filename = "landscape_ordination.png", plot = allplot)
# Run analyses on ordinated predictor variable
# create new, revised dataset that has coordinate 1 and coordinate 2 values together
# with the two uncorrelated response variables
# complicated code but it is just putting data all together

### Need to change long stocklet names to short stocklet names

longstocklet <- unique(stock.hucs$StockletID)

shortstocklet <- sapply(X = longstocklet, FUN = convertNames)

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
human_by_stocklet <- rename(human_by_stocklet, stocklet = stockletID)
human_by_stocklet$stocklet <- sapply(X = human_by_stocklet$stocklet,
                                     FUN = convertNames)
new.df <- inner_join(new.df,
                     human_by_stocklet)

# add imperviousness
imp_by_stocklet <- rename(imp_by_stocklet, stocklet = stockletID)
imp_by_stocklet$stocklet <- sapply(X = imp_by_stocklet$stocklet,
                                   FUN = convertNames)

new.df <- inner_join(new.df,
                     imp_by_stocklet)

# try some simple linear models, use weighted regression for the betas

# convert Yes / No in "lowest_biomass" to 1 and 0
new.df$lb_num <- NA
new.df$lb <- as.character(new.df$lb) # remove factors from low biomass
new.df$lb_num[new.df$lb=="no"] <- 0
new.df$lb_num[new.df$lb=="yes"] <- 1
new.df$wt <- 1/new.df$beta_se


# model selection on the probability that lowest biomass is in recent time period
AIClb <- rep(NA, times = 7)
fit.null <- glm(lb_num ~ 1 , family = binomial, data = new.df)
fit.lc <- glm(lb_num ~ ImpervPer, family = binomial, data = new.df)
fit.imp <- glm(lb_num ~ ImpervPer, family = binomial, data = new.df)
fit.human <- glm(lb_num ~ density, family = binomial, data = new.df)
fit.lc.change <- glm(lb_num ~  changecoord1 + changecoord2 + changecoord3, family = binomial, data = new.df)
fit.imp.change <- glm(lb_num ~ changeImperv, family = binomial, data = new.df)
fit.human.change <- glm(lb_num ~ change, family = binomial, data = new.df)


AIClb[1] <- calc_aicc(fit.null)
AIClb[2] <- calc_aicc(fit.lc)
AIClb[3] <- calc_aicc(fit.imp)
AIClb[4] <- calc_aicc(fit.human)
AIClb[5] <- calc_aicc(fit.lc.change)
AIClb[6] <- calc_aicc(fit.imp.change)
AIClb[7] <- calc_aicc(fit.human.change)

names(AIClb) <- c("null", "landcover", "impervious", "human density", "change landcover", "change impervious", "change human density")
DAIClb <- AIClb - min(AIClb)
print(DAIClb)

#summary(fit.human)
# model selection on estimated rate of population change

AICbeta <- rep(NA, times = 7)
names(AICbeta) <- c("null", "landcover", "impervious", "human density", "change landcover", "change impervious", "change human density")

fit.null <- glm(beta ~ 1, data = new.df, weights = wt)
fit.lc <- glm(beta ~ ImpervPer,  data = new.df, weights = wt)
fit.imp <- glm(beta~ ImpervPer, data = new.df, weights = wt)
fit.human <- glm(beta ~ density,  data = new.df, weights = wt)
fit.lc.change <- glm(beta ~  changecoord1 + changecoord2 + changecoord3,  data = new.df, weights = wt)
fit.imp.change <- glm(beta ~ changeImperv, data = new.df, weights = wt)
fit.human.change <- glm(beta ~ change,  data = new.df, weights = wt)

AICbeta[1] <- calc_aicc(fit.null)
AICbeta[2] <- calc_aicc(fit.lc)
AICbeta[3] <- calc_aicc(fit.imp)
AICbeta[4] <- calc_aicc(fit.human)
AICbeta[5] <- calc_aicc(fit.lc.change)
AICbeta[6] <- calc_aicc(fit.imp.change)
AICbeta[7] <- calc_aicc(fit.human.change)


DAICbeta <- AICbeta - min(AICbeta)
print(DAICbeta)
# make some plots
par(mfrow = c(1,2))
with(new.df, plot(density, lb_num,
                  pch = 21,
                  bg = "black",
                  xlab = "human density",
                  ylab = "is lowest biomass in recent period?",
                  las = 1)
)


with(new.df, plot(change, beta,
                  pch = 21,
                  bg = "black",
                  xlab = "rate of human popualtion change",
                  ylab = "rate of herring substock  change",
                  las = 1)
)
