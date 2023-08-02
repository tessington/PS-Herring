# import summary data
library(readxl)
library(dplyr)
library(vegan)

# get the regression.output
source("R/run_regressions.R")

# import rest of data
thedata <-read_excel(path = "data/HerringProjectR.xlsx")
names(thedata) <- tolower(names(thedata)) # just a trick to make all column names be spelled in lower case letters

# Create new column to hold t values and new slope
thedata$tval <- NA
thedata$logslope <- NA
# now add to "thedata" the t statistic of the linear regression, and the new unitless slope
for (i in 1:nrow(regression.output)) {
  stockid <- regression.output[i,1]
  index <- which(thedata$stocklet == stockid)
  thedata$tval[index] <- regression.output[i,4]
  thedata$logslope[index] <- regression.output[i,5]
}
ssModel <- readRDS("state_space_output.RDS")

# get correlation matrix of response variables
response_vars <- dplyr::select(thedata,
                        tval, percent_zeros, ave_to_highest, first_to_last, logslope)
print(cor(response_vars, use = "complete.obs"))


# This tells us that the regression and the average to highest biomass are telling us the same thing (very high correlation).  Also the 
#fancier regression gives us a metric that seems to cover all of our indices.  Good news!

# Now, let's look at the predictor variables
# First, pull out all of the predictor variables

# from here forward, remove the rows with NA
thedata <- thedata %>%
  select(-notes, -area, -huc12, - pic)
# remove rows with an NA
thedata <- thedata[complete.cases(thedata),]

predictor_vars <- dplyr::select(thedata, human_density, p_ag, p_residential_rural, p_residential_urban, p_urban_intensive, 
                          p_industrial_light, p_industrial_heavy)

# get correlations:
print(cor(predictor_vars))
# Super crazy correlated!  Consider ordination.  Try Non-Metric Dimensional Scaling (nMDS)


# Now run the MDS ordination to generate two dimensions

ord_predictor <- metaMDS(predictor_vars, distance = "bray", k = 3)

x <- ord_predictor$points[,1]
y <- ord_predictor$points[,2]
z <- ord_predictor$points[,3]
# view the sites in new ordination space
par(mfrow=c(2,2))
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n")
text(x, y, labels = thedata$stocklet, cex=.7)

plot(x, z, xlab="Coordinate 1", ylab="Coordinate 3",
     main="Metric MDS", type="n")
text(x, z, labels = thedata$stocklet, cex=.7)


plot(y, z, xlab="Coordinate 2", ylab="Coordinate 3",
     main="Metric MDS", type="n")
text(y, z, labels = thedata$stocklet, cex=.7)
par(mfrow= c(1,1))
# Run analyses on ordinated predictor variable
# create new, revised dataset that has coordinate 1 and coordinate 2 values together
# with the two uncorrelated response variables
# complicated code but it is just putting data all together
new.df <- inner_join(data.frame(stocklet = thedata$stocklet,
                     coord1 = x,
                     coord2 = y,
                     coord3 = z),
                     dplyr::select(thedata,
                                   stocklet,
                                   tval, 
                                   percent_zeros,
                                   first_to_last,
                                   low_biomass))


# try some simple linear models!

tval.fit <- lm(tval ~ coord1 + coord2+ coord3, data = new.df)
summary(tval.fit)

zeros.fit <- lm(percent_zeros ~ coord1 + coord2+ coord3, data = new.df)
summary(zeros.fit)

first.last.fit <- lm(first_to_last ~coord1 + coord2 + coord2, data = new.df)
summary(first.last.fit)

# convert Yes / No in "lowest_biomass" to 1 and 0
new.df$lb <- NA
new.df$lb[new.df$low_biomass=="No"] <- 0
new.df$lb[new.df$low_biomass=="Yes"] <- 1

lowest.biomass <- glm(lb ~ coord1 + coord2 + coord2 + coord3, family = binomial, data = new.df)
summary(lowest.biomass)

# make some plots
par(mfrow = c(2,2))
with(new.df, plot(coord1, tval,
                  pch = 21,
                  bg = "black",
                  xlab = "coordinate 1",
                  ylab = "t value of regression",
                  las = 1)
)
text(new.df$coord1, new.df$tval, labels = new.df$stocklet, cex = 0.7, pos = 3)


with(new.df, plot(coord2, tval,
                  pch = 21,
                  bg = "black",
                  xlab = "coordinate 2",
                  ylab = "t value of regression",
                  las = 1)
)
text(new.df$coord2, new.df$tval, labels = new.df$stocklet, cex = 0.7, pos = 3)

with(new.df, plot(coord3, tval,
                  pch = 21,
                  bg = "black",
                  xlab = "coordinate 3",
                  ylab = "t value of regression",
                  las = 1)
)
text(new.df$coord3, new.df$tval, labels = new.df$stocklet, cex = 0.7, pos = 3)


par(mfrow = c(2,2))
with(new.df, plot(coord1, percent_zeros,
                  pch = 21,
                  bg = "black",
                  xlab = "coordinate 1",
                  ylab = "Frequency of zeros",
                  las = 1)
)
text(new.df$coord1, new.df$percent_zeros, labels = new.df$stocklet, cex = 0.7, pos = 3)

    
with(new.df, plot(coord2, percent_zeros,
                  pch = 21,
                  bg = "black",
                  xlab = "coordinate 2",
                  ylab = "Frequency of zeros",
                  las = 1)
)
text(new.df$coord2, new.df$percent_zeros, labels = new.df$stocklet, cex = 0.7, pos = 3)

with(new.df, plot(coord3, percent_zeros,
                  pch = 21,
                  bg = "black",
                  xlab = "coordinate 3",
                  ylab = "Frequency of zeros",
                  las = 1)
)
text(new.df$coord3, new.df$percent_zeros, labels = new.df$stocklet, cex = 0.7, pos = 3)

