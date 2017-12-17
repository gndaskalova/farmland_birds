# Bird responses to AES in Northeastern Scotland
# Written by Gergana Daskalova
# 25th Nov 2017
# gndaskalova@gmail.com

# Libraries ----
library(MCMCglmm)
library(Rmisc)
library(dplyr)
library(ggplot2)
library(simr)
library(lme4)
library(spaMM)
library(rgdal)
library(maptools)
library(maps)
library(OpenStreetMap)
library(raster)
library(gridExtra)
library(fields)
library(stargazer)

# Functions ----
# Custom ggplot2 theme
farm_theme <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title.x=element_text(size=16, face="plain"),             
          axis.title.y=element_text(size=16, face="plain"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=20, vjust=1, hjust=0.5),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_blank(),
          legend.background = element_rect(color = "black", fill = "transparent", 
                                           size = 4, linetype="blank"))
}

# With angled x axis labels
farm_theme2 <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=14, angle = 45, hjust = 1, vjust =1),
          axis.text.y=element_text(size=14),
          axis.title.x=element_text(size=16, face="plain"),             
          axis.title.y=element_text(size=16, face="plain"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=20, vjust=1, hjust=0.5),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_blank())
}

# Function to extract MCMCglmm model summary outputs ----
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

# Load data ----
# Bird abundance on AES and control farms
birds_aes <- read.csv2("birds_aes.csv")
head(birds_aes)
tail(birds_aes)
summary(birds_aes)

# Legend - AES agri-environmental scheme; control - conventional farming
# LI linnet; TS reed bunting
# S skylark; TS tree sparrow; Y yellowhammer

# Land management options data - presence/absence
birds_lmo <- read.csv("birds_lmo.csv")
summary(birds_lmo)

# Spatial data - farm lat and long and abundance
spatial <- read.csv("spatial.csv")
summary(spatial)

# Data manipulation ----
LI <- filter(birds_aes, species == "LI")
LI$latcenter <- LI$latitude - mean(LI$latitude)
LI$visits_f <- as.factor(LI$visits)
LI$loncenter <- LI$longitude - mean(LI$longitude)
LI$areacenter <- log(LI$area) - mean(log(LI$area))
summary(LI)

RB <- filter(birds_aes, species == "RB")
RB$latcenter <- RB$latitude - mean(RB$latitude)
RB$visits_f <- as.factor(RB$visits)
RB$loncenter <- RB$longitude - mean(RB$longitude)
RB$areacenter <- log(RB$area) - mean(log(RB$area))
summary(RB)

TS <- filter(birds_aes, species == "TS")
TS$latcenter <- TS$latitude - mean(TS$latitude)
TS$visits_f <- as.factor(TS$visits)
TS$loncenter <- TS$longitude - mean(TS$longitude)
TS$areacenter <- log(TS$area) - mean(log(TS$area))
summary(TS)

S <- filter(birds_aes, species == "S")
S$latcenter <- S$latitude - mean(S$latitude)
S$visits_f <- as.factor(S$visits)
S$loncenter <- S$longitude - mean(S$longitude)
S$areacenter <- log(S$area) - mean(log(S$area))
summary(S)

Y <- filter(birds_aes, species == "Y")
Y$latcenter <- Y$latitude - mean(Y$latitude)
Y$visits_f <- as.factor(Y$visits)
Y$loncenter <- Y$longitude - mean(Y$longitude)
Y$areacenter <- log(Y$area) - mean(log(Y$area))
summary(Y)

# AES effect on abundance ----

# Model structure: 
# abundance ~ treatment  + areacenter + latcenter + loncenter + visits_f; 
# random = farm + year_f
# Area, latitude and longitude centered 
# Their intercept estimates correspond to the average value of the predictor
# Duration is a continuous variable
# visits is a categorical variable

# Sol contains the distribution for the mean
# VCV contains the distribution for the variance

# Defining parameter-expanded priors
prior1 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))

prior2 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))

# ** Linnet ----
s_LI <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                 random = ~farm + year, pr = TRUE, nitt = 200000, 
                 data = LI, family = "poisson", prior = prior2, burnin = 20000)
summary(s_LI)

# No significant effect of treatment type on LI

# Calculating predicted values
exp(mean(s_LI$Sol[,"(Intercept)"] + s_LI$Sol[,"visits_f3"]))
exp(HPDinterval(s_LI$Sol[,"(Intercept)"] + s_LI$Sol[,"visits_f3"]))
exp(mean(s_LI$Sol[,"(Intercept)"] + (s_LI$Sol[,"treatmentcontrol"]) + s_LI$Sol[,"visits_f3"]))
exp(HPDinterval(s_LI$Sol[,"(Intercept)"] + (s_LI$Sol[,"treatmentcontrol"]) + s_LI$Sol[,"visits_f3"]))

# Checking autocorrelation and model convergence
autocorr(s_LI$VCV) 
autocorr(s_LI$Sol) 

plot(s_LI$VCV)
plot(s_LI$Sol)

# ** Reed bunting ----


s_RB <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                 random = ~farm + year, pr = TRUE, nitt = 200000, 
                 data = RB, family = "poisson", prior = prior2, burnin = 20000)
summary(s_RB)

# No significant effect of treatment type on RB

# Calculating predicted values
exp(mean(s_RB$Sol[,"(Intercept)"] + s_RB$Sol[,"visits_f3"]))
exp(HPDinterval(s_RB$Sol[,"(Intercept)"] + s_RB$Sol[,"visits_f3"]))
exp(mean(s_RB$Sol[,"(Intercept)"] + (s_RB$Sol[,"treatmentcontrol"]) + s_RB$Sol[,"visits_f3"]))
exp(HPDinterval(s_RB$Sol[,"(Intercept)"] + (s_RB$Sol[,"treatmentcontrol"]) + s_RB$Sol[,"visits_f3"]))

# Checking autocorrelation and model convergence
autocorr(s_RB$VCV) 
autocorr(s_RB$Sol) 

plot(s_RB$VCV)
plot(s_RB$Sol)

# ** Tree sparrow ----
s_TS <- MCMCglmm(abundance~treatment+areacenter+latcenter+loncenter+visits_f, random = ~farm+year, pr=TRUE, nitt=200000, data=TS, family="poisson", prior=prior2,burnin=20000)
summary(s_TS)

# No significant effect of treatment type on TS

# Calculating predicted values
exp(mean(s_TS$Sol[,"(Intercept)"] + s_TS$Sol[,"visits_f3"]))
exp(HPDinterval(s_TS$Sol[,"(Intercept)"] + s_TS$Sol[,"visits_f3"]))
exp(mean(s_TS$Sol[,"(Intercept)"] + (s_TS$Sol[,"treatmentcontrol"]) + s_TS$Sol[,"visits_f3"]))
exp(HPDinterval(s_TS$Sol[,"(Intercept)"] + (s_TS$Sol[,"treatmentcontrol"]) + s_TS$Sol[,"visits_f3"]))

# Checking autocorrelation and model convergence
autocorr(s_TS$VCV) 
autocorr(s_TS$Sol) 

plot(s_TS$VCV)
plot(s_TS$Sol)

# ** Skylark ----
s_S <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                random = ~farm + year, pr = TRUE, nitt = 200000, 
                data = S, family = "poisson", prior = prior2, burnin = 20000)
summary(s_S)

# No significant effect of treatment type on S

# Calculating predicted values
exp(mean(s_S$Sol[,"(Intercept)"] + s_S$Sol[,"visits_f3"]))
exp(HPDinterval(s_S$Sol[,"(Intercept)"] + s_S$Sol[,"visits_f3"]))
exp(mean(s_S$Sol[,"(Intercept)"] + (s_S$Sol[,"treatmentcontrol"]) + s_S$Sol[,"visits_f3"]))
exp(HPDinterval(s_S$Sol[,"(Intercept)"] + (s_S$Sol[,"treatmentcontrol"]) + s_S$Sol[,"visits_f3"]))

# Checking autocorrelation and model convergence
autocorr(s_S$VCV) 
autocorr(s_S$Sol) 

plot(s_S$VCV)
plot(s_S$Sol)

# ** Yellowhammer ----
s_Y <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                random = ~farm + year, pr = TRUE, nitt = 200000, 
                data = Y, family = "poisson", prior = prior2, burnin = 20000)
summary(s_Y)
# No significant effect of treatment type on Y

# Calculating predicted values
exp(mean(s_Y$Sol[,"(Intercept)"] + s_Y$Sol[,"visits_f3"]))
exp(HPDinterval(s_Y$Sol[,"(Intercept)"] + s_Y$Sol[,"visits_f3"]))

exp(mean(s_Y$Sol[,"(Intercept)"] + (s_Y$Sol[,"treatmentcontrol"]) + s_Y$Sol[,"visits_f3"]))
exp(HPDinterval(s_Y$Sol[,"(Intercept)"] + (s_Y$Sol[,"treatmentcontrol"]) + s_Y$Sol[,"visits_f3"]))

# Checking autocorrelation and model convergence
autocorr(s_Y$VCV) 
autocorr(s_Y$Sol) 

plot(s_Y$VCV)
plot(s_Y$Sol)

# Reordering factors in the treatment category to get estimates for AES not control
LI$treatment <- relevel(LI$treatment, "control", "AES")
LI_ordered <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                       random = ~farm + year, pr = TRUE, nitt = 600000, 
                       data = LI, family = "poisson", prior = prior2, burnin = 20000)
summary(LI_ordered)

RB$treatment <- relevel(RB$treatment, "control", "AES")
RB_ordered <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                       random = ~farm + year, pr = TRUE, nitt = 600000, 
                       data = RB, family = "poisson", prior = prior2, burnin = 20000)
summary(RB_ordered)

TS$treatment <- relevel(TS$treatment, "control", "AES")
TS_ordered <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                       random = ~farm + year, pr = TRUE, nitt = 600000, 
                       data = TS, family = "poisson", prior = prior2, burnin = 20000)
summary(TS_ordered)

S$treatment <- relevel(S$treatment, "control", "AES")
S_ordered <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                      random = ~farm + year, pr = TRUE, nitt = 600000, 
                      data = S, family = "poisson", prior = prior2, burnin = 20000)
summary(S_ordered)

Y$treatment <- relevel(Y$treatment, "control", "AES")
Y_ordered <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                      random = ~farm + year, pr = TRUE, nitt = 200000, 
                      data = Y, family = "poisson", prior = prior2, burnin = 20000)
summary(Y_ordered)

# Summary tables
LI_ordered_sum <- clean.MCMC(LI_ordered)
stargazer(LI_ordered_sum, type = "html", summary = FALSE, digits = 2)
RB_ordered_sum <- clean.MCMC(RB_ordered)
stargazer(RB_ordered_sum, type = "html", summary = FALSE, digits = 2)
TS_ordered_sum <- clean.MCMC(TS_ordered)
stargazer(TS_ordered_sum, type = "html", summary = FALSE, digits = 2)
S_ordered_sum <- clean.MCMC(S_ordered)
stargazer(S_ordered_sum, type = "html", summary = FALSE, digits = 2)
Y_ordered_sum <- clean.MCMC(Y_ordered)
stargazer(Y_ordered_sum, type = "html", summary = FALSE, digits = 2)

# Population trend models ----
# Models with year mean centered

# Linnet
LI$yearcenter <- LI$year - mean(LI$year)

s_LI3 <- MCMCglmm(abundance ~ treatment + yearcenter + areacenter + latcenter + loncenter + visits_f, 
                  random = ~farm + year, pr = TRUE, nitt = 600000, 
                  data = LI, family = "poisson", prior = prior2, burnin = 20000)
summary(s_LI3)

# Summary table
LI_sum3 <- clean.MCMC(s_LI3)
stargazer(LI_sum3, type = "html", summary = FALSE, digits = 2)

# Reed bunting
RB$yearcenter <- RB$year - mean(RB$year)

s_RB3 <- MCMCglmm(abundance ~ treatment + yearcenter + areacenter + latcenter + loncenter + visits_f, 
                  random = ~farm + year, pr = TRUE, nitt = 601000, 
                  data = RB, family = "poisson", prior = prior2, burnin = 20000)
summary(s_RB3)

# Summary table
RB_sum3 <- clean.MCMC(s_RB3)
stargazer(RB_sum3, type = "html", summary = FALSE, digits = 2)

# Tree sparrow
TS$yearcenter <- TS$year - mean(TS$year)

s_TS3 <- MCMCglmm(abundance ~ treatment + yearcenter + areacenter + latcenter + loncenter + visits_f, 
                  random = ~farm + year, pr = TRUE, nitt = 600000, 
                  data = TS, family = "poisson", prior = prior2, burnin = 20000)
summary(s_TS3)

# Summary table
TS_sum3 <- clean.MCMC(s_TS3)
stargazer(TS_sum3, type = "html", summary = FALSE, digits = 2)

# Skylark
S$yearcenter <- S$year - mean(S$year)

s_S3 <- MCMCglmm(abundance ~ treatment + yearcenter + areacenter + latcenter + loncenter + visits_f, 
                 random = ~farm + year, pr = TRUE, nitt = 600000, 
                 data = S, family = "poisson", prior = prior2, burnin = 20000)
summary(s_S3)

# Summary table
S_sum3 <- clean.MCMC(s_S3)
stargazer(S_sum3, type = "html", summary = FALSE, digits = 2)

# Yellowhammer
Y$yearcenter <- Y$year - mean(Y$year)

s_Y3 <- MCMCglmm(abundance ~ treatment + yearcenter + areacenter + latcenter + loncenter + visits_f, 
                 random = ~farm + year, pr = TRUE, nitt = 600000, 
                 data = Y, family = "poisson", prior = prior2, burnin = 20000)
summary(s_Y3)

# Summary table
Y_sum3 <- clean.MCMC(s_Y3)
stargazer(Y_sum3, type = "html", summary = FALSE, digits = 2)


# Model with year mean centered and a yearcenter*treatmetn interaction
# Linnet
s_LI2 <- MCMCglmm(abundance ~ treatment*yearcenter + areacenter + latcenter + loncenter + visits_f, 
                  random = ~farm + year, pr = TRUE, nitt = 600000, 
                  data = LI, family = "poisson", prior = prior2, burnin = 20000)
summary(s_LI2)

# Summary table
LI_sum2 <- clean.MCMC(s_LI2)
stargazer(LI_sum2, type = "html", summary = FALSE, digits = 2)


# Reed bunting
s_RB2 <- MCMCglmm(abundance ~ treatment*yearcenter + areacenter + latcenter + loncenter + visits_f, 
                  random = ~farm + year, pr = TRUE, nitt = 600000, 
                  data = RB, family = "poisson", prior = prior2, burnin = 20000)
summary(s_RB2)

# Summary table
RB_sum2 <- clean.MCMC(s_RB2)
stargazer(RB_sum2, type = "html", summary = FALSE, digits = 2)

# Tree sparrow
s_TS2 <- MCMCglmm(abundance ~ treatment*yearcenter + areacenter + latcenter + loncenter + visits_f, 
                  random = ~farm + year, pr = TRUE, nitt = 600000, 
                  data = TS, family = "poisson", prior = prior2, burnin = 20000)
summary(s_TS2)

# Summary table
TS_sum2 <- clean.MCMC(s_TS2)
stargazer(TS_sum2, type = "html", summary = FALSE, digits = 2)

# Skylark
s_S2 <- MCMCglmm(abundance ~ treatment*yearcenter + areacenter + latcenter + loncenter + visits_f, 
                  random = ~farm + year, pr = TRUE, nitt = 600000, 
                  data = S, family = "poisson", prior = prior2, burnin = 20000)
summary(s_S2)

# Summary table
S_sum2 <- clean.MCMC(s_S2)
stargazer(S_sum2, type = "html", summary = FALSE, digits = 2)

# Yellowhammer
s_Y2 <- MCMCglmm(abundance ~ treatment*yearcenter + areacenter + latcenter + loncenter + visits_f, 
                 random = ~farm + year, pr = TRUE, nitt = 600000, 
                 data = Y, family = "poisson", prior = prior2, burnin = 20000)
summary(s_Y2)

# Summary table
Y_sum2 <- clean.MCMC(s_Y2)
stargazer(Y_sum2, type = "html", summary = FALSE, digits = 2)


# Years in AES models ----

# Linnet
LI$AESyears <- LI$duration
LI$AESyears[which(LI$treatment == "control")] <- 0
summary(LI$AESyears)

LI_years <- MCMCglmm(abundance ~ treatment + AESyears + areacenter + latcenter + loncenter + visits_f, 
                     random = ~farm + year, pr = TRUE, nitt = 600000, 
                     data = LI, family = "poisson", prior = prior2, burnin = 20000)
summary(LI_years)

# Reed bunting
RB$AESyears <- RB$duration
RB$AESyears[which(RB$treatment == "control")] <- 0
summary(RB$AESyears)

RB_years <- MCMCglmm(abundance ~ treatment + AESyears + areacenter + latcenter + loncenter + visits_f, 
                     random = ~farm + year, pr = TRUE, nitt = 600000, 
                     data = RB, family = "poisson", prior = prior2, burnin = 20000)
summary(RB_years)

# Tree sparrow
TS$AESyears <- TS$duration
TS$AESyears[which(TS$treatment == "control")] <- 0
summary(TS$AESyears)

TS_years <- MCMCglmm(abundance ~ treatment + AESyears + areacenter + latcenter + loncenter + visits_f, 
                     random = ~farm + year, pr = TRUE, nitt = 600000, 
                     data = TS, family = "poisson", prior = prior2, burnin = 20000)
summary(TS_years)

# Skylark
S$AESyears <- S$duration
S$AESyears[which(S$treatment == "control")] <- 0
summary(S$AESyears)

S_years <- MCMCglmm(abundance ~ treatment + AESyears + areacenter + latcenter + loncenter + visits_f, 
                    random = ~farm + year, pr = TRUE, nitt = 600000, 
                    data = TS, family = "poisson", prior = prior2, burnin = 20000)
summary(S_years)

# Yellowhammer
Y$AESyears <- Y$duration
Y$AESyears[which(Y$treatment == "control")] <- 0
summary(Y$AESyears)

Y_years <- MCMCglmm(abundance ~ treatment + AESyears + areacenter + latcenter + loncenter + visits_f, 
                    random = ~farm + year, pr = TRUE, nitt = 600000, 
                    data = Y, family = "poisson", prior = prior2, burnin = 20000)
summary(Y_years)

# Summary tables
LI_years_sum <- clean.MCMC(LI_years)
stargazer(LI_years_sum, type = "html", summary = FALSE, digits = 2)
RB_years_sum <- clean.MCMC(RB_years)
stargazer(RB_years_sum, type = "html", summary = FALSE, digits = 2)
TS_years_sum <- clean.MCMC(TS_years)
stargazer(TS_years_sum, type = "html", summary = FALSE, digits = 2)
S_years_sum <- clean.MCMC(S_years)
stargazer(S_years_sum, type = "html", summary = FALSE, digits = 2)
Y_years_sum <- clean.MCMC(Y_years)
stargazer(Y_years_sum, type = "html", summary = FALSE, digits = 2)

# Offset MCMC models ----
# Specify priors for the fixed effect, make it a slope of 1 for log area
prior3 <- list(B = list (mu = matrix(c(0,0,1,0,0,0),6),V = diag(6)*(10)), 
               R = list(V = 1, nu=0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))
diag(prior3$B$V)[3] <- 1e-9 

# Linnet
s_LI_off <- MCMCglmm(abundance ~ treatment + log(area) + latcenter + loncenter + visits_f, 
                     random = ~ farm + year, pr=TRUE, nitt = 600000, 
                     burnin = 20000, data = LI, family = "poisson", prior = prior3)
summary(s_LI_off)

# Reed bunting
s_RB_off <- MCMCglmm(abundance ~ treatment + log(area) + latcenter + loncenter + visits_f, 
                     random = ~ farm + year, pr = TRUE, nitt = 600000, 
                     burnin = 20000, data = RB, family = "poisson", prior = prior3)
summary(s_RB_off)

# Tree sparrow
s_TS_off <- MCMCglmm(abundance ~ treatment + log(area) + latcenter + loncenter + visits_f, 
                     random = ~ farm + year, pr = TRUE, nitt = 600000, 
                     burnin = 20000, data = TS, family = "poisson", prior = prior3)
summary(s_TS_off)

# Skylark
s_S_off <- MCMCglmm(abundance ~ treatment + log(area) + latcenter + loncenter + visits_f, 
                    random = ~ farm + year, pr = TRUE, nitt = 600000, 
                    burnin = 20000, data = S, family = "poisson", prior = prior3)
summary(s_S_off)

# Yellowhammer
s_Y_off <- MCMCglmm(abundance ~ treatment + log(area) + latcenter + loncenter + visits_f, 
                    random = ~ farm + year, pr = TRUE, nitt = 600000, 
                    burnin = 20000, data = Y, family = "poisson", prior = prior3)
summary(s_Y_off)

# Summary tables
s_LI_off_sum <- clean.MCMC(s_LI_off)
stargazer(s_LI_off_sum, type = "html", summary = FALSE, digits = 2)
s_RB_off_sum <- clean.MCMC(s_RB_off)
stargazer(s_RB_off_sum, type = "html", summary = FALSE, digits = 2)
s_TS_off_sum <- clean.MCMC(s_TS_off)
stargazer(s_TS_off_sum, type = "html", summary = FALSE, digits = 2)
s_S_off_sum <- clean.MCMC(s_S_off)
stargazer(s_S_off_sum, type = "html", summary = FALSE, digits = 2)
s_Y_off_sum <- clean.MCMC(s_Y_off)
stargazer(s_Y_off_sum, type = "html", summary = FALSE, digits = 2)


# Spatial modelling of bird abundance on AES and control farms ----

# Turning lat and long into easting and northing
# Coerce to a spatial object
coordinates(spatial) = c('longitude','latitude')
plot(spatial)

# Transformation
wgs84 = '+proj=longlat +datum=WGS84'
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
mrc = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'

# Specify coordinate system as WGS84 (typical latlons)
spatial@proj4string   # slot will be empty
spatial@proj4string = CRS(wgs84)

# reproject data to Mercator coordinate system
spatial_merc = spTransform(spatial, CRS(mrc))

# download Mercator tiles covering the data area
bbox = spatial@bbox
map = openmap(c(bbox[2,2],bbox[1,1]), c(bbox[2,1],bbox[1,2]), 7, 'osm')

plot(map)
plot(spatial_merc, add = T)

# conversion to British National Grid (OSGB36)
map_bng = projectRaster(raster(openproj(map)), crs = bng)
spatial_bng = spTransform(spatial, CRS(bng))

plotRGB(map_bng)
plot(spatial_bng, add=T)

spatial.east.north <- as.data.frame(spatial_bng)

names(spatial.east.north)[names(spatial.east.north) == 'latitude'] <- 'easting'
names(spatial.east.north)[names(spatial.east.north) == 'longitude'] <- 'northing'

spatial2 <- spatial.east.north

# Linnet
linnet <- filter(spatial2, species == "LI")
linnet$farm <- as.factor(linnet$farm)
linnet$treatment <- as.factor(linnet$treatment)
linnet$visits <- as.factor(linnet$visits)
linnet$year <- as.factor(linnet$year)
linnet$abundance <- as.numeric(linnet$abundance)
linnet$areacenter <- log(linnet$area) - mean(log(linnet$area))
linnet <- as.data.frame(linnet)

li.spamm <- corrHLfit(abundance ~ treatment + areacenter + visits + Matern(1|northing + easting), 
                      HLmethod = "HL(0, 1)", data = linnet, family = poisson(), 
                      ranFix = list(nu = 0.5))
summary(li.spamm)

# Reed bunting
reed.bunting <- filter(spatial2, species == "RB")
reed.bunting$farm <- as.factor(reed.bunting$farm)
reed.bunting$treatment <- as.factor(reed.bunting$treatment)
reed.bunting$visits <- as.factor(reed.bunting$visits)
reed.bunting$year <- as.factor(reed.bunting$year)
reed.bunting$abundance <- as.numeric(reed.bunting$abundance)
reed.bunting$areacenter <- log(reed.bunting$area) - mean(log(reed.bunting$area))
reed.bunting <- as.data.frame(reed.bunting)

rb.spamm <- corrHLfit(abundance ~ treatment + areacenter + visits + Matern(1|northing + easting), 
                      HLmethod = "HL(0, 1)", data = reed.bunting, family = poisson(), 
                      ranFix = list(nu = 0.5))
summary(rb.spamm)

# Tree sparrow
tree.sparrow <- filter(spatial2, species == "TS")
tree.sparrow$farm <- as.factor(tree.sparrow$farm)
tree.sparrow$treatment <- as.factor(tree.sparrow$treatment)
tree.sparrow$visits <- as.factor(tree.sparrow$visits)
tree.sparrow$year <- as.factor(tree.sparrow$year)
tree.sparrow$abundance <- as.numeric(tree.sparrow$abundance)
tree.sparrow$areacenter <- log(tree.sparrow$area) - mean(log(tree.sparrow$area))
tree.sparrow <- as.data.frame(tree.sparrow)

ts.spamm <- corrHLfit(abundance ~ treatment + areacenter + visits + Matern(1|northing + easting), 
                      HLmethod = "HL(0, 1)", data = tree.sparrow, family = poisson(), 
                      ranFix = list(nu = 0.5))
summary(ts.spamm)

# Skylark
skylark <- filter(spatial2, species == "S")
skylark$farm <- as.factor(skylark$farm)
skylark$treatment <- as.factor(skylark$treatment)
skylark$visits <- as.factor(skylark$visits)
skylark$year <- as.factor(skylark$year)
skylark$abundance <- as.numeric(skylark$abundance)
skylark$areacenter <- log(skylark$area) - mean(log(skylark$area))
skylark <- as.data.frame(skylark)

s.spamm <- corrHLfit(abundance ~ treatment + areacenter + visits + Matern(1|northing + easting), 
                     HLmethod = "HL(0, 1)", data = skylark, family = poisson(), 
                     ranFix = list(nu = 0.5))
summary(s.spamm)

# Yellowhammer
yellowhammer <- filter(spatial2, species == "Y")
yellowhammer$farm <- as.factor(yellowhammer$farm)
yellowhammer$treatment <- as.factor(yellowhammer$treatment)
yellowhammer$visits <- as.factor(yellowhammer$visits)
yellowhammer$year <- as.factor(yellowhammer$year)
yellowhammer$abundance <- as.numeric(yellowhammer$abundance)
yellowhammer$areacenter <- log(yellowhammer$area) - mean(log(yellowhammer$area))
yellowhammer <- as.data.frame(yellowhammer)

y.spamm <- corrHLfit(abundance ~ treatment + areacenter + visits + Matern(1|northing + easting), 
                     HLmethod = "HL(0, 1)", data = yellowhammer, family = poisson(), 
                     ranFix = list(nu = 0.5))
summary(y.spamm)

# Multi-membership LMO models ----
LI_lmo <- filter(birds_lmo, species == "LI", treatment == "AES")
LI_lmo$latcenter <- LI_lmo$latitude - mean (LI_lmo$latitude)
LI_lmo$visits_f <- as.factor(LI_lmo$visits)
LI_lmo$loncenter <- LI_lmo$longitude - mean(LI_lmo$longitude)
LI_lmo$areacenter <- log(LI_lmo$area) - mean(log(LI_lmo$area))

RB_lmo <- filter(birds_lmo, species == "RB", treatment == "AES")
RB_lmo$latcenter <- RB_lmo$latitude - mean (RB_lmo$latitude)
RB_lmo$visits_f <- as.factor(RB_lmo$visits)
RB_lmo$loncenter <- RB_lmo$longitude - mean(RB_lmo$longitude)
RB_lmo$areacenter <- log(RB_lmo$area) - mean(log(RB_lmo$area))

TS_lmo <- filter(birds_lmo, species == "TS", treatment == "AES")
TS_lmo$latcenter <- TS_lmo$latitude - mean (TS_lmo$latitude)
TS_lmo$visits_f <- as.factor(TS_lmo$visits)
TS_lmo$loncenter <- TS_lmo$longitude - mean(TS_lmo$longitude)
TS_lmo$areacenter <- log(TS_lmo$area) - mean(log(TS_lmo$area))

S_lmo <- filter(birds_lmo, species == "S", treatment == "AES")
S_lmo$latcenter <- S_lmo$latitude - mean (S_lmo$latitude)
S_lmo$visits_f <- as.factor(S_lmo$visits)
S_lmo$loncenter <- S_lmo$longitude - mean(S_lmo$longitude)
S_lmo$areacenter <- log(S_lmo$area) - mean(log(S_lmo$area))

Y_lmo <- filter(birds_lmo, species == "Y", treatment == "AES")
Y_lmo$latcenter <- Y_lmo$latitude - mean (Y_lmo$latitude)
Y_lmo$visits_f <- as.factor(Y_lmo$visits)
Y_lmo$loncenter <- Y_lmo$longitude - mean(Y_lmo$longitude)
Y_lmo$areacenter <- log(Y_lmo$area) - mean(log(Y_lmo$area))

# Defining parameter-expanded priors
prior3 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000),
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000),
                        G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))

LI_lmo$uc<-as.numeric(LI_lmo$uc)
LI_lmo$mg.all <-as.numeric(LI_lmo$mg.all)
LI_lmo$ogwgw <-as.numeric(LI_lmo$ogwgw)
LI_lmo$cmsrg <-as.numeric(LI_lmo$cmsrg)
LI_lmo$mhm <-as.numeric(LI_lmo$mhm)
LI_lmo$mwet.all <-as.numeric(LI_lmo$mwet.all)
LI_lmo$wmar <-as.numeric(LI_lmo$wmar)
LI_lmo$mhedge.all <-as.numeric(LI_lmo$mhedge.all)
LI_lmo$gmarbeet <-as.numeric(LI_lmo$gmarbeet)

# Linnet
LI_lmo_m2 <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                      random = ~farm + year + 
                        idv(uc+ mg.all+ ogwgw+ cmsrg+ mhm+ mwet.all+ wmar+ mhedge.all+ gmarbeet),
                      pr = TRUE, nitt = 600000, data = LI_lmo, family = "poisson", 
                      prior = prior3, burnin = 20000)

summary(LI_lmo_m2)
#the effect of treatments varies significantly. Next we should look at the blups.
colnames(LI_lmo_m2$Sol)
plot(LI_lmo_m2$Sol[,56:64])

# Calculating mean blup and 95% HPD for each treatment
ucTRUE <- mean(LI_lmo_m2$Sol[,56])
ucTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,56])
mg.allTRUE <- mean(LI_lmo_m2$Sol[,57])
mg.allTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,57])
ogwgwTRUE <- mean(LI_lmo_m2$Sol[,58])
ogwgwTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,58])
cmsrgTRUE <- mean(LI_lmo_m2$Sol[,59])
cmsrgTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,59])
mhmTRUE <- mean(LI_lmo_m2$Sol[,60])
mhmTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,60])
mwet.allTRUE <- mean(LI_lmo_m2$Sol[,61])
mwet.allTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,61])
wmarTRUE <- mean(LI_lmo_m2$Sol[,62])
wmarTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,62])
mhedge.allTRUE <- mean(LI_lmo_m2$Sol[,63])
mhedge.allTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,63])
gmarbeetTRUE <- mean(LI_lmo_m2$Sol[,64])
gmarbeetTRUE_CI <- HPDinterval(LI_lmo_m2$Sol[,64])

# Combining the results into one data frame
LI_BLUPS <- data.frame(matrix("NA", nrow = 9, ncol = 5))
colnames(LI_BLUPS) <- c("species", "LMO", "BLUP", "upper", "lower")
LI_BLUPS$species <- as.factor("LI")
LI_BLUPS$LMO <- as.factor(c("UC", "MG", "OGWGW",
                            "CMSRG", "MHM", "MWET", "WMAR",
                            "MHEDGE", "GMARBEET"))
LI_BLUPS$BLUP <- as.numeric(c(ucTRUE, mg.allTRUE, ogwgwTRUE, cmsrgTRUE,
                              mhmTRUE, mwet.allTRUE, wmarTRUE, mhedge.allTRUE,
                              gmarbeetTRUE))
LI_BLUPS$upper <- as.numeric(c(ucTRUE_CI[,2], mg.allTRUE_CI[,2],
                               ogwgwTRUE_CI[,2], cmsrgTRUE_CI[,2], mhmTRUE_CI[,2], 
                               mwet.allTRUE_CI[,2], wmarTRUE_CI [,2], 
                               mhedge.allTRUE_CI[,2], gmarbeetTRUE_CI[,2]))
LI_BLUPS$lower <- as.numeric(c(ucTRUE_CI[,1], mg.allTRUE_CI[,1],
                               ogwgwTRUE_CI[,1], cmsrgTRUE_CI[,1], mhmTRUE_CI[,1], 
                               mwet.allTRUE_CI[,1], wmarTRUE_CI [,1], 
                               mhedge.allTRUE_CI[,1], gmarbeetTRUE_CI[,1]))

# Reed bunting
RB_lmo$uc<-as.numeric(RB_lmo$uc)
RB_lmo$mg.all <-as.numeric(RB_lmo$mg.all)
RB_lmo$ogwgw <-as.numeric(RB_lmo$ogwgw)
RB_lmo$cmsrg <-as.numeric(RB_lmo$cmsrg)
RB_lmo$mhm <-as.numeric(RB_lmo$mhm)
RB_lmo$mwet.all <-as.numeric(RB_lmo$mwet.all)
RB_lmo$wmar <-as.numeric(RB_lmo$wmar)
RB_lmo$mhedge.all <-as.numeric(RB_lmo$mhedge.all)
RB_lmo$gmarbeet <-as.numeric(RB_lmo$gmarbeet)

RB_lmo_m2 <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                      random = ~farm + year + 
                        idv(uc+ mg.all+ ogwgw+ cmsrg+ mhm+ mwet.all+ wmar+ mhedge.all+ gmarbeet),
                      pr = TRUE, nitt = 600000, data = RB_lmo, family = "poisson", 
                      prior = prior3, burnin = 20000)

summary(RB_lmo_m2)

colnames(RB_lmo_m2$Sol)
plot(RB_lmo_m2$Sol[,56:64])

# Calculating mean blup and 95% HPD for each treatment
ucTRUE <- (mean(RB_lmo_m2$Sol[,56]))
ucTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,56]))
mg.allTRUE <- (mean(RB_lmo_m2$Sol[,57]))
mg.allTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,57]))
ogwgwTRUE <- (mean(RB_lmo_m2$Sol[,58]))
ogwgwTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,58]))
cmsrgTRUE <- (mean(RB_lmo_m2$Sol[,59]))
cmsrgTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,59]))
mhmTRUE <- (mean(RB_lmo_m2$Sol[,60]))
mhmTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,60]))
mwet.allTRUE <- (mean(RB_lmo_m2$Sol[,61]))
mwet.allTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,61]))
wmarTRUE <- (mean(RB_lmo_m2$Sol[,62]))
wmarTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,62]))
mhedge.allTRUE <- (mean(RB_lmo_m2$Sol[,63]))
mhedge.allTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,63]))
gmarbeetTRUE <- (mean(RB_lmo_m2$Sol[,64]))
gmarbeetTRUE_CI <- (HPDinterval(RB_lmo_m2$Sol[,64]))

# Combining the results into one data frame
RB_BLUPS <- data.frame(matrix("NA", nrow = 9, ncol = 5))
colnames(RB_BLUPS) <- c("species", "LMO", "BLUP", "upper", "lower")
RB_BLUPS$species <- as.factor("RB")
RB_BLUPS$LMO <- as.factor(c("UC", "MG", "OGWGW",
                            "CMSRG", "MHM", "MWET", "WMAR",
                            "MHEDGE", "GMARBEET"))
RB_BLUPS$BLUP <- as.numeric(c(ucTRUE, mg.allTRUE, ogwgwTRUE, cmsrgTRUE,
                              mhmTRUE, mwet.allTRUE, wmarTRUE, mhedge.allTRUE,
                              gmarbeetTRUE))
RB_BLUPS$upper <- as.numeric(c(ucTRUE_CI[,2], mg.allTRUE_CI[,2],
                               ogwgwTRUE_CI[,2], cmsrgTRUE_CI[,2], mhmTRUE_CI[,2], 
                               mwet.allTRUE_CI[,2], wmarTRUE_CI [,2], 
                               mhedge.allTRUE_CI[,2], gmarbeetTRUE_CI[,2]))
RB_BLUPS$lower <- as.numeric(c(ucTRUE_CI[,1], mg.allTRUE_CI[,1],
                               ogwgwTRUE_CI[,1], cmsrgTRUE_CI[,1], mhmTRUE_CI[,1], 
                               mwet.allTRUE_CI[,1], wmarTRUE_CI [,1], 
                               mhedge.allTRUE_CI[,1], gmarbeetTRUE_CI[,1]))

# Tree sparrow
TS_lmo$uc<-as.numeric(TS_lmo$uc)
TS_lmo$mg.all <-as.numeric(TS_lmo$mg.all)
TS_lmo$ogwgw <-as.numeric(TS_lmo$ogwgw)
TS_lmo$cmsrg <-as.numeric(TS_lmo$cmsrg)
TS_lmo$mhm <-as.numeric(TS_lmo$mhm)
TS_lmo$mwet.all <-as.numeric(TS_lmo$mwet.all)
TS_lmo$wmar <-as.numeric(TS_lmo$wmar)
TS_lmo$mhedge.all <-as.numeric(TS_lmo$mhedge.all)
TS_lmo$gmarbeet <-as.numeric(TS_lmo$gmarbeet)

TS_lmo_m2 <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                      random = ~farm + year + 
                        idv(uc+ mg.all+ ogwgw+ cmsrg+ mhm+ mwet.all+ wmar+ mhedge.all+ gmarbeet),
                      pr = TRUE, nitt = 600000, data = TS_lmo, family = "poisson", 
                      prior = prior3, burnin = 20000)

summary(TS_lmo_m2)

colnames(TS_lmo_m2$Sol)
plot(TS_lmo_m2$Sol[,56:64])

# Calculating mean blup and 95% HPD for each treatment
ucTRUE <- (mean(TS_lmo_m2$Sol[,56]))
ucTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,56]))
mg.allTRUE <- (mean(TS_lmo_m2$Sol[,57]))
mg.allTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,57]))
ogwgwTRUE <- (mean(TS_lmo_m2$Sol[,58]))
ogwgwTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,58]))
cmsrgTRUE <- (mean(TS_lmo_m2$Sol[,59]))
cmsrgTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,59]))
mhmTRUE <- (mean(TS_lmo_m2$Sol[,60]))
mhmTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,60]))
mwet.allTRUE <- (mean(TS_lmo_m2$Sol[,61]))
mwet.allTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,61]))
wmarTRUE <- (mean(TS_lmo_m2$Sol[,62]))
wmarTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,62]))
mhedge.allTRUE <- (mean(TS_lmo_m2$Sol[,63]))
mhedge.allTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,63]))
gmarbeetTRUE <- (mean(TS_lmo_m2$Sol[,64]))
gmarbeetTRUE_CI <- (HPDinterval(TS_lmo_m2$Sol[,64]))

# Combining the results into one data frame
TS_BLUPS <- data.frame(matrix("NA", nrow = 9, ncol = 5))
colnames(TS_BLUPS) <- c("species", "LMO", "BLUP", "upper", "lower")
TS_BLUPS$species <- as.factor("TS")
TS_BLUPS$LMO <- as.factor(c("UC", "MG", "OGWGW",
                            "CMSRG", "MHM", "MWET", "WMAR",
                            "MHEDGE", "GMARBEET"))
TS_BLUPS$BLUP <- as.numeric(c(ucTRUE, mg.allTRUE, ogwgwTRUE, cmsrgTRUE,
                              mhmTRUE, mwet.allTRUE, wmarTRUE, mhedge.allTRUE,
                              gmarbeetTRUE))
TS_BLUPS$upper <- as.numeric(c(ucTRUE_CI[,2], mg.allTRUE_CI[,2],
                               ogwgwTRUE_CI[,2], cmsrgTRUE_CI[,2], mhmTRUE_CI[,2], 
                               mwet.allTRUE_CI[,2], wmarTRUE_CI [,2], 
                               mhedge.allTRUE_CI[,2], gmarbeetTRUE_CI[,2]))
TS_BLUPS$lower <- as.numeric(c(ucTRUE_CI[,1], mg.allTRUE_CI[,1],
                               ogwgwTRUE_CI[,1], cmsrgTRUE_CI[,1], mhmTRUE_CI[,1], 
                               mwet.allTRUE_CI[,1], wmarTRUE_CI [,1], 
                               mhedge.allTRUE_CI[,1], gmarbeetTRUE_CI[,1]))

# Skylark
S_lmo$uc<-as.numeric(S_lmo$uc)
S_lmo$mg.all <-as.numeric(S_lmo$mg.all)
S_lmo$ogwgw <-as.numeric(S_lmo$ogwgw)
S_lmo$cmsrg <-as.numeric(S_lmo$cmsrg)
S_lmo$mhm <-as.numeric(S_lmo$mhm)
S_lmo$mwet.all <-as.numeric(S_lmo$mwet.all)
S_lmo$wmar <-as.numeric(S_lmo$wmar)
S_lmo$mhedge.all <-as.numeric(S_lmo$mhedge.all)
S_lmo$gmarbeet <-as.numeric(S_lmo$gmarbeet)

S_lmo_m2 <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                     random = ~farm + year + 
                       idv(uc+ mg.all+ ogwgw+ cmsrg+ mhm+ mwet.all+ wmar+ mhedge.all+ gmarbeet),
                     pr = TRUE, nitt = 600000, data = S_lmo, family = "poisson", 
                     prior = prior3, burnin = 20000)

summary(S_lmo_m2)

colnames(S_lmo_m2$Sol)
plot(S_lmo_m2$Sol[,56:64])

# Calculating mean blup and 95% HPD for each treatment
ucTRUE <- (mean(S_lmo_m2$Sol[,56]))
ucTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,56]))
mg.allTRUE <- (mean(S_lmo_m2$Sol[,57]))
mg.allTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,57]))
ogwgwTRUE <- (mean(S_lmo_m2$Sol[,58]))
ogwgwTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,58]))
cmsrgTRUE <- (mean(S_lmo_m2$Sol[,59]))
cmsrgTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,59]))
mhmTRUE <- (mean(S_lmo_m2$Sol[,60]))
mhmTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,60]))
mwet.allTRUE <- (mean(S_lmo_m2$Sol[,61]))
mwet.allTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,61]))
wmarTRUE <- (mean(S_lmo_m2$Sol[,62]))
wmarTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,62]))
mhedge.allTRUE <- (mean(S_lmo_m2$Sol[,63]))
mhedge.allTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,63]))
gmarbeetTRUE <- (mean(S_lmo_m2$Sol[,64]))
gmarbeetTRUE_CI <- (HPDinterval(S_lmo_m2$Sol[,64]))

# Combining the results into one data frame
S_BLUPS <- data.frame(matrix("NA", nrow = 9, ncol = 5))
colnames(S_BLUPS) <- c("species", "LMO", "BLUP", "upper", "lower")
S_BLUPS$species <- as.factor("S")
S_BLUPS$LMO <- as.factor(c("UC", "MG", "OGWGW",
                           "CMSRG", "MHM", "MWET", "WMAR",
                           "MHEDGE", "GMARBEET"))
S_BLUPS$BLUP <- as.numeric(c(ucTRUE, mg.allTRUE, ogwgwTRUE, cmsrgTRUE,
                             mhmTRUE, mwet.allTRUE, wmarTRUE, mhedge.allTRUE,
                             gmarbeetTRUE))
S_BLUPS$upper <- as.numeric(c(ucTRUE_CI[,2], mg.allTRUE_CI[,2],
                              ogwgwTRUE_CI[,2], cmsrgTRUE_CI[,2], mhmTRUE_CI[,2], 
                              mwet.allTRUE_CI[,2], wmarTRUE_CI [,2], 
                              mhedge.allTRUE_CI[,2], gmarbeetTRUE_CI[,2]))
S_BLUPS$lower <- as.numeric(c(ucTRUE_CI[,1], mg.allTRUE_CI[,1],
                              ogwgwTRUE_CI[,1], cmsrgTRUE_CI[,1], mhmTRUE_CI[,1], 
                              mwet.allTRUE_CI[,1], wmarTRUE_CI [,1], 
                              mhedge.allTRUE_CI[,1], gmarbeetTRUE_CI[,1]))

# Yellowhammer
Y_lmo$uc<-as.numeric(Y_lmo$uc)
Y_lmo$mg.all <-as.numeric(Y_lmo$mg.all)
Y_lmo$ogwgw <-as.numeric(Y_lmo$ogwgw)
Y_lmo$cmsrg <-as.numeric(Y_lmo$cmsrg)
Y_lmo$mhm <-as.numeric(Y_lmo$mhm)
Y_lmo$mwet.all <-as.numeric(Y_lmo$mwet.all)
Y_lmo$wmar <-as.numeric(Y_lmo$wmar)
Y_lmo$mhedge.all <-as.numeric(Y_lmo$mhedge.all)
Y_lmo$gmarbeet <-as.numeric(Y_lmo$gmarbeet)

Y_lmo_m2 <- MCMCglmm(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f, 
                     random = ~farm + year + 
                       idv(uc+ mg.all+ ogwgw+ cmsrg+ mhm+ mwet.all+ wmar+ mhedge.all+ gmarbeet),
                     pr = TRUE, nitt = 600000, data = Y_lmo, family = "poisson", 
                     prior = prior3, burnin = 20000)

summary(Y_lmo_m2)

colnames(Y_lmo_m2$Sol)
plot(Y_lmo_m2$Sol[,56:64])

# Calculating mean blup and 95% HPD for each treatment
ucTRUE <- (mean(Y_lmo_m2$Sol[,56]))
ucTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,56]))
mg.allTRUE <- (mean(Y_lmo_m2$Sol[,57]))
mg.allTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,57]))
ogwgwTRUE <- (mean(Y_lmo_m2$Sol[,58]))
ogwgwTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,58]))
cmsrgTRUE <- (mean(Y_lmo_m2$Sol[,59]))
cmsrgTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,59]))
mhmTRUE <- (mean(Y_lmo_m2$Sol[,60]))
mhmTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,60]))
mwet.allTRUE <- (mean(Y_lmo_m2$Sol[,61]))
mwet.allTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,61]))
wmarTRUE <- (mean(Y_lmo_m2$Sol[,62]))
wmarTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,62]))
mhedge.allTRUE <- (mean(Y_lmo_m2$Sol[,63]))
mhedge.allTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,63]))
gmarbeetTRUE <- (mean(Y_lmo_m2$Sol[,64]))
gmarbeetTRUE_CI <- (HPDinterval(Y_lmo_m2$Sol[,64]))

# Combining the results into one data frame
Y_BLUPS <- data.frame(matrix("NA", nrow = 9, ncol = 5))
colnames(Y_BLUPS) <- c("species", "LMO", "BLUP", "upper", "lower")
Y_BLUPS$species <- as.factor("Y")
Y_BLUPS$LMO <- as.factor(c("UC", "MG", "OGWGW",
                           "CMSRG", "MHM", "MWET", "WMAR",
                           "MHEDGE", "GMARBEET"))
Y_BLUPS$BLUP <- as.numeric(c(ucTRUE, mg.allTRUE, ogwgwTRUE, cmsrgTRUE,
                             mhmTRUE, mwet.allTRUE, wmarTRUE, mhedge.allTRUE,
                             gmarbeetTRUE))
Y_BLUPS$upper <- as.numeric(c(ucTRUE_CI[,2], mg.allTRUE_CI[,2],
                              ogwgwTRUE_CI[,2], cmsrgTRUE_CI[,2], mhmTRUE_CI[,2], 
                              mwet.allTRUE_CI[,2], wmarTRUE_CI [,2], 
                              mhedge.allTRUE_CI[,2], gmarbeetTRUE_CI[,2]))
Y_BLUPS$lower <- as.numeric(c(ucTRUE_CI[,1], mg.allTRUE_CI[,1],
                              ogwgwTRUE_CI[,1], cmsrgTRUE_CI[,1], mhmTRUE_CI[,1], 
                              mwet.allTRUE_CI[,1], wmarTRUE_CI [,1], 
                              mhedge.allTRUE_CI[,1], gmarbeetTRUE_CI[,1]))

LI_lmo_m2_sum <- clean.MCMC(LI_lmo_m2)
stargazer(LI_lmo_m2_sum, type = "html", summary = FALSE, digits = 2)
RB_lmo_m2_sum <- clean.MCMC(RB_lmo_m2)
stargazer(RB_lmo_m2_sum, type = "html", summary = FALSE, digits = 2)
TS_lmo_m2_sum <- clean.MCMC(TS_lmo_m2)
stargazer(TS_lmo_m2_sum, type = "html", summary = FALSE, digits = 2)
S_lmo_m2_sum <- clean.MCMC(S_lmo_m2)
stargazer(S_lmo_m2_sum, type = "html", summary = FALSE, digits = 2)
Y_lmo_m2_sum <- clean.MCMC(Y_lmo_m2)
stargazer(Y_lmo_m2_sum, type = "html", summary = FALSE, digits = 2)

# Figures ----

# ** Predictions for bird abundance on AES and control and effect sizes ----
# Importing the predicted values data
prediction <- read.csv2("aes_predictions.csv")

# Plotting predicted values for each species across treatments
(pred_plot <- ggplot(prediction, aes(x = species, y = mean, 
                                     colour = treatment, group = treatment)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, 
                  colour = treatment), width = .3, size = 1.3, position = position_dodge(.3)) + 
    geom_point(position = position_dodge(.3), size = 4, alpha = 0.8) +
    geom_point(shape = 1, size = 4,colour = "black", position = position_dodge(.3)) +
    farm_theme() +
    labs(x = "", y="Abundance") +
    scale_colour_manual(values = c("#EE7600", "#00868B")) + 
    theme(legend.text = element_text(size = 14, face = "bold"),    
          legend.title = element_blank(),                        
          legend.position = c(0.9, 0.9),
          legend.background = element_rect(fill = "transparent", size = .5),
          legend.key = element_rect(fill = "transparent", colour = "transparent")))

# Plotting effect sizes and CIs with a line for zero
eff_size4 <- read.csv("ordered_AESeffects.csv")

(eff_plot4 <- ggplot(eff_size4, aes(x = species, y = post.mean)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "#EE7600", width = .3, size = 1.3, 
                  position = position_dodge(.3)) + 
    geom_point(position = position_dodge(.3), size = 4, alpha = 0.8, colour = "#EE7600") +
    geom_point(shape = 1, size = 4,colour = "black", position = position_dodge(.3)) +
    farm_theme() +
    labs(x = "", y = "Coefficient") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(legend.text = element_text(size = 12, face = "bold"),    
          legend.title = element_blank(),                        
          legend.position = c(0.9, 0.9),
          legend.background = element_rect(fill = "transparent", size = .5),
          legend.key = element_rect(fill = "transparent", colour = "transparent")))


# ** AES years effects ----
AESyears_eff <- read.csv("AESyears_effects.csv")

(AES_years_plot <- ggplot(AESyears_eff, aes(x = species, y = slope)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), colour = "#EE7600", width = .3, size = 1.3, 
                  position = position_dodge(.3)) + 
    geom_point(position = position_dodge(.3), size = 4, alpha = 0.8, colour = "#EE7600") +
    geom_point(shape = 1, size = 4,colour = "black", position = position_dodge(.3)) +
    farm_theme() +
    labs(x = "", y = "Coefficient") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(legend.text = element_text(size = 12, face = "bold"),    
          legend.title = element_blank(),                        
          legend.position = c(0.9, 0.9),
          legend.background = element_rect(fill = "transparent", size = .5),
          legend.key = element_rect(fill = "transparent", colour = "transparent")))

panel <- grid.arrange(pred_plot, eff_plot4, AES_years_plot, ncol = 1)
ggsave(panel, file="eff_size_panel3.png", width = 8, height = 12)

# ** Plotting raw abundance data across years ----
Ssum <- summarySE(S, measurevar = "abundance", groupvars = c("treatment", "year"))
TSsum <- summarySE(TS, measurevar = "abundance", groupvars = c("treatment", "year"))
LIsum <- summarySE(LI, measurevar = "abundance", groupvars = c("treatment", "year"))
RBsum <- summarySE(RB, measurevar = "abundance", groupvars = c("treatment", "year"))
Ysum <- summarySE(Y, measurevar = "abundance", groupvars = c("treatment", "year"))

(LIraw <- ggplot(LIsum, aes(x = year, y = abundance, 
                            colour = treatment))+
    geom_errorbar(aes(ymin = abundance - se, ymax = abundance + se, colour = treatment), width = .6, size = 1.3) + 
    geom_line(size = 1.5) + geom_point(size = 4, alpha = 0.8) +
    geom_point(shape = 1, size = 4,colour = "black") +
    farm_theme() +
    labs(y = "Abundance", title = "Linnet") +
    scale_x_continuous(breaks = seq(2003, 2015, 3)) +
    scale_y_continuous(limits = c(3, 11), breaks = c(4, 6, 8, 10)) +
    scale_colour_manual(values = c("#EE7600", "#00868B")) + 
    theme(axis.title.x = element_blank(),
          legend.position = c(0.8, 0.84)))

(RBraw <- ggplot(RBsum, aes(x = year, y = abundance, 
                            colour = treatment))+
    geom_errorbar(aes(ymin = abundance - se, ymax = abundance + se, colour = treatment), width = .6, size = 1.3) + 
    geom_line(size = 1.5) + geom_point(size = 4, alpha = 0.8) +
    geom_point(shape = 1, size = 4,colour = "black") +
    farm_theme() +
    labs(y = "Abundance", title = "Reed bunting") +
    scale_x_continuous(breaks = seq(2003, 2015, 3)) +
    scale_y_continuous(limits = c(0, 7), breaks = c(2, 4, 6)) +
    scale_colour_manual(values = c("#EE7600", "#00868B")) + 
    theme(axis.title.x = element_blank(), legend.position = "none"))

(TSraw <- ggplot(TSsum, aes(x = year, y = abundance, 
                            colour = treatment))+
    geom_errorbar(aes(ymin = abundance - se, ymax = abundance + se, colour = treatment), width = .6, size = 1.3) + 
    geom_line(size = 1.5) + geom_point(size = 4, alpha = 0.8) +
    geom_point(shape = 1, size = 4,colour = "black") +
    farm_theme() +
    labs(y = "Abundance", title = "Tree sparrow") +
    scale_x_continuous(breaks = seq(2003, 2015, 3)) +
    scale_y_continuous(limits = c(0, 10), breaks = c(2, 4, 6, 8)) +
    scale_colour_manual(values = c("#EE7600", "#00868B")) + 
    theme(axis.title.x = element_blank(), legend.position = "none"))

(Sraw <- ggplot(Ssum, aes(x = year, y = abundance, 
                          colour = treatment)) +
    geom_errorbar(aes(ymin = abundance - se, ymax = abundance + se, colour = treatment), width = .6, size = 1.3) + 
    geom_line(size = 1.5) + geom_point(size = 4, alpha = 0.8) +
    geom_point(shape = 1, size = 4,colour = "black") +
    farm_theme() +
    labs(y = "Abundance", title = "Skylark") +
    scale_x_continuous(breaks = seq(2003, 2015, 3)) +
    scale_y_continuous(limits = c(10, 29), breaks = c(15, 20, 25)) +
    scale_colour_manual(values = c("#EE7600", "#00868B")) + 
    theme(axis.title.x = element_blank(), legend.position = "none"))

(Yraw <- ggplot(Ysum, aes(x = year, y = abundance, 
                          colour = treatment))+
    geom_errorbar(aes(ymin = abundance - se, ymax = abundance + se, colour = treatment), width = .6, size = 1.3) + 
    geom_line(size = 1.5) + geom_point(size = 4, alpha = 0.8) +
    geom_point(shape = 1, size = 4,colour = "black") +
    farm_theme() +
    labs(y = "Abundance", title = "Yellowhammer") +
    scale_x_continuous(breaks = seq(2003, 2015, 3)) +
    scale_y_continuous(limits = c(5, 23), breaks = c(5, 10, 15, 20)) +
    scale_colour_manual(values = c("#EE7600", "#00868B")) + 
    theme(axis.title.x = element_blank(), legend.position = "none"))

yearly_panel_raw <- grid.arrange(LIraw, Sraw, RBraw, Yraw, TSraw, ncol = 2)
ggsave(yearly_panel_raw, file="yearly_panel4.png", width = 10, height = 10)

# ** Spatial autocorrelation ----
# Based on rho values from spaMM models
# Linnet
rho_LI <- 0.0002109955
whendoesitreach0.1_LI <- seq(1, 50000, 1)[which(abs(Exponential(seq(1, 50000, 1),
                                                                range = 1/rho_LI) - 0.1) == min(abs(Exponential(seq(1, 50000, 1), 
                                                                                                                range = 1/rho_LI) - 0.1)))]

# Reed bunting 
rho_RB <- 0.0003079227
whendoesitreach0.1_RB <- seq(1, 50000, 1)[which(abs(Exponential(seq(1, 50000, 1),
                                                                range = 1/rho_RB) - 0.1) == min(abs(Exponential(seq(1, 50000, 1), 
                                                                                                                range = 1/rho_RB) - 0.1)))]

# Tree sparrow 
rho_TS <- 0.0005441002
whendoesitreach0.1_TS <- seq(1, 50000, 1)[which(abs(Exponential(seq(1, 50000, 1),
                                                                range = 1/rho_TS) - 0.1) == min(abs(Exponential(seq(1, 50000, 1),
                                                                                                                range = 1/rho_TS) - 0.1)))]

# Skylark 
rho_S <- 0.00123716
whendoesitreach0.1_S <- seq(1, 50000, 1)[which(abs(Exponential(seq(1, 50000, 1),
                                                               range = 1/rho_S) - 0.1) == min(abs(Exponential(seq(1, 50000, 1), 
                                                                                                              range = 1/rho_S) - 0.1)))]

# Yellowhammer 
rho_Y <- 0.0003745226 
whendoesitreach0.1_Y <- seq(1, 50000, 1)[which(abs(Exponential(seq(1, 50000, 1),
                                                               range = 1/rho_Y) - 0.1) == min(abs(Exponential(seq(1, 50000, 1), 
                                                                                                              range = 1/rho_Y) - 0.1)))]

cor_LI <- ggplot() +
  geom_point(aes(x = seq(1,50000,1), y = Exponential(seq(1,50000,1), range = 1/rho_LI)), 
             size = 0.8) +
  farm_theme() +
  labs(y = "Correlation", x = "Distance (m)", title = "Linnet") +
  geom_vline(xintercept = whendoesitreach0.1_LI)

cor_RB <- ggplot() +
  geom_point(aes(x = seq(1,50000,1), y = Exponential(seq(1,50000,1), range = 1/rho_RB)), 
             size = 0.8) +
  farm_theme() +
  labs(y = "Correlation", x = "Distance (m)", title = "Reed bunting") +
  geom_vline(xintercept = whendoesitreach0.1_RB)

cor_TS <- ggplot() +
  geom_point(aes(x = seq(1,50000,1), y = Exponential(seq(1,50000,1), range = 1/rho_TS)), 
             size = 0.8) +
  farm_theme() +
  labs(y = "Correlation", x = "Distance (m)", title = "Tree sparrow") +
  geom_vline(xintercept = whendoesitreach0.1_TS)

cor_S <- ggplot() +
  geom_point(aes(x = seq(1,50000,1), y = Exponential(seq(1,50000,1), range = 1/rho_S)), 
             size = 0.8) +
  farm_theme() +
  labs(y = "Correlation", x = "Distance (m)", title = "Skylark") +
  geom_vline(xintercept = whendoesitreach0.1_S)

cor_Y <- ggplot() +
  geom_point(aes(x = seq(1,50000,1), y = Exponential(seq(1,50000,1), range = 1/rho_Y)), 
             size = 0.8) +
  farm_theme() +
  labs(y = "Correlation", x = "Distance (m)", title = "Yellowhammer") +
  geom_vline(xintercept = whendoesitreach0.1_Y)

panel_cor <- grid.arrange(cor_LI, cor_RB, cor_TS, cor_S, cor_Y, ncol = 1)
ggsave(panel_cor, file = "cor_panel.png", width = 5, height = 15)

# ** LMO BLUPs ----

(LI_BLUP_p <- ggplot(LI_BLUPS, aes(x = LMO, y = BLUP)) +
   geom_errorbar(aes(ymin = lower, ymax = upper), 
                 colour = "#EE7600", width = .3, size = 1.3) + 
   geom_point(size = 4, alpha = 0.8, colour = "#EE7600") +
   geom_point(shape = 1, size = 4,colour = "black") +
   geom_hline(yintercept = 0, linetype = "dashed") +
   farm_theme2() +
   labs(y = "BLUP", title = "Linnet", x = ""))

(RB_BLUP_p <- ggplot(RB_BLUPS, aes(x = LMO, y = BLUP)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  colour = "#EE7600", width = .3, size = 1.3) + 
    geom_point(size = 4, alpha = 0.8, colour = "#EE7600") +
    geom_point(shape = 1, size = 4,colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    farm_theme2() +
    labs(y = "BLUP", title = "Reed bunting", x = ""))

(TS_BLUP_p <- ggplot(TS_BLUPS, aes(x = LMO, y = BLUP)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  colour = "#EE7600", width = .3, size = 1.3) + 
    geom_point(size = 4, alpha = 0.8, colour = "#EE7600") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(shape = 1, size = 4,colour = "black") +
    farm_theme2() +
    labs(y = "BLUP", title = "Tree sparrow", x = ""))

(S_BLUP_p <- ggplot(S_BLUPS, aes(x = LMO, y = BLUP)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  colour = "#EE7600", width = .3, size = 1.3) + 
    geom_point(size = 4, alpha = 0.8, colour = "#EE7600") +
    geom_point(shape = 1, size = 4,colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    farm_theme2() +
    labs(y = "BLUP", title = "Skylark", x = ""))

(Y_BLUP_p <- ggplot(Y_BLUPS, aes(x = LMO, y = BLUP)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                  colour = "#EE7600", width = .3, size = 1.3) + 
    geom_point(size = 4, alpha = 0.8, colour = "#EE7600") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(shape = 1, size = 4,colour = "black") +
    farm_theme2() +
    labs(y = "BLUP", title = "Yellowhammer", x = ""))

# Making a panel with all 5 graphs
lmo_panel3 <- grid.arrange(LI_BLUP_p, RB_BLUP_p, TS_BLUP_p, S_BLUP_p, Y_BLUP_p, ncol = 2)
ggsave(lmo_panel3, file = "lmoBLUPs2.png", width = 10, height = 12)

# Farmland power analysis ----
# Linnet
LI <- filter(birds_aes, species == "LI")
LI$areacenter <- log(LI$area) - mean(log(LI$area))

# Reed bunting
RB <- filter(birds_aes, species == "RB")
RB$areacenter <- log(RB$area) - mean(log(RB$area))

# Tree sparrow
TS <- filter(birds_aes, species == "TS")
TS$areacenter <- log(TS$area) - mean(log(TS$area))

# Skylark
S <- filter(birds_aes, species == "S")
S$areacenter <- log(S$area) - mean(log(S$area))

# Yellowhammer
Y <- filter(birds_aes, species == "Y")
Y$areacenter <- log(Y$area) - mean(log(Y$area))

# Power analysis for AES being 25% and 10% better than control
# I.e. an effect size of log(1.25) = 0.2231436 and log(1.10) = 0.09531018

# Linnet
LI_pow <- glmer(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f + (1|farm) + (1|year), family = "poisson", data = LI)
summary(LI_pow)

fixef(LI_pow)["treatmentAES"]  # 0.007739446 
powerSim(LI_pow)
fixef(LI_pow)["treatmentAES"] <- 0.2231436
powerSim(LI_pow)

fixef(LI_pow)["treatmentAES"] <- 0.09531018
powerSim(LI_pow)

# Reed bunting - doesn't converge with  + visits_f  + latcenter + loncenter + as a fixed effect, simplified the model
RB_pow <- glmer(abundance ~ treatment + areacenter + (1|farm) + (1|year), family = "poisson", data = RB)
summary(RB_pow)

fixef(RB_pow)["treatmentAES"]  # 0.10 original
powerSim(RB_pow)

fixef(RB_pow)["treatmentAES"] <- 0.2231436
powerSim(RB_pow)

fixef(RB_pow)["treatmentAES"] <- 0.09531018
powerSim(RB_pow)

# Tree sparrow
TS_pow <- glmer(abundance ~ treatment + areacenter  + latcenter + loncenter + visits_f + (1|farm) + (1|year), family = "poisson", data = TS)
summary(TS_pow)

fixef(TS_pow)["treatmentAES"]  #0.19
powerSim(TS_pow)

fixef(TS_pow)["treatmentAES"] <- 0.2231436
powerSim(TS_pow)

fixef(TS_pow)["treatmentAES"] <- 0.09531018
powerSim(TS_pow)

# Skylark
S_pow <- glmer(abundance ~ treatment + areacenter + latcenter + loncenter + visits_f + (1|farm) + (1|year), family = "poisson", data = S)
summary(S_pow)

fixef(S_pow)["treatmentAES"]  # -0.12
powerSim(S_pow)


fixef(S_pow)["treatmentAES"] <- 0.2231436
powerSim(S_pow)

fixef(S_pow)["treatmentAES"] <- 0.09531018
powerSim(S_pow)

# Yellowhammer
Y_pow <- glmer(abundance ~ treatment + areacenter  + latcenter + loncenter + visits_f + (1|farm) + (1|year), family = "poisson", data = Y)
summary(Y_pow)

fixef(Y_pow)["treatmentAES"]  # 0.11
powerSim(Y_pow)

fixef(Y_pow)["treatmentAES"] <- 0.2231436
powerSim(Y_pow)

fixef(Y_pow)["treatmentAES"] <- 0.09531018
powerSim(Y_pow) 