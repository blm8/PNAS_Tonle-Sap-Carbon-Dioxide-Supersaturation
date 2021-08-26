

##===============================  Statistical Analyses


##  Author(s):  B. Miller

##  Required dataframes:  Master_TSL-Corrected.csv

##  Description:  Checks for normality and homogeneity of variance;
##                Runs ANOVA tests for significance; 
##                Runs Tukey HSD tests for significance;  
##                Calculates Cohen's d

#install.packages("readr")
library(readr) 
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("mgcv")
library(mgcv)
#install.packages("AICcmodavg")
library(AICcmodavg)
#install.packages("rLakeAnalyzer")
library(rLakeAnalyzer) # Produced by GLEON
#install.packages("lubridate")
library(lubridate) # Includes functions for time
#install.packages("LakeMetabolizer")
library(LakeMetabolizer) # Produced by GLEON
#install.packages("plotrix")
library(plotrix)
#install.packages("rcompanion")
library(rcompanion)

setwd("~/Desktop/Working")

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

High <- subset(newData, STAGE=="High")
mean(na.omit(High$WATER_T))
Falling <- subset(newData, STAGE=="Falling")
mean(na.omit(Falling$WATER_T))

cGas <- newData %>%
  dplyr::select(SITE, CLASS_Z, DATE, PCO2, PCH4, D13C_CO2, D13C_CH4, ENVIRON, STAGE, CO2_molL, CH4_molL) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCO2 = mean(PCO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.PCH4 = mean(PCH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.CO2_molL = mean(CO2_molL, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.CH4_molL = mean(CH4_molL, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CO2 = mean(D13C_CO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CH4 = mean(D13C_CH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(eC = avg.D13C_CO2 - avg.D13C_CH4) %>% # Apply the function mean() to every column
  distinct(avg.PCO2, .keep_all=TRUE) %>%
  ungroup() 
cGas <- cGas[!cGas$SITE == "KgPr", ]
cGas <- cGas[!cGas$SITE == "StSe", ]
cGas <- cGas[!cGas$SITE == "PtSd.O", ] 
cGas <- cGas[!cGas$SITE == "PtSd.E", ] 
cGas <- cGas[!cGas$SITE == "PtSd.F", ] 
cGas <- cGas[!cGas$SITE == "KgKl.O", ] 
cGas <- cGas[!cGas$SITE == "KgKl.E", ] 
cGas <- cGas[!cGas$STAGE == "Rising", ] 
cGas <- cGas[complete.cases(cGas), ]
#View(cGas)

High <- subset(cGas, STAGE=="High")
Falling <- subset(cGas, STAGE=="Falling")

Open <- subset(Falling, ENVIRON=="Pelagic")
Edge <- subset(Falling, ENVIRON=="Edge")
Floodplain <- subset(Falling, ENVIRON=="Floodplain")
Table(Edge$avg.PCH4)


##===============================  Functions for Descriptive Statistics and Effect Sizes  

Table <- function(x) { # Calculate means and standard deviations
  n <- length(na.omit(x))
  Mean <- mean(na.omit(x))
  SD <- sd(na.omit(x))
  Results <- list(n, Mean, SD)
  return(Results)
}

d <- function(x,y) { # Calculate Cohen's d for effect size
  xM <- mean(na.omit(x))
  yM <- mean(na.omit(y))
  sdx <- sd(na.omit(x))
  sdy <- sd(na.omit(y))
  SD <- sqrt(((sdx^2)+(sdy^2))/2)
  d <- (xM-yM)/SD
  return(d)
}


##=============================== Normality Tests & Statistical Comparisons for PCO2 

# Determine whether data follow normal or non-normal distributions in order to choose appropriate statistical test

# By flood stage

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(High$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(High$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(High$avg.PCO2) # p-value = 0.05783 # Distribution of data is not significantly different from normal distribution 

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Falling$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Falling$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Falling$avg.PCO2) # p-value = 0.6506 # Distribution of data is not significantly different from normal distribution 
qqnorm(log10(Falling$avg.PCO2), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(Falling$avg.PCO2), col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(log10(Falling$avg.PCO2)) # p-value = 0.6506 # Distribution of data is not significantly different from normal distribution 

df <- group_by(cGas, STAGE) %>%
  dplyr::select(avg.PCO2, STAGE)
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ STAGE, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) # p-value = 5e-07 # High avg.PCO2 is significantly different than Falling avg.PCO2
boxplot(High$avg.PCO2, Falling$avg.PCO2)

Table(High$avg.PCO2) # 13275.8
Table(High$PCO2) # 13299.92
range(High$avg.PCO2) # 3163.21 22502.92
Table(Falling$avg.PCO2) # 25490.84
Table(Falling$PCO2) # 26297.79
range(Falling$avg.PCO2) # 15336.24 41891.69
d(High$avg.PCO2, Falling$avg.PCO2) # d= 1.835836

View(cGas)
# By depth

Surface <- subset(cGas, CLASS_Z=="Surface")
Bottom <- subset(cGas, CLASS_Z=="Bottom")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Surface$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Surface$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Surface$avg.PCO2) # p-value = 0.3565 # Distribution of data is not significantly different from normal distribution 

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Bottom$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Bottom$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Bottom$avg.PCO2) # p-value = 0.02027 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(Bottom$avg.PCO2)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Bottom$avg.PCO2)), col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(log10(abs(Bottom$avg.PCO2))) # p-value = 0.2204 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(cGas, CLASS_Z) %>%
  dplyr::select(avg.PCO2, CLASS_Z) %>%
  mutate(avg.PCO2_log = log(abs(avg.PCO2)))
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2_log ~ CLASS_Z, data = df)
# Summary of the analysis
summary(res.aov)
TukeyHSD(res.aov) 

# By lake environment
  
High <- subset(cGas, STAGE=="High")
Falling <- subset(cGas, STAGE=="Falling")

  # High-water stage

Open <- subset(High, ENVIRON=="Pelagic")
Edge <- subset(High, ENVIRON=="Edge")
Floodplain <- subset(High, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.PCO2) # p-value = 0.156 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Edge$avg.PCO2) # p-value = 0.4802 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Floodplain$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Floodplain$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Floodplain$avg.PCO2) # p-value = 0.07363 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(High, ENVIRON) %>%
  dplyr::select(avg.PCO2, ENVIRON)
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

Table(Open$avg.PCO2)
Table(Edge$avg.PCO2)
Table(Floodplain$avg.PCO2)

  # Falling-water stage

Open <- subset(Falling, ENVIRON=="Pelagic")
Edge <- subset(Falling, ENVIRON=="Edge")
Floodplain <- subset(Falling, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.PCO2) # p-value = 0.5626 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.PCO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.PCO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Edge$avg.PCO2) # p-value = 0.889 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(Falling, ENVIRON) %>%
  dplyr::select(avg.PCO2, ENVIRON)
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCO2 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

Table(Open$avg.PCO2)
Table(Edge$avg.PCO2)
Table(Floodplain$avg.PCO2)


##=============================== Normality Tests and Statistical Comparisons for PCH4 

# Determine whether data follow normal or non-normal distributions in order to choose appropriate statistical test

# By flood stage

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(High$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(High$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(High$avg.PCH4) # p-value = 2.264e-05 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(High$avg.PCH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(High$avg.PCH4)), col = "gray40", lwd=3, lty=1) # Data are skewed 
shapiro.test(log10(abs(High$avg.PCH4))) # p-value = 0.004785 # Distribution of data is significantly different from normal distribution # Use Kruskal Wallis

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Falling$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Falling$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Falling$avg.PCH4) # p-value = 8.951e-05 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(Falling$avg.PCH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Falling$avg.PCH4)), col = "gray40", lwd=3, lty=1) # Data are skewed 
shapiro.test(log10(abs(Falling$avg.PCH4))) # p-value = 0.1618 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(cGas, STAGE) %>%
  dplyr::select(avg.PCH4, STAGE) %>%
  mutate(avg.PCH4_log = log10(abs(avg.PCH4)))
View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCH4_log ~ STAGE, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) # p-value = 0.0011549 # High avg.PCH4_log is significantly different than Falling avg.PCH4_log

Table(High$avg.PCH4)
Table(Falling$avg.PCH4)
d(log10(abs(High$avg.PCH4)), log10(abs(Falling$avg.PCH4))) # d = 1.144518

# By depth

Surface <- subset(cGas, CLASS_Z=="Surface")
Bottom <- subset(cGas, CLASS_Z=="Bottom")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Surface$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Surface$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Transform
shapiro.test(Surface$avg.PCH4) # p-value = 5.641e-07 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(Surface$avg.PCH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Surface$avg.PCH4)), col = "gray40", lwd=3, lty=1) # Data are skewed 
shapiro.test(log10(abs(Surface$avg.PCH4))) # p-value = 0.2351 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Bottom$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Bottom$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Bottom$avg.PCH4) # p-value = 0.004178 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(Bottom$avg.PCH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Bottom$avg.PCH4)), col = "gray40", lwd=3, lty=1) # Data are skewed 
shapiro.test(log10(abs(Bottom$avg.PCH4))) # p-value = 0.004178 # Distribution of data is significantly different from normal distribution # Krusal Wallis test

df <- group_by(cGas, CLASS_Z) %>%
  dplyr::select(avg.PCH4, CLASS_Z) %>%
  mutate(avg.PCH4_log = log10(abs(avg.PCH4)))
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCH4_log ~ CLASS_Z, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

# By lake environment

High <- subset(cGas, STAGE=="High")
Falling <- subset(cGas, STAGE=="Falling")
  
  # High-water stage

Open <- subset(High, ENVIRON=="Pelagic")
Edge <- subset(High, ENVIRON=="Edge")
Floodplain <- subset(High, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.PCH4) # p-value = 0.08656 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Edge$avg.PCH4) # p-value = 0.04739 # Distribution of data is significantly different from normal distribution
qqnorm(log10(abs(Edge$avg.PCH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Edge$avg.PCH4)), col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(log10(abs(Edge$avg.PCH4))) # p-value = 0.1081 # Distribution of data is not significantly different from normal distribution # ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Floodplain$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Floodplain$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Floodplain$avg.PCH4) # p-value = 0.007613 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(Floodplain$avg.PCH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Floodplain$avg.PCH4)), col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(log10(abs(Floodplain$avg.PCH4))) # p-value = 0.03501 # Distribution of data is significantly different from normal distribution # Krusal Wallis 

df <- group_by(High, ENVIRON) %>%
  dplyr::select(avg.PCH4, ENVIRON) %>%
  mutate(avg.PCH4_log = log10(abs(avg.PCH4)))
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCH4_log ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) # High-water Pelagic avg.PCH4_log is signficantly different than High-water Floodplain avg.PCH4_log (p-value = 0.0073185) 

Table(Open$avg.PCH4)
Table(Edge$avg.PCH4)
Table(Floodplain$avg.PCH4)
d(log10(abs(Open$avg.PCH4)), log10(abs(Floodplain$avg.PCH4))) # d = 1.259184

  # Falling-water stage

Open <- subset(Falling, ENVIRON=="Pelagic")
Edge <- subset(Falling, ENVIRON=="Edge")
Floodplain <- subset(Falling, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.PCH4) # p-value = 0.08656 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.PCH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.PCH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Edge$avg.PCH4) # p-value = 0.04739 # Distribution of data is significantly different from normal distribution
qqnorm(log10(abs(Edge$avg.PCH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Edge$avg.PCH4)), col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(log10(abs(Edge$avg.PCH4))) # p-value = 0.1081 # Distribution of data is not significantly different from normal distribution # ANOVA

df <- group_by(Falling, ENVIRON) %>%
  dplyr::select(avg.PCH4, ENVIRON) %>%
  mutate(avg.PCH4_log = log10(abs(avg.PCH4)))
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.PCH4_log ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) # Falling-water Pelagic avg.PCH4_log is signficantly different than High-water Floodplain avg.PCH4_log (p-value = 0.0085569) 

Table(Open$avg.PCH4)
Table(Edge$avg.PCH4)
d(log10(abs(Open$avg.PCH4)), log10(abs(Edge$avg.PCH4))) # d = 1.882848


##=============================== Normality Tests and Statistical Comparisons for D13C_CO2 

# Determine whether data follow normal or non-normal distributions in order to choose appropriate statistical test

# By flood stage

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(High$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(High$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(High$avg.D13C_CO2) # p-value = 0.001664 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(High$avg.D13C_CO2)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(High$avg.D13C_CO2)), col = "gray40", lwd=3, lty=1) # Data are skewed 
shapiro.test(log10(abs(High$avg.D13C_CO2))) # p-value = 0.04055 # Distribution of data is significantly different from normal distribution # Use Kruskal Wallis

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Falling$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Falling$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Falling$avg.D13C_CO2) # p-value = 0.3666 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(cGas, STAGE) %>%
  dplyr::select(avg.D13C_CO2, STAGE) %>%
  mutate(avg.D13C_CO2_log = log10(abs(avg.D13C_CO2)))
#View(df)
#df <- df[-c(3, 5, 8, 9), ]
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CO2_log ~ STAGE, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

Table(High$avg.D13C_CO2)
Table(Falling$avg.D13C_CO2)
d(log10(abs(High$avg.D13C_CO2)), log10(abs(Falling$avg.D13C_CO2))) 

# By depth

Surface <- subset(cGas, CLASS_Z=="Surface")
Bottom <- subset(cGas, CLASS_Z=="Bottom")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
#View(Surface)
#df <- df[-3, ]
qqnorm(Surface$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Surface$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Surface$avg.D13C_CO2) # p-value = 0.001642 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Bottom$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Bottom$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Bottom$avg.D13C_CO2) # p-value = 0.1479 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(cGas, CLASS_Z) %>%
  dplyr::select(avg.D13C_CO2, CLASS_Z) 
#View(df)
#df <- df[-c(3, 5, 8, 9), ]
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CO2 ~ CLASS_Z, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

# By lake environment

High <- subset(cGas, STAGE=="High")
Falling <- subset(cGas, STAGE=="Falling")

  # High-water stage

Open <- subset(High, ENVIRON=="Pelagic")
Edge <- subset(High, ENVIRON=="Edge")
Floodplain <- subset(High, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.D13C_CO2) # p-value = 0.7724 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Edge$avg.D13C_CO2) # p-value = 0.6789 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
#View(Floodplain)
#Floodplain <- Floodplain[-c(1, 4), ]
qqnorm(Floodplain$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Floodplain$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Floodplain$avg.D13C_CO2) # p-value = 0.001439 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(Floodplain$avg.D13C_CO2)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Floodplain$avg.D13C_CO2)), col = "gray40", lwd=3, lty=1) # Data are skewed 
shapiro.test(log10(abs(Floodplain$avg.D13C_CO2))) # p-value = 0.001439 # Distribution of data is significantly different from normal distribution # Use Kruskal Wallis

df <- group_by(High, ENVIRON) %>%
  dplyr::select(avg.D13C_CO2, ENVIRON) 
#View(df)
#df <- df[-c(3, 5, 8, 9),]
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CO2 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

Table(Open$avg.D13C_CO2)
Table(Edge$avg.D13C_CO2)
Table(Floodplain$avg.D13C_CO2)

  # Falling-water stage

Open <- subset(Falling, ENVIRON=="Pelagic")
Edge <- subset(Falling, ENVIRON=="Edge")
Floodplain <- subset(Falling, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.D13C_CO2) # p-value = 0.3719 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.D13C_CO2, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.D13C_CO2, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Edge$avg.D13C_CO2) # p-value = 0.8289 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(Falling, ENVIRON) %>%
  dplyr::select(avg.D13C_CO2, ENVIRON) 
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CO2 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

Table(Open$avg.D13C_CO2)
Table(Edge$avg.D13C_CO2)


##=============================== Normality Tests and Statistical Comparisons for D13C_CH4

# Determine whether data follow normal or non-normal distributions in order to choose appropriate statistical test

# By flood stage

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(High$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(High$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(High$avg.D13C_CH4) # p-value = 0.1153 # Distribution of data is not significantly different from normal distribution # Use ANOVA 

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Falling$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Falling$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Falling$avg.D13C_CH4) # p-value = 0.9617 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(cGas, STAGE) %>%
  dplyr::select(avg.D13C_CH4, STAGE) 
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CH4 ~ STAGE, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) # p-value = 8.9e-06 # High avg.D13C_CH4 is significantly different than Falling avg.D13C_CH4

Table(High$avg.D13C_CH4)
Table(Falling$avg.D13C_CH4)
d(High$avg.D13C_CH4, Falling$avg.D13C_CH4) # d = 1.663052

# By depth

Surface <- subset(cGas, CLASS_Z=="Surface")
Bottom <- subset(cGas, CLASS_Z=="Bottom")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Surface$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Surface$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Surface$avg.D13C_CH4) # p-value = 0.2068 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Bottom$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Bottom$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Bottom$avg.D13C_CH4) # p-value = 0.1322 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(cGas, CLASS_Z) %>%
  dplyr::select(avg.D13C_CH4, CLASS_Z) 
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CH4 ~ CLASS_Z, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

# By lake environment

High <- subset(cGas, STAGE=="High")
#View(High)
High <- High[-c(9, 19, 20, 27), ]
Falling <- subset(cGas, STAGE=="Falling")

  # High-water stage

Open <- subset(High, ENVIRON=="Pelagic")
Edge <- subset(High, ENVIRON=="Edge")
Floodplain <- subset(High, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are skewed # Log transform
shapiro.test(Open$avg.D13C_CH4) # p-value = 0.03421 # Distribution of data is significantly different from normal distribution 
qqnorm(log10(abs(Open$avg.D13C_CH4)), pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(log10(abs(Open$avg.D13C_CH4)), col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(log10(abs(Open$avg.D13C_CH4))) # p-value = 0.09062 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Edge$avg.D13C_CH4) # p-value = 0.5607 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Floodplain$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Floodplain$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Floodplain$avg.D13C_CH4) # p-value = 0.4169 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(High, ENVIRON) %>%
  dplyr::select(avg.D13C_CH4, ENVIRON) 
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CH4 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

Table(Open$avg.D13C_CH4)
Table(Edge$avg.D13C_CH4)
Table(Floodplain$avg.D13C_CH4)

  # Falling-water stage

Open <- subset(Falling, ENVIRON=="Pelagic")
Edge <- subset(Falling, ENVIRON=="Edge")
Floodplain <- subset(Falling, ENVIRON=="Floodplain")

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Open$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Open$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Open$avg.D13C_CH4) # p-value = 0.9442 # Distribution of data is not significantly different from normal distribution # Use ANOVA

par(mfrow=c(1, 1)) # mfrow=c(nrows, ncols)
par(oma=c(5, 5, 5, 5))
par(mar=c(5, 5, 3, 3)) # Bottom, left, top, right
par(xpd=FALSE)
qqnorm(Edge$avg.D13C_CH4, pch=1, cex=3, lwd=2, cex.axis=2, cex.lab=2, frame=FALSE, main=NULL)
qqline(Edge$avg.D13C_CH4, col = "gray40", lwd=3, lty=1) # Data are not skewed 
shapiro.test(Edge$avg.D13C_CH4) # p-value = 0.5368 # Distribution of data is not significantly different from normal distribution # Use ANOVA

df <- group_by(Falling, ENVIRON) %>%
  dplyr::select(avg.D13C_CH4, ENVIRON) 
#View(df)
# Compute the analysis of variance
res.aov <- aov(avg.D13C_CH4 ~ ENVIRON, data = df)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

Table(Open$avg.D13C_CH4)
Table(Edge$avg.D13C_CH4)


##===============================  Figure 1

setwd("~/Desktop/Working")

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

hypoxic_O2 <- ((3/1000)/32)*1000000 # 93.75 umol L-1 # Threshold O2 concentration for hypoxia in limnology and oceanography
hypoxic_O2 <- hypoxic_O2 - ((209460/1000000)*0.00118078*1000000)
hypoxic_O2 # -153.5762 umol L-1 # # Threshold O2 deficit for hypoxia in limnology and oceanography

cGas <- newData %>% 
  dplyr::select(DATE, STAGE, SITE, ENVIRON, CLASS_Z, O2, KH_O2, KH_PCO2, gradConcCO2, gradConcCH4, PCO2, ACTUAL_Z, O2_PERCENT) %>%
  mutate(CO2_conc = (PCO2/1000000)*KH_PCO2*1000000) %>% # CO2 in umol^L-1
  mutate(O2_conc = ((O2/1000)/32)*1000000) # O2 in umol^L-1
  #mutate(iCO2_umolL = (gradConcCO2/44.01)*1000000) %>% # Excess CO2 in umol^L-1
  #mutate(iCH4_umolL = (gradConcCH4/16.04)*1000000) %>% # Excess CH4 in umol^L-1
  #mutate(iO2_umolL = ((O2/1000)/32) - (0.00118078*(209460/1000000))*1000000) # O2 deficit in umol^L-1
#cGas <- cGas[cGas$CLASS_Z == "Surface", ]
cGas <- cGas[!cGas$SITE == "StSe", ] 
cGas <- cGas[!cGas$SITE == "KgPr", ] 
cGas <- cGas[!cGas$SITE == "PtSd.O", ]
cGas <- cGas[!cGas$SITE == "PtSd.E", ] 
cGas <- cGas[!cGas$SITE == "PtSd.F", ] 
cGas <- cGas[!cGas$SITE == "KgKl.O", ] 
cGas <- cGas[!cGas$SITE == "KgKl.E", ] 
cGas <- cGas[!cGas$STAGE == "Rising", ] 
cGas <- cGas[complete.cases(cGas), ]
#View(cGas)
cGas <- cGas %>% 
  dplyr::select(DATE, STAGE, SITE, ENVIRON, CLASS_Z, ACTUAL_Z, CO2_conc, O2_conc, PCO2, O2) %>%
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(CO2_umolL = mean(CO2_conc)-((400/1000000)*0.029098*1000000)) %>%
  mutate(O2_umolL = mean(O2_conc)-((209460/1000000)*0.00118078*1000000)) %>% 
  distinct(CO2_umolL, .keep_all=TRUE) %>%
  ungroup() 
#View(cGas)

High <- cGas[cGas$STAGE == "High", ]
Falling <- cGas[cGas$STAGE == "Falling", ]

new_High <- data.frame(CO2_umolL = seq(min(High$CO2_umolL), max(High$CO2_umolL), length.out = 100)) 
lm_y <- lm(O2_umolL ~ CO2_umolL, data=High)
summary(lm_y) # Multiple R-squared:  0.2112 # p-value: 0.005483
pred_High <- predict(lm_y, newdata=new_High, se.fit=TRUE, type="response")
new_High$pred <- predict(lm_y, newdata=new_High)
#plot(High$CO2_mmolL, resid(lm_y))

new_Falling <- data.frame(CO2_umolL = seq(min(Falling$CO2_umolL), max(Falling$CO2_umolL), length.out = 100)) 
lm_y <- lm(O2_umolL ~ CO2_umolL, data=Falling)
summary(lm_y) # Multiple R-squared:  0.09874 # p-value: 0.2957
pred_Falling <- predict(lm_y, newdata=new_Falling, se.fit=TRUE, type="response")
new_Falling$pred <- predict(lm_y, newdata=new_Falling)
#plot(Falling$CO2_mmolL, resid(lm_y))

par(mfrow=c(3, 2))
par(oma=c(4, 4, 4, 4))
par(mar=c(4, 8, 5, 2)) # Bottom, left, top, right
range(cGas$CO2_umolL) #  51.66945 963.43589
range(cGas$O2_umolL) # -247.3262 -247.3259
par(xpd=FALSE)
plot(cGas$CO2_umolL, cGas$O2_umolL, xlim=c(-100, 900), ylim=c(-300, 100), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="blue3")
box(lty=1, lwd=1, col="black")

abline(0, -1, col="black", lty=2, lwd=2)
abline(h=0, col="darkorange", lwd=2)
abline(v=0, col="darkorange", lwd=2)

polygon(c(seq(min(High$CO2_umolL), max(High$CO2_umolL), length.out = 100), rev(seq(min(High$CO2_umolL), max(High$CO2_umolL), length.out = 100))), 
        c(pred_High$fit - (1.96 * pred_High$se.fit), rev(pred_High$fit + (1.96 * pred_High$se.fit))), 
        col=adjustcolor("cornflowerblue", alpha.f=0.4), border = NA) # High
lines(new_High$CO2_umolL, new_High$pred, lwd=2, lty=2, col="cornflowerblue") 

mtext(expression(O[2] ~ Deficit ~ (mu*mol ~ L^{-1})), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(CO[2] ~ "Supersaturation" ~ (mu*mol ~ L^{-1})), side=1, line=4.5, cex=1.8, font=1)
axis(2, ylim=c(-300, 100), col="black", las=0, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-100, 900), las=1, cex=2, cex.axis=2, cex.lab=2, font=1)  #las=1 makes horizontal labels

levels(factor(High$ENVIRON)) # Open pch=21 # Edge pch=22 # Floodplain pch=25
points(jitter(High$CO2_umolL, factor=1), High$O2_umolL, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 25, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)

legend("topright", legend=expression("Aerobic" ~ "ER" ~ italic(m) ~ "=" ~ "-1.0"),
       text.col="black", text.font=3, cex=1.8, bty="n")
legend("bottomright", legend=c(expression(italic(R^2) ~ "=" ~ "0.21"),
                               expression(italic(p) ~ "=" ~ "0.005"),
                               expression(italic(df) ~ "=" ~ "33")),
       text.col="cornflowerblue", text.font=3, cex=1.8, bty="n")
legend("right", legend=expression(italic(m) ~ "=" ~ "-0.1"), 
       text.col="cornflowerblue", text.font=3, cex=1.8, bty="n")

par(xpd=TRUE)
legend("top", inset=-0.25, legend=c("Open", "Edge", "Floodplain"), 
       text.col=c("black"), horiz=TRUE, text.font=c(1, 1, 1),  pch=c(21, 22, 25), col="black", bg="white", pt.cex=c(2.5, 2.5, 2.5), cex=1.8, bty="n")
legend("topleft", inset=-0.21, legend=expression("a)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")

#arrows(-275, -153.5762, -275, -300, lty=1, lwd=2, length=0.1, xpd=TRUE, col="gray40")
#text(-275, -90, labels="Hypoxia", col="gray40", cex=1.8, font=1, srt=90)

#par(mfrow=c(3, 2))
#par(oma=c(4, 4, 4, 4))
par(mar=c(4, 3, 5, 7)) # Bottom, left, top, right
range(cGas$CO2_umolL) #  51.66945 963.43589
range(cGas$O2_umolL) # -247.3262 -247.3259
par(xpd=FALSE)
plot(cGas$CO2_umolL, cGas$O2_umolL, xlim=c(-100, 900), ylim=c(-300, 100), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="blue3")
box(lty=1, lwd=1, col="black")

abline(0, -1, col="black", lty=2, lwd=2)
abline(h=0, col="darkorange", lwd=2)
abline(v=0, col="darkorange", lwd=2)

#polygon(c(seq(min(Falling$CO2_mmolL), max(Falling$CO2_mmolL), length.out = 100), rev(seq(min(Falling$CO2_mmolL), max(Falling$CO2_mmolL), length.out = 100))), 
        #c(pred_Falling$fit - (1.96 * pred_Falling$se.fit), rev(pred_Falling$fit + (1.96 * pred_Falling$se.fit))), 
        #col=adjustcolor("blue3", alpha.f=0.4), border = NA) # Falling
lines(new_Falling$CO2_umolL, new_Falling$pred, lwd=2, lty=2, col="blue3") 

#mtext(expression(O[2] ~ Deficit ~ (mu*mol ~ L^{-1})), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(CO[2] ~ "Supersaturation" ~ (mu*mol ~ L^{-1})), side=1, line=4.5, cex=1.8, font=1)
#axis(2, ylim=c(-350, 50), col="black", las=0, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-100, 900), las=1, cex=2, cex.axis=2, cex.lab=2, font=1)  #las=1 makes horizontal labels

levels(factor(Falling$ENVIRON)) # Open pch=21 # Edge pch=22 # Floodplain pch=25
points(jitter(Falling$CO2_umolL, factor=1), Falling$O2_umolL, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 25, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)
View(Falling)
legend("topright", legend=expression("Atm." ~ CO[2] ~ "&" ~ O[2]), 
       text.col="darkorange", text.font=3, cex=1.8, bty="n")
legend("right", legend=expression(italic(m) ~ "=" ~ "-0.1"), 
       text.col="blue3", text.font=3, cex=1.8, bty="n")

par(xpd=TRUE)
legend("top", inset=-0.25, legend=c("High-water", "Falling-water"), 
       text.col=c("cornflowerblue", "blue3"), horiz=TRUE, text.font=c(1, 1, 1),  cex=1.8, bty="n")
legend("topleft", inset=-0.2, legend=expression("b)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")

#arrows(-220, -153.5762, -220, -300, lty=1, lwd=2, length=0.1, xpd=TRUE, col="gray40")
#text(-220, -90, labels="Hypoxia", col="gray40", cex=1.8, font=1, srt=90)


##===============================  Figure 3a

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

cGas <- newData %>% 
  dplyr::select(DATE, CLASS_Z, STAGE, SITE, ENVIRON, PCH4, D13C_CH4, PCO2, D13C_CO2) 
#mutate(PCH4_x = (1/PCH4)*1000) # Corresponds to 1/ppb of CH4
#cGas <- cGas[cGas$CLASS_Z == "Surface", ]
cGas <- cGas[!cGas$SITE == "StSe", ] # Omit Stueng Saen tributary site
cGas <- cGas[!cGas$SITE == "KgPr", ] # Omit Kampong Prak tributary site
#cGas$PCH4 <- log10(abs(cGas$PCH4)) 
#cGas <- cGas[!cGas$PCH4 == "-Inf", ]
cGas <- cGas[complete.cases(cGas), ]
#names(cGas)
#View(cGas)
cGas <- cGas %>% 
  dplyr::select(DATE, STAGE, SITE, ENVIRON, CLASS_Z, PCH4, D13C_CH4, PCO2, D13C_CO2) %>%
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCH4 = mean(PCH4)) %>%
  mutate(avg.D13C_CH4 = mean(D13C_CH4)) %>% 
  mutate(avg.PCO2 = mean(PCO2)) %>% 
  mutate(avg.D13C_CO2 = mean(D13C_CO2)) %>% 
  distinct(avg.PCH4, .keep_all=TRUE) %>%
  ungroup() 
#View(cGas)

# High

Data <- subset(cGas, STAGE=="High") 
#Data <- Data[-c(3, 5, 8, 9), ]
range(Data$D13C_CO2)
range(Data$D13C_CH4)
names(Data)
#View(Data)

par(mfrow=c(3, 2))
par(oma=c(3, 3, 3, 3))
par(mar=c(4, 8, 5, 2)) # Bottom, left, top, right par(mar=c(4, 2, 4, 8)) # Bottom, left, top, right
par(xpd=FALSE)
range(Data$D13C_CH4) # -54.325 -14.747
range(cGas$D13C_CO2) # -58.843 -29.194
plot(Data$avg.D13C_CH4, Data$avg.D13C_CO2, xlim=c(-85, -10), ylim=c(-70, -10), xaxt="n", yaxt="n", xlab="", ylab="", pch="", 
     font=2, las=1, cex=3, col="gray40", bg="gray")
box(lty=1, lwd=1, col="black")

abline(5, 1, lty=2, lwd=1, col="black")
legend(-30, -15, legend=expression(epsilon[C] ~ "=" ~ "5"), text.col=c("black"), text.font=4, cex=1.8, bty="n")
abline(20, 1, lty=2, lwd=1, col="black")
#legend(-50, -20, legend=expression(epsilon[C] ~ "=" ~ "20"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(30, 1, lty=2, lwd=1, col="black")
legend(-55, -15, legend=expression(epsilon[C] ~ "=" ~ "30"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(40, 1, lty=2, lwd=1, col="black")
legend(-70, -20, legend=expression(epsilon[C] ~ "=" ~ "40"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(55, 1, lty=2, lwd=1, col="black")
legend(-80, -15, legend=expression(epsilon[C] ~ "=" ~ "55"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
polygon
polygon(x=c(-72, -50, -50, -72), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray50", alpha.f=0.2), border="gray40")
#text(-59, -16, expression(italic("Acetate \n Ferm.")), cex=1.55, font=3)

polygon(x=c(-100, -69.5, -69.5, -100), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray50", alpha.f=0.2), border="gray40")
text(-78, -16, expression(italic("Carbonate")), cex=1.8, font=3)

polygon(x=c(-65, -10, -10, -65), y=c(-35, -35, -75, -75), col=adjustcolor("gray50", alpha.f=0.2), border="gray40")
text(-35, -57, expression(italic("Oxid.")), cex=1.8, font=3)

points(jitter(Data$avg.D13C_CH4, factor=4), Data$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(Data$ENVIRON))], cex=3)
#points(jitter(Edge$avg.D13C_CH4, factor=4), Edge$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=22, cex=3)
#points(jitter(Floodplain$avg.D13C_CH4, factor=4), Floodplain$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=25, cex=3)
points(-48, -10, col="darkorange", pch=19, cex=1.5)
legend(-54, -7, legend=expression("Atm." ~ CO[2] ~ "&" ~ CH[4]), text.col=c("darkorange"), text.font=4, cex=1.8, bty="n")
text(-59, -16, expression(italic("Acetate")), cex=1.8, font=3)
legend(-50, -20, legend=expression(epsilon[C] ~ "=" ~ "20"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
mtext(expression(paste(delta^13 ~ C-CO[2] ~ "(\u2030)")), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(paste(delta^13 ~ C-CH[4] ~ "(\u2030)")), side=1, line=4.5, cex=1.8, font=1)
axis(2, ylim=c(-70, -10), at=c(-10, -25, -40, -55, -70), 
     labels=c("-10", "-25", "-40", "-55", "-70"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-85, -10), at=c(-10, -25, -40, -55, -70, -85), 
     labels=c("-10", "-25", "-40", "-65", "-70", "-85"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
par(xpd=TRUE)
legend("top", inset=-0.25, legend=c("Open", "Edge", "Floodplain"), 
       text.col=c("black"), horiz=TRUE, text.font=c(1, 1, 1),  pch=c(21, 22, 25), col="black", bg="white", pt.cex=c(2.5, 2.5, 2.5), cex=1.8, bty="n")
legend("topleft", inset=-0.21, legend=expression("a)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")


##===============================  Figure 3b

# Falling

Data <- subset(cGas, STAGE=="Falling") 
#names(Data)
#View(Data)

#par(mfrow=c(3, 2))
#par(oma=c(3, 3, 3, 3))
par(mar=c(4, 3, 5, 7)) # Bottom, left, top, right
par(xpd=FALSE)
range(Data$D13C_CH4) # -54.325 -14.747
range(cGas$D13C_CO2) # -58.843 -29.194
plot(Data$avg.D13C_CH4, Data$avg.D13C_CO2, xlim=c(-85, -10), ylim=c(-70, -10), xaxt="n", yaxt="n", xlab="", ylab="", pch="", 
     font=2, las=1, cex=3, col="gray40", bg="gray")
box(lty=1, lwd=1, col="black")

abline(5, 1, lty=2, lwd=1, col="black")
legend(-30, -15, legend=expression(epsilon[C] ~ "=" ~ "5"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(20, 1, lty=2, lwd=1, col="black")
legend(-50, -20, legend=expression(epsilon[C] ~ "=" ~ "20"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(30, 1, lty=2, lwd=1, col="black")
legend(-55, -15, legend=expression(epsilon[C] ~ "=" ~ "30"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(40, 1, lty=2, lwd=1, col="black")
legend(-70, -20, legend=expression(epsilon[C] ~ "=" ~ "40"), text.col=c("black"), text.font=4, cex=1.75, bty="n")
abline(55, 1, lty=2, lwd=1, col="black")
legend(-80, -15, legend=expression(epsilon[C] ~ "=" ~ "55"), text.col=c("black"), text.font=4, cex=1.75, bty="n")

polygon(x=c(-72, -50, -50, -72), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray50", alpha.f=0.2), border="gray40")
text(-59, -16, expression(italic("Acetate")), cex=1.8, font=3)

polygon(x=c(-100, -69.5, -69.5, -100), y=c(-31.5, -31.5, -5, -5), col=adjustcolor("gray50", alpha.f=0.2), border="gray40")
text(-78, -16, expression(italic("Carbonate")), cex=1.8, font=3)

polygon(x=c(-65, -10, -10, -65), y=c(-35, -35, -75, -75), col=adjustcolor("gray50", alpha.f=0.2), border="gray40")
text(-35, -57, expression(italic("Oxid.")), cex=1.8, font=3)

points(jitter(Data$avg.D13C_CH4, factor=2), Data$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 24, 21)[as.numeric(as.factor(Data$ENVIRON))], cex=3)
#points(jitter(Edge$avg.D13C_CH4, factor=2), Edge$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=22, cex=3)
#points(jitter(Floodplain$D13C_CH4, factor=2), Floodplain$D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=25, cex=3)
points(-48, -10, col="darkorange", pch=19, cex=1.5)
legend(-54, -7, legend=expression("Atm." ~ CO[2] ~ "&" ~ CH[4]), text.col=c("darkorange"), text.font=4, cex=1.8, bty="n")

#mtext(expression(paste(delta^13 ~ C-CO[2] ~ "(\u2030)")), side=2, line=4.5, cex=1.8, font=1)
mtext(expression(paste(delta^13 ~ C-CH[4] ~ "(\u2030)")), side=1, line=4.5, cex=1.8, font=1)
#axis(2, ylim=c(-70, -10), at=c(-10, -25, -40, -55, -70), 
     #labels=c("-10", "-25", "-40", "-55", "-70"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(-85, -10), at=c(-10, -25, -40, -55, -70, -85), 
     labels=c("-10", "-25", "-40", "-65", "-70", "-85"), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
par(xpd=TRUE)
legend("top", inset=-0.25, legend=c("High-water", "Falling-water"), 
       text.col=c("cornflowerblue", "blue3"), horiz=TRUE, text.font=c(1, 1, 1),  cex=1.8, bty="n")
legend("topleft", inset=-0.2, legend=expression("b)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")


##===============================  Figure 3c

setwd("~/Desktop/Working")

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

cGas <- newData %>%
  dplyr::select(SITE, CLASS_Z, DATE, PCO2, PCH4, D13C_CO2, D13C_CH4, ENVIRON, STAGE, O2, MOX, ENVIRON) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.MOX = mean(MOX, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.PCO2 = mean(PCO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.O2 = mean(O2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  #mutate(avg.ORP = mean(ORP, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.PCH4 = mean(PCH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CO2 = mean(D13C_CO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CH4 = mean(D13C_CH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(eC = avg.D13C_CO2 - avg.D13C_CH4) %>% # Apply the function mean() to every column
  distinct(avg.PCO2, .keep_all=TRUE) %>%
  ungroup() 
cGas <- cGas[!cGas$SITE == "KgPr", ]
cGas <- cGas[!cGas$SITE == "StSe", ]
cGas <- cGas[!cGas$SITE == "PtSd.O", ] 
cGas <- cGas[!cGas$SITE == "PtSd.E", ] 
cGas <- cGas[!cGas$SITE == "PtSd.F", ] 
cGas <- cGas[!cGas$SITE == "KgKl.O", ] 
cGas <- cGas[!cGas$SITE == "KgKl.E", ] 
cGas <- cGas[!cGas$STAGE == "Rising", ] 
cGas <- cGas[complete.cases(cGas), ]
#View(cGas)

cGas <- cGas %>%  
  dplyr::select(avg.D13C_CO2, avg.MOX, STAGE, ENVIRON, CLASS_Z) %>%
  mutate(avg.MOX_cbrt = sign(avg.MOX*24) * (abs(avg.MOX*24))^(1/2) * -1) # Convert from nmol L-1 h-1 to nmol L-1 d-1 and take cubed root
cGas <- cGas[-c(3, 5, 8, 9), ]
#View(cGas)

High <- subset(cGas, STAGE=="High")
Falling <- subset(cGas, STAGE=="Falling")

summary(lm(avg.D13C_CO2 ~ avg.MOX_cbrt, data=High))
summary(lm(avg.D13C_CO2 ~ avg.MOX, data=High))
range(High$avg.MOX)

new_High <- data.frame(avg.MOX_cbrt = seq(min(High$avg.MOX_cbrt), max(High$avg.MOX_cbrt), length.out = 100)) 
lm_y <- lm(avg.D13C_CO2 ~ avg.MOX_cbrt, data=High)
poly_y <- lm(avg.D13C_CO2 ~ avg.MOX_cbrt + I(avg.MOX_cbrt^2), data=High)
gam_y <- gam(avg.D13C_CO2 ~ s(avg.MOX_cbrt, k=5), method="REML", data=High)
AICc1 <- AICc(lm_y, k=2, REML=NULL)
AICc2 <- AICc(poly_y, k=2, REML=NULL)
AICc3 <- AICc(gam_y, k=2, REML=NULL)
delAICc1 <- AICc1 - min(c(AICc1, AICc2, AICc3)) # Delta AICc
delAICc1 
delAICc2 <- AICc2 - min(c(AICc1, AICc2, AICc3))
delAICc2 
delAICc3 <- AICc3 - min(c(AICc1, AICc3))
delAICc3 
summary(lm_y) # Multiple R-squared:  0.3372 # p-value: 0.0007668
pred_High <- predict(lm_y, newdata=new_High, se.fit=TRUE, type="response")
new_High$pred <- predict(lm_y, newdata=new_High)
#plot(High$avg.MOX_cbrt, resid(lm_y))

new_Falling <- data.frame(avg.MOX_cbrt = seq(min(Falling$avg.MOX_cbrt), max(Falling$avg.MOX_cbrt), length.out = 100)) 
lm_y <- lm(avg.D13C_CO2 ~ avg.MOX_cbrt, data=Falling)
poly_y <- lm(avg.D13C_CO2 ~ avg.MOX_cbrt + I(avg.MOX_cbrt^2), data=Falling)
gam_y <- gam(avg.D13C_CO2 ~ s(avg.MOX_cbrt, k=5), method="REML", data=Falling)
AICc1 <- AICc(lm_y, k=2, REML=NULL)
AICc2 <- AICc(poly_y, k=2, REML=NULL)
AICc3 <- AICc(gam_y, k=2, REML=NULL)
delAICc1 <- AICc1 - min(c(AICc1, AICc2, AICc3)) # Delta AICc
delAICc1 
delAICc2 <- AICc2 - min(c(AICc1, AICc2, AICc3))
delAICc2 
delAICc3 <- AICc3 - min(c(AICc1, AICc3))
delAICc3 
summary(lm_y) 
pred_Falling <- predict(poly_y, newdata=new_Falling, se.fit=TRUE, type="response")
new_Falling$pred <- predict(poly, newdata=new_Falling)
#plot(Falling$avg.MOX_cbrt, resid(lm_y))

#par(mfrow=c(3, 2))
#par(oma=c(4, 4, 4, 4))
par(mar=c(5, 8, 4, 2)) # Bottom, left, top, right
range(High$avg.MOX_cbrt) #  0.6214465 9.7225042
range(High$avg.D13C_CO2) # -46.78400 -30.97733
par(xpd=FALSE)
plot(cGas$avg.MOX_cbrt, cGas$avg.D13C_CO2, xlim=c(0, 40), ylim=c(-50, -30), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="cornflowerblue")
box(lty=1, lwd=1, col="black")

polygon(c(seq(min(High$avg.MOX_cbrt), max(High$avg.MOX_cbrt), length.out = 100), rev(seq(min(High$avg.MOX_cbrt), max(High$avg.MOX_cbrt), length.out = 100))), 
        c(pred_High$fit - (1.96 * pred_High$se.fit), rev(pred_High$fit + (1.96 * pred_High$se.fit))), 
        col=adjustcolor("cornflowerblue", alpha.f=0.4), border = NA) # High
lines(new_High$avg.MOX_cbrt, new_High$pred, lwd=2, lty=2, col="cornflowerblue") 

mtext(expression("Net" ~ CH[4] ~ "Oxid." ~ (nmol ~ L^{-1} ~ d^{-1})), side=1, line=4.5, cex=1.8, font=1)
mtext(expression(paste(delta^13 ~ C-CO[2] ~ "(\u2030)")), side=2, line=4.5, cex=1.8, font=1)
axis(2, ylim=c(-50, -30), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(0, 40), at=c(0, 10, 20, 30, 40), 
     labels=c(expression(0), expression(10), expression(20), expression(30), expression(40)), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 

levels(factor(High$ENVIRON)) # Open pch=21 # Edge pch=22 # Floodplain pch=25
points(jitter(High$avg.MOX_cbrt, factor=1), High$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 25, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)

legend("topright", legend=c(expression(italic(R^2) ~ "=" ~ "0.33"),
                               expression(italic(p) ~ "<" ~ "0.001"),
                               expression(italic(df) ~ "=" ~ "28")),
       text.col="cornflowerblue", text.font=3, cex=1.8, bty="n")

par(xpd=TRUE)
legend("topleft", inset=-0.21, legend=expression("c)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")


##===============================  Figure 3d

#par(mfrow=c(3, 2))
#par(oma=c(4, 4, 4, 4))
par(mar=c(5, 3, 4, 7)) # Bottom, left, top, right
range(Falling$avg.MOX_cbrt) #  -2.186975  0.664437
range(Falling$avg.D13C_CO2) # -49.41767 -34.37233
par(xpd=FALSE)
plot(cGas$avg.MOX_cbrt, cGas$avg.D13C_CO2, xlim=c(0, 40), ylim=c(-50, -30), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="blue3")
box(lty=1, lwd=1, col="black")

#polygon(c(seq(min(Falling$avg.MOX_cbrt), max(Falling$avg.MOX_cbrt), length.out = 100), rev(seq(min(Falling$avg.MOX_cbrt), max(Falling$avg.MOX_cbrt), length.out = 100))), 
        #c(pred_Falling$fit - (1.96 * pred_Falling$se.fit), rev(pred_Falling$fit + (1.96 * pred_Falling$se.fit))), 
        #col=adjustcolor("blue3", alpha.f=0.4), border = NA) # High
#lines(new_Falling$avg.MOX_cbrt, new_Falling$pred, lwd=2, lty=2, col="blue3") 

mtext(expression("Net" ~ CH[4] ~ "Oxid." ~ (nmol ~ L^{-1} ~ d^{-1})), side=1, line=4.5, cex=1.8, font=1)
#mtext(expression(paste(delta^13 ~ C-CO[2] ~ "(\u2030)")), side=2, line=4.5, cex=1.8, font=1)
#axis(2, ylim=c(-50, -30), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(0, 40), at=c(0, 10, 20, 30, 40), 
     labels=c(expression(0), expression(10), expression(20), expression(30), expression(40)), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1)

levels(factor(Falling$ENVIRON)) # Open pch=21 # Edge pch=22 # Floodplain pch=25
points(jitter(Falling$avg.MOX_cbrt, factor=1), Falling$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 21)[as.numeric(as.factor(High$ENVIRON))], cex=3)

#legend("bottomright", legend=c(expression(italic(R^2) ~ "=" ~ "0.34"),
                               #expression(italic(p) ~ "<" ~ "0.001"),
                               #expression(italic(df) ~ "=" ~ "28")),
       #text.col="blue3", text.font=3, cex=1.8, bty="n")

par(xpd=TRUE)
legend("topleft", inset=-0.2, legend=expression("d)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")


##===============================  Figure 3e

setwd("~/Desktop/Working")

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

MProd <- subset(newData, CLASS_Z=="Sediment") %>%
  dplyr::select(SITE, DATE, STAGE, MPROD, ENVIRON) %>% # Select columns of interest 
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE) %>% # Group data by these columns
  mutate(avg.MPROD = mean(MPROD, na.rm=TRUE)) %>%
  distinct(avg.MPROD, .keep_all=TRUE) %>%
  ungroup()
MProd <- MProd[!MProd$SITE == "KgPr", ]
MProd <- MProd[!MProd$SITE == "StSe", ]
MProd <- MProd[!MProd$SITE == "PtSd.O", ] 
MProd <- MProd[!MProd$SITE == "PtSd.E", ] 
MProd <- MProd[!MProd$SITE == "PtSd.F", ] 
MProd <- MProd[!MProd$SITE == "KgKl.O", ] 
MProd <- MProd[!MProd$SITE == "KgKl.E", ] 
MProd <- MProd[complete.cases(MProd), ]
#View(MProd)

del13C <- subset(newData, CLASS_Z=="Bottom") %>%
  dplyr::select(SITE, DATE, STAGE, D13C_CO2, ENVIRON) %>% # Select columns of interest 
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE) %>% # Group data by these columns
  mutate(avg.D13C_CO2 = mean(D13C_CO2, na.rm=TRUE)) %>%
  distinct(avg.D13C_CO2, .keep_all=TRUE) %>%
  ungroup()
del13C <- del13C[!del13C$SITE == "KgPr", ]
del13C <- del13C[!del13C$SITE == "StSe", ]
del13C <- del13C[!del13C$SITE == "PtSd.O", ] 
del13C <- del13C[!del13C$SITE == "PtSd.E", ] 
del13C <- del13C[!del13C$SITE == "PtSd.F", ] 
del13C <- del13C[!del13C$SITE == "KgKl.O", ] 
del13C <- del13C[!del13C$SITE == "KgKl.E", ] 
del13C <- del13C[complete.cases(del13C), ]
#View(del13C)

cGas <- left_join(MProd, del13C, by=c("SITE", "DATE"))
#View(cGas)

cGas <- cGas %>%  
  dplyr::select(avg.D13C_CO2, avg.MPROD, STAGE.x, ENVIRON.x) %>%
  mutate(avg.MPROD_cbrt = sign(avg.MPROD*24) * (abs(avg.MPROD*24))^(1/2)) # Convert from nmol L-1 h-1 to nmol L-1 d-1 and take cubed root
cGas <- cGas[-c(1, 4), ]
#View(cGas)

High <- subset(cGas, STAGE.x=="High")
Falling <- subset(cGas, STAGE.x=="Falling")

summary(lm(avg.D13C_CO2 ~ avg.MPROD_cbrt, data=High))
summary(lm(avg.D13C_CO2 ~ avg.MPROD, data=High))
range(High$avg.MPROD)

new_High <- data.frame(avg.MPROD_cbrt = seq(min(High$avg.MPROD_cbrt), max(High$avg.MPROD_cbrt), length.out = 100)) 
lm_y <- lm(avg.D13C_CO2 ~ avg.MPROD_cbrt, data=High)
poly_y <- lm(avg.D13C_CO2 ~ avg.MPROD_cbrt + I(avg.MPROD_cbrt^2), data=High)
gam_y <- gam(avg.D13C_CO2 ~ s(avg.MPROD_cbrt, k=5), method="REML", data=High)
AICc1 <- AICc(lm_y, k=2, REML=NULL)
AICc2 <- AICc(poly_y, k=2, REML=NULL)
AICc3 <- AICc(gam_y, k=2, REML=NULL)
delAICc1 <- AICc1 - min(c(AICc1, AICc2, AICc3)) # Delta AICc
delAICc1 
delAICc2 <- AICc2 - min(c(AICc1, AICc2, AICc3))
delAICc2 
delAICc3 <- AICc3 - min(c(AICc1, AICc3))
delAICc3 
summary(lm_y) # Multiple R-squared:  0.5268 # p-value: 0.001458
pred_High <- predict(lm_y, newdata=new_High, se.fit=TRUE, type="response")
new_High$pred <- predict(lm_y, newdata=new_High)
#plot(High$avg.MPROD_cbrt, resid(lm_y))

new_Falling <- data.frame(avg.MPROD_cbrt = seq(min(Falling$avg.MPROD_cbrt), max(Falling$avg.MPROD_cbrt), length.out = 100)) 
lm_y <- lm(avg.D13C_CO2 ~ avg.MPROD_cbrt, data=Falling)
poly_y <- lm(avg.D13C_CO2 ~ avg.MPROD_cbrt + I(avg.MPROD_cbrt^2), data=Falling)
gam_y <- gam(avg.D13C_CO2 ~ s(avg.MPROD_cbrt, k=5), method="REML", data=Falling)
AICc1 <- AICc(lm_y, k=2, REML=NULL)
AICc2 <- AICc(poly_y, k=2, REML=NULL)
AICc3 <- AICc(gam_y, k=2, REML=NULL)
delAICc1 <- AICc1 - min(c(AICc1, AICc2, AICc3)) # Delta AICc
delAICc1 
delAICc2 <- AICc2 - min(c(AICc1, AICc2, AICc3))
delAICc2 
delAICc3 <- AICc3 - min(c(AICc1, AICc3))
delAICc3 
summary(lm_y) 
pred_Falling <- predict(poly_y, newdata=new_Falling, se.fit=TRUE, type="response")
new_Falling$pred <- predict(poly, newdata=new_Falling)
#plot(Falling$avg.MPROD_cbrt, resid(lm_y))

#par(mfrow=c(3, 2))
#par(oma=c(4, 4, 4, 4))
par(mar=c(5, 8, 4, 2)) # Bottom, left, top, right
range(High$avg.MPROD_cbrt) #  1.149779 9.661469
range(High$avg.D13C_CO2) # -46.78400 -36.55967
par(xpd=FALSE)
plot(cGas$avg.MPROD_cbrt, cGas$avg.D13C_CO2, xlim=c(0, 40), ylim=c(-50, -30), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="cornflowerblue")
box(lty=1, lwd=1, col="black")

polygon(c(seq(min(High$avg.MPROD_cbrt), max(High$avg.MPROD_cbrt), length.out = 100), rev(seq(min(High$avg.MPROD_cbrt), max(High$avg.MPROD_cbrt), length.out = 100))), 
        c(pred_High$fit - (1.96 * pred_High$se.fit), rev(pred_High$fit + (1.96 * pred_High$se.fit))), 
        col=adjustcolor("cornflowerblue", alpha.f=0.4), border = NA) # High
lines(new_High$avg.MPROD_cbrt, new_High$pred, lwd=2, lty=2, col="cornflowerblue") 

mtext(expression(CH[4] ~ "Prod." ~ (nmol ~ cm^{-3} ~ d^{-1})), side=1, line=4.5, cex=1.8, font=1)
mtext(expression(paste(delta^13 ~ C-CO[2] ~ "(\u2030)")), side=2, line=4.5, cex=1.8, font=1)
axis(2, ylim=c(-50, -30), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(0, 40), at=c(0, 10, 20, 30, 40), 
     labels=c(expression(0), expression(10), expression(20), expression(30), expression(40)), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1)

levels(factor(High$ENVIRON.x)) # Open pch=21 # Edge pch=22 # Floodplain pch=25
points(jitter(High$avg.MPROD_cbrt, factor=1), High$avg.D13C_CO2, bg=adjustcolor("cornflowerblue", alpha.f=0.8), col="black", pch=c(22, 25, 21)[as.numeric(as.factor(High$ENVIRON.x))], cex=3)

legend("topright", legend=c(expression(italic(R^2) ~ "=" ~ "0.44"),
                            expression(italic(p) ~ "=" ~ "0.005"),
                            expression(italic(df) ~ "=" ~ "14")),
       text.col="cornflowerblue", text.font=3, cex=1.8, bty="n")

par(xpd=TRUE)
legend("topleft", inset=-0.21, legend=expression("e)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")


##===============================  Figure 3f

#par(mfrow=c(3, 2))
#par(oma=c(4, 4, 4, 4))
par(mar=c(5, 3, 4, 7)) # Bottom, left, top, right
range(Falling$avg.MPROD_cbrt) #  1.841033 11.654009
range(Falling$avg.D13C_CO2) # -49.41767 -34.37233
par(xpd=FALSE)
plot(cGas$avg.MPROD_cbrt, cGas$avg.D13C_CO2, xlim=c(0, 40), ylim=c(-50, -30), xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch="", 
     font=2, las=0, cex=3, col="black", bg="blue3")
box(lty=1, lwd=1, col="black")

#polygon(c(seq(min(Falling$avg.MPROD_cbrt), max(Falling$avg.MPROD_cbrt), length.out = 100), rev(seq(min(Falling$avg.MPROD_cbrt), max(Falling$avg.MPROD_cbrt), length.out = 100))), 
#c(pred_Falling$fit - (1.96 * pred_Falling$se.fit), rev(pred_Falling$fit + (1.96 * pred_Falling$se.fit))), 
#col=adjustcolor("blue3", alpha.f=0.4), border = NA) # High
#lines(new_Falling$avg.MPROD_cbrt, new_Falling$pred, lwd=2, lty=2, col="blue3") 

mtext(expression(CH[4] ~ "Prod." ~ (nmol ~ cm^{-3} ~ d^{-1})), side=1, line=4.5, cex=1.8, font=1)
#mtext(expression(paste(delta^13 ~ C-CO[2] ~ "(\u2030)")), side=2, line=4.5, cex=1.8, font=1)
#axis(2, ylim=c(-50, -30), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 
axis(1, xlim=c(0, 40), at=c(0, 10, 20, 30, 40), 
     labels=c(expression(0), expression(10), expression(20), expression(30), expression(40)), col="black", las=1, cex=2, cex.axis=2, cex.lab=2, font=1)

levels(factor(Falling$ENVIRON.x)) # Open pch=21 # Edge pch=22 # Floodplain pch=25
points(jitter(Falling$avg.MPROD_cbrt, factor=1), Falling$avg.D13C_CO2, bg=adjustcolor("blue3", alpha.f=0.8), col="black", pch=c(22, 21)[as.numeric(as.factor(High$ENVIRON.x))], cex=3)

#legend("bottomright", legend=c(expression(italic(R^2) ~ "=" ~ "0.34"),
#expression(italic(p) ~ "<" ~ "0.001"),
#expression(italic(df) ~ "=" ~ "28")),
#text.col="blue3", text.font=3, cex=1.8, bty="n")

par(xpd=TRUE)
legend("topleft", inset=-0.2, legend=expression("f)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")



