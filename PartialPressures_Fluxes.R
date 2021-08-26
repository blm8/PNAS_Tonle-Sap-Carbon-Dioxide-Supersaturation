

##===============================  Partial Pressures and Gas Fluxes 


##  Author(s):  B. Miller

##  Required dataframes:  Master_TSL.csv, metPNH.csv

##  Description:  Interpolates inundation times;
##                Calculates temperature-dependent schmidt numbers 
##                Organizes data for Monte Carlo simulations for k and flux calculations  
##                Runs Monte Carlo simulations for k & calculates fluxes


rm(list=ls())

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

setwd("~/Desktop/Working")


##===============================  Inundation Time Interpolation

data <- read_csv("Master_TSL.csv")
#names(data)
#View(data)
range(na.omit(data$DIST_FPL)) # -4859  5150

check <- data %>%
  dplyr::select(DIST_FPL, INUN_T) %>%
  na.omit(INUN_T)
#plot(check$DIST_FPL, check$INUN_T)
lm_y <- lm(INUN_T ~ DIST_FPL, data=data)
summary(lm_y) # Adjusted R-squared:  0.8814 # p-value: < 2.2e-16
poly_y <- lm(INUN_T ~ DIST_FPL + I(DIST_FPL^2) + I(DIST_FPL^3), data=data)
summary(poly_y) # Adjusted R-squared:  0.9461 # p-value: < 2.2e-16
AICc1 <- AICc(lm_y, k=2, REML=NULL)
AICc2 <- AICc(poly_y, k=2, REML=NULL)
delAICc1 <- AICc1 - min(c(AICc1, AICc2)) # Delta AICc
delAICc1 
delAICc2 <- AICc2 - min(c(AICc1, AICc2))
delAICc2 # This suggests that the cubic polynomial relationship best fits the data; use a cubic (or Hermite) spline function for interpolation

#?splinefun()
SplineFun <- splinefun(x=check$DIST_FPL, y=check$INUN_T, method="natural", ties=mean) # Create the spline function for interpolation
range(na.omit(data$DIST_FPL)) # -4859  5150
DIST_FPL <- seq(-4859, 5150, by = 1) # Range of x-values for interpolation
interpT <- SplineFun(DIST_FPL)
length(interpT)
toMerge <- as.data.frame(cbind(DIST_FPL, interpT))
length(toMerge$interpT)
#names(toMerge)
#View(toMerge) # This returns unrealistic values for the most negative INUN_T values; use linear interpolation

#?approxfun()
ApproxFun <- approxfun(x=check$DIST_FPL, y=check$INUN_T, method="linear", rule=2) # Create the linear interpolation function
DIST_FPL <- seq(-4859, 5150, by = 1)
length(DIST_FPL)
interpT <- ApproxFun(DIST_FPL)
length(interpT)
toMerge <- as.data.frame(cbind(DIST_FPL, interpT))
length(toMerge$interpT)
#names(toMerge)
#View(toMerge)

length(data$DATE) # Check
newData <- left_join(data, toMerge, by="DIST_FPL")
newData <- subset(newData, select=-INUN_T)
length(newData$DATE) # Check
#names(newData)
#View(newData)


##===============================  Schmidt Numbers

newData$Sc_CH4 <- 1897.8-(114.28*newData$WATER_T)+(3.2902*(newData$WATER_T^2))-(0.0339061*(newData$WATER_T^3)) # Calculate Schmidt number
range(na.omit(newData$Sc_CH4)) # Should be higher than Sc.CO2
newData$CH4_molL <- newData$KH_PCH4*(newData$PCH4/1000000)
newData$gradConcCH4 <- (newData$CH4_molL-(newData$KH_PCH4*(1.8/1000000)))*16.04 # g cm^-3
#names(newData)
#View(newData)

newData$Sc_CO2 <- 1911.1-(118.11*newData$WATER_T)+(3.4527*(newData$WATER_T^2))-(0.04132*(newData$WATER_T^3)) # Calculate Schmidt number
range(na.omit(newData$Sc_CO2))
newData$CO2_molL <- newData$KH_PCO2*(newData$PCO2/1000000)
newData$gradConcCO2 <- (newData$CO2_molL-(newData$KH_PCO2*(400/1000000)))*44.01 # g L-1
#names(newData)
#View(newData)


##===============================  Monte Carlo Simulations

metPNH <- read_csv("metPNH.csv")

#names(metPNH)
#View(metPNH)
metPNH.mcDraws <- metPNH %>% # Calculate daily mean and standard deviation for U10 windspeed
  dplyr::select(DATE, wind_speed_ms) %>% # Select columns of interest
  #mutate(logwind_speed_ms = log10(wind_speed_ms)) %>%
  group_by(DATE) %>% # Group data by these columns
  summarise(avg_wind_speed_ms = mean(wind_speed_ms, na.rm=TRUE), 
            sd_wind_speed_ms = sd(wind_speed_ms)) %>%
  ungroup()
#View(metPNH.mcDraws) 

#names(newData)
#View(newData)
flux.mcDraws <- newData[newData$CLASS_Z == "Surface", ]
flux.mcDraws <- flux.mcDraws[!flux.mcDraws$SITE == "StSe", ] # Omit Stueng Saen tributary site
flux.mcDraws <- flux.mcDraws[!flux.mcDraws$SITE == "KgPr", ] # Omit Kampong Prak tributary site

flux.mcDraws <- flux.mcDraws %>% 
  dplyr::select(SITE, DATE, Sc_CH4, Sc_CO2)  %>%
  group_by(SITE, DATE) %>%
  summarise(ScCH4 = mean(Sc_CH4, na.rm=TRUE), ScCO2 = mean(Sc_CO2, na.rm=TRUE)) 
flux.mcDraws <- flux.mcDraws[complete.cases(flux.mcDraws), ]
#names(flux.mcDraws)
#View(flux.mcDraws)

mcDraws <- left_join(flux.mcDraws, metPNH.mcDraws, by="DATE") 
class(flux.mcDraws$DATE)
class(metPNH.mcDraws$DATE)
#names(mcDraws)
#View(mcDraws)

##===========================================================================================##
##===========================================================================================##

sim <- mcDraws[62,] # Row of the dataframe containing mean and sd for parameters of interest
result <- vector("numeric")

for (i in 1:10000) {
  U10 <- rnorm(sim$avg_wind_speed_ms, sim$sd_wind_speed_ms)
  k600 <- 2.07 + (0.215*(U10^1.7))
  k <- k600*((sim$ScCH4/600)^-0.66) # Cole and Caracao (1998) model # cm h^-1
  result[i] <- k
}

length(result)
mean(na.omit(result))
sd(na.omit(result)) # Enter into Excel spreadsheet

##===========================================================================================##
##===========================================================================================##

mcDraws.CH4 <- read_csv("mcDraws.CH4.csv")
#names(mcDraws.CH4)
#View(mcDraws.CH4)

mcDraws.CH4 <- cbind.data.frame(mcDraws, mcDraws.CH4)
newData <- left_join(newData, mcDraws.CH4, by=c("DATE", "SITE"))
#names(newData)
#View(newData)

newData$diffCH4 <- (newData$gradConcCH4*newData$avg_mcCH4_k)*100*100*24*(12.01/16.04) # mg C-CH4 m^-2 d^-1
newData$sdCH4 <- (newData$gradConcCH4*newData$sd_mcCH4_k)*100*100*24*(12.01/16.04) # mg C-CH4 m^-2 d^-1
#names(newData)
#View(newData)

##===========================================================================================##
##===========================================================================================##

sim <- mcDraws[62,] # Row of the dataframe containing mean and sd for parameters of interest
result <- vector("numeric")

for (i in 1:10000) {
  U10 <- rnorm(sim$avg_wind_speed_ms, sim$sd_wind_speed_ms)
  k600 <- 2.07 + (0.215*(U10^1.7))
  k <- k600*((sim$ScCO2/600)^-0.66) # Cole and Caracao (1998) model # cm h^-1
  result[i] <- k600
}

length(result)
mean(na.omit(result))
sd(na.omit(result)) # Enter into Excel spreadsheet

##===========================================================================================##
##===========================================================================================##

mcDraws.CO2 <- read_csv("mcDraws.CO2.csv")
#names(mcDraws.CO2)
#View(mcDraws.CO2)

mcDraws.CO2 <- cbind.data.frame(mcDraws, mcDraws.CO2)
mcDraws.CO2 <- subset(mcDraws.CO2, select=c("SITE", "DATE", "avg_mcCO2_k", "sd_mcCO2_k", "avg_k600", "sd_k600"))
newData <- left_join(newData, mcDraws.CO2, by=c("DATE", "SITE"))

newData$diffCO2 <- (newData$gradConcCO2*newData$avg_mcCO2_k)*100*100*24*(12.01/16.04) # mg C-CO2 m^-2 d^-1
newData$sdCO2 <- (newData$gradConcCO2*newData$sd_mcCO2_k)*100*100*24*(12.01/16.04) # mg C-CO2 m^-2 d^-1

names(newData)
#View(newData)

write.csv(newData, file="Master_TSL-Corrected.csv")


##===============================  Table 1

Table2 <- function(x) { # Calculate means and standard errors
  n <- length(na.omit(x))
  Mean <- mean(na.omit(x))
  SE <- sd(na.omit(x))/sqrt(length(na.omit(x))) 
  Results <- list(n, Mean, SE)
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

names(newData)
View(newData)

cGas <- newData %>%
  dplyr::select(SITE, CLASS_Z, DATE, PCO2, PCH4, diffCO2, diffCH4, TRANSECT_POINT, LOCATION, ENVIRON, STAGE, avg_k600, sd_k600, CO2_molL, CH4_molL) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCO2 = mean(PCO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.PCH4 = mean(PCH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.CO2_molL = mean(CO2_molL, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.CH4_molL = mean(CH4_molL, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.diffCO2 = mean(diffCO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.diffCH4 = mean(diffCH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  distinct(avg.PCO2, .keep_all=TRUE) %>%
  ungroup() 
cGas <- cGas[!cGas$SITE == "KgPr", ]
cGas <- cGas[!cGas$SITE == "StSe", ]
#cGas <- fp[complete.cases(cGas), ]
View(cGas)


Table(cGas$avg.PCO2)

High <- subset(cGas, STAGE=="High")
Open <- subset(High, ENVIRON=="Pelagic")
Edge <- subset(High, ENVIRON=="Edge")
Floodplain <- subset(High, ENVIRON=="Floodplain")

Table(Floodplain$avg.PCH4)

Falling <- subset(cGas, STAGE=="Falling")
Open <- subset(Falling, ENVIRON=="Pelagic")
Edge <- subset(Falling, ENVIRON=="Edge")
Floodplain <- subset(Falling, ENVIRON=="Floodplain")

Table(Edge$avg.PCH4)



