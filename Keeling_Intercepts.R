

##===============================  Keeling Intercepts


##  Author(s):  B. Miller

##  Required dataframes:  Master_TSL-Corrected.csv

##  Description:  Calculates Keeling intercepts from 1/PCO2 and d13C_CO2

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

cGas <- newData %>%
  dplyr::select(SITE, CLASS_Z, DATE, PCO2, D13C_CO2, ENVIRON, STAGE) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCO2 = mean(PCO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(inv.PCO2 = 1/avg.PCO2) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CO2 = mean(D13C_CO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  distinct(avg.PCO2, .keep_all=TRUE) %>%
  dplyr::select(SITE, CLASS_Z, DATE, avg.PCO2, inv.PCO2, avg.D13C_CO2, ENVIRON, STAGE) 
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

plot(High$inv.PCO2, High$avg.D13C_CO2)  
summary(lm(High$avg.D13C_CO2~High$inv.PCO2))
abline(lm(High$avg.D13C_CO2~High$inv.PCO2)) # Intercept:  -41.341 # Multiple R-squared:  0.7192 # p-value:  0.03287 # df:  5

plot(Open$inv.PCO2, Open$avg.D13C_CO2)  
summary(lm(Open$avg.D13C_CO2~Open$inv.PCO2))
abline(lm(Open$avg.D13C_CO2~Open$inv.PCO2)) # Intercept:  -50.940 # Multiple R-squared:  0.7192 # p-value:  0.03287 # df:  5

Stage <- subset(cGas, STAGE=="High") 

Open <- subset(Stage, ENVIRON=="Pelagic")
Edge <- subset(Stage, ENVIRON=="Edge")
Floodplain <- subset(Stage, ENVIRON=="Floodplain")

plot(Open$inv.PCO2, Open$avg.D13C_CO2)  
summary(lm(Open$avg.D13C_CO2~Open$inv.PCO2))
abline(lm(Open$avg.D13C_CO2~Open$inv.PCO2)) # Intercept:  -50.940 # Multiple R-squared:  0.7192 # p-value:  0.03287 # df:  5

#View(Edge)
Edge <- Edge[-2, ]
plot(Edge$inv.PCO2, Edge$avg.D13C_CO2)  
summary(lm(Edge$avg.D13C_CO2~Edge$inv.PCO2))
abline(lm(Edge$avg.D13C_CO2~Edge$inv.PCO2)) # Intercept:  -43.579 # Multiple R-squared:  0.8593 # p-value:  0.02342 # df:  4

plot(Floodplain$inv.PCO2, Floodplain$avg.D13C_CO2)  
summary(lm(Floodplain$avg.D13C_CO2~Floodplain$inv.PCO2))
abline(lm(Floodplain$avg.D13C_CO2~Floodplain$inv.PCO2)) # Intercept:  -40.810 # Multiple R-squared:  0.002994 # p-value:  0.8042 # df:  21  

Stage <- subset(cGas, STAGE=="Falling") 

Open <- subset(Stage, ENVIRON=="Pelagic")
Edge <- subset(Stage, ENVIRON=="Edge")
Floodplain <- subset(Stage, ENVIRON=="Floodplain")

plot(Open$inv.PCO2, Open$avg.D13C_CO2)  
summary(lm(Open$avg.D13C_CO2~Open$inv.PCO2))
abline(lm(Open$avg.D13C_CO2~Open$inv.PCO2)) # Intercept:  -43.145 # Multiple R-squared:  0.9058 # p-value:  0.00344 # df:  4

#View(Edge)
#Edge <- Edge[-2, ]
plot(Edge$inv.PCO2, Edge$avg.D13C_CO2)  
summary(lm(Edge$avg.D13C_CO2~Edge$inv.PCO2))
abline(lm(Edge$avg.D13C_CO2~Edge$inv.PCO2)) # Intercept:  -1.112e+01 # Multiple R-squared:  0.8492 # p-value:  0.008997 # df:  4 



