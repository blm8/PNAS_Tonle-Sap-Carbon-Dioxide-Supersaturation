

##===============================  CO2 Mass Balance


##  Author(s):  B. Miller

##  Required dataframes:  metab.csv, Master_DIC-Corrected.csv, Master_TSL-Corrected.csv

##  Description:  Runs CO2 mass balance

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


# GPP and ER

setwd("~/Desktop/Diel TSL O2")
metab <- read.csv("metab.csv", header=T)
names(metab)
#View(metab)

GPP_ER.1 <- metab %>%
  dplyr::select(SITE, DATE, STAGE, DESCRIP, GPP, R) %>% # Select columns of interest 
  mutate(DATE = as.Date(DATE)) %>% # Add Julian DOY column
  mutate(ENVIRON = DESCRIP) %>%
  mutate(mb_GPP = (GPP/1000)/32*1000*1000) %>% # mmol CO2 m^-3 d^-1
  mutate(mb_ER = -1*(R/1000)/32*1000*1000) %>% # mmol CO2 m^-3 d^-1
  dplyr::select(SITE, DATE, STAGE, ENVIRON, mb_GPP, mb_ER)
#View(GPP_ER.1)

GPP_ER.2 <- GPP_ER.1 %>%
  dplyr::select(SITE, DATE, STAGE, ENVIRON, mb_GPP, mb_ER) %>%
  group_by(ENVIRON, STAGE) %>% # Group data by these columns
  mutate(avg_GPP = mean(mb_GPP)) %>%
  mutate(sd_GPP = sd(mb_GPP)) %>%
  mutate(n_GPP = length(mb_GPP)) %>%
  mutate(avg_ER = mean(mb_ER)) %>%
  mutate(sd_ER = sd(mb_ER)) %>%
  mutate(n_ER = length(mb_ER)) %>%
  dplyr::select(ENVIRON, STAGE, avg_GPP, sd_GPP, n_GPP, avg_ER, sd_ER, n_ER) %>%
  distinct(avg_ER, .keep_all=TRUE) 
GPP_ER.2 <- GPP_ER.2[!GPP_ER.2$STAGE == "Rising", ] 
GPP_ER.2 <- GPP_ER.2[-5, ]
GPP_ER.2$ENVIRON <- as.character(GPP_ER.2$ENVIRON)
GPP_ER.2$STAGE <- as.character(GPP_ER.2$STAGE)
#View(GPP_ER.2)

# DIC

setwd("~/Desktop/Working")
cGas <- read_csv("Master_DIC-Corrected.csv")
names(cGas)
#View(cGas)
Table(na.omit(cGas$D13C_DIC))

DIC.1 <- cGas %>%
  dplyr::select(SITE, DATE, ENVIRON, CLASS_Z, STAGE, DIC, PCO2, KH_PCO2, FP_Z) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  mutate(DIC_conv = (DIC/1000000)*KH_PCO2*1000*1000*(FP_Z*0.5)) %>% # mmol DIC m^-2
  mutate(PCO2_conv = (PCO2/1000000)*KH_PCO2*1000*1000*(FP_Z*0.5)) %>% # mmol CO2 m^-2
  mutate(DIC_comb = DIC_conv + PCO2_conv) %>% # mmol DIC m^-3
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(mb_DIC = mean(DIC_comb, na.rm=TRUE)) %>% # Apply the function mean() to every column
  distinct(mb_DIC, .keep_all=TRUE) %>%
  dplyr::select(SITE, DATE, ENVIRON, CLASS_Z, STAGE, mb_DIC, FP_Z)
DIC.1 <- DIC[!DIC.1$STAGE == "Rising", ] 
DIC.1 <- DIC[!DIC.1$ENVIRON == "Tributary", ] 
DIC.1 <- DIC[complete.cases(DIC.1), ]
#View(DIC.1)

DIC.2 <- DIC.1 %>%
  dplyr::select(SITE, DATE, ENVIRON, CLASS_Z, STAGE, mb_DIC, FP_Z) %>%
  group_by(ENVIRON, STAGE) %>% # Group data by these columns
  mutate(avg_DIC = mean(mb_DIC)) %>%
  mutate(sd_DIC = sd(mb_DIC)) %>%
  mutate(n_DIC = length(mb_DIC)) %>%
  dplyr::select(ENVIRON, STAGE, avg_DIC, sd_DIC, n_DIC, FP_Z) %>%
  distinct(avg_DIC, .keep_all=TRUE) 
#View(DIC.2)

TERMS <- left_join(DIC.2, GPP_ER.2, by=c("STAGE", "ENVIRON")) %>%
  dplyr::select(STAGE, ENVIRON, avg_DIC, sd_DIC, n_DIC, FP_Z, avg_GPP, sd_GPP, n_GPP, avg_ER, sd_ER, n_ER) %>%
  mutate(avg_GPP = avg_GPP*(FP_Z*0.5)) %>%
  mutate(sd_GPP = sd_GPP*(FP_Z*0.5)) %>%
  mutate(avg_ER = avg_ER*(FP_Z*0.5)) %>%
  mutate(sd_ER = sd_ER*(FP_Z*0.5)) 
#View(TERMS)

# MOX and CO2

setwd("~/Desktop/Working")
cGas <- read_csv("Master_TSL-Corrected.csv")
#View(cGas)

MOX_CO2.1 <- cGas %>%
  dplyr::select(SITE, DATE, ENVIRON, CLASS_Z, STAGE, PCO2, KH_PCO2, MOX, FP_Z) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  mutate(CO2_conv = (PCO2/1000000)*KH_PCO2*1000*1000) %>% # mmol CO2 m^-3
  mutate(MOX_conv = -1*(MOX/1000000)*24*1000) %>% # mmol CO2 m^-3 d^-1
  mutate(MOX_Barbosa = -1*(MOX/1000000)*24*1000*16.04*(12.01/16.04)) %>% # mg C-CH4 m^-3 d^-1
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(mb_MOX = mean(MOX_conv, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(mb_CO2 = mean(CO2_conv, na.rm=TRUE)) %>% # Apply the function mean() to every column
  distinct(mb_CO2, .keep_all=TRUE) %>%
  dplyr::select(SITE, DATE, ENVIRON, CLASS_Z, STAGE, mb_MOX, mb_CO2, FP_Z, MOX_Barbosa)
MOX_CO2.1 <- MOX_CO2.1[complete.cases(MOX_CO2.1), ]
#View(MOX_CO2.1)

range(MOX_CO2.1$MOX_Barbosa)

Surface <- subset(MOX_CO2.1, CLASS_Z=="Surface")
Bottom <- subset(MOX_CO2.1, CLASS_Z=="Bottom")

MOX_CO2.2 <- left_join(Surface, Bottom, by=c("SITE", "DATE", "ENVIRON")) %>% 
  replace(is.na(.), 0)
#View(MOX_CO2.2)

MOX_CO2.3 <- MOX_CO2.2 %>%
  dplyr::select(SITE, DATE, ENVIRON, STAGE.x, mb_MOX.x, mb_MOX.y, mb_CO2.x, mb_CO2.y, FP_Z.x) %>% # Select columns of interest 
  mutate(STAGE = STAGE.x) %>% 
  mutate(CLASS_Z = CLASS_Z.x) %>% 
  mutate(FP_Z = FP_Z.x) %>%
  mutate(mb_MOX = (mb_MOX.x*(FP_Z*0.5))+(mb_MOX.y*(FP_Z*0.5))) %>% # MOX in mmol m^-2 d^-1
  mutate(mb_CO2 = (mb_CO2.x*(FP_Z*0.5))+(mb_CO2.y*(FP_Z*0.5))) # MOX in mmol m^-2 d^-1
#View(MOX_CO2.3)


MOX_CO2.4 <- MOX_CO2.3 %>%
  dplyr::select(SITE, DATE, STAGE, ENVIRON, DATE, mb_MOX, mb_CO2) %>%
  group_by(STAGE, ENVIRON) %>%
  mutate(avg_MOX = mean(mb_MOX)) %>%
  mutate(sd_MOX = sd(mb_MOX)) %>%
  mutate(n_MOX = length(mb_MOX)) %>%
  mutate(avg_CO2 = mean(mb_CO2)) %>%
  mutate(sd_CO2 = sd(mb_CO2)) %>%
  mutate(n_CO2 = length(mb_CO2)) %>%
  distinct(avg_MOX, .keep_all=TRUE) 
#View(MOX_CO2.4) # Merge this df to TERMS

# mb_diffCO2 and mb_sdCO2

setwd("~/Desktop/Working")
diff.1 <- read_csv("Master_TSL-Corrected.csv") %>%
  dplyr::select(SITE, DATE, ENVIRON, CLASS_Z, STAGE, diffCO2, sdCO2) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  mutate(diffCO2_conv = (diffCO2/1000)/44.01*1000) %>% # diffCO2 in mmol m^-2 d^-1
  mutate(sdCO2_conv = (sdCO2/1000)/44.01*1000) %>% # diffCO2 in mmol m^-2 d^-1
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(mb_diffCO2 = mean(diffCO2_conv, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(mb_sdCO2 = mean(sdCO2_conv, na.rm=TRUE)) %>% # Apply the function mean() to every column
  distinct(mb_diffCO2, .keep_all=TRUE) %>%
  dplyr::select(SITE, DATE, ENVIRON, STAGE, CLASS_Z, mb_diffCO2, mb_sdCO2)
diff.1 <- diff.1[!diff.1$SITE == "KgPr", ]
diff.1 <- diff.1[!diff.1$SITE == "StSe", ]
diff.1 <- diff.1[!diff.1$SITE == "PtSd.O", ] 
diff.1 <- diff.1[!diff.1$SITE == "PtSd.E", ] 
diff.1 <- diff.1[!diff.1$SITE == "PtSd.F", ] 
diff.1 <- diff.1[!diff.1$SITE == "KgKl.O", ] 
diff.1 <- diff.1[!diff.1$SITE == "KgKl.E", ] 
diff.1 <- diff.1[!diff.1$STAGE == "Rising", ] 
diff.1 <- diff.1[!diff.1$CLASS_Z == "Bottom", ]
diff.1 <- diff.1[complete.cases(diff.1), ]
#View(diff.1) 

diff.2 <- diff.1 %>%
  dplyr::select(SITE, DATE, ENVIRON, STAGE, CLASS_Z, mb_diffCO2, mb_sdCO2) %>%
  group_by(STAGE, ENVIRON) %>%
  mutate(avg_diff = mean(mb_diffCO2)) %>%
  mutate(sd_diff = mean(mb_sdCO2)) %>%
  mutate(n_diff = length(mb_diffCO2)) %>%
  distinct(avg_diff, .keep_all=TRUE) 
diff.2 <- diff.2[-6, ]
#View(diff.2) # Merge this df to TERMS

# mb_MPROD

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

MPROD.1 <- subset(newData, CLASS_Z=="Sediment") %>%
  dplyr::select(SITE, DATE, STAGE, MPROD, ENVIRON) %>% # Select columns of interest 
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE) %>% # Group data by these columns
  mutate(mb_MPROD = (mean(MPROD, na.rm=TRUE)*24*1000000*0.1)/1000000) %>% # Measured rates valid over 10 cm, in mmol m-2 d-1
  distinct(mb_MPROD, .keep_all=TRUE) 
MPROD.1 <- MPROD.1[complete.cases(MPROD.1), ]
View(MPROD.1)

MPROD.2 <- MPROD.1 %>%
  dplyr::select(SITE, DATE, ENVIRON, STAGE, mb_MPROD) %>%
  group_by(STAGE, ENVIRON) %>%
  mutate(avg_MPROD = mean(mb_MPROD)) %>%
  mutate(sd_MPROD = sd(mb_MPROD)) %>%
  mutate(n_MPROD = length(mb_MPROD)) %>%
  distinct(avg_MPROD, .keep_all=TRUE) 
#View(MPROD.2) # Merge this df to TERMS

# Merged df

df.1 <- left_join(TERMS, MOX_CO2.4, by=c("STAGE", "ENVIRON"))
#View(df.1)
df.2 <- left_join(df.1, diff.2, by=c("STAGE", "ENVIRON"))
#View(df.2)
df.3 <- left_join(df.2, MPROD.2, by=c("STAGE", "ENVIRON"))
#View(df.3)
names(df.3)
df.4 <- df.3 %>%
  dplyr::select(STAGE, ENVIRON, avg_CO2, sd_CO2, n_CO2, avg_ER, sd_ER, n_ER, avg_GPP, sd_GPP, n_GPP,
                avg_MPROD, sd_MPROD, n_MPROD, avg_MOX, sd_MOX, n_MOX, avg_DIC, sd_DIC, n_DIC, avg_diff, sd_diff, n_diff) # Select columns of interest 
#View(df.4) 

# Monte Carlo simulations to solve for ER.exSitu

sim <- subset(df.4, STAGE=="High")
sim <- subset(sim, ENVIRON=="Pelagic")

result <- vector("numeric")

for (i in 1:10000) {
  CO2 <- rnorm(n=sim$n_CO2, mean=sim$avg_CO2, sd=sim$sd_CO2)
  ER.inSitu <- 210.5282
  GPP <- 50.33331
  #ER.inSitu <- rnorm(n=sim$n_ER, mean=sim$avg_ER, sd=sim$sd_ER)
  #GPP <- rnorm(n=sim$n_GPP, mean=sim$avg_GPP, sd=sim$sd_GPP)
  MPROD <- rnorm(n=sim$n_MPROD, mean=sim$avg_MPROD, sd=sim$sd_MPROD)
  MOX <- rnorm(n=sim$n_MOX, mean=sim$avg_MOX, sd=sim$sd_MOX)
  #DIC <- rnorm(n=sim$n_DIC, mean=sim$avg_DIC, sd=sim$sd_DIC)
  #diff <- rnorm(n=sim$n_diff, mean=sim$avg_DIC, sd=sim$sd_DIC)
  ER.exSitu = CO2 - ER.inSitu + GPP - MPROD - MOX #- DIC + diff
  result[i] <- ER.exSitu
}

length(result)
mean(na.omit(result))
sd(na.omit(result)) # Enter into Excel spreadsheet

#View(sim)


# Error propagation for Monte Carlo simulations 

a <- 430

b <- 0.03

d_a <- 20

d_b <- 0.02

perc <- b/a
perc*100
error <- (sqrt(((d_a/a)^2) + ((d_b/b)^2)))*perc
error*100



