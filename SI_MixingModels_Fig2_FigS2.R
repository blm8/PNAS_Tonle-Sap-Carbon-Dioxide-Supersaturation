

##===============================  Stable Isotope Mixing Models 

##  Author:  B. Miller

##  Required:  Master_TSL-Corrected.csv, fractions_0.050.csv

##  Description:  Determines all possible fractional contributions by sources of organic and inorganic C to dissolved carbon dioxide given either a normal distribution or continuous uniform distribution of values for each source;
##                Functions as IsoSource by U.S. EPA and later SIAR, a precursor to MixSIAR

rm(list=ls())

#install.packages("plyr")
require(plyr)
#install.packages("readr")
require(readr)
#install.packages("dplyr")
require(dplyr)
#install.packages("tidyr")
require(tidyr)
#install.packages("gtools")
require(gtools)
#install.packages("Rfast")
require(Rfast)


##===============================  Dataframe of Possible Fractional Contributions by Sources

# Create a dataframe with a single column of fractions to be applied across n end members in the stable isotope mixing model
# Increase recursion limit in R
options(expressions=500000)
options("expressions")==500000
# Minimum fraction for each end member is either 0 % or 2.5 % # Memory exceeded when minimum is decreased to 1 %
#df <- as.data.frame(seq(from=0.00, to=1.00, by=0.01)) 
#df <- as.data.frame(seq(from=0.00, to=1.00, by=0.025)) 
df <- as.data.frame(seq(from=0.00, to=1.00, by=0.05)) 
colnames(df) <- "fractions" 
df

# List all possible combinations of these n numbers, allowing numbers to repeat in combinations
#?combinations() 
combos <- combinations(length(df$fractions), 8, df$fractions, repeats.allowed = TRUE)
unique(combos[,1]) # Check that all expected values appear in a given column
subset_vctr <- which(Rfast::rowsums(combos[,1:8]) == 1) # Subset rows that sum to 1.00
combos2 <- combos[subset_vctr, ]
head(combos2)
tail(combos2)
result <- as.data.frame(combos2)
colnames(result) <- c("1", "2", "3", "4", "5", "6", "7", "8")

# Re-order to allow each end member each combination
col_order <- c("8", "1", "2", "3", "4", "5", "6", "7") 
result2 <- result[, col_order]
col_order <- c("7", "8", "1", "2", "3", "4", "5", "6")
result3 <- result[, col_order]
col_order <- c("6", "7", "8", "1", "2", "3", "4", "5")
result4 <- result[, col_order]
col_order <- c("5", "6", "7", "8", "1", "2", "3", "4")
result5 <- result[, col_order]
col_order <- c("4", "5", "6", "7", "8", "1", "2", "3")
result6 <- result[, col_order]
col_order <- c("3", "4", "5", "6", "7", "8", "1", "2")
result7 <- result[, col_order]
col_order <- c("2", "3", "4", "5", "6", "7", "8", "1")
result8 <- result[, col_order]
colnames(result) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
colnames(result2) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
colnames(result3) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
colnames(result4) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
colnames(result5) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
colnames(result6) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
colnames(result7) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
colnames(result8) <- c("f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8")
fractions <- rbind(result, result2, result3, result4, result5, result6, result7, result8)
fractions <- distinct(fractions, .keep_all = FALSE)
unique(fractions[,1]) # Check that all expected values appear in a given column
View(fractions)

setwd("~/Desktop/Working")
write.csv(fractions,'fractions_0.050.csv')


##===============================  Stable Isotope Mixing Models A and B

setwd("~/Desktop/Working")
fractions <- read_csv("fractions_0.050.csv")
fractions <- fractions[,-1]
#View(fractions)
newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 

cGas <- newData %>%
  dplyr::select(SITE, CLASS_Z, DATE, PCO2, PCH4, D13C_CO2, D13C_CH4, ENVIRON, STAGE) %>% # Select columns of interest 
  na.omit(CLASS_Z) %>% # Remove NAs
  mutate(date = as.Date(DATE)) %>% # Add Julian DOY column
  group_by(SITE, DATE, CLASS_Z) %>% # Group data by these columns
  mutate(avg.PCO2 = mean(PCO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CO2 = mean(D13C_CO2, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.PCH4 = mean(PCH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.D13C_CH4 = mean(D13C_CH4, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(eC = avg.D13C_CO2-avg.D13C_CH4, na.rm=TRUE) %>% # Apply the function mean() to every column
  distinct(avg.PCO2, .keep_all=TRUE) %>%
  dplyr::select(SITE, CLASS_Z, DATE, avg.PCO2, avg.D13C_CO2, avg.PCH4, avg.D13C_CH4, ENVIRON, STAGE, eC) 
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

unique(cGas$STAGE)
stage <- subset(cGas, STAGE=="High") 
unique(cGas$ENVIRON)
environ <- subset(stage, ENVIRON=="Pelagic")

# High Pelagic f_GPP=0.027 # Falling Pelagic f_GPP=0.03
# High Edge f_GPP=0.07 # Falling Edge f_GPP=0.014
# High Floodplain f_GPP=0.11

# High Pelagic f_diff=0.147 # Falling Pelagic f_diff=0.192
# High Edge f_diff=0.167 # Falling Edge f_diff=0.744
# High Floodplain f_diff=0.288

hist(rnorm(n=length(fractions$f2), mean=-24, sd=4)) # Run repeatedly to see different continuous uniform distributions possible # Use Monte Carlo sensitivity analysis

hist(rnorm(n=length(fractions$f3), mean=-28, sd=4)) # Run repeatedly to see different continuous uniform distributions possible # Use Monte Carlo sensitivity analysis

for (i in 1:dim(fractions)[1]) { # v.2 using n, mean, and sd in a normal distribution for each end member
  atm <- rnorm(n=length(fractions$f1), mean=-8.1, sd=0.4)
  f_atm <- fractions$f1
  #phyto <- -24
  phyto <- rnorm(n=length(fractions$f2), mean=-24, sd=4)
  f_phyto <- fractions$f2
  #peri <- -28
  peri <- rnorm(n=length(fractions$f3), mean=-28, sd=4)
  f_peri <- fractions$f3
  #c3plant <- -29
  c3plant <- rnorm(n=length(fractions$f4), mean=-29, sd=2)
  f_c3plant <- fractions$f4
  #macro <- -27.6
  macro <- rnorm(n=length(fractions$f5), mean=-27.6, sd=0.9)
  f_macro <- fractions$f5
  #c4plant <- -12.2
  c4plant <- rnorm(n=length(fractions$f6), mean=-12.2, sd=0.3)
  f_c4plant <- fractions$f6
  CH4 <- runif(length(fractions$f8), min=-110, max=-50) # Re-populate from either acetate fermentation, or both acetate fermentation and carbonate reduction
  f_CH4 <- fractions$f8
  CO2_obs <- mean(environ$avg.D13C_CO2) 
  f_GPP <- 0.027 # Re-populate from CO2 mass balances
  eps_GPP <- -19
  f_diff <- 0.147 # Re-populate from CO2 mass balances
  eps_diff <- -1.1
  CO2_calc = (atm*f_atm) + (phyto*f_phyto) + (peri*f_peri) + (c3plant*f_c3plant) + (macro*f_macro) + (c4plant*f_c4plant) + (CH4*f_CH4) - ((CO2_obs+eps_GPP)*f_GPP) - ((CO2_obs+eps_diff)*f_diff)
  outcome1 <- data.frame(CO2_calc, f_atm, f_phyto, f_peri, f_c3plant, f_macro, f_c4plant, f_CH4)
  criteria1 <- mean(environ$avg.D13C_CO2 + 1) # Use +/- 1 per mill, se, or sd
  criteria2 <- mean(environ$avg.D13C_CO2 - 1)
  outcome2 <- outcome1[!(outcome1$CO2_calc>criteria1),]
  outcome3 <- outcome2[!(outcome2$CO2_calc<criteria2),]
  outcome4 <- mean(outcome3$f_CH4)
  return(outcome4)
}

outcome4 # Run repeatedly to see different continuous uniform distributions possible # Use Monte Carlo sensitivity analysis

unique(cGas$STAGE)
stage <- subset(cGas, STAGE=="Falling") 
unique(cGas$ENVIRON)
environ <- subset(stage, ENVIRON=="Edge")

# High Pelagic f_GPP=0.027 # Falling Pelagic f_GPP=0.03
# High Edge f_GPP=0.07 # Falling Edge f_GPP=0.014
# High Floodplain f_GPP=0.11

# High Pelagic f_diff=0.147 # Falling Pelagic f_diff=0.192
# High Edge f_diff=0.167 # Falling Edge f_diff=0.744
# High Floodplain f_diff=0.288

x <- replicate(10000, { # Using Monte Carlo sensitivity analysis
  for (i in 1:dim(fractions)[1]) { # v.2 using n, mean, and sd in a normal distribution for each end member
    atm <- rnorm(n=length(fractions$f1), mean=-8.1, sd=0.4)
    f_atm <- fractions$f1
    #phyto <- -24
    phyto <- rnorm(n=length(fractions$f2), mean=-24, sd=4)
    f_phyto <- fractions$f2
    #peri <- -28
    peri <- rnorm(n=length(fractions$f3), mean=-28, sd=4)
    f_peri <- fractions$f3
    #c3plant <- -29
    c3plant <- rnorm(n=length(fractions$f4), mean=-29, sd=2)
    f_c3plant <- fractions$f4
    #macro <- -27.6
    macro <- rnorm(n=length(fractions$f5), mean=-27.6, sd=0.9)
    f_macro <- fractions$f5
    #c4plant <- -12.2
    c4plant <- rnorm(n=length(fractions$f6), mean=-12.2, sd=0.3)
    f_c4plant <- fractions$f6
    CH4 <- runif(length(fractions$f8), min=-110, max=-50) # Re-populate from either acetate fermentation, or both acetate fermentation and carbonate reduction
    f_CH4 <- fractions$f8
    CO2_obs <- mean(environ$avg.D13C_CO2) 
    f_GPP <- 0.014 # Re-populate from CO2 mass balances
    eps_GPP <- -19
    f_diff <- 0.744 # Re-populate from CO2 mass balances
    eps_diff <- -1.1
    CO2_calc = (atm*f_atm) + (phyto*f_phyto) + (peri*f_peri) + (c3plant*f_c3plant) + (macro*f_macro) + (c4plant*f_c4plant) + (CH4*f_CH4) - ((CO2_obs+eps_GPP)*f_GPP) - ((CO2_obs+eps_diff)*f_diff)
    outcome1 <- data.frame(CO2_calc, f_atm, f_phyto, f_peri, f_c3plant, f_macro, f_c4plant, f_CH4)
    criteria1 <- mean(environ$avg.D13C_CO2 + 1) # Use +/- 1 per mill, se, or sd
    criteria2 <- mean(environ$avg.D13C_CO2 - 1)
    outcome2 <- outcome1[!(outcome1$CO2_calc>criteria1),]
    outcome3 <- outcome2[!(outcome2$CO2_calc<criteria2),]
    outcome4 <- mean(outcome3$f_CH4)
    return(outcome4)
  }
})

length(x)
y <- na.omit(x)
length(y)
mean(y)
sd(y)


##===============================  Stable Isotope Mixing Models C and D

unique(cGas$STAGE)
stage <- subset(cGas, STAGE=="High") 
unique(cGas$ENVIRON)
environ <- subset(stage, ENVIRON=="Pelagic")

# High Pelagic f_GPP=0.027 # Falling Pelagic f_GPP=0.03
# High Edge f_GPP=0.07 # Falling Edge f_GPP=0.014
# High Floodplain f_GPP=0.11

# High Pelagic f_diff=0.147 # Falling Pelagic f_diff=0.192
# High Edge f_diff=0.167 # Falling Edge f_diff=0.744
# High Floodplain f_diff=0.288

hist(runif(length(fractions$f1), min=-33, max=-7.7)) # Run repeatedly to see different continuous uniform distributions possible # Use Monte Carlo sensitivity analysis

hist(runif(length(fractions$f2), min=-110, max=-50)) # Run repeatedly to see different continuous uniform distributions possible # Use Monte Carlo sensitivity analysis

for (i in 1:dim(fractions)[1]) { # v.3 using n, mean, and sd in a continuous uniform distribution encompassing all end members
  other <- runif(length(fractions$f1), min=-37, max=-7.7)
  f_other <- fractions$f1
  CH4 <- runif(length(fractions$f2), min=-110, max=-50) # Re-populate from either acetate fermentation, or both acetate fermentation and carbonate reduction
  f_CH4 <- fractions$f2
  CO2_obs <- mean(environ$avg.D13C_CO2) 
  f_GPP <- 0.027 # Re-populate from CO2 mass balances
  eps_GPP <- -19
  f_diff <- 0.147 # Re-populate from CO2 mass balances
  eps_diff <- -1.1
  CO2_calc = (other*f_other) + (CH4*f_CH4) - ((CO2_obs+eps_GPP)*f_GPP) - ((CO2_obs+eps_diff)*f_diff)
  outcome1 <- data.frame(CO2_calc, f_other, f_CH4)
  criteria1 <- mean(environ$avg.D13C_CO2 + 1) # Use +/- 1 per mill, se, or sd
  criteria2 <- mean(environ$avg.D13C_CO2 - 1)
  outcome2 <- outcome1[!(outcome1$CO2_calc>criteria1),]
  outcome3 <- outcome2[!(outcome2$CO2_calc<criteria2),]
  outcome4 <- outcome3[!(outcome3$f_other<0.15),] # Allow for all other end members to an input of at least 20%
  #outcome4 <- outcome3[which.max(outcome3$f_other),] # Allow maximum input for all other end members
  outcome5 <- mean(outcome4$f_CH4)
  return(outcome5)
}

outcome5 # Run repeatedly to see different continuous uniform distributions possible # Use Monte Carlo sensitivity analysis

unique(cGas$STAGE)
stage <- subset(cGas, STAGE=="High") 
unique(cGas$ENVIRON)
environ <- subset(stage, ENVIRON=="Pelagic")

# High Pelagic f_GPP=0.027 # Falling Pelagic f_GPP=0.03
# High Edge f_GPP=0.07 # Falling Edge f_GPP=0.014
# High Floodplain f_GPP=0.11

# High Pelagic f_diff=0.147 # Falling Pelagic f_diff=0.192
# High Edge f_diff=0.167 # Falling Edge f_diff=0.744
# High Floodplain f_diff=0.288

x <- replicate(10000, { # Using Monte Carlo sensitivity analysis
  for (i in 1:dim(fractions)[1]) { # v.3 using n, mean, and sd in a continuous uniform distribution encompassing all end members
  other <- runif(length(fractions$f1), min=-37, max=-7.7)
  f_other <- fractions$f1
  CH4 <- runif(length(fractions$f2), min=-110, max=-50) # Re-populate from either acetate fermentation, or both acetate fermentation and carbonate reduction
  f_CH4 <- fractions$f2
  CO2_obs <- mean(environ$avg.D13C_CO2) 
  f_GPP <- 0.027 # Re-populate from CO2 mass balances
  eps_GPP <- -19
  f_diff <- 0.147 # Re-populate from CO2 mass balances
  eps_diff <- -1.1
  CO2_calc = (other*f_other) + (CH4*f_CH4) - ((CO2_obs+eps_GPP)*f_GPP) - ((CO2_obs+eps_diff)*f_diff)
  outcome1 <- data.frame(CO2_calc, f_other, f_CH4)
  criteria1 <- mean(environ$avg.D13C_CO2 + 1) # Use +/- 1 per mill, se, or sd
  criteria2 <- mean(environ$avg.D13C_CO2 - 1)
  outcome2 <- outcome1[!(outcome1$CO2_calc>criteria1),]
  outcome3 <- outcome2[!(outcome2$CO2_calc<criteria2),]
  outcome4 <- outcome3[!(outcome3$f_other<0.15),] # Allow for all other end members to an input of at least 20%
  #outcome4 <- outcome3[which.max(outcome3$f_other),] # Allow maximum input for all other end members
  outcome5 <- mean(outcome4$f_CH4)
  return(outcome5)
  }
})

length(x)
y <- na.omit(x)
length(y)
mean(y)
sd(y)


##===============================  Figure 2a

atm <- -8.1
phyto <- c(-29.9851, -24.9366, -23.0619, -22.3132, -21.5092, -19.4513) # rnorm(n=6, mean=-24, sd=4)
mean(phyto) # -23.54288
peri <- c(-20.5708, -22.5733, -23.4248, -23.44, -25.701, -25.7054, -26.6462, -27.9043, -28.1134, -28.3895, -29.7284, -31.3265, -31.7701, -32.1064, -32.1889, -32.9124) #rnorm(n=18, mean=-28, sd=4)
mean(peri) # -27.65634
c3plant <-c(-31.4654, -30.77, -28.6339, -28.3252, -28.1488, -27.1574, -26.1541)  #rnorm(n=7, mean=-29, sd=2)
mean(c3plant) # -28.66497
macro <- c(-29.5429, -29.9328, -35.4164, -36.5418) #rnorm(n=3, mean=-27.6, sd=0.9)
mean(macro) # -32.85847
sd(macro)
c4plant <- rnorm(n=2, mean=-12.2, sd=0.3)
DIC <- rnorm(n=98, mean=-13.8, sd=0.4)
acetate <- runif(98, min=-65, max=-50)
carbonate <- runif(98, min=-110, max=10)
High <- subset(cGas, STAGE=="High")
high <- High$avg.D13C_CO2
mean(high)
Falling <- subset(cGas, STAGE=="Falling")
falling <- Falling$avg.D13C_CO2
mean(falling)

res.aov <- aov(phyto ~ peri)
# Summary of the analysis
summary(res.aov) 
TukeyHSD(res.aov) 

par(mfrow=c(3, 2))
par(oma=c(2, 2, 2, 2))
par(mar=c(0.25, 0.25, 0.25, 0.25)) # Bottom, left, top, right
boxplot(acetate, carbonate, high, falling, macro, c3plant, peri, phyto, DIC, c4plant, atm,
        names=c(expression(""), 
                expression(""), 
                expression("High-water" ~ C-CO[2]), 
                expression("Falling-water" ~ C-CO[2]),
                expression("C3" ~ "plants"),
                expression("Periphyton"),
                expression("Macrophytes"),
                expression("Phytoplankton"),
                expression("DIC"),
                expression("C4" ~ "plants"),
                expression("Atmospheric" ~ C-CO[2])),
        col=c("white",
              "white",
              adjustcolor("cornflowerblue", alpha.f=0.8), 
              adjustcolor("blue3", alpha.f=0.8),
              "white", 
              "white", 
              "white", 
              "white",
              "white",
              adjustcolor("gray50", alpha.f=0.2), 
              adjustcolor("gray50", alpha.f=0.2)),
        border=c("white",
                 "white",
                 "black", 
                 "black", 
                 "black", 
                 "black", 
                 "black", 
                 "black", 
                 "black",
                 "gray40", 
                 "gray40"),
        lwd=1.5, las=2, horizontal=T, range=0, xaxt="", yaxt="")

polygon(x=c(-8.1, -8.1, -37, -37), y=c(0.75, 1.5, 1.5, 0.75), col="white", border="black", lwd=1.5)
polygon(x=c(-50, -50, -110, -110), y=c(0.75, 1.5, 1.5, 0.75), col=adjustcolor("gray50", alpha.f=0.2), border="gray40", lwd=1.5)

polygon(x=c(-50, -50, -65, -65), y=c(1.75, 2.5, 2.5, 1.75), col=adjustcolor("gray50", alpha.f=0.2), border="gray40", lwd=1.5)
polygon(x=c(-60, -60, -110, -110), y=c(1.75, 2.5, 2.5, 1.75), col=adjustcolor("gray50", alpha.f=0.2), border="gray40", lwd=1.5)

axis(1, xlim=c(-115, 0), at=c(-5, -20, -35, -50, -65, -80, -95, -110), 
     labels=c("-5", "-20", "-35", "-50", "-65", "-80", "-95", "-110"), las=1, cex=2, cex.axis=2, cex.lab=2, font=1) 

mtext(expression(paste(delta^13*C ~ "(\u2030)")), side=1, line=4.5, cex=1.75, font=1)

text(2.5, 11, expression("Atm." ~ CO[2]), cex=2, font=3, col="gray40")
text(0, 10, expression("C4" ~ "Plants"), cex=2, font=3, col="gray40")
text(-7.5, 9, expression("DIC"), cex=2, font=3, col="black")
text(-11, 8, expression("Phyto."), cex=2, font=3, col="black")
text(-13, 7, expression("Peri."), cex=2, font=3, col="black")
text(-14, 6, expression("C3" ~ "Plants"), cex=2, font=3, col="black")
text(-15, 5, expression("Macrophytes"), cex=2, font=3, col="black")

text(-70, 4, expression("Falling-water" ~ CO[2]), cex=2, font=3, col="blue3")
text(-81, 3, expression("High-water" ~ CO[2]), cex=2, font=3, col="cornflowerblue")

text(-23.1, 1.1, expression(italic("Other")), cex=2, font=3, col="black")
text(-80, 1.1, expression(italic("Methane")), cex=2, font=3, col="gray40")

text(-58, 2.1, expression(italic("Acetate")), cex=2, font=3, col="gray40")
text(-83, 2.1, expression(italic("Carbonate")), cex=2, font=3, col="gray40")

text(-105, 11, expression("a)"), cex=2.75, font=3, col="black")

newData <- read_csv("Master_TSL-Corrected.csv")
names(newData)
#View(newData) 


##===============================  Figure 2b

cGas <- newData %>%
  dplyr::select(SITE, CLASS_Z, DATE, PCO2, PCH4, D13C_CO2, D13C_CH4, ENVIRON, STAGE, O2, MOX, ENVIRON, CO2_molL, CH4_molL) %>% # Select columns of interest 
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
  mutate(app = D13C_CO2 / D13C_CH4) %>% # Apply the function mean() to every column
  mutate(avg.CO2_molL = mean(CO2_molL, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.CH4_molL = mean(CH4_molL, na.rm=TRUE)) %>% # Apply the function mean() to every column
  mutate(avg.app = mean(app, na.rm=TRUE)) %>% # Apply the function mean() to every column
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
mean(cGas$avg.PCO2)
sd(cGas$avg.PCO2)

High <- subset(cGas, STAGE=="High")
mean(High$avg.PCO2)
sd(High$avg.PCO2)

mean(High$avg.CO2_molL)*1000000
sd(High$avg.CO2_molL)*1000000
Falling <- subset(cGas, STAGE=="Falling")
mean(Falling$avg.CO2_molL)*1000000
sd(Falling$avg.CO2_molL)*1000000

range(High$avg.app)
mean(High$avg.app)
sd(High$avg.app) 

#par(mfrow=c(3, 2))
#par(oma=c(4, 4, 4, 4))
par(mar=c(5, 10, 5, 10)) # Bottom, left, top, right
boxplot(cGas$avg.app~cGas$STAGE, ylim=c(0.4, 2.1), xaxt="n", xlab="", yaxt="n", ylab="", bty="n", varwidth=FALSE, border=c("black", "black"), col=c(adjustcolor("cornflowerblue", alpha.f=0.8), adjustcolor("blue3", alpha.f=0.8)), cex.axis=2, cex.lab=2, lwd=1.5, boxwex=c(0.6, 0.6, 0.6))
box(lty=1, lwd=1, col="black")
mtext(expression(paste(alpha)[app]), side=2, line=4.5, cex=2, font=1)
laby <- c(expression(0.5), expression(1), expression(1.5), expression(2))
axis(2, ylim=c(0.4, 2.1), at=c(0.5, 1, 1.5, 2), labels=laby, col="black", las=0, cex=2, cex.axis=2, cex.lab=2, font=1) 
par(xpd=TRUE)
#legend("topleft", inset=-0.2, legend=expression("b)"), text.col=c("black"), text.font=2, cex=2.75, bty="n")
text(0.3, 2.25, expression("b)"), cex=2.75, font=3, col="black")
par(xpd=FALSE)
abline(h=1.055, col="gray40", lwd=2, lty=2)

text(0.7, 0.95, expression(italic("Acetate")), cex=2, font=3, col="gray40")
text(0.7, 1.15, expression(italic("Carbonate")), cex=2, font=3, col="gray40")
text(1, 2, expression("High-water"), cex=2, font=3, col="cornflowerblue")
text(2, 2, expression("Falling-water"), cex=2, font=3, col="blue3")


##===============================  Figure S2

par(mfrow=c(1, 1))
par(oma=c(2, 2, 2, 2))
par(mar=c(5, 1, 4, 1)) # Bottom, left, top, right
xPhyto <- seq(-3, 3, length=1000)*4 + (-24)
yPhyto <- dnorm(xPhyto, mean=-24, sd=4)
plot(xPhyto, yPhyto, xlim=c(-130, 10), pch="", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
mtext(expression(paste(delta^13 ~ C ~ "(\u2030)")), side=1, line=4.5, cex=2.5, font=1)
axis(1, xlim=c(-120, 0), at=c(-120, -100, -80, -60, -40, -20, 0), 
     labels=c("-120", "-100", "-80", "-60", "-40", "-20", "0"), col="black", las=0, cex=2, cex.axis=2, cex.lab=2, font=1)
rect(-39, -0.17, -6.1, 0.10, border="black", lwd=2, col="white")
rect(-110, -0.17, -50, 0.10, border="gray40", lwd=2, col=adjustcolor("gray50", alpha.f=0.2))
par(new=TRUE)

plot(xPhyto, yPhyto, type="l", lwd=2, xlim=c(-130, 10),  xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="chartreuse4")
par(new=TRUE)

xPeri <- seq(-3, 3, length=1000)*4 + (-28)
yPeri <- dnorm(xPeri, mean=-28, sd=4)
plot(xPeri, yPeri, type="l", lwd=2, xlim=c(-130, 10), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="chartreuse3")
par(new=TRUE)

xC3plant <- seq(-3, -3, length=1000)*2 + (-29)
yC3plant <- dnorm(xC3plant, mean=-29, sd=2)
plot(xC3plant, yC3plant, type="l", lwd=2, xlim=c(-130, 10), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="chartreuse2")
par(new=TRUE)

xMacro <- seq(-3, 3, length=1000)*0.9 + (-27.6)
yMacro <- dnorm(xMacro, mean=-27.6, sd=0.9)
plot(xMacro, yMacro, type="l", lwd=2, xlim=c(-130, 10), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="chartreuse1")
par(new=TRUE)

xC4plant <- seq(-3, 3, length=1000)*0.3 + (-12.2)
yC4plant <- dnorm(xC4plant, mean=-12.2, sd=0.3)
plot(xC4plant, yC4plant, type="l", lwd=2, xlim=c(-130, 10), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="chartreuse")
par(new=TRUE)

xDIC <- seq(-3, 3, length=1000)*0.4 + (-13.8)
yDIC <- dnorm(xDIC, mean=-13.8, sd=0.4)
plot(xC4plant, yC4plant, type="l", lwd=2, xlim=c(-130, 10), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", col="cadetblue")

legend("right", inset=0, legend=c("Phyto.", "Peri.", "C3 Plants", "Macrophytes", "C4 Plants", "DIC"), 
       text.col=c("chartreuse4", "chartreuse3", "chartreuse2", "chartreuse1", "chartreuse", "cadetblue"), text.font=1, cex=1.75, bty="n", horiz=FALSE)

text(-80, 1, expression(italic("Methane")), cex=2, font=3, col="gray40")
text(-25, 1, expression(italic("Other")), cex=2, font=3, col="black")
points(-8.1, 1.33, col="darkorange", pch=19, cex=1.5)
legend(-15.5, 1.33, legend=expression("Atm." ~ CO[2]), text.col="darkorange", text.font=4, cex=1.75, bty="n")


##===============================  Delta Notation for C Isotopes

r13C_sample <- 0.011*1
r12C_sample <- 0.989*1

ratio_standard <- 0.0112372

del13C <- (((r13C_sample/r12C_sample)/(ratio_standard))-1)*1000

del13C



