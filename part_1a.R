########################################################################################
# Project: Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison
# Prepared by: S. Arcusa
# Version: 1
# Date: 07-10-2019
# Part 1a: Script to calculate multiple linear regression worldwide and the Southwest
#######################################################################################

library('ncdf4')
library(plyr)
library(dplyr)
library(relaimpo)

# Run this part for each model ensemble
setwd("/projects/pd_lab/sha59/LME/run003") # Change to folder containing all variables of one model run

ncname <- list.files(pattern=".nc")
dname <- vector()
for(i in 1:length(ncname)){
  dname[i] <- unlist(strsplit(ncname[i], split = ".", fixed = T))[[1]]
}

ncin <- nc_open(ncname[1])
DSTFLXT <- ncvar_get(ncin,dname[1])
ncin <- nc_open(ncname[2])
FSNO_012 <- ncvar_get(ncin,dname[2])
ncin <- nc_open(ncname[3])
SOILICE <- ncvar_get(ncin,dname[3])
liq_012 <- SOILICE[,,1,]
SOILICE_top <- SOILICE[,,1,]
ncin <- nc_open(ncname[4])
SOILLIQ <- ncvar_get(ncin,dname[3])
soilliq_012_SW <- SOILLIQ[xmin:xmax,ymin:ymax,1,]
SOILLIQ_top <- SOILLIQ[,,1,]
ncin <- nc_open(ncname[5])
SOILWATER <- ncvar_get(ncin,dname[5])
ncin <- nc_open(ncname[6])
TLAI <- ncvar_get(ncin,dname[6])
ncin <- nc_open(ncname[7])
TSAI <- ncvar_get(ncin,dname[7])
ncin <- nc_open(ncname[8])
U10 <- ncvar_get(ncin,dname[8])
ncin  <- nc_open(ncname)
H2OSOI_all  <- ncvar_get(ncin,dname)
H2OSOI <- H2OSOI_all[,,1,]

# Calculate f_v
f_v <- array(data = NA, dim = c(144,96,13872))
for(i in 1:144){
  for(j in 1:96){
    f_v[i,j,] <- (TLAI[i,j,]+TSAI[i,j,])/0.3
  }
}

f_v[f_v >1] = 1
f_v[f_v <0] = 0

# Calculate liquid water ratio
liq_water_ratio <- array(data = NA, dim = c(144,96,13872))
for(i in 1:144){
  for(j in 1:96){
    liq_water_ratio[i,j,] <- SOILLIQ_top[i,j,]/(SOILLIQ_top[i,j,] + SOILICE_top[i,j,])
  }
}

# Calculate f_m
f_m  <- array(data = NA, dim = c(144,96,13872))
for(i in 1:144){
  for(j in 1:96){
    f_m[i,j,] <- (100- LAKE[i,j] -WETLAND[i,j])/100 * (1-FSNO[i,j,])*(1-f_v[i,j,])*liq_water_ratio[i,j,]
  }
}

# Calculate Wind speed cubed
U10cubed <- array(data = NA, dim = c(144,96,13872)) 
for(i in 1:144){
  for(j in 1:96){
    U10cubed[i,j,]= U10[i,j,]^3 
  }
}

# Calculate DSTFLXT in g/m2/s
DSTFLXTg  <- array(data = NA, dim = c(144,96,13872)) 
for(i in 1:144){
  for(j in 1:96){
    DSTFLXTg[i,j,]= DSTFLXT[i,j,]*1000 
  }
}

# Calculate annuals

DSTFLXT_ann <- array(data= NA, dim= c(144,96,1156))
SOILWATER_ann <- array(data= NA, dim= c(144,96,1156))
f_m_ann <- array(data= NA, dim= c(144,96,1156))
U10cubed_ann <- array(data= NA, dim= c(144,96,1156))
H2OSOI_ann  <- array(data= NA, dim= c(144,96,1156))

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    D <- data.frame(FLUX = DSTFLXTg[i,j,])
    D$Years  <- Years
    D_ann <- ddply(D, .(Years), summarise,
                   Emissions = mean(FLUX, na.rm = T))
    DSTFLXT_ann[i,j,] <- D_ann$Emissions
    
    S <- data.frame(WATER = SOILWATER[i,j,])
    S$Years  <-  Years
    S_ann <- ddply(S, .(Years), summarise,
                   mean_moisture = mean(WATER, na.rm = T))
    SOILWATER_ann[i,j,] <- S_ann$mean_moisture
    
    F_M  <- data.frame(BARE = f_m[i,j,])
    F_M$Years  <- Years
    F_M_ann <- ddply(F_M, .(Years), summarise,
                     Potential = mean(BARE, na.rm = T))
    f_m_ann[i,j,] <- F_M_ann$Potential
    
    U <- data.frame(WIND = U10cubed[i,j,])
    U$Years  <- Years
    U_ann <- ddply(U, .(Years), summarise,
                   meanU10cubed = mean(WIND, na.rm = T))
    U10cubed_ann[i,j,] <- U_ann$meanU10cubed
    
  }
}

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    F_M  <- data.frame(BARE = f_m[i,j,])
    F_M$Years  <- Years
    F_M_ann <- ddply(F_M, .(Years), summarise,
                     Potential = mean(BARE, na.rm = T))
    f_m_ann[i,j,] <- F_M_ann$Potential
  }
}

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    H  <- data.frame(VOL = H2OSOI[i,j,])
    H$Years  <- Years
    H_ann <- ddply(H, .(Years), summarise,
                     Vol_water = mean(VOL, na.rm = T))
    H2OSOI_ann[i,j,] <- H_ann$Vol_water
  }
}

SOILWATER_ann[SOILWATER_ann == "NaN"] <- NA
U10cubed_ann[U10cubed_ann == "NaN"] <- NA
f_m_ann[f_m_ann == "NaN"]  <- NA
DSTFLXT_ann[DSTFLXT_ann == "NaN"]  <- NA
H2OSOI_ann[H2OSOI_ann == "NaN"] <- NA

save.image(file="LME_003_ann.RData") # Saves R 

# Regression analysis

# Do this only to test H2OSOI
#SOILWATER <- H2OSOI
#SOILWATER_ann  <- H2OSOI_ann

R2 <- array(data = NA, dim = c(144,96))
lmg_f_m <- array(data = NA, dim = c(144,96))
lmg_U10cubed <- array(data = NA, dim = c(144,96))
lmg_SOILWATER <- array(data = NA, dim = c(144,96))
coef_f_m<- array(data = NA, dim = c(144,96))
coef_U10cubed <- array(data = NA, dim = c(144,96))
coef_SOILWATER <- array(data = NA, dim = c(144,96))

R2.2 <- array(data = NA, dim = c(144,96))
lmg_f_m.2 <- array(data = NA, dim = c(144,96))
lmg_U10cubed.2 <- array(data = NA, dim = c(144,96))
lmg_SOILWATER.2 <- array(data = NA, dim = c(144,96))
coef_f_m.2 <- array(data = NA, dim = c(144,96))
coef_U10cubed.2 <- array(data = NA, dim = c(144,96))
coef_SOILWATER.2 <- array(data = NA, dim = c(144,96))

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    if(is.na(SOILWATER[i,j,]) == TRUE | is.na(f_m[i,j,]) == TRUE | (sum(DSTFLXTg[i,j,]) == 0) == TRUE | (abs(max(SOILWATER[i,j,]) - min(SOILWATER[i,j,])) < .Machine$double.eps ^ 0.5) == TRUE){
      R2[i,j] <- NA
      lmg_f_m[i,j]  <- NA
      lmg_U10cubed[i,j]  <- NA
      lmg_SOILWATER[i,j]  <- NA
      coef_f_m[i,j]  <- NA
      coef_U10cubed[i,j]  <- NA
      coef_SOILWATER[i,j]  <- NA
      
      R2.2[i,j] <- NA
      lmg_f_m.2[i,j]  <- NA
      lmg_U10cubed.2[i,j]  <- NA
      lmg_SOILWATER.2[i,j]  <- NA
      coef_f_m.2[i,j]  <- NA
      coef_U10cubed.2[i,j]  <- NA
      coef_SOILWATER.2[i,j]  <- NA
      
    } else {
      
      MLR <- lm(DSTFLXTg[i,j,] ~ f_m[i,j,] + U10cubed[i,j,] + SOILWATER[i,j,], na.action = "na.omit")
      RELAIMPO <- try(calc.relimp(MLR,type = c("lmg"), rela = T, na.action = "na.omit"))
      
      if(!class(RELAIMPO)=="try-error"){
        R2[i,j] <- RELAIMPO@R2
        lmg_f_m[i,j]  <- RELAIMPO@lmg[[1]]
        lmg_U10cubed[i,j]  <- RELAIMPO@lmg[[2]]
        lmg_SOILWATER[i,j] <- RELAIMPO@lmg[[3]]
        coef_f_m[i,j] <- RELAIMPO@ave.coeffs[[1,3]]
        coef_U10cubed[i,j] <- RELAIMPO@ave.coeffs[[2,3]]
        coef_SOILWATER[i,j] <- RELAIMPO@ave.coeffs[[3,3]]    
      } else {
        
        MLR <- lm(DSTFLXTg[i,j,] ~ U10cubed[i,j,] + SOILWATER[i,j,], na.action = "na.omit")
        RELAIMPO <- try(calc.relimp(MLR,type = c("lmg"), rela = T, na.action = "na.omit"))
        
        R2.2[i,j] <- RELAIMPO@R2
        lmg_U10cubed.2[i,j]  <- RELAIMPO@lmg[[1]]
        lmg_SOILWATER.2[i,j] <- RELAIMPO@lmg[[2]]
        coef_U10cubed.2[i,j] <- RELAIMPO@ave.coeffs[[1,2]]
        coef_SOILWATER.2[i,j] <- RELAIMPO@ave.coeffs[[2,2]] 
        
        
      }
    }
  }
}

save.image(file="LME_003_coef_H2OSOI.RData")

R2 <- array(data = NA, dim = c(144,96))
lmg_f_m <- array(data = NA, dim = c(144,96))
lmg_U10cubed <- array(data = NA, dim = c(144,96))
lmg_SOILWATER <- array(data = NA, dim = c(144,96))
coef_f_m<- array(data = NA, dim = c(144,96))
coef_U10cubed <- array(data = NA, dim = c(144,96))
coef_SOILWATER <- array(data = NA, dim = c(144,96))

R2.2 <- array(data = NA, dim = c(144,96))
lmg_f_m.2 <- array(data = NA, dim = c(144,96))
lmg_U10cubed.2 <- array(data = NA, dim = c(144,96))
lmg_SOILWATER.2 <- array(data = NA, dim = c(144,96))
coef_f_m.2 <- array(data = NA, dim = c(144,96))
coef_U10cubed.2 <- array(data = NA, dim = c(144,96))
coef_SOILWATER.2 <- array(data = NA, dim = c(144,96))

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    if(is.na(SOILWATER_ann[i,j,]) == TRUE | is.na(f_m_ann[i,j,]) == TRUE | (sum(DSTFLXT_ann[i,j,]) == 0) == TRUE | (abs(max(SOILWATER_ann[i,j,]) - min(SOILWATER_ann[i,j,])) < .Machine$double.eps ^ 0.5) == TRUE){
      R2[i,j] <- NA
      lmg_f_m[i,j]  <- NA
      lmg_U10cubed[i,j]  <- NA
      lmg_SOILWATER[i,j]  <- NA
      coef_f_m[i,j]  <- NA
      coef_U10cubed[i,j]  <- NA
      coef_SOILWATER[i,j]  <- NA
      
      R2.2[i,j] <- NA
      lmg_f_m.2[i,j]  <- NA
      lmg_U10cubed.2[i,j]  <- NA
      lmg_SOILWATER.2[i,j]  <- NA
      coef_f_m.2[i,j]  <- NA
      coef_U10cubed.2[i,j]  <- NA
      coef_SOILWATER.2[i,j]  <- NA
      
    } else {
      
      MLR <- lm(DSTFLXT_ann[i,j,] ~ f_m_ann[i,j,] + U10cubed_ann[i,j,] + SOILWATER_ann[i,j,], na.action = "na.omit")
      RELAIMPO <- try(calc.relimp(MLR,type = c("lmg"), rela = T, na.action = "na.omit"))
      
      if(!class(RELAIMPO)=="try-error"){
        R2[i,j] <- RELAIMPO@R2
        lmg_f_m[i,j]  <- RELAIMPO@lmg[[1]]
        lmg_U10cubed[i,j]  <- RELAIMPO@lmg[[2]]
        lmg_SOILWATER[i,j] <- RELAIMPO@lmg[[3]]
        coef_f_m[i,j] <- RELAIMPO@ave.coeffs[[1,3]]
        coef_U10cubed[i,j] <- RELAIMPO@ave.coeffs[[2,3]]
        coef_SOILWATER[i,j] <- RELAIMPO@ave.coeffs[[3,3]]    
      } else {
        
        MLR <- lm(DSTFLXT_ann[i,j,] ~ U10cubed_ann[i,j,] + SOILWATER_ann[i,j,], na.action = "na.omit")
        RELAIMPO <- try(calc.relimp(MLR,type = c("lmg"), rela = T, na.action = "na.omit"))
        
        R2.2[i,j] <- RELAIMPO@R2
        lmg_U10cubed.2[i,j]  <- RELAIMPO@lmg[[1]]
        lmg_SOILWATER.2[i,j] <- RELAIMPO@lmg[[2]]
        coef_U10cubed.2[i,j] <- RELAIMPO@ave.coeffs[[1,2]]
        coef_SOILWATER.2[i,j] <- RELAIMPO@ave.coeffs[[2,2]] 
        
        
      }
    }
  }
}

save.image(file="LME_003_ann_coef_H2OSOI.RData")

# Calculate decades

Time <- seq(as.Date("850-01-01"), by = "month", length.out = 13872)
Years <- as.numeric(format(as.Date(Time), format = "%Y"))
Months <- as.numeric(format(as.Date(Time), format = "%m"))

Years_ann <- H_ann$Years
Decades <- as.data.frame(Years_ann) %>% mutate(Decades = Years_ann - (Years_ann %% 10))
Decade  <- Decades$Decades

DSTFLXT_dec <- array(data= NA, dim= c(144,96,116))
SOILWATER_dec <- array(data= NA, dim= c(144,96,116))
f_m_dec <- array(data= NA, dim= c(144,96,116))
U10cubed_dec <- array(data= NA, dim= c(144,96,116))
H2OSOI_dec <- array(data= NA, dim= c(144,96,116))

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    D <- data.frame(FLUX = DSTFLXTg[i,j,])
    D$Decade  <- Decade
    D_dec <- ddply(D, .(Decade), summarise,
                   Emissions = mean(FLUX, na.rm = T))
    DSTFLXT_dec[i,j,] <- D_dec$Emissions
    
    S <- data.frame(WATER = SOILWATER[i,j,])
    S$Decade  <-  Decade
    S_dec <- ddply(S, .(Decade), summarise,
                   mean_moisture = mean(WATER, na.rm = T))
    SOILWATER_dec[i,j,] <- S_dec$mean_moisture
    
    F_M  <- data.frame(BARE = f_m[i,j,])
    F_M$Decade  <- Decade
    F_M_dec <- ddply(F_M, .(Decade), summarise,
                     Potential = mean(BARE, na.rm = T))
    f_m_dec[i,j,] <- F_M_dec$Potential
    
    U <- data.frame(WIND = U10cubed[i,j,])
    U$Decade  <- Decade
    U_dec <- ddply(U, .(Decade), summarise,
                   meanU10cubed = mean(WIND, na.rm = T))
    U10cubed_dec[i,j,] <- U_dec$meanU10cubed
    
  }
}

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    H <- data.frame(H2O = H2OSOI[i,j,])
    H$Decade  <- Decade
    H_dec <- ddply(H, .(Decade), summarise,
                   meanH2O = mean(H2O, na.rm = T))
    H2OSOI_dec[i,j,] <- H_dec$meanH2O
    
  }
}


SOILWATER_dec[SOILWATER_dec == "NaN"] <- NA
U10cubed_dec[U10cubed_dec == "NaN"] <- NA
f_m_dec[f_m_dec == "NaN"]  <- NA
DSTFLXT_dec[DSTFLXT_dec == "NaN"]  <- NA
H2OSOI_dec[H2OSOI_dec == "NaN"] <- NA

#Do this only for H2OSOI
SOILWATER_dec <- H2OSOI_dec

R2 <- array(data = NA, dim = c(144,96))
lmg_f_m <- array(data = NA, dim = c(144,96))
lmg_U10cubed <- array(data = NA, dim = c(144,96))
lmg_SOILWATER <- array(data = NA, dim = c(144,96))
coef_f_m<- array(data = NA, dim = c(144,96))
coef_U10cubed <- array(data = NA, dim = c(144,96))
coef_SOILWATER <- array(data = NA, dim = c(144,96))

R2.2 <- array(data = NA, dim = c(144,96))
lmg_f_m.2 <- array(data = NA, dim = c(144,96))
lmg_U10cubed.2 <- array(data = NA, dim = c(144,96))
lmg_SOILWATER.2 <- array(data = NA, dim = c(144,96))
coef_f_m.2 <- array(data = NA, dim = c(144,96))
coef_U10cubed.2 <- array(data = NA, dim = c(144,96))
coef_SOILWATER.2 <- array(data = NA, dim = c(144,96))

for(i in 1:144){
  print(i)
  for(j in 1:96){
    print(j)
    if(is.na(SOILWATER_dec[i,j,]) == TRUE | is.na(f_m_dec[i,j,]) == TRUE | (sum(DSTFLXT_dec[i,j,]) == 0) == TRUE | (abs(max(SOILWATER_dec[i,j,]) - min(SOILWATER_dec[i,j,])) < .Machine$double.eps ^ 0.5) == TRUE | is.na(DSTFLXT_dec[i,j,]) == TRUE){
      R2[i,j] <- NA
      lmg_f_m[i,j]  <- NA
      lmg_U10cubed[i,j]  <- NA
      lmg_SOILWATER[i,j]  <- NA
      coef_f_m[i,j]  <- NA
      coef_U10cubed[i,j]  <- NA
      coef_SOILWATER[i,j]  <- NA
      
      R2.2[i,j] <- NA
      lmg_f_m.2[i,j]  <- NA
      lmg_U10cubed.2[i,j]  <- NA
      lmg_SOILWATER.2[i,j]  <- NA
      coef_f_m.2[i,j]  <- NA
      coef_U10cubed.2[i,j]  <- NA
      coef_SOILWATER.2[i,j]  <- NA
      
    } else {
      
      MLR <- lm(DSTFLXT_dec[i,j,] ~ f_m_dec[i,j,] + U10cubed_dec[i,j,] + SOILWATER_dec[i,j,], na.action = "na.omit")
      RELAIMPO <- try(calc.relimp(MLR,type = c("lmg"), rela = T, na.action = "na.omit"))
      
      if(!class(RELAIMPO)=="try-error"){
        R2[i,j] <- RELAIMPO@R2
        lmg_f_m[i,j]  <- RELAIMPO@lmg[[1]]
        lmg_U10cubed[i,j]  <- RELAIMPO@lmg[[2]]
        lmg_SOILWATER[i,j] <- RELAIMPO@lmg[[3]]
        coef_f_m[i,j] <- RELAIMPO@ave.coeffs[[1,3]]
        coef_U10cubed[i,j] <- RELAIMPO@ave.coeffs[[2,3]]
        coef_SOILWATER[i,j] <- RELAIMPO@ave.coeffs[[3,3]]    
      } else {
        
        MLR <- lm(DSTFLXT_dec[i,j,] ~ U10cubed_dec[i,j,] + SOILWATER_dec[i,j,], na.action = "na.omit")
        RELAIMPO <- try(calc.relimp(MLR,type = c("lmg"), rela = T, na.action = "na.omit"))
        
        R2.2[i,j] <- RELAIMPO@R2
        lmg_U10cubed.2[i,j]  <- RELAIMPO@lmg[[1]]
        lmg_SOILWATER.2[i,j] <- RELAIMPO@lmg[[2]]
        coef_U10cubed.2[i,j] <- RELAIMPO@ave.coeffs[[1,2]]
        coef_SOILWATER.2[i,j] <- RELAIMPO@ave.coeffs[[2,2]] 
        
        
      }
    }
  }
}

save.image(file="LME_003_dec_coef_H2OSOI.RData")

# END OF SPECIFIC MODEL ENSEMBLE MEMBER ANALYSIS

####################################################
# Relative importance of variables to bare ground, only for the SW. Same analysis for the Southwest
library('ncdf4')
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggsignif)
library(relaimpo)


ncname <- list.files(pattern="nc")

dname <- vector()
for(i in 1:length(ncname)){
  dname[i] <- unlist(strsplit(ncname[i], split = ".085", fixed = T))[[1]]
}

newname <- paste0("SOILWATER_", seq(002,013,001))
SOILWATER_SW <- list()
for(i in 1:length(ncname)){
  ncin  <- nc_open(ncname[i])
  var_array <- ncvar_get(ncin,"SOILWATER_10CM")
  SOILWATER_SW[[i]] <- var_array[xmin:xmax,ymin:ymax,]
  print(i)
}

SOILWATER_ann <- array(data = NA, dim = c(5,6,1156,11))

for(k in 1:11){
  print(k)
  for(i in 1:5){
    #print(i)
    for(j in 1:6){
      #print(j)
      S  <- data.frame(SW = SOILWATER_SW[[k]][i,j,])
      S$Years  <- Years
      S_ann <- ddply(S, .(Years), summarise,
                     Soilwater = mean(SW, na.rm = T))
      SOILWATER_ann[i,j,,k] <- S_ann$Soilwater
    }
  }
}

# Prepare f_v

ncname <- list.files(pattern="nc")

dname <- vector()
for(i in 1:length(ncname)){
  dname[i] <- unlist(strsplit(ncname[i], split = ".085", fixed = T))[[1]]
}

newname <- paste0("TSAI", seq(002,013,001))
TSAI_SW <- list()
for(i in 1:length(ncname)){
  ncin  <- nc_open(ncname[i])
  var_array <- ncvar_get(ncin,"TSAI")
  TSAI_SW[[i]] <- var_array[xmin:xmax,ymin:ymax,]
  print(i)
}

ncname <- list.files(pattern="nc")

dname <- vector()
for(i in 1:length(ncname)){
  dname[i] <- unlist(strsplit(ncname[i], split = ".085", fixed = T))[[1]]
}

newname <- paste0("TLAI", seq(002,013,001))
TLAI_SW <- list()
for(i in 1:length(ncname)){
  ncin  <- nc_open(ncname[i])
  var_array <- ncvar_get(ncin,"TLAI")
  TLAI_SW[[i]] <- var_array[xmin:xmax,ymin:ymax,]
  print(i)
}

f_v <- array(data = NA, dim = c(5,6,13872,11))
for(k in 1:11){
  for(i in 1:5){
    for(j in 1:6){
      f_v[i,j,,k] <- (TLAI_SW[[k]][i,j,]+TSAI_SW[[k]][i,j,])/0.3
    }
  }
}
f_v[f_v >1] = 1
f_v[f_v <0] = 0

f_v_ann <- array(data = NA, dim = c(5,6,1156,11))

for(k in 1:11){
  print(k)
  for(i in 1:5){
    #print(i)
    for(j in 1:6){
      #print(j)
      S  <- data.frame(FV = f_v[i,j,,k])
      S$Years  <- Years
      S_ann <- ddply(S, .(Years), summarise,
                     Vegcover = mean(FV, na.rm = T))
      f_v_ann[i,j,,k] <- S_ann$Vegcover
    }
  }
}

# Prepare FM

# for 003
LAKE_SW <- LAKE[xmin:xmax,ymin:ymax]
WETLAND_SW <- WETLAND[xmin:xmax,ymin:ymax]
f_v_003_SW <- f_v[xmin:xmax,ymin:ymax,]
f_m_003_SW  <- array(data = NA, dim = c(5,6,13872))
for(i in 1:5){
  for(j in 1:6){
      f_m_003_SW[i,j,] <- (100- LAKE_SW[i,j] -WETLAND_SW[i,j])/100 * (1-FSNO_003_SW[i,j,])*(1-f_v_003_SW[i,j,])*liq_003_SW[i,j,]
  }
}
f_m_ann_003 <- array(data = NA, dim = c(5,6,1156))
for(i in 1:5){
  #print(i)
  for(j in 1:6){
    #print(j)
    S  <- data.frame(FM = f_m_003_SW[i,j,])
    S$Years  <- Years
    S_ann <- ddply(S, .(Years), summarise,
                   Bare = mean(FM, na.rm = T))
    f_m_ann_003[i,j,] <- S_ann$Bare
  }
}
f_m_ann_003[f_m_ann_003 == "NaN"] <- NA

FM_ann_all <- array(c(f_m_ann_002[xmin:xmax,ymin:ymax,], f_m_ann_003,f_m_ann_004[xmin:xmax,ymin:ymax,],f_m_ann_005[xmin:xmax,ymin:ymax,],f_m_ann_006[xmin:xmax,ymin:ymax,],f_m_ann_007[xmin:xmax,ymin:ymax,],f_m_ann_008[xmin:xmax,ymin:ymax,],f_m_ann_009[xmin:xmax,ymin:ymax,],f_m_ann_010[xmin:xmax,ymin:ymax,],f_m_ann_011[xmin:xmax,ymin:ymax,],f_m_ann_012[xmin:xmax,ymin:ymax,]), dim = c(5,6,1156,11))

#FM_SW <- FM_ann_all[xmin:xmax,ymin:ymax,,]

# Prepare U10cubed

U10cubed_ann_all <- array(c(U10cubed_ann_002, U10cubed_ann_003,U10cubed_ann_004,U10cubed_ann_005,U10cubed_ann_006,U10cubed_ann_007,U10cubed_ann_008,U10cubed_ann_009,U10cubed_ann_010,U10cubed_ann_011,U10cubed_ann_012), dim = c(144,96,1156,11))

U10cubed_SW <- U10cubed_ann_all[xmin:xmax,ymin:ymax,,]

# Preapre FSNO

FSNO_all <- array(c(FSNO_002_SW,FSNO_003_SW,FSNO_004_SW,FSNO_005_SW,FSNO_006_SW,FSNO_007_SW,FSNO_008_SW,FSNO_009_SW,FSNO_010_SW,FSNO_011_SW,FSNO_012_SW), dim = c(5,6,13872,11))

FSNO_ann_SW <- array(data = NA, dim = c(5,6,1156,11))

for(k in 1:11){
  print(k)
  for(i in 1:5){
    #print(i)
    for(j in 1:6){
      #print(j)
      S  <- data.frame(SNO = FSNO_all[i,j,,k])
      S$Years  <- Years
      S_ann <- ddply(S, .(Years), summarise,
                     Snow = mean(SNO, na.rm = T))
      FSNO_ann_SW[i,j,,k] <- S_ann$Snow
    }
  }
}

FSNO_ann_SW[FSNO_ann_SW == "NaN"] <- NA

# Prepare liq water ratio for 002, 009, 012

soilliq_all <- array(c(soilliq_002_SW, soilliq_009_SW, soilliq_012_SW), dim = c(5,6,13872,3))
soilice_all <- array(c(soilice_002_SW, soilice_009_SW, soilice_012_SW), dim = c(5,6,13872,3))

# Calculate liquid water ratio
liq_water_ratio <- array(data = NA, dim = c(5,6,13872,3))
for(k in 1:3){
  print(k)
  for(i in 1:5){
    for(j in 1:6){
      liq_water_ratio[i,j,,k] <- soilliq_all[i,j,,k]/(soilliq_all[i,j,,k] + soilice_all[i,j,,k])
    }
  }
}

liq_all <- array(c(liq_water_ratio[,,,1],liq_003_SW,liq_004_SW, liq_005_SW, liq_006_SW, liq_007_SW, liq_008_SW, liq_water_ratio[,,,2], liq_010_SW, liq_011_SW, liq_water_ratio[,,,1]), dim = c(5,6,13872,11))

liq_ann_SW <- array(data = NA, dim = c(5,6,1156,11))

for(k in 1:11){
  print(k)
  for(i in 1:5){
    #print(i)
    for(j in 1:6){
      #print(j)
      S  <- data.frame(liq = liq_all[i,j,,k])
      S$Years  <- Years
      S_ann <- ddply(S, .(Years), summarise,
                     liq_water = mean(liq, na.rm = T))
      liq_ann_SW[i,j,,k] <- S_ann$liq_water
    }
  }
}

liq_ann_SW[liq_ann_SW == "NaN"] <- NA
FM_SW[FM_SW == "NaN"] <- NA
f_v_SW[f_v_SW == "NaN"] <- NA

SW_lmg_fv <- array(data = NA, dim = c(5,6,11))
SW_lmg_fsno <- array(data = NA, dim = c(5,6,11))
SW_lmg_liq <- array(data = NA, dim = c(5,6,11))
SW_coef_fv <- array(data = NA, dim = c(5,6,11))
SW_coef_fsno <- array(data = NA, dim = c(5,6,11))
SW_coef_liq <- array(data = NA, dim = c(5,6,11))

SW_lmg_fv.2 <- array(data = NA, dim = c(5,6,11))
SW_lmg_fsno.2 <- array(data = NA, dim = c(5,6,11))
SW_lmg_liq.2 <- array(data = NA, dim = c(5,6,11))
SW_coef_fv.2 <- array(data = NA, dim = c(5,6,11))
SW_coef_fsno.2 <- array(data = NA, dim = c(5,6,11))
SW_coef_liq.2 <- array(data = NA, dim = c(5,6,11))

for(k in 1:11){
  for(i in 1:5){
    print(i)
    for(j in 1:6){
      print(j)
      if(is.na(FM_SW[i,j,,k]) == TRUE | is.na(f_v_ann_SW[i,j,,k]) == TRUE | is.na(FSNO_ann_SW[i,j,,k]) == TRUE | (abs(max(liq_ann_SW[i,j,,k]) - min(liq_ann_SW[i,j,,k])) < .Machine$double.eps ^ 0.5) == TRUE | all(FM_SW[i,j,,k] == 0) == TRUE | all(f_v_ann_SW[i,j,,k] == 1)){
        
        SW_lmg_fv[i,j,k]  <- NA
        SW_lmg_fsno[i,j,k]  <- NA
        SW_lmg_liq[i,j,k]  <- NA
        SW_coef_fv[i,j,k]  <- NA
        SW_coef_fsno[i,j,k]  <- NA
        SW_coef_liq[i,j,k]  <- NA
        
        SW_lmg_fv.2[i,j,k]  <- NA
        SW_lmg_fsno.2[i,j,k]  <- NA
        SW_lmg_liq.2[i,j,k]  <- NA
        SW_coef_fv.2[i,j,k]  <- NA
        SW_coef_fsno.2[i,j,k]  <- NA
        SW_coef_liq.2[i,j,k]  <- NA
        
      } else {
        
        MLR <- lm(FM_SW[i,j,,k] ~  FSNO_ann_SW[i,j,,k] + f_v_ann_SW[i,j,,k] + liq_ann_SW[i,j,,k], na.action = "na.omit")
        RELAIMPO <- try(calc.relimp(MLR,type = c("lmg"), rela = T, na.action = "na.omit"))
        
        if(!class(RELAIMPO)=="try-error"){
          SW_lmg_fsno[i,j,k]  <- RELAIMPO@lmg[[1]]
          SW_lmg_fv[i,j,k] <- RELAIMPO@lmg[[2]]
          SW_lmg_liq[i,j,k] <- RELAIMPO@lmg[[3]]
          SW_coef_fsno[i,j,k] <- RELAIMPO@ave.coeffs[[1,3]]
          SW_coef_fv[i,j,k] <- RELAIMPO@ave.coeffs[[2,3]]
          SW_coef_liq[i,j,k] <- RELAIMPO@ave.coeffs[[3,3]]
        } else {
          
          SW_lmg_fsno.2[i,j,k]  <- NA
          SW_lmg_fv.2[i,j,k] <- NA
          SW_lmg_liq.2[i,j,k] <- NA
          SW_coef_fsno.2[i,j,k] <- NA
          SW_coef_fv.2[i,j,k] <- NA
          SW_coef_liq.2[i,j,k] <- NA
          
          
        }
      }
    }
  }
}

mean(SW_lmg_fsno, na.rm = T)*100
mean(SW_lmg_fv, na.rm = T)*100
mean(SW_lmg_liq, na.rm = T)*100
sd(SW_lmg_fsno, na.rm = T)*100
sd(SW_lmg_fv, na.rm = T)*100
sd(SW_lmg_liq, na.rm = T)*100
