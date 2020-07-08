########################################################################################
# Project: Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison
# Prepared by: S. Arcusa
# Version: 1
# Date: 07-10-2019
# Part II: Creating map objects and plotting for the world and the Southwest
#######################################################################################

library(tools)
library(raster)
library(rgdal)
library(maps)
library(rasterVis)
library(viridis)
library(ggplot2)
library(ncdf4)
library(RColorBrewer)
library(leaflet)
library(maptools)

setwd("/scratch/sha59/LME/NetCDF") #Change to directory of surfdata file

ncin  <- nc_open("surfdata_1.9x2.5_simyr1350_c131018.nc")
LAKE <- ncvar_get(ncin,"PCT_LAKE")
WETLAND <- ncvar_get(ncin,"PCT_WETLAND")
lon <- ncvar_get(ncin,"LONGXY")[,1]
lat <- ncvar_get(ncin,"LATIXY")[1,]

folder.dir = "/projects/pd_lab/sha59/LME" # change to directory of parent folder from where files from part_1 were saved

all.ens  <- sprintf('%0.3d',2:11)
all.ens.names <- paste0(sprintf('%0.3d',2:11),"/")
all.wd.names <- paste0("run", all.ens.names)

data  <- list()
data_ <- list()

for(i in 1:length(all.ens)){
  
  setwd(file.path(folder.dir,all.wd.names[i]))
  
  dname <- vector()
  file.names  <- list.files(pattern=paste0(all.ens[i],".csv"))
  for(k in 1: length(file.names)){
    dname[k]  <- file_path_sans_ext(file.names[k])
  }
  
  for(j in 1:length(file.names)){
    data[[j]]  <- read.csv(file.names[j])[,-1]
  }
  names(data) <- dname
  
  data_[[i]] <- data
  
}

index = grep("R2",names(data))[1] # must be selecting R2
R2.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("R2.ann",names(data)) # must be selecting R2.ann
R2.ann.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("coef_SOILWATER",names(data))[1] # must be selecting coef_SOILWATER
coef_SOILWATER.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("coef_SOILWATER.ann",names(data)) # must be selecting coef_SOILWATER.ann
coef_SOILWATER.ann.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("coef_U10cubed",names(data))[1] # must be selecting coef_U10cubed
coef_U10cubed.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("coef_U10cubed.ann",names(data)) # must be selecting coef_U10cubed.ann
coef_U10cubed.ann.all <-array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("coef_f_m",names(data))[1] # must be selecting coef_f_m
coef_f_m.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("coef_f_m.ann",names(data)) # must be selecting coef_f_m.ann
coef_f_m.ann.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("lmg_SOILWATER",names(data))[1] # must be selecting lmg_SOILWATER
lmg_SOILWATER.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("lmg_SOILWATER.ann",names(data)) # must be selecting lmg_SOILWATER
lmg_SOILWATER.ann.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("lmg_f_m",names(data))[1] # must be selecting lmg_f_m
lmg_f_m.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("lmg_f_m.ann",names(data)) # must be selecting lmg_f_m
lmg_f_m.ann.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("lmg_U10cubed",names(data))[1] # must be selecting lmg_U10cubed
lmg_U10cubed.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("lmg_U10cubed.ann",names(data)) # must be selecting lmg_U10cubed.ann
lmg_U10cubed.ann.all <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(144,96,10))

index = grep("importance",names(data)) # must be selecting importance_to_f_m
importance_to_fm.all  <- array(c(unlist(data_[[1]][index]), unlist(data_[[2]][index]), unlist(data_[[3]][index]), unlist(data_[[4]][index]), unlist(data_[[5]][index]), unlist(data_[[6]][index]), unlist(data_[[7]][index]), unlist(data_[[8]][index]), unlist(data_[[9]][index]), unlist(data_[[10]][index])), dim = c(3,3,10))

R2.all[R2.all == "NaN"] <- NA
R2.ann.all[R2.ann.all == "NaN"] <- NA
coef_SOILWATER.ann.all[coef_SOILWATER.ann.all == "NaN"] <- NA
coef_U10cubed.all[coef_U10cubed.all == "NaN"] <- NA
coef_SOILWATER.all[coef_SOILWATER.all == "NaN"] <- NA
coef_U10cubed.ann.all[coef_U10cubed.ann.all == "NaN"] <- NA
coef_f_m.all[coef_f_m.all == "NaN"] <- NA
coef_f_m.ann.all[coef_f_m.ann.all == "NaN"] <- NA
lmg_SOILWATER.all[lmg_SOILWATER.all == "NaN"] <- NA
lmg_SOILWATER.ann.all[lmg_SOILWATER.ann.all == "NaN"] <- NA
lmg_f_m.all[lmg_f_m.all == "NaN"] <- NA
lmg_f_m.ann.all[lmg_f_m.ann.all == "NaN"] <- NA
lmg_U10cubed.all[lmg_U10cubed.all == "NaN"] <- NA
lmg_U10cubed.ann.all[lmg_U10cubed.ann.all == "NaN"] <- NA
importance_to_fm.all[importance_to_fm.all == "NaN"] <- NA

print("all arrays created")

R2 <- matrix(data = NA, nrow = 144, ncol = 96)
coef_SOILWATER <- matrix(data = NA, nrow = 144, ncol = 96)
coef_U10cubed <- matrix(data = NA, nrow = 144, ncol = 96)
coef_f_m <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_SOILWATER <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_U10cubed <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_f_m <- matrix(data = NA, nrow = 144, ncol = 96)

R2.sd <- matrix(data = NA, nrow = 144, ncol = 96)
coef_SOILWATER.sd <- matrix(data = NA, nrow = 144, ncol = 96)
coef_U10cubed.sd <- matrix(data = NA, nrow = 144, ncol = 96)
coef_f_m.sd <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_SOILWATER.sd <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_U10cubed.sd <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_f_m.sd <- matrix(data = NA, nrow = 144, ncol = 96)

R2.ann <- matrix(data = NA, nrow = 144, ncol = 96)
coef_SOILWATER.ann <- matrix(data = NA, nrow = 144, ncol = 96)
coef_U10cubed.ann <- matrix(data = NA, nrow = 144, ncol = 96)
coef_f_m.ann <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_SOILWATER.ann <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_U10cubed.ann <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_f_m.ann <- matrix(data = NA, nrow = 144, ncol = 96)

R2.sd.ann <- matrix(data = NA, nrow = 144, ncol = 96)
coef_SOILWATER.sd.ann <- matrix(data = NA, nrow = 144, ncol = 96)
coef_U10cubed.sd.ann <- matrix(data = NA, nrow = 144, ncol = 96)
coef_f_m.sd.ann <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_SOILWATER.sd.ann <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_U10cubed.sd.ann <- matrix(data = NA, nrow = 144, ncol = 96)
lmg_f_m.sd.ann <- matrix(data = NA, nrow = 144, ncol = 96)

for(i in 1:144){
  for(j in 1:96){
    R2[i,j] <- mean(R2.all[i,j,], na.rm = T)
    coef_SOILWATER[i,j] <- mean(coef_SOILWATER.all[i,j,], na.rm = T)
    coef_U10cubed[i,j] <- mean(coef_U10cubed.all[i,j,], na.rm = T)
    coef_f_m[i,j] <- mean(coef_f_m.all[i,j,], na.rm = T)
    lmg_SOILWATER[i,j] <- mean(lmg_SOILWATER.all[i,j,], na.rm = T)
    lmg_U10cubed[i,j] <- mean(lmg_U10cubed.all[i,j,], na.rm = T)
    lmg_f_m[i,j] <- mean(lmg_f_m.all[i,j,], na.rm = T)
    
    R2.sd[i,j] <- sd(R2.all[i,j,], na.rm = T)
    coef_SOILWATER.sd[i,j] <- sd(coef_SOILWATER.all[i,j,], na.rm = T)
    coef_U10cubed.sd[i,j] <- sd(coef_U10cubed.all[i,j,], na.rm = T)
    coef_f_m.sd[i,j] <- sd(coef_f_m.all[i,j,], na.rm = T)
    lmg_SOILWATER.sd[i,j] <- sd(lmg_SOILWATER.all[i,j,], na.rm = T)
    lmg_U10cubed.sd[i,j] <- sd(lmg_U10cubed.all[i,j,], na.rm = T)
    lmg_f_m.sd[i,j] <- sd(lmg_f_m.all[i,j,], na.rm = T)   
    
    R2.ann[i,j] <- mean(R2.all[i,j,], na.rm = T)
    coef_SOILWATER.ann[i,j] <- mean(coef_SOILWATER.all[i,j,], na.rm = T)
    coef_U10cubed.ann[i,j] <- mean(coef_U10cubed.all[i,j,], na.rm = T)
    coef_f_m.ann[i,j] <- mean(coef_f_m.all[i,j,], na.rm = T)
    lmg_SOILWATER.ann[i,j] <- mean(lmg_SOILWATER.all[i,j,], na.rm = T)
    lmg_U10cubed.ann[i,j] <- mean(lmg_U10cubed.all[i,j,], na.rm = T)
    lmg_f_m.ann[i,j] <- mean(lmg_f_m.all[i,j,], na.rm = T)
    
    R2.sd.ann[i,j] <- sd(R2.all[i,j,], na.rm = T)
    coef_SOILWATER.sd.ann[i,j] <- sd(coef_SOILWATER.all[i,j,], na.rm = T)
    coef_U10cubed.sd.ann[i,j] <- sd(coef_U10cubed.all[i,j,], na.rm = T)
    coef_f_m.sd.ann[i,j] <- sd(coef_f_m.all[i,j,], na.rm = T)
    lmg_SOILWATER.sd.ann[i,j] <- sd(lmg_SOILWATER.all[i,j,], na.rm = T)
    lmg_U10cubed.sd.ann[i,j] <- sd(lmg_U10cubed.all[i,j,], na.rm = T)
    lmg_f_m.sd.ann[i,j] <- sd(lmg_f_m.all[i,j,], na.rm = T)    
    
  }
}

importance <- matrix(data = NA, nrow = 3, ncol = 3)

for(i in 1:3){
  for(j in 1:3){
 importance[i,j] <- mean(importance_to_fm.all[i,j,], na.rm = T)
  }
}
importance[,1] <- c("FSNO", "fv", "LWR")
colnames(importance) <- c("Parameter", "Means", "SD")

print(importance)

R2[R2 == "NaN"] <- NA
R2.ann[R2.ann == "NaN"] <- NA
coef_SOILWATER.ann[coef_SOILWATER.ann == "NaN"] <- NA
coef_U10cubed[coef_U10cubed == "NaN"] <- NA
coef_SOILWATER[coef_SOILWATER == "NaN"] <- NA
coef_U10cubed.ann[coef_U10cubed.ann == "NaN"] <- NA
coef_f_m[coef_f_m == "NaN"] <- NA
coef_f_m.ann[coef_f_m.ann == "NaN"] <- NA
lmg_SOILWATER[lmg_SOILWATER == "NaN"] <- NA
lmg_SOILWATER.ann[lmg_SOILWATER.ann == "NaN"] <- NA
lmg_f_m[lmg_f_m == "NaN"] <- NA
lmg_f_m.ann[lmg_f_m.ann == "NaN"] <- NA
lmg_U10cubed[lmg_U10cubed == "NaN"] <- NA
lmg_U10cubed.ann[lmg_U10cubed.ann == "NaN"] <- NA

print("calculated all means and sd for each variable")

# Preparing the maps

Mybrks <- seq(0,1,by=0.1)
nb <- length(Mybrks)-1

my.reds <- colorRampPalette(c("white", "red", "red4"), interpolate = ("spline"), bias = 1, space = "rgb")
my.blues <- colorRampPalette(c("white","dodgerblue", "blue3"), interpolate = ("spline"), bias = 1, space = "rgb")
my.greens <- colorRampPalette(c("white", "green", "green4"), interpolate = ("spline"), bias = 1, space = "rgb")

XFrm <- -180
XTo <- 180
YFrm <- -90
YTo <- 90
XStp <- 60
YStp <- 30

R2.map <- raster(t(R2), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
R2.map <- flip(R2.map, direction = "y")
R2.map <- rotate(R2.map)
R2.map <- trim(R2.map)

R2.ann.map <- raster(t(R2.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
R2.ann.map <- flip(R2.ann.map, direction = "y")
R2.ann.map <- rotate(R2.ann.map)
R2.ann.map <- trim(R2.ann.map)

SOILWATER.map <- raster(t(lmg_SOILWATER), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
SOILWATER.map <- flip(SOILWATER.map, direction = "y")
SOILWATER.map <- rotate(SOILWATER.map)

SOILWATER.ann.map <- raster(t(lmg_SOILWATER.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
SOILWATER.ann.map <- flip(SOILWATER.ann.map, direction = "y")
SOILWATER.ann.map <- rotate(SOILWATER.ann.map)

U10cubed.map <- raster(t(lmg_U10cubed), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
U10cubed.map <- flip(U10cubed.map, direction = "y")
U10cubed.map <- rotate(U10cubed.map)

U10cubed.ann.map <- raster(t(lmg_U10cubed.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
U10cubed.ann.map <- flip(U10cubed.ann.map, direction = "y")
U10cubed.ann.map <- rotate(U10cubed.ann.map)

f_m.map <- raster(t(lmg_f_m), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
f_m.map <- flip(f_m.map, direction = "y")
f_m.map <- rotate(f_m.map)

f_m.ann.map <- raster(t(lmg_f_m.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
f_m.ann.map <- flip(f_m.ann.map, direction = "y")
f_m.ann.map <- rotate(f_m.ann.map)

R2sd.map <- raster(t(R2.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
R2sd.map <- flip(R2sd.map, direction = "y")
R2sd.map <- rotate(R2sd.map)

R2sd.ann.map <- raster(t(R2.sd.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
R2sd.map <- flip(R2sd.map, direction = "y")
R2sd.map <- rotate(R2sd.map)

SOILWATERsd.map <- raster(t(lmg_SOILWATER.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
SOILWATERsd.map <- flip(SOILWATERsd.map, direction = "y")
SOILWATERsd.map <- rotate(SOILWATERsd.map)

SOILWATERsd.ann.map <- raster(t(lmg_SOILWATER.sd.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
SOILWATERsd.ann.map <- flip(SOILWATERsd.ann.map, direction = "y")
SOILWATERsd.ann.map <- rotate(SOILWATERsd.ann.map)

U10cubedsd.map <- raster(t(lmg_U10cubed.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
U10cubedsd.map <- flip(U10cubedsd.map, direction = "y")
U10cubedsd.map <- rotate(U10cubedsd.map)

U10cubedsd.ann.map <- raster(t(lmg_U10cubed.sd.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
U10cubedsd.ann.map <- flip(U10cubedsd.ann.map, direction = "y")
U10cubedsd.ann.map <- rotate(U10cubedsd.ann.map)

f_msd.map <- raster(t(lmg_f_m.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
f_msd.map <- flip(f_msd.map, direction = "y")
f_msd.map <- rotate(f_msd.map)

f_msd.ann.map <- raster(t(lmg_f_m.sd.ann), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
f_msd.ann.map <- flip(f_msd.ann.map, direction = "y")
f_msd.ann.map <- rotate(f_msd.ann.map)

print("created all map objects. Now will start plotting")

# Plotting

setwd("/projects/pd_lab/sha59/LME")

pdf("Ens_R2.pdf")
plot(R2.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_R2.png", width = 350, height = 350)
plot(R2.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_R2_ann.pdf")
plot(R2.ann.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_R2_ann.png", width = 350, height = 350)
plot(R2.ann.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_SOILWATER.pdf")
plot(SOILWATER.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_SOILWATER.png", width = 350, height = 350)
plot(SOILWATER.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_SOILWATER_ann.pdf")
plot(SOILWATER.ann.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_SOILWATER_ann.png", width = 350, height = 350)
plot(SOILWATER.ann.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_U10cubed.pdf")
plot(U10cubed.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_U10cubed.png", width = 350, height = 350)
plot(U10cubed.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_U10cubed_ann.pdf")
plot(U10cubed.ann.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_U10cubed_ann.png", width = 350, height = 350)
plot(U10cubed.ann.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_f_m.pdf")
plot(f_m.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_f_m.png", width = 350, height = 350)
plot(f_m.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_f_m_ann.pdf")
plot(f_m.ann.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1),colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_f_m_ann.png", width = 350, height = 350)
plot(f_m.ann.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_corr_map.pdf")
s <- stack(SOILWATER.map, f_m.map, U10cubed.map)
cor.map <- plotRGB(s, stretch = 'lin', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=F, labels = F, pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

png("Ens_corr_map.png", width = 350, height = 350)
s <- stack(SOILWATER.map, f_m.map, U10cubed.map)
cor.map <- plotRGB(s, stretch = 'hist', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=F, labels = F, pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

pdf("Ens_corr_map_ann.pdf")
s <- stack(SOILWATER.ann.map, f_m.ann.map, U10cubed.ann.map)
cor.map <- plotRGB(s, stretch = 'lin', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=F, labels = F, pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

png("Ens_corr_map_ann.png", width = 350, height = 350)
s <- stack(SOILWATER.ann.map, f_m.ann.map, U10cubed.ann.map)
cor.map <- plotRGB(s, stretch = 'hist', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=F, labels = F, pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

pdf("Ens_corr_R2.pdf")
old.par <- par(mfrow=c(2,1))
plot(R2.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
cor.map <- plotRGB(s, stretch = 'hist', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

png("Ens_corr_R2.png")
old.par <- par(mfrow=c(2,1))
plot(R2.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
cor.map <- plotRGB(s, stretch = 'hist', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

pdf("Ens_corr_R2_ann.pdf")
old.par <- par(mfrow=c(2,1))
plot(R2.ann.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
cor.map <- plotRGB(s, stretch = 'hist', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

png("Ens_corr_R2_ann.png")
old.par <- par(mfrow=c(2,1))
plot(R2.ann.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
cor.map <- plotRGB(s, stretch = 'hist', axes = T)
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
#axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
dev.off()

print("finished creating world plots")

#### SD maps

png("Ens_R2_sd.png", width = 350, height = 350)
plot(R2sd.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_R2_sd.pdf")
plot(R2sd.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_SOILWATER_sd.png", width = 350, height = 350)
plot(SOILWATERsd.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_SOILWATER_sd.pdf")
plot(SOILWATERsd.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_U10cubed_sd.png", width = 350, height = 350)
plot(U10cubedsd.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_U10cubed_sd.pdf")
plot(U10cubedsd.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_f_m_sd.png", width = 350, height = 350)
plot(f_msd.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_f_m_sd.pdf")
plot(f_msd.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_R2_sd_ann.png", width = 350, height = 350)
plot(R2sd.ann.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_R2_sd_ann.pdf")
plot(R2sd.ann.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_SOILWATER_sd_ann.png", width = 350, height = 350)
plot(SOILWATERsd.ann.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_SOILWATER_sd_ann.pdf")
plot(SOILWATERsd.ann.map, col = my.reds(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_U10cubed_sd_ann.png", width = 350, height = 350)
plot(U10cubedsd.ann.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_U10cubed_sd_ann.pdf")
plot(U10cubedsd.ann.map, col = my.blues(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

png("Ens_f_m_sd_ann.png", width = 350, height = 350)
plot(f_msd.ann.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_f_m_sd_ann.pdf")
plot(f_msd.ann.map, col = my.greens(nb), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

print("finished plotting SD maps")

Highs = c(high <- (length(which(R2 >= 0.5))/(length(which(is.na(R2) == F))))*100)
print(paste0("The highest R2 is ",Highs))
Highs_ann = c(high <- (length(which(R2.ann >= 0.5))/(length(which(is.na(R2.ann) == F))))*100)
print(paste0("The highest R2 for annual data is ", Highs_ann))

Means = mean(R2, na.rm = T)*100
SDs = sd(R2, na.rm = T)*100
Mins = min(R2, na.rm = T)*100
Maxs = max(R2, na.rm = T)*100
Median = median(R2, na.rm = T)*100
print(paste0("The mean variance and SD explained worldwide is ", Means, " and ", SDs))
hist(R2, xlab = parse(text = "R^2"))
IQR(R2, na.rm = T)*100

getmode <- function(v){
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v,uniqv)))]
}

getmode(R2)*100

print(paste0("SD for SOILWATER worldwide is ", mean(SOILWATERsd.ann.map@data@values, na.rm = T)))
print(paste0("SD for U10cubed worldwide is ", mean(U10cubedsd.ann.map@data@values, na.rm = T)))
print(paste0("SD for fm worldwide is ", mean(f_msd.ann.map@data@values, na.rm = T)))

# North Africa

NAfri <- R2[which.min(abs(350-LON)):which.min(abs(40-LON)),which.min(abs(0-LAT)):which.min(abs(41-LAT))]
length(which(!is.na(NAfri))) / length(which(!is.na(R2))) *100
min(NAfri, na.rm = T) *100
hist(NAfri)

###############################################################################
# Legend triangle

png("lengend_triangle.PNG", width = 350, height = 350)
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),asp=1,bty="n",axes=F,xlab="",ylab="")
sm <- 500
x <- do.call(c, sapply(1:(sm*sqrt(3)/2)/sm, function(i)(i*sm/sqrt(3)):(sm-i*sm/sqrt(3))/sm))
y <- do.call(c, sapply(1:(sm*sqrt(3)/2)/sm, function(i)rep(i, length((i*sm/sqrt(3)):(sm-i*sm/sqrt(3))))))
d.red = y
d.green = abs(sqrt(3)* x -y)/sqrt(3+1)
d.blue = abs(- sqrt(3)*x-y + sqrt(3))/ sqrt(3+1)
points(x,y,col=rgb(1-d.red, 1-d.green, 1-d.blue), pch=19)
dev.off()


############################################
# Analysis for the Southwest

setwd("/scratch/sha59/LME/NetCDF") #Change to directory of surfdata file

ncin  <- nc_open("surfdata_1.9x2.5_simyr1350_c131018.nc")
LON <- ncvar_get(ncin,"LONGXY")[,1]
LAT <- ncvar_get(ncin,"LATIXY")[1,]

xmin  <- which.min(abs(240-LON))
xmax  <- which.min(abs(252-LON))
ymin  <- which.min(abs(31-LAT))
ymax  <- which.min(abs(41-LAT))

SW_soilwater <- lmg_SOILWATER[xmin:xmax,ymin:ymax]
SW_f_m <- lmg_f_m[xmin:xmax,ymin:ymax]
SW_U10cubed <- lmg_U10cubed[xmin:xmax,ymin:ymax]
SW_R2 <- R2[xmin:xmax,ymin:ymax]

SW_soilwater_ann <- lmg_SOILWATER.ann[xmin:xmax,ymin:ymax]
SW_f_m_ann <- lmg_f_m.ann[xmin:xmax,ymin:ymax]
SW_U10cubed_ann <- lmg_U10cubed.ann[xmin:xmax,ymin:ymax]
SW_R2_ann <- R2.ann[xmin:xmax,ymin:ymax]

SW  <- data.frame(Means = c(mean(SW_soilwater, na.rm = T),
                            mean(SW_f_m, na.rm = T),
                            mean(SW_U10cubed, na.rm = T),
                            mean(SW_R2, na.rm = T),
                            mean(SW_soilwater_ann, na.rm = T),
                            mean(SW_f_m_ann, na.rm = T),
                            mean(SW_U10cubed_ann, na.rm = T),
                            mean(SW_R2_ann, na.rm = T))*100,
                  SDs = c(sd(SW_soilwater, na.rm = T),
                          sd(SW_f_m, na.rm = T),
                          sd(SW_U10cubed, na.rm = T),
                          sd(SW_R2, na.rm = T),
                          sd(SW_soilwater_ann, na.rm = T),
                          sd(SW_f_m_ann, na.rm = T),
                          sd(SW_U10cubed_ann, na.rm = T),
                          sd(SW_R2_ann, na.rm = T))*100,
                   Parameter = c("Soilwater", "fm", "U10cubed", "R2",
                                "Soilwater_ann", "fm_ann", "U10cubed_ann", "R2_ann"))

write.csv(SW, "SW_parameter.csv")

#########
# Comparison to the max soilwater input

highest <- tail(sort(lmg_SOILWATER),10)*100
find <- highest[10]/100

lmg_SOILWATER[which.max(lmg_SOILWATER[,ceiling(which(lmg_SOILWATER == find)/length(LON))]),ceiling(which(lmg_SOILWATER == find)/length(LON))]

print(paste0("the highest soilwater input in the world is ", round(max(lmg_SOILWATER, na.rm = T)*100), "% located at 32.8125N 93.0816E"))

# Plot the Southwest

SW.region = extent(-121,-105,30,40)
SW.SOILWATER <- raster::crop(SOILWATER.map, SW.region)
SW.R2 <- raster::crop(R2.map, SW.region)
SW.fm <- raster::crop(f_m.map, SW.region)
SW.U10 <- raster::crop(U10cubed.map, SW.region)
usa <- map("state", add =F)

setwd("/projects/pd_lab/sha59/LME") #Change to directory of surfdata file

pdf("SW_LME_soil_moisture.pdf")
plot(SW.SOILWATER, col = my.reds(nb), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
map(usa, add = T)
map("world", add = T, fill = F)
dev.off()

png("SW_LME_soil_moisture.png", width = 350, height = 350)
plot(SW.SOILWATER, col = my.reds(nb), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
map(usa, add = T)
map("world", add = T, fill = F)
dev.off()

pdf("SW_LME_fm.pdf")
plot(SW.fm, col = my.greens(nb), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
map(usa, add = T)
map("world", add = T, fill = F)
dev.off()

png("SW_LME_fm.png", width = 350, height = 350)
plot(SW.fm, col = my.greens(nb), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
map(usa, add = T)
map("world", add = T, fill = F)
dev.off()

pdf("SW_LME_U10cubed.pdf")
plot(SW.U10, col = my.blues(nb), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
map(usa, add = T)
map("world", add = T, fill = F) 
dev.off()

png("SW_LME_U10cubed.png", width = 350, height = 350)
plot(SW.U10, col = my.blues(nb), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1), colNA = "light grey")
map(usa, add = T)
map("world", add = T, fill = F)
dev.off()

pdf("SW_LME_R2.pdf")
plot(SW.R2,col = rev(viridis(nb)),breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
map(usa, add = T)
map("world", add = T, fill = F) 
dev.off()

png("SW_LME_R2.png", width = 350, height = 350)
plot(SW.R2,col = rev(viridis(nb)), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
map(usa, add = T)
map("world", add = T, fill = F) 
dev.off()

png("SW_ens_corr_map.png", width = 350, height = 350)
a <- stack(SW.SOILWATER, SW.fm, SW.U10)
cor.map <- plotRGB(a, stretch = 'hist', axes = T)
axis(1,tick=T,pos=30, las =1, at=seq(-120,-108,2))
axis(3, tick =T, labels= F, pos=42, las =1, at = seq(-120,-109,2))
axis(4,tick =T, labels= T, pos=-106, at=seq(30,41,2))
map("world", add=T, fill=F)
map(usa, add = T)
dev.off()

pdf("SW_ens_corr_map.pdf", width = 350, height = 350)
a <- stack(SW.SOILWATER, SW.fm, SW.U10)
cor.map <- plotRGB(a, stretch = 'lin', axes = T)
axis(1,tick=T,pos=30, las =1, at=seq(-120,-108,2))
axis(3, tick =T, labels= F, pos=42, las =1, at = seq(-120,-109,2))
axis(4,tick =T, labels= T, pos=-106, at=seq(30,41,2))
map("world", add=T, fill=F)
map(usa, add = T)
dev.off()

print("this is the end of part 2")
