###
# Script for Fig. 1 in Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison Approach. Written by S.Arcusa

library(ncdf4)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(tidyr)
library(tools)
library(raster)
library(rgdal)
library(maps)
library(rasterVis)
library(viridis)
library(autoimage)

load_filenames <- function(where){
  
  setwd(where)
  ncname <- list.files(pattern=".nc")
  filenames <- list(ncname)
  
}

setwd("/scratch/sha59/LME/NetCDF") #Change to directory of surfdata file

ncin  <- nc_open("surfdata_1.9x2.5_simyr1350_c131018.nc")
LON <- ncvar_get(ncin,"LONGXY")[,1]
LAT <- ncvar_get(ncin,"LATIXY")[1,]
Time <- seq(as.Date("850-01-01"), by = "month", length.out = 13872)
Years <- as.numeric(format(as.Date(Time), format = "%Y"))
Months <- as.numeric(format(as.Date(Time), format = "%m"))

xmin  <- which.min(abs(240-LON))
xmax  <- which.min(abs(251-LON))
ymin  <- which.min(abs(31-LAT))
ymax  <- which.min(abs(41-LAT))

nsteps <- 13872

#path to parent folder containing data
folder.dir = "/projects/pd_lab/sha59/LME" 

all.ens  <- sprintf('%0.3d',2:11)
all.ens.names <- paste0(sprintf('%0.3d',2:11),"/")
all.wd.names <- paste0("run", all.ens.names)
all.data <- list()
all.data_ <- list()
SW_av <- vector()
spatial_ <- list()
spatial  <- list()
full_ <- list()
full  <- list()
var_av <- matrix(NA, ncol = 96, nrow = 144)

for(i in 1:length(all.ens)){
  print(i)
  setwd(file.path(folder.dir,all.wd.names[i]))
  where = getwd()
  run_files <- unlist(load_filenames(where))
  
  dname <- vector()
  for(k in 1:length(run_files)){
    dname[k] <- unlist(strsplit(run_files[k], split = ".", fixed = T))[[1]]
  }
    
  get_vars <- dname[c(grep(pattern = "DSTFLXT", dname))]
  get_files <- run_files[c(grep(pattern = "DSTFLXT", run_files))]
  
  for(j in 1:length(get_vars)){
    print(j)
    ncin <- nc_open(get_files[j])
    var <- ncvar_get(ncin,get_vars[j])
    SW_array <- var[xmin:xmax,ymin:ymax,]
    
    for(n in 1:144){
      for(m in 1:96){
        var_av[n,m] <- mean(var[n,m,], na.rm = T)
      }
    }
    
    for(m in 1:nsteps){
      SW_av[m] <- mean(SW_array[,,m], na.rm = T)
    }
    
    full_ <- var_av
    all.data[[j]] <- SW_av
    spatial_[[j]]  <- SW_array
    
  }
  
  full[[i]] <- full_
  names(all.data) <- get_vars
  all.data_[[i]] <- all.data
  spatial[[i]] <- spatial_

}
names(all.data_) <- all.ens

f <- array(c(full[[1]], full[[2]], full[[3]], full[[4]], full[[5]], full[[6]], full[[7]], full[[8]], full[[9]], full[[10]]), dim = c(144,96,10))

f[f == "NaN"]  <- NA

spatial.mean_ <- list()
ens <- matrix(NA, ncol = 6, nrow = 5)
spatial.mean  <- matrix(NA, ncol = 6, nrow = 5)
spatial.sd  <- matrix(NA, ncol = 6, nrow = 5)

sp <- c(spatial[[1]], spatial[[2]], spatial[[3]], spatial[[4]], spatial[[5]],spatial[[6]],spatial[[7]], spatial[[8]], spatial[[9]], spatial[[10]])

sp[sp == "NaN"] <- NA
remove(spatial, spatial_, var, ncin)

for(k in 1:10){
  for(i in 1:5){
    for(j in 1:6){
      
      cell_mean <- mean(sp[[k]][i,j,], na.rm = T)
      ens[i,j] <- cell_mean
      
    }
  }
  
  spatial.mean_[[k]] <- ens
}

data <- array(as.numeric(unlist(spatial.mean_)), dim = c(5,6,10))

for(i in 1:5){
  for(j in 1:6){
    spatial.mean[i,j] <- mean(data[i,j,], na.rm = T)
    spatial.sd[i,j] <- sd(data[i,j,], na.rm = T)
    
  }
}

spatial.mean[spatial.mean == "NaN"] <- NA
spatial.sd[spatial.sd == "NaN"]  <- NA

spatial.mean <- round(spatial.mean*31536000000, digits = 2) # convert from kg/m2/s to g/m2/yr

##### Plotting

Mybrks  <- exp(seq(log(0.01), log(65), length.out = 5))

nb <- length(Mybrks)-1

my.yellows <- colorRampPalette(c("white", "yellow", "brown"), interpolate = ("spline"), bias = 1, space = "rgb")

pal = colorRampPalette(c("white", "yellow", "brown"))

lon = c(-121,-105)
lat = c(30,40)

spatial.map <- raster(t(spatial.mean), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
spatial.map <- flip(spatial.map, direction = "y")
#spatial.map <- rotate(spatial.map)
#spatial.map <- trim(spatial.map)

SW.region = extent(-121,-105,30,40)
usa <- map("state", add =F)

setwd("/projects/pd_lab/sha59/LME")

pdf("SW_DSTFLXT_mean.pdf")
plot(spatial.map, col = my.yellows(nb), breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,65), colNA = "light grey")
map(usa, add = T)
map("world", add = T, fill = F)
dev.off()
