########################################################################################
# Project: Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison
# Prepared by: S. Arcusa
# Version: 1
# Date: 07-10-2019
# Part III: Investigating megadroughts
#######################################################################################

library(ncdf4)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(readxl)
#devtools::install_github("nickmckay/lipd-utilities",subdir = "R")
library(lipdR)
#devtools::install_github("nickmckay/nuspectral")
#install.packages("rbacon")
library(BiocManager)
#devtools::install_github("nickmckay/geoChronR")
library(geoChronR)

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

folder.dir = "/projects/pd_lab/sha59/LME" #path to parent folder

all.ens  <- sprintf('%0.3d',2:11)
all.ens.names <- paste0(sprintf('%0.3d',2:11),"/")
all.wd.names <- paste0("run", all.ens.names)
all.data <- list()
all.data_ <- list()
SW_av <- vector()

for(i in 1:length(all.ens)){
  print(i)
  setwd(file.path(folder.dir,all.wd.names[i]))
  where = getwd()
  run_files <- unlist(load_filenames(where))
  
  dname <- vector()
  for(k in 1:length(run_files)){
    dname[k] <- unlist(strsplit(run_files[k], split = ".", fixed = T))[[1]]
  }
  dname = gsub("SOILWATER","SOILWATER_10CM",dname)
  
  get_vars <- dname[c(grep(pattern = "DSTFLXT", dname),grep(pattern = "DSTDEP", dname),grep(pattern = "SOILWATER", dname))]
  get_files <- run_files[c(grep(pattern = "DSTFLXT", run_files),grep(pattern = "DSTDEP", run_files),grep(pattern = "SOILWATER", run_files))]
  
  for(j in 1:length(get_vars)){
    print(j)
    ncin <- nc_open(get_files[j])
    var <- ncvar_get(ncin,get_vars[j])
    SW_array <- var[xmin:xmax,ymin:ymax,]
    
    for(m in 1:nsteps){
      SW_av[m] <- mean(SW_array[,,m], na.rm = T)
    }
    
    all.data[[j]] <- SW_av
    
  }
  
  names(all.data) <- get_vars
  all.data_[[i]] <- all.data
  
}
names(all.data_) <- all.ens

ann.data <- list()
for(i in 1:length(all.ens)){
  
  all.data_[[i]]$Years  <- Years
  ann.data[[i]] <- aggregate(. ~ Years, FUN= sum, all.data_[[i]][-grep("SOILWATER_10CM",names(all.data_[[i]]))])  
  ann.data[[i]]$SOILWATER_10CM <- unlist(aggregate(. ~ Years, FUN= mean, all.data_[[i]][grep(c("SOILWATER_10CM","Years"),names(all.data_[[i]]))])[2])
  
}  
names(ann.data) <- all.ens

DSTDEP <- matrix(data = NA, nrow = nrow(ann.data[[1]]), ncol = length(all.ens))
SOILW <- matrix(data = NA, nrow = nrow(ann.data[[1]]), ncol = length(all.ens))
for(i in 1:length(all.ens)){
  
  DSTDEP[,i] <- ann.data[[i]]$DSTDEP
  SOILW[,i] <- ann.data[[i]]$SOILWATER_10CM
  
}

write.csv(DSTDEP, "LME_ann_dep_sum.csv")
write.csv(SOILW, "LME_ann_soilw.csv")

##############################################
## Defining megadroughts

# Routson et al (2016): PDSI averaged and smoothed with 50 yr spline with intervals that exceed 0.2 PDSI units below the regional mean. Will use soil moisture instead of computing PDSI.
# Ault and St George (2018): soil moisture index 0.5 standard deviation below the long term mean for 35 years or more 

# Method: using a 35 year moving average on soil moisture (from both sides), find periods of 35 years 
# Second method: using a 11 year moving average on soil moisture from both sides, then find periods of 35 years or longer

duration = 35
mav <- function(x,n=duration){stats::filter(x,rep(1/n,n), sides=2)}

sm_dat <- list()
ann.data.35 <- list()

for(i in 1:length(all.ens)){
  for(j in 2:length(names(ann.data[[1]]))){
    
    sm_dat[[j-1]] <- as.vector(mav(ann.data[[i]][j]))
  }
  names(sm_dat) <- names(ann.data[[1]])[-1]
  ann.data.35[[i]] <- sm_dat
  ann.data.35[[i]]$Years  <- seq(850,2005,1)
  
}
names(ann.data.35) <- all.ens

####
# Calculating soil stats

soil.stats <- matrix(data = NA, nrow = 4, ncol = length(all.ens))
for(i in 1:length(all.ens)){
  soil.stats[1,i] <- mean(ann.data.35[[i]]$SOILWATER_10CM,na.rm = T)
  soil.stats[2,i] <- sd(ann.data.35[[i]]$SOILWATER_10CM, na.rm = T)/2
  soil.stats[3,i] <- mean(ann.data[[i]]$SOILWATER_10CM,na.rm = T)
  soil.stats[4,i] <- sd(ann.data[[i]]$SOILWATER_10CM, na.rm = T)/2
}

# Calculating wet/dry periods

index.mega <- list()
index.pluvial <- list()
for(i in 1:length(all.ens)){
  index.mega[[i]] <- which(ann.data.35[[i]]$SOILWATER_10CM <= soil.stats[1,i] - soil.stats[2,i])
  index.pluvial[[i]] <- which(ann.data.35[[i]]$SOILWATER_10CM >= soil.stats[1,i] + soil.stats[2,i]) 
}

cumsum.mega <- vector()
flx_during_mega <- list()
dep_during_mega <- list()
soilw_during_mega <- list()

for(i in 1:length(all.ens)){
  cumsum.mega <- cumsum(c(1,abs(index.mega[[i]][-length(index.mega[[i]])] - index.mega[[i]][-1]) >1))
  result <- table(cumsum.mega)
  group <- as.integer(which(result >=35)) 
  want.index <- which(cumsum.mega %in% group)
  inlist.index <- index.mega[[i]][want.index]
  
  flx_during_mega[[i]] <- ann.data[[i]]$DSTFLXT[inlist.index]
  dep_during_mega[[i]] <- ann.data[[i]]$DSTDEP[inlist.index]
  soilw_during_mega[[i]] <- ann.data[[i]]$SOILWATER_10CM[inlist.index]
  
}

flx_during_mega.list <- flx_during_mega
soilw_during_mega.list <- soilw_during_mega
dep_during_mega.list <- dep_during_mega

flx_during_mega <- unlist(flx_during_mega)
soilw_during_mega <- unlist(soilw_during_mega)
dep_during_mega <- unlist(dep_during_mega)

flx_during_pluvial <- list()
dep_during_pluvial <- list()
cumsum.mega <- vector()

for(i in 1:length(all.ens)){
  cumsum.mega <- cumsum(c(1,abs(index.pluvial[[i]][-length(index.pluvial[[i]])] - index.pluvial[[i]][-1]) >1))
  result <- table(cumsum.mega)
  group <- as.integer(which(result >=35)) 
  want.index <- which(cumsum.mega %in% group)
  inlist.index <- index.pluvial[[i]][want.index]
  flx_during_pluvial[[i]] <- ann.data[[i]]$DSTFLXT[inlist.index]
  dep_during_pluvial[[i]] <- ann.data[[i]]$DSTDEP[inlist.index]
}

flx_during_pluvial.list <- flx_during_pluvial
dep_during_pluvial.list <- dep_during_pluvial

flx_during_pluvial <- unlist(flx_during_pluvial)
dep_during_pluvial <- unlist(dep_during_pluvial)

flx_drought_years <- list()
dep_drought_years <- list()
flx_wet_years <- list()
dep_wet_years <- list()

for(i in 1:length(all.ens)){
  
  flx_drought_years[[i]] <- ann.data[[i]]$DSTFLXT[which(ann.data[[i]]$SOILWATER_10CM <= soil.stats[3,i]- soil.stats[4,i])]
  dep_drought_years[[i]] <- ann.data[[i]]$DSTDEP[which(ann.data[[i]]$SOILWATER_10CM <= soil.stats[3,i] - soil.stats[4,i])]
  
  flx_wet_years[[i]] <- ann.data[[i]]$DSTFLXT[which(ann.data[[i]]$SOILWATER_10CM >= soil.stats[3,i] + soil.stats[4,i])]
  dep_wet_years[[i]] <- ann.data[[i]]$DSTDEP[which(ann.data[[i]]$SOILWATER_10CM >= soil.stats[3,i] + soil.stats[4,i])]
}

flx_drought_years.list <- flx_drought_years
dep_drought_years.list <- dep_drought_years
flx_wet_years.list <- flx_wet_years
dep_wet_years.list <- dep_wet_years

flx_drought_years <- unlist(flx_drought_years)
dep_drought_years <- unlist(dep_drought_years)

flx_wet_years <- unlist(flx_wet_years)
dep_wet_years <- unlist(dep_wet_years)

############
# Plotting density plots

df <- data.frame(Dust_emission = c(flx_during_mega*31536000000,flx_during_pluvial*31536000000), Dataset = c(rep("Megadrought",length(flx_during_mega)), rep("Megapluvial", length(flx_during_pluvial)))) # converted to g/m2/yr from kg/m2/s
# to convert to g/m2/yr ... 31536000 s in a yr (used to be *1000 to get mg/m2/s)

df.dep <- data.frame(Dust_deposition = c(dep_during_mega*31536000000,dep_during_pluvial*31536000000), Dataset = c(rep("Megadrought",length(dep_during_mega)), rep("Megapluvial", length(dep_during_pluvial))))
# converted to g/m2/yr

df.indiv.yr <- data.frame(Dust_emission = c(flx_drought_years*31536000000, flx_wet_years*31536000000), Dataset = c(rep("Drought",length(flx_drought_years)),rep("Wet", length(flx_wet_years))))# converted to g/m2/yr

df.dep.indiv.yr <- data.frame(Dust_deposition = c(dep_drought_years*31536000000, dep_wet_years*31536000000), Dataset = c(rep("Drought",length(dep_drought_years)),rep("Wet", length(dep_wet_years))))# converted to g/m2/yr

#####
# Plots

Emission.mega.plot <- ggplot(df, aes(x = Dust_emission))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  scale_x_continuous(limits = c(0,200))+
  #scale_y_continuous(limits = c(0,0.05))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Emission.mega.plot

Emission.yr.plot <- ggplot(df.indiv.yr, aes(x = Dust_emission))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual emissions (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,200))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Emission.yr.plot

Emissions.plot <- cowplot::plot_grid(Emission.mega.plot, Emission.yr.plot, ncol = 2, labels = c("(a)", "(b)"))

setwd("/projects/pd_lab/sha59/LME")
ggsave(filename = paste0("LME_emissions_gm2yr_wet_dry_only_",Sys.Date(),".pdf"), plot = Emissions.plot)

Depo.mega.plot <- ggplot(df.dep, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Depo.mega.plot

Depo.yr.plot <- ggplot(df.dep.indiv.yr, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Depo.yr.plot

Depo.plot <- cowplot::plot_grid(Depo.mega.plot, Depo.yr.plot, ncol = 2, labels = c("(c)", "(d)"))

ggsave(filename = paste0("LME_deposition_plot_gm2yr_wet_dry_only_",Sys.Date(),".pdf"), plot = Depo.plot)

######
# Plot separately

rep.times.dep_during_mega <- vector()
rep.times.dep_during_pluvial <- vector()
rep.times.dep_drought_years <- vector()
rep.times.dep_wet_years  <- vector()

rep.times.flx_during_mega <- vector()
rep.times.flx_during_pluvial <- vector()
rep.times.flx_drought_years <- vector()
rep.times.flx_wet_years  <- vector()

for(i in 1: length(all.ens.names)){
  rep.times.dep_during_mega[i] <- length(dep_during_mega.list[[i]])
  rep.times.dep_during_pluvial[i] <- length(dep_during_pluvial.list[[i]])
  rep.times.dep_drought_years[i]  <- length(dep_drought_years.list[[i]])
  rep.times.dep_wet_years[i]  <- length(dep_wet_years.list[[i]])
  
  rep.times.flx_during_mega[i] <- length(flx_during_mega.list[[i]])
  rep.times.flx_during_pluvial[i] <- length(flx_during_pluvial.list[[i]])
  rep.times.flx_drought_years[i]  <- length(flx_drought_years.list[[i]])
  rep.times.flx_wet_years[i]  <- length(flx_wet_years.list[[i]])
}


dep.ens.mega <- data.frame(Dust_deposition = c(unlist(dep_during_mega.list)*31536000000, unlist(dep_during_pluvial.list)*31536000000), Ensemble = c(rep("1", rep.times.dep_during_mega[1]), rep("2", rep.times.dep_during_mega[2]), rep("3", rep.times.dep_during_mega[3]), rep("4", rep.times.dep_during_mega[4]), rep("5", rep.times.dep_during_mega[5]), rep("6", rep.times.dep_during_mega[6]), rep("7", rep.times.dep_during_mega[7]), rep("8",rep.times.dep_during_mega[8]), rep("9", rep.times.dep_during_mega[9]), rep("10", rep.times.dep_during_mega[10]),rep("1", rep.times.dep_during_pluvial[1]), rep("2", rep.times.dep_during_pluvial[2]), rep("3", rep.times.dep_during_pluvial[3]), rep("4", rep.times.dep_during_pluvial[4]), rep("5", rep.times.dep_during_pluvial[5]), rep("6", rep.times.dep_during_pluvial[6]), rep("7", rep.times.dep_during_pluvial[7]), rep("8",rep.times.dep_during_pluvial[8]), rep("9", rep.times.dep_during_pluvial[9]), rep("10", rep.times.dep_during_pluvial[10])), Dataset = c(rep("Megadrought",length(dep_during_mega)), rep("Megapluvial", length(dep_during_pluvial))))

dep.ens.indiv <- data.frame(Dust_deposition = c(unlist(dep_drought_years.list)*31536000000, unlist(dep_wet_years.list)*31536000000), Ensemble = c(rep("1", rep.times.dep_drought_years[1]), rep("2", rep.times.dep_drought_years[2]), rep("3", rep.times.dep_drought_years[3]), rep("4", rep.times.dep_drought_years[4]), rep("5", rep.times.dep_drought_years[5]), rep("6", rep.times.dep_drought_years[6]), rep("7", rep.times.dep_drought_years[7]), rep("8",rep.times.dep_drought_years[8]), rep("9", rep.times.dep_drought_years[9]), rep("10", rep.times.dep_drought_years[10]),rep("1", rep.times.dep_wet_years[1]), rep("2", rep.times.dep_wet_years[2]), rep("3", rep.times.dep_wet_years[3]), rep("4", rep.times.dep_wet_years[4]), rep("5", rep.times.dep_wet_years[5]), rep("6", rep.times.dep_wet_years[6]), rep("7", rep.times.dep_wet_years[7]), rep("8",rep.times.dep_wet_years[8]), rep("9", rep.times.dep_wet_years[9]), rep("10", rep.times.dep_wet_years[10])), Dataset = c(rep("Drought",length(dep_drought_years)), rep("Wet", length(dep_wet_years))))

Depo.mega.ens.plot <- ggplot(dep.ens.mega, aes(Dust_deposition, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.65),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Depo.mega.ens.plot

Depo.indiv.ens.plot <- ggplot(dep.ens.indiv, aes(Dust_deposition, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.65),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Depo.indiv.ens.plot

flx.ens.mega <- data.frame(Dust_emission = c(unlist(flx_during_mega.list)*31536000000, unlist(flx_during_pluvial.list)*31536000000), Ensemble = c(rep("1", rep.times.flx_during_mega[1]), rep("2", rep.times.flx_during_mega[2]), rep("3", rep.times.flx_during_mega[3]), rep("4", rep.times.flx_during_mega[4]), rep("5", rep.times.flx_during_mega[5]), rep("6", rep.times.flx_during_mega[6]), rep("7", rep.times.flx_during_mega[7]), rep("8",rep.times.flx_during_mega[8]), rep("9", rep.times.flx_during_mega[9]), rep("10", rep.times.flx_during_mega[10]),rep("1", rep.times.flx_during_pluvial[1]), rep("2", rep.times.flx_during_pluvial[2]), rep("3", rep.times.flx_during_pluvial[3]), rep("4", rep.times.flx_during_pluvial[4]), rep("5", rep.times.flx_during_pluvial[5]), rep("6", rep.times.flx_during_pluvial[6]), rep("7", rep.times.flx_during_pluvial[7]), rep("8",rep.times.flx_during_pluvial[8]), rep("9", rep.times.flx_during_pluvial[9]), rep("10", rep.times.flx_during_pluvial[10])), Dataset = c(rep("Megadrought",length(flx_during_mega)), rep("Megapluvial", length(flx_during_pluvial))))

flx.ens.indiv <- data.frame(Dust_emission = c(unlist(flx_drought_years.list)*31536000000, unlist(flx_wet_years.list)*31536000000), Ensemble = c(rep("1", rep.times.flx_drought_years[1]), rep("2", rep.times.flx_drought_years[2]), rep("3", rep.times.flx_drought_years[3]), rep("4", rep.times.flx_drought_years[4]), rep("5", rep.times.flx_drought_years[5]), rep("6", rep.times.flx_drought_years[6]), rep("7", rep.times.flx_drought_years[7]), rep("8",rep.times.flx_drought_years[8]), rep("9", rep.times.flx_drought_years[9]), rep("10", rep.times.flx_drought_years[10]),rep("1", rep.times.flx_wet_years[1]), rep("2", rep.times.flx_wet_years[2]), rep("3", rep.times.flx_wet_years[3]), rep("4", rep.times.flx_wet_years[4]), rep("5", rep.times.flx_wet_years[5]), rep("6", rep.times.flx_wet_years[6]), rep("7", rep.times.flx_wet_years[7]), rep("8",rep.times.flx_wet_years[8]), rep("9", rep.times.flx_wet_years[9]), rep("10", rep.times.flx_wet_years[10])), Dataset = c(rep("Drought",length(flx_drought_years)), rep("Wet", length(flx_wet_years))))

Flx.indiv.ens.plot <- ggplot(flx.ens.indiv, aes(Dust_emission, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual emission (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,200))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.65),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Flx.indiv.ens.plot

Flx.mega.ens.plot <- ggplot(flx.ens.mega, aes(Dust_emission, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual emission (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,200))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.65),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Flx.mega.ens.plot

All.ens.plot <- cowplot::plot_grid(Depo.mega.ens.plot, Depo.indiv.ens.plot,Flx.mega.ens.plot,Flx.indiv.ens.plot, ncol = 2, labels = c("(a)", "(b)", "(c)", "(d)"))

setwd("/projects/pd_lab/sha59/LME")

ggsave(filename = paste0("LME_all_ens_plot_",Sys.Date(),".pdf"), plot = All.ens.plot)


#######
# Stat analysis

mega.stats <- df %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_emission),
            grp.IQR = IQR(Dust_emission),
            grp.median = median(Dust_emission))
mega.stats

indiv.stats <- df.indiv.yr %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_emission),
            grp.IQR = IQR(Dust_emission),
            grp.median = median(Dust_emission))
indiv.stats

var.test(flx_during_mega, flx_during_pluvial) # if p > 0.05 then there is no significant difference in variance
var.test(dep_during_mega, dep_during_pluvial)

t.test(flx_during_mega, flx_during_pluvial, var.equal = T, alternative = "greater")
t.test(flx_drought_years, flx_wet_years, var.equal = T, alternative = "greater")

t.test(dep_during_mega, dep_during_pluvial, var.equal = T, alternative = "greater")
t.test(dep_drought_years, dep_wet_years, var.equal = T, alternative = "greater")

t.test(flx_during_mega*31536000, flx_during_pluvial*31536000, var.equal = F, alternative = "greater") # if p <0.05 then the median is significantly different
#t.test(drought_years*31536000,nondrought_years*31536000, alternative = "greater")

flx.indiv.ens.stats <- matrix(data = NA, ncol = 2, nrow = 10)
for(i in 1:10){
p <- ks.test(flx_drought_years.list[[i]], flx_wet_years.list[[i]])$p.value
D <- ks.test(flx_drought_years.list[[i]], flx_wet_years.list[[i]])$statistic

flx.indiv.ens.stats [i,1] <- D
flx.indiv.ens.stats [i,2] <- p
}

flx.mega.ens.stats <- matrix(data = NA, ncol = 2, nrow = 10)
for(i in c(1,2,3,5,6,7,8,9,10)){
  p <- ks.test(flx_during_mega.list[[i]], flx_during_pluvial.list[[i]])$p.value
  D <- ks.test(flx_during_mega.list[[i]], flx_during_pluvial.list[[i]])$statistic
  
  flx.mega.ens.stats [i,1] <- D
  flx.mega.ens.stats [i,2] <- p
}

dep.mega.ens.stats <- matrix(data = NA, ncol = 2, nrow = 10)
for(i in c(1,2,3,5,6,7,8,9,10)){
  p <- ks.test(dep_during_mega.list[[i]], dep_during_pluvial.list[[i]])$p.value
  D <- ks.test(dep_during_mega.list[[i]], dep_during_pluvial.list[[i]])$statistic
  
  dep.mega.ens.stats [i,1] <- D
  dep.mega.ens.stats [i,2] <- p
}

one <- hist(flx_during_mega.list[[2]])
two <- hist(flx_during_pluvial.list[[2]])
plot(one, col=rgb(0,0,1,1/4), add = T)
plot(two, col=rgb(1,0,0,1/4))


#####
# Correlations between soilwater_10cm and dust emission in each ensemble member

correlation <- matrix(data = NA, ncol = length(all.ens), nrow = 4)
for(i in 1:length(all.ens)){
  
  correlation[1,i]  <- as.numeric(round(cor.test(x = ann.data.35[[i]]$DSTFLXT,y = ann.data.35[[i]]$SOILWATER_10CM, method = "kendall")$estimate, 3))
  correlation[2,i]  <- cor.test(x = ann.data.35[[i]]$DSTFLXT,y = ann.data.35[[i]]$SOILWATER_10CM, method = "kendall")$p.value
  correlation[3,i]  <- as.numeric(round(cor.test(x = ann.data.35[[i]]$DSTDEP,y = ann.data.35[[i]]$SOILWATER_10CM, method = "kendall")$estimate, 3))
  correlation[4,i]  <- cor.test(x = ann.data.35[[i]]$DSTDEP,y = ann.data.35[[i]]$SOILWATER_10CM, method = "kendall")$p.value
  
}
colnames(correlation) <- all.ens
correlation <- as.data.frame(correlation)
correlation$Parameters <- c("Kendall corr", "pvalue")
correlation$Variable <- c("FLXT","FLXT","DEP","DEP")
correlation

range(correlation[1,1:10])
range(correlation[3,1:10])

############################################################3
## Dust in the paleo-records

# Data needed include: FourCornersPMDI.xlsx and lipd file BlueLake.Routson.2018.lpd 

# Dust data 

B <- readLipd(path = "/home/sha59")
core <- B

core <- runBacon(core, labIDVar = NULL, age14CVar = "age14c", ageVar = "age", ageUncertaintyVar = "ageUncertainty", depthVar = "depth", reservoirAge14CVar = NULL, reservoirAge14CUncertaintyVar = NULL, rejectedAgesVar = NULL, baconThick = 3)

start = -70
end = 4500
max.ens = 1000

core <- mapAgeEnsembleToPaleoData(core)

ae  <- selectData(L = core, varName = "ageensemble")
DBD  <- selectData(L = core, varName = "dbd")
DBDm  <- matrix(DBD$values, nrow = length(DBD$values), ncol = max.ens)
depth <- selectData(L = core, varName = "depth")
df  <- selectData(L = core, varName = "dust")
ae.CE <- convertBP2AD(ae$values)

ageDiff  <- diff(ae$values)
depthDiff  <- diff(depth$values)
replicated.depth.matrix = matrix(data = depthDiff, nrow = length(depthDiff), ncol = ncol(ageDiff))
srEns = replicated.depth.matrix/ageDiff
empty.values  <- as.matrix(t(rep(NA, max.ens)))
srEns  <- rbind(srEns, empty.values)

MAR  <- DBDm * srEns
DMAR <- df$values * MAR # this contains the dust measurements from Blue Lake

# Hydroclimate data

PMDI  <- as.data.frame(read_xlsx("/home/sha59/FourCornersPMDI.xlsx"))
PMDI.df <- data.frame(Time = PMDI$'Time (CE)', PMDI = PMDI$'Four Courners PMDI (35°N to 37°N, -107°E to -112°E')

# Analysis

PMDI.35 <- mav(PMDI.df$PMDI)
DMAR.35 <- matrix(data = NA, nrow = nrow(DMAR), ncol = max.ens)
for(i in 1: max.ens){
  DMAR.35[,i] <- mav(DMAR[,i])
}

PMDI.stats <- vector()
PMDI.stats[1] <- mean(PMDI.35, na.rm = T)
PMDI.stats[2] <- sd(PMDI.35, na.rm = T)/2
PMDI.stats[3] <- mean(PMDI.df$PMDI, na.rm = T)
PMDI.stats[4] <- sd(PMDI.df$PMDI, na.rm = T)/2

index.mega <- which(PMDI.35 <= PMDI.stats[1] - PMDI.stats[2])
cumsum.mega <- cumsum(c(1,abs(index.mega[-length(index.mega)] - index.mega[-1]) >1))
result  <- table(cumsum.mega)
group  <- as.integer(which(result >= 35))
want.index <- which(cumsum.mega %in% group)
inlist.index <- index.mega[want.index]
years.mega <- PMDI.df$Time[inlist.index]

age.index.mega <- list()
dust.mega  <- list()

for(i in 1:max.ens){
  age.index.mega[[i]] = match(years.mega, sort(round(ae.CE[complete.cases(ae.CE),i])))
  c  <- age.index.mega[[i]]
  c  <- c[!is.na(c)]
  dust.mega[[i]] <- DMAR.35[c,i]
  
}

index.pluvial <- which(PMDI.35 >= PMDI.stats[1] + PMDI.stats[2])
cumsum.pluvial <- cumsum(c(1,abs(index.pluvial[-length(index.pluvial)] - index.pluvial[-1]) >1))
result  <- table(cumsum.pluvial)
group  <- as.integer(which(result >= 35))
want.index <- which(cumsum.pluvial %in% group)
inlist.index <- index.pluvial[want.index]
years.pluvial <- PMDI.df$Time[inlist.index]

age.index.pluvial <- list()
dust.pluvial  <- list()

for(i in 1:max.ens){
  age.index.pluvial[[i]] = match(years.pluvial, sort(round(ae.CE[complete.cases(ae.CE),i])))
  c  <- age.index.pluvial[[i]]
  c  <- c[!is.na(c)]
  dust.pluvial[[i]] <- DMAR.35[c,i]
  
}

index.dry <- which(PMDI.df$PMDI <= PMDI.stats[3] - PMDI.stats[4])
index.wet <- which(PMDI.df$PMDI >= PMDI.stats[3] + PMDI.stats[4])
years.dry <- PMDI.df$Time[index.dry]
years.wet <- PMDI.df$Time[index.wet]

age.index.dry <- list()
age.index.wet <- list()
dust.dry <- list()
dust.wet <- list()

for(i in 1:max.ens){
  age.index.dry[[i]] <- match(years.dry, sort(round(ae.CE[complete.cases(ae.CE),i])))
  a <- age.index.dry[[i]]
  a[a == NaN]  <- NA
  a <- a[!is.na(a)]
  dust.dry[[i]] <- DMAR[a,i]
}

for(i in 1:max.ens){
  age.index.wet[[i]] <- match(years.wet, sort(round(ae.CE[complete.cases(ae.CE),i])))
  a <- age.index.wet[[i]]
  a[a == NaN]  <- NA
  a <- a[!is.na(a)]
  dust.wet[[i]] <- DMAR[a,i]
}

Depo.mega <- data.frame(Dust_deposition = c(unlist(dust.mega), unlist(dust.pluvial))*1000, Dataset = c(rep("Megadrought", length(unlist(dust.mega))), rep("Megapluvial", length(unlist(dust.pluvial)))))

Depo.yr <- data.frame(Dust_deposition = c(unlist(dust.dry), unlist(dust.wet))*1000, Dataset = c(rep("Dry", length(unlist(dust.dry))), rep("Wet", length(unlist(dust.wet)))))

unloadNamespace("ggplot2")
library(ggplot2)

Blue.depo.mega.plot <- ggplot(Depo.mega, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Blue.depo.mega.plot

Blue.depo.yr.plot <- ggplot(Depo.yr, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Blue.depo.yr.plot

Blue.depo.plot <- cowplot::plot_grid(Blue.depo.mega.plot, Blue.depo.yr.plot, ncol = 2, labels = c("(e)", "(f)"))

setwd("/projects/pd_lab/sha59/LME")
ggsave(filename = paste0("Blue_deposition_plot_gm2yr_wet_dry_only_",Sys.Date(),".pdf"), plot = Blue.depo.plot)

Final_figure <- cowplot::plot_grid(Emissions.plot, Depo.plot, Blue.depo.plot, ncol = 1, nrow = 3)
ggsave(filename = paste0("Fig5",Sys.Date(),".pdf"), plot = Final_figure, width = 7)

# Stat analysis

Blue.mega.stats <- Depo.mega %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition),
            grp.IQR = IQR(Dust_deposition),
            grp.median = median(Dust_deposition))
Blue.mega.stats

Blue.indiv.stats <- Depo.yr %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition, na.rm = T),
            grp.IQR = IQR(Dust_deposition, na.rm = T),
            grp.median = median(Dust_deposition, na.rm = T))
Blue.indiv.stats

var.test(unlist(dust.mega), unlist(dust.pluvial)) # if p > 0.05 then there is no significant difference in variance

t.test(unlist(dust.mega), unlist(dust.pluvial), var.equal = T, alternative = "greater")
# if p <0.05 then the median is significantly different

dep.unc.mega.KS <- matrix(data = NA, ncol = 2, nrow = 10)
for(i in c(1,2,3,4,5,6,7,8,9,10)){
  p <- ks.test(dust.mega[[i]], dust.pluvial[[i]])$p.value
  D <- ks.test(dust.mega[[i]], dust.pluvial[[i]])$statistic
  
  dep.unc.mega.KS [i,1] <- D
  dep.unc.mega.KS [i,2] <- p
}

dep.unc.indiv.KS <- matrix(data = NA, ncol = 2, nrow = 10)
for(i in c(1,2,3,4,5,6,7,8,9,10)){
  p <- ks.test(dust.dry[[i]], dust.wet[[i]])$p.value
  D <- ks.test(dust.dry[[i]], dust.wet[[i]])$statistic
  
  dep.unc.indiv.KS [i,1] <- D
  dep.unc.indiv.KS [i,2] <- p
}

one <- hist(flx_during_mega.list[[2]])
two <- hist(flx_during_pluvial.list[[2]])
plot(one, col=rgb(0,0,1,1/4), add = T)
plot(two, col=rgb(1,0,0,1/4))

flx_during_mega.list, flx_during_pluvial.list, dep_drought_years.list, dep_wet_years.list, dep_during_mega.list, dep_during_pluvial.list
#########################################
## Testing age uncertainty

# Data needed include: tree_PMDI.csv, Blue_DMAR_ens.csv, Blue_age_ens.csv, LME_ann_dep_sum.csv, LME_ann_soilw.csv

setwd("/home/sha59")
PMDI <- read.csv("tree_PMDI.csv")
Blue_DMAR <- read.csv("Blue_DMAR_ens.csv")
Blue_age <- read.csv("Blue_age_ens.csv")
LME_dust <- read.csv("LME_ann_dep_sum.csv")
LME_soilw <- read.csv("LME_ann_soilw.csv")

LME <- data.frame(Year = seq(850,2005,1), LME_dust[,-1], LME_soilw[,-1])
LME$mean.dust <- rowMeans(LME[,2:11])

mean.age <- round(rowMeans(Blue_age[,-1]))
Blue_age[,1] <- mean.age
names(Blue_age)[1] <- "mean.age"
Blue_age_abs <- round(Blue_age)
Blue_age_abs_com <- Blue_age_abs[min(which(Blue_age_abs$mean.age <= 2005)):max(which(Blue_age_abs$mean.age >= 850)),]
Blue_DMAR_com <- Blue_DMAR[min(which(Blue_age_abs$mean.age <= 2005)):max(which(Blue_age_abs$mean.age >= 850)),]
Blue_DMAR_com_sorted <- Blue_DMAR_com[order(-Blue_DMAR_com$X),]
Blue_age_abs_com_sorted <- Blue_age_abs_com[order(Blue_age_abs_com$mean.age),]

LME_binned <- list(list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)), list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)),list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)),list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)),list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)),list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)),list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)), list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)), list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)), list(Age =matrix(data = NA, nrow = 581, ncol = 100), Dust = matrix(data = NA, nrow = 581, ncol = 100)))

names(LME_binned) <- c("LME2","LME3", "LME4","LME5", "LME6", "LME7", "LME8", "LME9", "LME10", "LME11")

for(j in 1:99){
  
  calc <- binEns(time = LME$Year, values = LME[,2], binvec = Blue_age_abs_com_sorted[,j])
  LME_binned[[1]]$Age[,j] = calc$time
  LME_binned[[1]]$Dust[,j] = calc$matrix
}

a <- c(1,99,100,199,200,299,300,399,400,499,500,599,600,699,700,799,800,899,900,999,1000)

for(j in 1:100){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,2], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[1]]$Age[,i] = calc$time
    LME_binned[[1]]$Dust[,i] = calc$matrix
  }
}

for(j in 101:200){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,3], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[2]]$Age[,i] = calc$time
    LME_binned[[2]]$Dust[,i] = calc$matrix
  }
}

for(j in 201:300){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,4], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[3]]$Age[,i] = calc$time
    LME_binned[[3]]$Dust[,i] = calc$matrix
  }
}

for(j in 301:400){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,5], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[4]]$Age[,i] = calc$time
    LME_binned[[4]]$Dust[,i] = calc$matrix
  }
}

for(j in 401:500){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,6], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[5]]$Age[,i] = calc$time
    LME_binned[[5]]$Dust[,i] = calc$matrix
  }
}

for(j in 501:600){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,7], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[6]]$Age[,i] = calc$time
    LME_binned[[6]]$Dust[,i] = calc$matrix
  }
}

for(j in 601:700){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,8], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[7]]$Age[,i] = calc$time
    LME_binned[[7]]$Dust[,i] = calc$matrix
  }
}

for(j in 701:800){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,9], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[8]]$Age[,i] = calc$time
    LME_binned[[8]]$Dust[,i] = calc$matrix
  }
}
for(j in 801:900){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,10], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[9]]$Age[,i] = calc$time
    LME_binned[[9]]$Dust[,i] = calc$matrix
  }
}
for(j in 901:1000){
  for(i in 1:100){
    print(j)
    calc <- binEns(time = LME$Year, values = LME[,11], binvec = Blue_age_abs_com_sorted[,j])
    LME_binned[[10]]$Age[,i] = calc$time
    LME_binned[[10]]$Dust[,i] = calc$matrix
  }
}

newname <- paste0("SOILWATER_", all.ens)

ann.soil.35 <- matrix(data = NA, nrow = nrow(LME), ncol = length(all.ens))
for(i in 1:length(all.ens)){
  ann.soil.35[,i] <- mav(LME[,i+11])
}
ann.soil.35 <- as.data.frame(ann.soil.35)
names(ann.soil.35)[1:length(all.ens)] <- newname[1:length(all.ens)]
ann.soil.35$Years <- seq(850,2005,1)

ann.dep.35 <- LME_binned
for(i in 1:length(all.ens)){
  for(j in 1:100){
    ann.dep.35[[i]]$Dust[,j] <- mav(LME_binned[[i]]$Dust[,j])
  }
}

ann.soil.stats <- matrix(data = NA, nrow = 4, ncol = length(all.ens))
for(i in 1:length(all.ens)){
  ann.soil.stats[1,i] <- mean(ann.soil.35[,i],na.rm = T)
  ann.soil.stats[2,i] <- sd(ann.soil.35[,i], na.rm = T)/2
  ann.soil.stats[3,i] <- mean(LME[,i+11],na.rm = T)
  ann.soil.stats[4,i] <- sd(LME[,i+11], na.rm = T)/2
  
}

index.mega <- list()
index.pluvial <- list()
for(i in 1:length(all.ens)){
  index.mega[[i]] <- ann.soil.35[which(ann.soil.35[,i] <= ann.soil.stats[1,i]- ann.soil.stats[2,i]),11]
  index.pluvial[[i]] <- ann.soil.35[which(ann.soil.35[,i] >= ann.soil.stats[1,i] + ann.soil.stats[2,i]),11] 
} # now spits out the years rather than the index, necessary for later

cumsum.mega <- vector()
dep_during_mega_unc <- list()
soilw_during_mega_unc <- list()

for(i in 1:length(all.ens)){
  for(j in 1:100){
    print(i)
    cumsum.mega <- cumsum(c(1,abs(index.mega[[i]][-length(index.mega[[i]])] - index.mega[[i]][-1]) >1))
    result <- table(cumsum.mega)
    group <- as.integer(which(result >=35)) 
    want.index <- which(cumsum.mega %in% group)
    inlist.index <- index.mega[[i]][want.index]
    
    dep_during_mega_unc[[i]] <- na.omit(ann.dep.35[[i]]$Dust[which(inlist.index %in% ann.dep.35[[i]]$Age[,j]),j])
    soilw_during_mega_unc[[i]] <- na.omit(ann.soil.35[inlist.index,i])
    
  }
}
soilw_during_mega_unc <- unlist(soilw_during_mega_unc)
dep_during_mega_unc <- unlist(dep_during_mega_unc)

dep_during_pluvial_unc <- list()
cumsum.pluvial <- vector()

for(i in 1:length(all.ens)){
  cumsum.pluvial <- cumsum(c(1,abs(index.pluvial[[i]][-length(index.pluvial[[i]])] - index.pluvial[[i]][-1]) >1))
  result <- table(cumsum.pluvial)
  group <- as.integer(which(result >=35)) 
  want.index <- which(cumsum.pluvial %in% group)
  inlist.index <- index.pluvial[[i]][want.index]
  
  dep_during_pluvial_unc[[i]] <- na.omit(ann.dep.35[[i]]$Dust[which(inlist.index %in% ann.dep.35[[i]]$Age[,j]),j])
}
dep_during_pluvial_unc <- unlist(dep_during_pluvial_unc)

index.drought <- list()
index.wet <- list()
for(i in 1:length(all.ens)){
  index.drought[[i]] <- LME[which(LME[,i+11] <= ann.soil.stats[3,i]- ann.soil.stats[4,i]),1]
  index.wet[[i]] <- LME[which(LME[,i+11] >= ann.soil.stats[3,i] + ann.soil.stats[4,i]),1] 
}

dep_drought_years_unc <- list()
dep_wet_years_unc <- list()

for(i in 1:length(all.ens)){
  for(j in 1:100){
    dep_drought_years_unc[[i]] <- na.omit(LME_binned[[i]]$Dust[which(index.drought[[i]] %in% LME_binned[[i]]$Age[,j]),j])
    dep_wet_years_unc[[i]] <- na.omit(LME_binned[[i]]$Dust[which(index.wet[[i]] %in% LME_binned[[i]]$Age[,j]),j])
    
  }
}

dep_drought_years_unc <- unlist(dep_drought_years_unc)
dep_wet_years_unc <- unlist(dep_wet_years_unc)

df.dep.unc <- data.frame(Dust_deposition = c(dep_during_mega_unc*31536000000,dep_during_pluvial_unc*31536000000), Dataset = c(rep("Megadrought",length(dep_during_mega_unc)), rep("Megapluvial", length(dep_during_pluvial_unc))))

df.dep.indiv.yr.unc <- data.frame(Dust_deposition = c(dep_drought_years_unc*31536000000, dep_wet_years_unc*31536000000), Dataset = c(rep("Drought",length(dep_drought_years_unc)),rep("Wet", length(dep_wet_years_unc))))

Unc.mega.plot <- ggplot(df.dep.unc, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  xlab("")+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Unc.mega.plot

Unc.yr.plot <- ggplot(df.dep.indiv.yr.unc, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Unc.yr.plot

Unc.dens.plot <- cowplot::plot_grid(Unc.mega.plot, Unc.yr.plot, nrow = 2, labels = c("(a)", "(b)"))
setwd("/projects/pd_lab/sha59/LME")
ggsave(filename = "LME_dep_density_age_uncertainty.pdf", plot = Unc.dens.plot)
