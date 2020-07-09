########################################################################################
# Project: Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison
# Prepared by: S. Arcusa
# Version: 2
# Date: 07-09-2020
# Part IV: Investigating individual members, Blue Lake, and LME with age uncertainty
#######################################################################################

library(ncdf4)
library(tidyr)
library(plyr)
library(dplyr)
library(abind)
library(ggplot2)
library(geoChronR)
library(zoo)
library(ggpubr)
library(dabestr)
library(purrr)


# Functions

density.graphs <- function(df, na.rm = T){
  for(i in 1:length(all.ens)){
    plot  <- ggplot(subset(df[[i]]), aes(x = Value))+
      geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
      facet_wrap(Variable ~., scales = "free")+
      scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
      xlab("")+
      theme_bw()+
      theme(legend.title = element_blank(), legend.position = c(0.75,0.10),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      ggtitle(paste0("Ensemble member ",all.ens[i]))
    
    ggsave(plot, filename = paste0(results,"LME_",all.ens[i],"_",duration,"_",Sys.Date(),".pdf"))
  }
}

ts.graphs <- function(df, na.rm = T){
  for(i in 1:length(all.ens)){
    
    plot  <- ggplot(subset(df[[i]]), aes(x = Years, y = value))+
      geom_line(col = "darkgrey")+
      geom_path(aes(color = factor(Period),group = 1))+
      scale_color_manual(values = c( "MD" = "red", "MP" = "blue"))+
      facet_grid(key ~., scales = "free")+
      theme_bw()+
      theme(legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      ggtitle(paste0("Ensemble member ",all.ens[i]))
    
    ggsave(plot, filename = paste0(results,"LME_ts_ann_",all.ens[i],"_",duration,"_",Sys.Date(),".pdf"))
  }
}

#### Data

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

load_filenames <- function(where){
  
  setwd(where)
  ncname <- list.files(pattern=".nc")
  filenames <- list(ncname)
  
}

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
  
  get_vars <- dname[c(grep(pattern = "DSTFLXT", dname), grep(pattern = "DSTDEP", dname), grep(pattern = "TLAI",dname), grep(pattern = "TSAI",dname), grep(pattern = "SOILWATER_10CM",dname), grep(pattern = "FSNO",dname), grep(pattern =  "TREFHT",dname), grep(pattern =  "TOTPREC",dname), grep(pattern =  "U10", dname))]
  get_files <- run_files[c(grep(pattern = "DSTFLXT", run_files),grep(pattern = "DSTDEP", run_files) ,grep(pattern = "TLAI", run_files),grep(pattern = "TSAI", run_files),grep(pattern = "SOILWATER", run_files), grep(pattern = "FSNO", run_files), grep(pattern = "TREFHT", run_files), grep(pattern = "TOTPREC", run_files), grep(pattern = "U10", run_files))]
  
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

for(i in 1:length(all.ens)){
  all.data_[[i]]$fv <- all.data_[[i]]$TLAI + all.data_[[i]]$TSAI
}

ann.data <- list()
dust <- c("DSTFLXT","DSTDEP")
dustyr <- c("DSTFLXT","DSTDEP","Years")
for(i in 1:length(all.ens)){
  
  all.data_[[i]]$Years  <- Years
  ann.data[[i]] <- aggregate(. ~ Years, FUN= mean, all.data_[[i]][-grep(paste(dust,collapse="|"),names(all.data_[[i]]))])  
  ann.data[[i]]$DSTFLXT <- unlist(aggregate(. ~ Years, FUN= sum, all.data_[[i]][grep(c("DSTFLXT","Years"),names(all.data_[[i]]))])[2])
  ann.data[[i]]$DSTDEP <- unlist(aggregate(. ~ Years, FUN= sum, all.data_[[i]][grep(c("DSTDEP","Years"),names(all.data_[[i]]))])[2])
  
}
names(ann.data) <- all.ens

##### Multi-decadal drought

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

soil.stats <- matrix(data = NA, nrow = 4, ncol = length(all.ens))
for(i in 1:length(all.ens)){
  soil.stats[1,i] <- mean(ann.data.35[[i]]$SOILWATER_10CM,na.rm = T)
  soil.stats[2,i] <- sd(ann.data.35[[i]]$SOILWATER_10CM, na.rm = T)/2
  soil.stats[3,i] <- mean(ann.data[[i]]$SOILWATER_10CM,na.rm = T)
  soil.stats[4,i] <- sd(ann.data[[i]]$SOILWATER_10CM, na.rm = T)/2
}

index.mega <- list()
index.pluvial <- list()

for(i in 1:length(all.ens)){
  index.mega[[i]] <- which(ann.data.35[[i]]$SOILWATER_10CM <= soil.stats[1,i] - soil.stats[2,i])
  index.pluvial[[i]] <- which(ann.data.35[[i]]$SOILWATER_10CM >= soil.stats[1,i] + soil.stats[2,i]) 
 
}

cumsum.mega <- vector()
flx_during_mega <- list()
dep_during_mega <- list()
TLAI_during_mega <- list()
TSAI_during_mega <- list()
fv_during_mega <- list()
T_during_mega  <- list()
P_during_mega <- list()
SW_during_mega <- list()
FSNO_during_mega <- list()
U10_during_mega <- list()
mega.list <- list()
non.mega.list <- list()
dep_during_non_mega <- list()

for(i in 1:length(all.ens)){
  cumsum.mega <- cumsum(c(1,abs(index.mega[[i]][-length(index.mega[[i]])] - index.mega[[i]][-1]) >1))
  result <- table(cumsum.mega)
  group <- as.integer(which(result >=35)) 
  want.index <- which(cumsum.mega %in% group)
  inlist.index <- index.mega[[i]][want.index]
  mega.list[[i]] <- inlist.index
  non.mega.list[[i]] <- seq(1,1156,1)[-mega.list[[i]]]
  
  flx_during_mega[[i]] <- ann.data[[i]]$DSTFLXT[inlist.index]
  dep_during_mega[[i]] <- ann.data[[i]]$DSTDEP[inlist.index]
  TLAI_during_mega[[i]] <- ann.data[[i]]$TLAI[inlist.index]
  TSAI_during_mega[[i]] <- ann.data[[i]]$TSAI[inlist.index]
  fv_during_mega[[i]] <- ann.data[[i]]$fv[inlist.index]
  T_during_mega[[i]] <- ann.data[[i]]$TREFHT[inlist.index]
  P_during_mega[[i]] <- ann.data[[i]]$TOTPREC[inlist.index]
  SW_during_mega[[i]] <- ann.data[[i]]$SOILWATER_10CM[inlist.index]
  FSNO_during_mega[[i]] <- ann.data[[i]]$FSNO[inlist.index]
  U10_during_mega[[i]] <- ann.data[[i]]$U10[inlist.index]
  
  dep_during_non_mega[[i]] <- ann.data[[i]]$DSTDEP[-inlist.index]
  
}

flx_during_mega.list <- flx_during_mega
dep_during_mega.list <- dep_during_mega
SW_during_mega.list <- SW_during_mega
TLAI_during_mega.list <- TLAI_during_mega
TSAI_during_mega.list <- TSAI_during_mega
fv_during_mega.list <- fv_during_mega
T_during_mega.list <- T_during_mega
P_during_mega.list <- P_during_mega
FSNO_during_mega.list <- FSNO_during_mega
U10_during_mega.list <- U10_during_mega
dep_during_non_mega.list <- dep_during_non_mega

flx_during_mega <- unlist(flx_during_mega)
dep_during_mega <- unlist(dep_during_mega)
SW_during_mega <- unlist(SW_during_mega)
TLAI_during_mega <- unlist(TLAI_during_mega)
TSAI_during_mega <- unlist(TSAI_during_mega)
fv_during_mega <- unlist(fv_during_mega)
T_during_mega <- unlist(T_during_mega)
P_during_mega <- unlist(P_during_mega)
FSNO_during_mega <- unlist(FSNO_during_mega)
U10_during_mega <- unlist(U10_during_mega)
dep_during_non_mega <- unlist(dep_during_non_mega)

flx_during_pluvial <- list()
dep_during_pluvial <- list()
TLAI_during_pluvial <- list()
TSAI_during_pluvial <- list()
fv_during_pluvial <- list()
T_during_pluvial  <- list()
P_during_pluvial <- list()
SW_during_pluvial <- list()
FSNO_during_pluvial <- list()
U10_during_pluvial <- list()
cumsum.pluvial <- vector()
pluvial.list <- list()
non.pluvial.list <- list()
dep_during_non_pluvial <- list()

for(i in 1:length(all.ens)){
  cumsum.pluvial <- cumsum(c(1,abs(index.pluvial[[i]][-length(index.pluvial[[i]])] - index.pluvial[[i]][-1]) >1))
  result <- table(cumsum.pluvial)
  group <- as.integer(which(result >=35)) 
  want.index <- which(cumsum.pluvial %in% group)
  inlist.index <- index.pluvial[[i]][want.index]
  pluvial.list[[i]] <- inlist.index
  non.pluvial.list[[i]] <- seq(1,1156,1)[-pluvial.list[[i]]]
  
  flx_during_pluvial[[i]] <- ann.data[[i]]$DSTFLXT[inlist.index]
  dep_during_pluvial[[i]] <- ann.data[[i]]$DSTDEP[inlist.index]
  TLAI_during_pluvial[[i]] <- ann.data[[i]]$TLAI[inlist.index]
  TSAI_during_pluvial[[i]] <- ann.data[[i]]$TSAI[inlist.index]
  fv_during_pluvial[[i]] <- ann.data[[i]]$fv[inlist.index]
  T_during_pluvial[[i]] <- ann.data[[i]]$TREFHT[inlist.index]
  P_during_pluvial[[i]] <- ann.data[[i]]$TOTPREC[inlist.index]
  SW_during_pluvial[[i]] <- ann.data[[i]]$SOILWATER_10CM[inlist.index]
  FSNO_during_pluvial[[i]] <- ann.data[[i]]$FSNO[inlist.index]
  U10_during_pluvial[[i]] <- ann.data[[i]]$U10[inlist.index]
  
  dep_during_non_pluvial[[i]] <- ann.data[[i]]$DSTDEP[-inlist.index]
}

flx_during_pluvial.list <- flx_during_pluvial
dep_during_pluvial.list <- dep_during_pluvial
SW_during_pluvial.list <- SW_during_pluvial
TLAI_during_pluvial.list <- TLAI_during_pluvial
TSAI_during_pluvial.list <- TSAI_during_pluvial
fv_during_pluvial.list <- fv_during_pluvial
T_during_pluvial.list <- T_during_pluvial
P_during_pluvial.list <- P_during_pluvial
FSNO_during_pluvial.list <- FSNO_during_pluvial
U10_during_pluvial.list <- U10_during_pluvial
dep_during_non_pluvial.list <- dep_during_non_pluvial

flx_during_pluvial <- unlist(flx_during_pluvial)
dep_during_pluvial <- unlist(dep_during_pluvial)
SW_during_pluvial <- unlist(SW_during_pluvial)
TLAI_during_pluvial <- unlist(TLAI_during_pluvial)
TSAI_during_pluvial <- unlist(TSAI_during_pluvial)
fv_during_pluvial <- unlist(fv_during_pluvial)
T_during_pluvial <- unlist(T_during_pluvial)
P_during_pluvial <- unlist(P_during_pluvial)
FSNO_during_pluvial <- unlist(FSNO_during_pluvial)
U10_during_pluvial <- unlist(U10_during_pluvial)
dep_during_non_pluvial <- unlist(dep_during_non_pluvial)

dens.df <- list()
for(i in 1:length(all.ens)){
  dens.df[[i]] <- data.frame(
    Value = c(flx_during_mega.list[[i]]*31536000000,flx_during_pluvial.list[[i]]*31536000000,dep_during_mega.list[[i]]*31536000000,dep_during_pluvial.list[[i]] *31536000000,TLAI_during_mega.list[[i]], TLAI_during_pluvial.list[[i]], TSAI_during_mega.list[[i]], TSAI_during_pluvial.list[[i]], T_during_mega.list[[i]],T_during_pluvial.list[[i]], P_during_mega.list[[i]],P_during_pluvial.list[[i]], SW_during_mega.list[[i]], SW_during_pluvial.list[[i]], FSNO_during_mega.list[[i]], FSNO_during_pluvial.list[[i]], U10_during_mega.list[[i]], U10_during_pluvial.list[[i]]), 
    Dataset = c(rep("Megadrought",length(flx_during_mega.list[[i]])), rep("Megapluvial", length(flx_during_pluvial.list[[i]])),rep("Megadrought",length(dep_during_mega.list[[i]])), rep("Megapluvial", length(dep_during_pluvial.list[[i]])), rep("Megadrought",length(TLAI_during_mega.list[[i]])), rep("Megapluvial", length(TLAI_during_pluvial.list[[i]])), rep("Megadrought",length(TSAI_during_mega.list[[i]])), rep("Megapluvial", length(TSAI_during_pluvial.list[[i]])), rep("Megadrought",length(T_during_mega.list[[i]])), rep("Megapluvial", length(T_during_pluvial.list[[i]])),rep("Megadrought",length(P_during_mega.list[[i]])), rep("Megapluvial", length(P_during_pluvial.list[[i]])), rep("Megadrought",length(SW_during_mega.list[[i]])), rep("Megapluvial", length(SW_during_pluvial.list[[i]])),rep("Megadrought",length(FSNO_during_mega.list[[i]])), rep("Megapluvial", length(FSNO_during_pluvial.list[[i]])), rep("Megadrought", length(U10_during_mega.list[[i]])), rep("Megapluvial", length(U10_during_pluvial.list[[i]]))),
    Variable = c(rep("FLXT (g/m2/yr)", length(flx_during_mega.list[[i]])), rep("FLXT (g/m2/yr)", length(flx_during_pluvial.list[[i]])),rep("DEP (g/m2/yr)", length(dep_during_mega.list[[i]])), rep("DEP (g/m2/yr)", length(dep_during_pluvial.list[[i]])), rep("TLAI", length(TLAI_during_mega.list[[i]])), rep("TLAI", length(TLAI_during_pluvial.list[[i]])), rep("TSAI", length(TSAI_during_mega.list[[i]])), rep("TSAI", length(TSAI_during_pluvial.list[[i]])), rep("T (K)", length(T_during_mega.list[[i]])), rep("T (K)", length(T_during_pluvial.list[[i]])), rep("PPT (m/s)", length(P_during_mega.list[[i]])), rep("PPT (m/s)", length(P_during_pluvial.list[[i]])), rep("SW (kg/m2)", length(SW_during_mega.list[[i]])), rep("SW (kg/m2)", length(SW_during_pluvial.list[[i]])), rep("FSNO", length(FSNO_during_mega.list[[i]])), rep("FSNO", length(FSNO_during_pluvial.list[[i]])), rep("U10", length(U10_during_mega.list[[i]])), rep("U10", length(U10_during_pluvial.list[[i]]))))
  
}
names(dens.df) <- all.ens

results <- "/projects/pd_lab/sha59/LME/"

density.graphs(dens.df)

r <- seq(1,1156,1)

for(i in 1:length(all.ens)){
  ann.data.35[[i]]$Period <- rep(NA,1156)
  ann.data.35[[i]]$Period[mega.list[[i]]] <- "MD"
  ann.data.35[[i]]$Period[pluvial.list[[i]]] <- "MP"
}

ann.data.m <- list()
for(i in 1:length(all.ens)){
  ann.data.m[[i]] <- gather(as.data.frame(ann.data.35[[i]]), key, value, -Years, -Period) 
}

plot(ann.data.35[[10]]$SOILWATER_10CM, type = "l")
lines(y = ann.data.35[[10]]$SOILWATER_10CM[mega.list[[10]]], x = mega.list[[10]], col = "red")
lines(y = ann.data.35[[10]]$SOILWATER_10CM[pluvial.list[[10]]], x = pluvial.list[[10]], col = "blue")

ts.graphs(ann.data.m)

#### Comparing soil moisture from LME and tree rings



rep.times.dep_during_mega <- vector()
rep.times.dep_during_pluvial <- vector()
rep.times.dep_during_nonmega <- vector()
rep.times.dep_during_nonpluvial <- vector()


for(i in 1: length(all.ens.names)){
  rep.times.dep_during_mega[i] <- length(dep_during_mega.list[[i]])
  rep.times.dep_during_pluvial[i] <- length(dep_during_pluvial.list[[i]])
  
  rep.times.dep_during_nonmega[i] <- length(dep_during_non_mega.list[[i]])
  rep.times.dep_during_nonpluvial[i] <- length(dep_during_non_pluvial.list[[i]])
  
}

Nounc.dep.ens.mega <- data.frame(Dust_deposition = c(unlist(dep_during_mega.list)*31536000000, unlist(dep_during_pluvial.list)*31536000000), Ensemble = c(rep("2", rep.times.dep_during_mega[1]), rep("3", rep.times.dep_during_mega[2]), rep("4", rep.times.dep_during_mega[3]), rep("5", rep.times.dep_during_mega[4]), rep("6", rep.times.dep_during_mega[5]), rep("7", rep.times.dep_during_mega[6]), rep("8", rep.times.dep_during_mega[7]), rep("9",rep.times.dep_during_mega[8]), rep("10", rep.times.dep_during_mega[9]), rep("11", rep.times.dep_during_mega[10]),rep("2", rep.times.dep_during_pluvial[1]), rep("3", rep.times.dep_during_pluvial[2]), rep("4", rep.times.dep_during_pluvial[3]), rep("5", rep.times.dep_during_pluvial[4]), rep("6", rep.times.dep_during_pluvial[5]), rep("7", rep.times.dep_during_pluvial[6]), rep("8", rep.times.dep_during_pluvial[7]), rep("9",rep.times.dep_during_pluvial[8]), rep("10", rep.times.dep_during_pluvial[9]), rep("11", rep.times.dep_during_pluvial[10]) ), Dataset = c(rep("Megadrought",length(dep_during_mega)), rep("Megapluvial", length(dep_during_pluvial))))

NonMD <- data.frame(Dust_deposition = c(unlist(dep_during_mega.list)*31536000000,unlist(dep_during_non_mega.list)*31536000000), Ensemble = c(rep("2", rep.times.dep_during_mega[1]), rep("3", rep.times.dep_during_mega[2]), rep("4", rep.times.dep_during_mega[3]), rep("5", rep.times.dep_during_mega[4]), rep("6", rep.times.dep_during_mega[5]), rep("7", rep.times.dep_during_mega[6]), rep("8", rep.times.dep_during_mega[7]), rep("9",rep.times.dep_during_mega[8]), rep("10", rep.times.dep_during_mega[9]), rep("11", rep.times.dep_during_mega[10]),rep("2", rep.times.dep_during_nonmega[1]), rep("3", rep.times.dep_during_nonmega[2]), rep("4", rep.times.dep_during_nonmega[3]), rep("5", rep.times.dep_during_nonmega[4]), rep("6", rep.times.dep_during_nonmega[5]), rep("7", rep.times.dep_during_nonmega[6]), rep("8", rep.times.dep_during_nonmega[7]), rep("9",rep.times.dep_during_nonmega[8]), rep("10", rep.times.dep_during_nonmega[9]), rep("11", rep.times.dep_during_nonmega[10])), Dataset = c(rep("Megadrought",length(dep_during_mega)),rep("Other", length(dep_during_non_mega))))

LME.nonMD <- data.frame(Dust_deposition = c(unlist(dep_during_mega.list)*31536000000, unlist(dep_during_non_mega.list)*31536000000), Dataset = c(rep("Megadrought", length(unlist(dep_during_mega.list))), rep("Other", length(unlist(dep_during_non_mega.list)))))

NonMP <- data.frame(Dust_deposition = c(unlist(dep_during_pluvial.list)*31536000000,unlist(dep_during_non_pluvial.list)*31536000000), Ensemble = c(rep("2", rep.times.dep_during_pluvial[1]), rep("3", rep.times.dep_during_pluvial[2]), rep("4", rep.times.dep_during_pluvial[3]), rep("5", rep.times.dep_during_pluvial[4]), rep("6", rep.times.dep_during_pluvial[5]), rep("7", rep.times.dep_during_pluvial[6]), rep("8", rep.times.dep_during_pluvial[7]), rep("9",rep.times.dep_during_pluvial[8]), rep("10", rep.times.dep_during_pluvial[9]), rep("11", rep.times.dep_during_pluvial[10]),rep("2", rep.times.dep_during_nonpluvial[1]), rep("3", rep.times.dep_during_nonpluvial[2]), rep("4", rep.times.dep_during_nonpluvial[3]), rep("5", rep.times.dep_during_nonpluvial[4]), rep("6", rep.times.dep_during_nonpluvial[5]), rep("7", rep.times.dep_during_nonpluvial[6]), rep("8", rep.times.dep_during_nonpluvial[7]), rep("9",rep.times.dep_during_nonpluvial[8]), rep("10", rep.times.dep_during_nonpluvial[9]), rep("11", rep.times.dep_during_nonpluvial[10])), Dataset = c(rep("Megapluvial",length(dep_during_pluvial)),rep("Other", length(dep_during_non_pluvial))))

Nounc.mega.plot <- ggplot(Nounc.dep.ens.mega, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  xlab("")+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Nounc.mega.plot

ggsave(Nounc.mega.plot, filename = paste0(results,"LME_nounc_all_",duration,"_",Sys.Date(),".pdf"))

LME.ens.plot  <- ggplot(Nounc.dep.ens.mega, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
LME.ens.plot

ggsave(LME.ens.plot, filename = paste0(results,"LME_no_unc_ens_",duration,"_",Sys.Date(),".pdf"))

LME.nonMD.plot <- ggplot(NonMD, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("#EFC000FF", "darkgrey"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.9,0.20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
LME.nonMD.plot

ggsave(LME.nonMD.plot, filename = paste0(results,"LME_noUnc_nonMD_",duration,"_",Sys.Date(),".pdf"))

Nounc.nonMD.plot <- ggplot(LME.nonMD, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  scale_fill_manual(values = c("#EFC000FF", "darkgrey"))+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  xlab("")+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Nounc.nonMD.plot

ggsave(Nounc.nonMD.plot, filename = paste0(results,"LME_nounc_all_nonMD_",duration,"_",Sys.Date(),".pdf"))

LME.nonMP.plot <- ggplot(NonMP, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("#0073C2FF", "darkgrey"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
LME.nonMP.plot

ggsave(LME.nonMP.plot, filename = paste0(results,"LME_noUnc_nonMP_",duration,"_",Sys.Date(),".pdf"))

MD.nounc.stats <- NonMD %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition),
            grp.IQR = IQR(Dust_deposition),
            grp.median = median(Dust_deposition))
MD.nounc.stats

trues <- vector()
trues.t <- vector()

for(i in c(2:4,6:11)){
  tmpD  <- Nounc.dep.ens.mega %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megadrought")
  tmpP  <- Nounc.dep.ens.mega %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megapluvial")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  t <- t.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "greater")
  print(k$p.value)
  print(t$p.value)
  
  if(k$p.value <= 0.05){
    trues[i-1] <- T 
  }else{
    trues[i-1] <- F
  }
  if(t$p.value <= 0.05){
    trues.t[i-1] <- T 
  }else{
    trues.t[i-1] <- F
  }
 
}
sum(trues, na.rm = T)/(length(trues)-1)*100
sum(trues.t, na.rm = T)/(length(trues.t)-1)*100

trues <- vector()
trues.t <- vector()
for(i in c(2:4,6:11)){
  tmpD  <- NonMD %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megadrought")
  tmpP  <- NonMD %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Other")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  t <- t.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "greater")
  print(k$p.value)
  print(t$p.value)
  
  if(k$p.value <= 0.05){
    trues[i-1] <- T 
  }else{
    trues[i-1] <- F
  }
  if(t$p.value <= 0.05){
    trues.t[i-1] <- T 
  }else{
    trues.t[i-1] <- F
  }
  
}
sum(trues, na.rm = T)/(length(trues)-1)*100
sum(trues.t, na.rm = T)/(length(trues.t)-1)*100

#### Blue Lake

max.ens  <- 1000

setwd("/home/sha59")
PMDI <- read.csv("tree_PMDI.csv")
Blue_DMAR <- read.csv("Blue_DMAR_ens.csv")
Blue_age <- read.csv("Blue_age_ens.csv")
LME_dust <- read.csv("LME_ann_dep_sum.csv")
LME_soilw <- read.csv("LME_ann_soilw.csv")

PMDI.stats <- vector()
PMDI.stats[1] <- mean(PMDI$PMDI.35,na.rm = T)
PMDI.stats[2] <- sd(PMDI$PMDI.35, na.rm = T)/2

index.m <- list()
index.p <- list()

index.m<- which(PMDI$PMDI.35 <= PMDI.stats[1] - PMDI.stats[2])
index.p <- which(PMDI$PMDI.35 >= PMDI.stats[1] + PMDI.stats[2]) 

cumsum.m <- vector()
cumsum.p <- vector()
soilw_during_m <- list()

cumsum.m <- cumsum(c(1,abs(index.m[-length(index.m)] - index.m[-1]) >1))
result <- table(cumsum.m)
group <- as.integer(which(result >=35)) 
want.index <- which(cumsum.m %in% group)
inlist.index <- index.m[want.index]
mega.index  <- inlist.index
years.mega <- PMDI$Age[mega.index]
years.nonMD <- PMDI$Age[-mega.index]

cumsum.p <- cumsum(c(1,abs(index.p[-length(index.p)] - index.p[-1]) >1))
result <- table(cumsum.p)
group <- as.integer(which(result >=35)) 
want.index <- which(cumsum.p %in% group)
inlist.index <- index.p[want.index]
pluvial.index  <- inlist.index
soilw_during_m <- PMDI$PMDI.35[mega.index]
years.pluvial <- PMDI$Age[pluvial.index]
years.nonMP <- PMDI$Age[-pluvial.index]

plot(PMDI$Age, PMDI$PMDI.35, type = "l")
lines(y = PMDI$PMDI.35[mega.index], x = PMDI$Age[mega.index], col = "red")
lines(y = PMDI$PMDI.35[pluvial.index], x = PMDI$Age[pluvial.index], col = "blue")

mean.age <- round(rowMeans(Blue_age[,-1]))
Blue_age[,1] <- mean.age
names(Blue_age)[1] <- "mean.age"
Blue_age_abs <- round(Blue_age)
Blue_age_abs_com <- Blue_age_abs[min(which(Blue_age_abs$mean.age <= 2005)):max(which(Blue_age_abs$mean.age >= 850)),]
Blue_DMAR_com <- Blue_DMAR[min(which(Blue_age_abs$mean.age <= 2005)):max(which(Blue_age_abs$mean.age >= 850)),]
Blue_DMAR_com_sorted <- Blue_DMAR_com[order(-Blue_DMAR_com$X),]
Blue_DMAR_sorted <- Blue_DMAR[order(-Blue_DMAR$X),]
Blue_age_abs_com_sorted <- Blue_age_abs_com[order(Blue_age_abs_com$mean.age),]
Blue_age_abs_com_sorted_df <- Blue_age_abs_com_sorted
Blue_age_abs_com_sorted <- Blue_age_abs_com_sorted[,-1]
Blue_DMAR_com_sorted <- Blue_DMAR_com_sorted[,-1]

age.index.mega <- list()
dust.mega  <- list()
dust.nonMD <- list()
age.index.pluvial <- list()
dust.pluvial  <- list()
dust.nonMP <- list()

for(i in 1:max.ens){
  print(i)
  age.index.mega[[i]] <- which(sort(round(Blue_age_abs[complete.cases(Blue_age_abs),i+1])) %in% years.mega)
  c  <- age.index.mega[[i]]
  c  <- c[!is.na(c)]
  dust.mega[[i]] <- Blue_DMAR_sorted[c,i+1]
  dust.nonMD[[i]] <- Blue_DMAR_sorted[-c,i+1]
  
  age.index.pluvial[[i]] = which(sort(round(Blue_age_abs[complete.cases(Blue_age_abs),i+1])) %in% years.pluvial)
  d  <- age.index.pluvial[[i]]
  d  <- d[!is.na(d)]
  dust.pluvial[[i]] <- Blue_DMAR_sorted[d,i+1]
  dust.nonMP[[i]] <- Blue_DMAR_sorted[-d,i+1]
  
}

Blue.depo.mega <- data.frame(Dust_deposition = c(unlist(dust.mega), unlist(dust.pluvial))*1000, Dataset = c(rep("Megadrought", length(unlist(dust.mega))), rep("Megapluvial", length(unlist(dust.pluvial)))))

Blue.nonMD <- data.frame(Dust_deposition = c(unlist(dust.mega), unlist(dust.nonMD))*1000, Dataset = c(rep("Megadrought", length(unlist(dust.mega))), rep("Other", length(unlist(dust.nonMD)))))

Blue.nonMP <- data.frame(Dust_deposition = c(unlist(dust.pluvial), unlist(dust.nonMP))*1000, Dataset = c(rep("Megapluvial", length(unlist(dust.pluvial))), rep("Other", length(unlist(dust.nonMP)))))

Blue.depo.mega$Dust_deposition[Blue.depo.mega$Dust_deposition < 0]  <- 0

Blue.depo.mega.plot <- ggplot(Blue.depo.mega, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  scale_x_continuous(limits = c(0,20))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Blue.depo.mega.plot

ggsave(Blue.depo.mega.plot, filename = paste0(results,"Blue_lake_",duration,"_",Sys.Date(),".pdf"))

Blue.depo.nonMD.plot <- ggplot(Blue.nonMD, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "darkgrey"))+
  xlab("")+
  theme_bw()+
  scale_x_continuous(limits = c(0,20))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Blue.depo.nonMD.plot

ggsave(Blue.depo.nonMD.plot, filename = paste0(results,"Blue_lake_nonMD_",duration,"_",Sys.Date(),".pdf"))

Blue.depo.nonMP.plot <- ggplot(Blue.nonMP, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#0073C2FF", "darkgrey"))+
  xlab("")+
  theme_bw()+
  scale_x_continuous(limits = c(0,20))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Blue.depo.nonMP.plot

ggsave(Blue.depo.nonMP.plot, filename = paste0(results,"Blue_lake_nonMP_",duration,"_",Sys.Date(),".pdf"))


Blue.depo.mega.all <- data.frame(Dust_deposition = c(unlist(dust.mega)*1000, unlist(dust.pluvial)*1000), Ensemble = as.factor(c(rep(1:1000,rapply(dust.mega, length, how = "unlist")[1:1000]),rep(1:1000,rapply(dust.pluvial, length, how = "unlist")[1:1000]))), Dataset = c(rep("Megadrought",length(unlist(dust.mega))), rep("Megapluvial", length(unlist(dust.pluvial)))))

Blue.nonMD.mega.all <- data.frame(Dust_deposition = c(unlist(dust.mega)*1000, unlist(dust.nonMD)*1000), Ensemble = as.factor(c(rep(1:1000,rapply(dust.mega, length, how = "unlist")[1:1000]),rep(1:1000,rapply(dust.nonMD, length, how = "unlist")[1:1000]))), Dataset = c(rep("Megadrought",length(unlist(dust.mega))), rep("Other", length(unlist(dust.nonMD)))))

Blue.depo.mega.all.plot <- ggplot(Blue.depo.mega.all, aes(Dust_deposition, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,40))+
  theme(legend.position = "none")

Blue.depo.mega.all$Dust_deposition[Blue.depo.mega.all$Dust_deposition < 0]  <- 0

unpaired.plot <- na.omit(Blue.depo.mega.all) %>% 
  filter(Ensemble %in% as.character(sample(seq(1,1000,1),8))) %>%
  dabest(Dataset, Dust_deposition,
         idx = c("Megadrought", "Megapluvial"), paired = F)

unpaired.plot.viz <- plot(unpaired.plot, color.column = Ensemble)

ggsave(unpaired.plot.viz, filename = paste0(results,"Blue_lake_unpaired_",duration,"_",Sys.Date(),".pdf"))

MD.Blue.stats <- na.omit(Blue.nonMD.mega.all) %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition),
            grp.IQR = IQR(Dust_deposition),
            grp.median = median(Dust_deposition))
MD.Blue.stats

#plot(ecdf(tmpD$Dust_deposition))
#plot(ecdf(tmpP$Dust_deposition), col = "blue", add = T)


trues <- vector()
trues.t <- vector()

for(i in 1:1000){
  tmpD  <- Blue.depo.mega.all %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megadrought")
  tmpP  <- Blue.depo.mega.all %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megapluvial")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  t <- t.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "greater")
    
  if(k$p.value <= 0.05){
    trues[i] <- T 
  }else{
    trues[i] <- F
  }
  
  if(t$p.value <= 0.05){
    trues.t[i] <- T 
  }else{
    trues.t[i] <- F
  }
  
}
sum(trues)/length(trues)*100
sum(trues.t)/length(trues.t)*100

trues <- vector()
trues.t <- vector()

for(i in 1:1000){
  tmpD  <- Blue.nonMD.mega.all %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megadrought")
  tmpP  <- Blue.nonMD.mega.all %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Other")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  t <- t.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "greater")
  
  if(k$p.value <= 0.05){
    trues[i] <- T 
  }else{
    trues[i] <- F
  }
  
  if(t$p.value <= 0.05){
    trues.t[i] <- T 
  }else{
    trues.t[i] <- F
  }
  
}
sum(trues)/length(trues)*100
sum(trues.t)/length(trues.t)*100

#########
#### LME with uncertainty

Blue_age_abs_com_sorted_edge <- matrix(data = NA,nrow = nrow(Blue_age_abs_com_sorted), ncol=ncol(Blue_age_abs_com_sorted))
for(j in 1:ncol(Blue_age_abs_com_sorted)){
  print(j)
  Blue_age_abs_com_sorted_edge[1,j]  <- Blue_age_abs_com_sorted[1,j]-round(mean(diff(as.matrix(Blue_age_abs_com_sorted[1:2,]))))
}

for(i in 1:nrow(Blue_age_abs_com_sorted)-1){
  for(j in 1:ncol(Blue_age_abs_com_sorted)){
    print(i)
    print(j)
    
    Blue_age_abs_com_sorted_edge[i+1,j] <- mean(Blue_age_abs_com_sorted[i:i+1,j])
  }
}
save.image("/scratch/sha59/LME/all_ensemble_part4_v2.RData")

LME <- data.frame(Year = seq(850,2005,1), LME_dust[,-1], LME_soilw[,-1])

LMEbins <- list()
for(i in 1:length(all.ens)){
  d <- bin(time = LME$Year, values = LME[,i+1], binvec = Blue_age_abs_com_sorted_df$mean.age)
  LMEbins[[i]]  <- data.frame(Year = c(d$x,NA), Dust = c(d$y,NA))
}
names(LMEbins) <- c("LME2","LME3", "LME4","LME5", "LME6", "LME7", "LME8", "LME9", "LME10", "LME11")

index.mega <- list()
index.pluvial <- list()
for(i in 1:length(all.ens)){
  index.mega[[i]] <- LME$Year[which(ann.data.35[[i]]$SOILWATER_10CM <= soil.stats[1,i] - soil.stats[2,i])]
  index.pluvial[[i]] <- LME$Year[which(ann.data.35[[i]]$SOILWATER_10CM >= soil.stats[1,i] + soil.stats[2,i])]
}

cumsum.mega <- vector()
dep_during_mega_unc <- list()
dep_during_non_mega_unc <- list()
values <- list()
others <- list()

for(i in 1:length(all.ens)){
  for(j in 1:max.ens){
    print(i)
    cumsum.mega <- cumsum(c(1,abs(index.mega[[i]][-length(index.mega[[i]])] - index.mega[[i]][-1]) >1))
    result <- table(cumsum.mega)
    group <- as.integer(which(result >=35)) 
    want.index <- which(cumsum.mega %in% group)
    inlist.index <- index.mega[[i]][want.index]
    
    dataIWant <- which(Blue_age_abs_com_sorted[,j] %in% inlist.index)
    
    yearCheck <- Blue_age_abs_com_sorted[dataIWant,j]
    values[[j]] <- LMEbins[[i]]$Dust[dataIWant]
    others[[j]] <- LMEbins[[i]]$Dust[-dataIWant]
    
  }
  dep_during_mega_unc[[i]] <- values 
  dep_during_non_mega_unc[[i]] <- others
}
dep_during_mega_unc.list <- list()
dep_during_non_mega_unc.list <- list()
for(i in 1:length(all.ens)){
  dep_during_mega_unc.list[[i]]  <- unlist(dep_during_mega_unc[[i]])
  dep_during_non_mega_unc.list[[i]]  <- unlist(dep_during_non_mega_unc[[i]])
}

dep_during_mega_unc.sep <- dep_during_mega_unc
dep_during_non_mega_unc.sep <- dep_during_non_mega_unc

dep_during_mega_unc <- unlist(dep_during_mega_unc)
dep_during_non_mega_unc <- unlist(dep_during_non_mega_unc)

dep_during_pluvial_unc <- list()
dep_during_non_pluvial_unc <- list()
cumsum.pluvial <- vector()
values <- list()
others <- list()

for(i in 1:length(all.ens)){
  for(j in 1:max.ens){
    print(i)
    cumsum.pluvial <- cumsum(c(1,abs(index.pluvial[[i]][-length(index.pluvial[[i]])] - index.pluvial[[i]][-1]) >1))
    result <- table(cumsum.pluvial)
    group <- as.integer(which(result >=35)) 
    want.index <- which(cumsum.pluvial %in% group)
    inlist.index <- index.pluvial[[i]][want.index]
    
    dataIWant <- which(Blue_age_abs_com_sorted[,j] %in% inlist.index)
    
    yearCheck <- Blue_age_abs_com_sorted[dataIWant,j]
    values[[j]] <- LMEbins[[i]]$Dust[dataIWant]
    others[[j]] <- LMEbins[[i]]$Dust[-dataIWant]
  }
  dep_during_pluvial_unc[[i]] <- values
  dep_during_non_pluvial_unc[[i]] <- others
}
dep_during_pluvial_unc.list <- list()
dep_during_non_pluvial_unc.list <- list()
for(i in 1:length(all.ens)){
  dep_during_pluvial_unc.list[[i]]  <- unlist(dep_during_pluvial_unc[[i]])
  dep_during_non_pluvial_unc.list[[i]]  <- unlist(dep_during_non_pluvial_unc[[i]])
}

dep_during_pluvial_unc.sep <- dep_during_pluvial_unc
dep_during_non_pluvial_unc.sep <- dep_during_non_pluvial_unc

dep_during_pluvial_unc <- unlist(dep_during_pluvial_unc)
dep_during_non_pluvial_unc <- unlist(dep_during_non_pluvial_unc)

rep.times.dep_during_mega <- vector()
rep.times.dep_during_pluvial <- vector()
rep.times.dep_during_nonmega <- vector()
rep.times.dep_during_nonpluvial <- vector()


for(i in 1: length(all.ens.names)){
  rep.times.dep_during_mega[i] <- length(dep_during_mega_unc.list[[i]])
  rep.times.dep_during_pluvial[i] <- length(dep_during_pluvial_unc.list[[i]])
  
  rep.times.dep_during_nonmega[i] <- length(dep_during_non_mega_unc.list[[i]])
  rep.times.dep_during_nonpluvial[i] <- length(dep_during_non_pluvial_unc.list[[i]])
}

dep.ens.mega <- data.frame(Dust_deposition = c(unlist(dep_during_mega_unc.list)*31536000000, unlist(dep_during_pluvial_unc.list)*31536000000), Ensemble = c(rep("2", rep.times.dep_during_mega[1]), rep("3", rep.times.dep_during_mega[2]), rep("4", rep.times.dep_during_mega[3]), rep("5", rep.times.dep_during_mega[4]), rep("6", rep.times.dep_during_mega[5]), rep("7", rep.times.dep_during_mega[6]), rep("8", rep.times.dep_during_mega[7]), rep("9",rep.times.dep_during_mega[8]), rep("10", rep.times.dep_during_mega[9]), rep("11", rep.times.dep_during_mega[10]),rep("2", rep.times.dep_during_pluvial[1]), rep("3", rep.times.dep_during_pluvial[2]), rep("4", rep.times.dep_during_pluvial[3]), rep("5", rep.times.dep_during_pluvial[4]), rep("6", rep.times.dep_during_pluvial[5]), rep("7", rep.times.dep_during_pluvial[6]), rep("8", rep.times.dep_during_pluvial[7]), rep("9",rep.times.dep_during_pluvial[8]), rep("10", rep.times.dep_during_pluvial[9]), rep("11", rep.times.dep_during_pluvial[10])), Dataset = c(rep("Megadrought",length(dep_during_mega_unc)), rep("Megapluvial", length(dep_during_pluvial_unc))))

ens.nonMD <- data.frame(Dust_deposition = c(unlist(dep_during_mega_unc.list)*31536000000, unlist(dep_during_non_mega_unc.list)*31536000000), Ensemble = c(rep("2", rep.times.dep_during_mega[1]), rep("3", rep.times.dep_during_mega[2]), rep("4", rep.times.dep_during_mega[3]), rep("5", rep.times.dep_during_mega[4]), rep("6", rep.times.dep_during_mega[5]), rep("7", rep.times.dep_during_mega[6]), rep("8", rep.times.dep_during_mega[7]), rep("9",rep.times.dep_during_mega[8]), rep("10", rep.times.dep_during_mega[9]), rep("11", rep.times.dep_during_mega[10]), rep("2", rep.times.dep_during_nonmega[1]), rep("3", rep.times.dep_during_nonmega[2]), rep("4", rep.times.dep_during_nonmega[3]), rep("5", rep.times.dep_during_nonmega[4]), rep("6", rep.times.dep_during_nonmega[5]), rep("7", rep.times.dep_during_nonmega[6]), rep("8", rep.times.dep_during_nonmega[7]), rep("9",rep.times.dep_during_nonmega[8]), rep("10", rep.times.dep_during_nonmega[9]), rep("11", rep.times.dep_during_nonmega[10])), Dataset = c(rep("Megadrought",length(dep_during_mega_unc)), rep("Other", length(dep_during_non_mega_unc))))

ens.nonMP <- data.frame(Dust_deposition = c(unlist(dep_during_pluvial_unc.list)*31536000000, unlist(dep_during_non_pluvial_unc.list)*31536000000), Ensemble = c(rep("2", rep.times.dep_during_pluvial[1]), rep("3", rep.times.dep_during_pluvial[2]), rep("4", rep.times.dep_during_pluvial[3]), rep("5", rep.times.dep_during_pluvial[4]), rep("6", rep.times.dep_during_pluvial[5]), rep("7", rep.times.dep_during_pluvial[6]), rep("8", rep.times.dep_during_pluvial[7]), rep("9",rep.times.dep_during_pluvial[8]), rep("10", rep.times.dep_during_pluvial[9]), rep("11", rep.times.dep_during_pluvial[10]), rep("2", rep.times.dep_during_nonpluvial[1]), rep("3", rep.times.dep_during_nonpluvial[2]), rep("4", rep.times.dep_during_nonpluvial[3]), rep("5", rep.times.dep_during_nonpluvial[4]), rep("6", rep.times.dep_during_nonpluvial[5]), rep("7", rep.times.dep_during_nonpluvial[6]), rep("8", rep.times.dep_during_nonpluvial[7]), rep("9",rep.times.dep_during_nonpluvial[8]), rep("10", rep.times.dep_during_nonpluvial[9]), rep("11", rep.times.dep_during_nonpluvial[10])), Dataset = c(rep("Megapluvial",length(dep_during_pluvial_unc)), rep("Other", length(dep_during_non_pluvial_unc))))

Depo.mega.ens.plot <- ggplot(dep.ens.mega, aes(Dust_deposition, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,20))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.65),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

uncertainty.plot  <- ggplot(dep.ens.mega, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
uncertainty.plot

ggsave(uncertainty.plot, filename = paste0(results,"LME_unc_ens_",duration,"_",Sys.Date(),".pdf"))

uncMD.plot  <- ggplot(ens.nonMD, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("#EFC000FF", "darkgrey"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.9,0.2),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
uncMD.plot

ggsave(uncMD.plot, filename = paste0(results,"LME_unc_ens_nonMD_",duration,"_",Sys.Date(),".pdf"))

uncMP.plot  <- ggplot(ens.nonMP, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("darkgrey", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
uncMP.plot

ggsave(uncMP.plot, filename = paste0(results,"LME_unc_ens_nonMP_",duration,"_",Sys.Date(),".pdf"))

df.dep.unc <- data.frame(Dust_deposition = c(dep_during_mega_unc*31536000000,dep_during_pluvial_unc*31536000000), Dataset = c(rep("Megadrought",length(dep_during_mega_unc)), rep("Megapluvial", length(dep_during_pluvial_unc))))

Unc.mega.plot <- ggplot(df.dep.unc, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  xlab("")+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Unc.mega.plot

ggsave(Unc.mega.plot, filename = paste0(results,"LME_unc_",duration,"_",Sys.Date(),".pdf"))

LME.unc.nonMD <- data.frame(Dust_deposition = c(unlist(dep_during_mega_unc.list)*31536000000, unlist(dep_during_non_mega_unc.list)*31536000000), Dataset = c(rep("Megadrought", length(unlist(dep_during_mega_unc.list))), rep("Other", length(unlist(dep_during_non_mega_unc.list)))))

Unc.nonMD.plot <- ggplot(LME.unc.nonMD, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  scale_fill_manual(values = c("#EFC000FF", "darkgrey"))+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  xlab("")+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Unc.nonMD.plot

ggsave(Unc.nonMD.plot, filename = paste0(results,"LME_unc_all_nonMD_",duration,"_",Sys.Date(),".pdf"))

MD.unc.stats <- na.omit(LME.unc.nonMD) %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition),
            grp.IQR = IQR(Dust_deposition),
            grp.median = median(Dust_deposition))
MD.unc.stats

trues <- vector()
trues.t <- vector()

for(i in c(2:4,6:11)){
  tmpD  <- dep.ens.mega %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megadrought")
  tmpP  <- dep.ens.mega %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megapluvial")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  t  <- t.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "greater")
  print(k$p.value)
  print(t$p.value)
  if(k$p.value <= 0.05){
    trues[i-1] <- T 
  }else{
    trues[i-1] <- F
  }
  
  if(t$p.value <= 0.05){
    trues.t[i-1] <- T 
  }else{
    trues.t[i-1] <- F
  }
}
sum(trues, na.rm = T)/(length(trues)-1)*100
sum(trues.t, na.rm = T)/(length(trues.t)-1)*100

trues <- vector()
trues.t <- vector()

for(i in c(2:4,6:11)){
  tmpD  <- ens.nonMD %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Megadrought")
  tmpP  <- ens.nonMD %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Other")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  t  <- t.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "greater")
  print(k$p.value)
  print(t$p.value)
  if(k$p.value <= 0.05){
    trues[i-1] <- T 
  }else{
    trues[i-1] <- F
  }
  
  if(t$p.value <= 0.05){
    trues.t[i-1] <- T 
  }else{
    trues.t[i-1] <- F
  }
}
sum(trues, na.rm = T)/(length(trues)-1)*100
sum(trues.t, na.rm = T)/(length(trues.t)-1)*100

# For each age ensemble member in each model ensemble
trues <- vector()
trues.t <- vector()
all.trues <- list()
all.trues.t <- list()

for(j in c(1:10)){
  for(i in 1:1000){
    tryCatch({
    print(j)
    tmpD  <- unlist(dep_during_mega_unc.sep[[j]][i])
    tmpP  <- unlist(dep_during_pluvial_unc.sep[[j]][i])
    k  <- ks.test(tmpD,tmpP,alternative = "less")
    t <- t.test(tmpD,tmpP,alternative = "greater")
    
    if(k$p.value <= 0.05){
      trues[i] <- T 
    }else{
      trues[i] <- F
    }
    
    if(t$p.value <= 0.05){
      trues.t[i] <- T 
    }else{
      trues.t[i] <- F
    }} ,error = function(e){})
  }
  all.trues[[j]] <- trues
  all.trues.t[[j]] <- trues.t
}
sum(unlist(all.trues))/length(unlist(all.trues))*100
sum(unlist(all.trues.t))/length(unlist(all.trues.t))*100

for(i in 1:10){
  a[i] = sum(unlist(all.trues[[i]]))/length(unlist(all.trues[[i]]))*100
  print(sum(unlist(all.trues[[i]]))/length(unlist(all.trues[[i]]))*100)
}
for(i in 1:10){
  b[i] = sum(unlist(all.trues.t[[i]]))/length(unlist(all.trues.t[[i]]))*100
  print(sum(unlist(all.trues.t[[i]]))/length(unlist(all.trues.t[[i]]))*100)
}

# For MD and other years
trues <- vector()
trues.t <- vector()
all.trues <- list()
all.trues.t <- list()

for(j in c(1:10)){
  for(i in 1:1000){
    tryCatch({
      print(j)
      tmpD  <- unlist(dep_during_mega_unc.sep[[j]][i])
      tmpP  <- unlist(dep_during_non_pluvial_unc.sep[[j]][i])
      k  <- ks.test(tmpD,tmpP,alternative = "less")
      t <- t.test(tmpD,tmpP,alternative = "greater")
      
      if(k$p.value <= 0.05){
        trues[i] <- T 
      }else{
        trues[i] <- F
      }
      
      if(t$p.value <= 0.05){
        trues.t[i] <- T 
      }else{
        trues.t[i] <- F
      }} ,error = function(e){})
  }
  all.trues[[j]] <- trues
  all.trues.t[[j]] <- trues.t
}
sum(unlist(all.trues))/length(unlist(all.trues))*100
sum(unlist(all.trues.t))/length(unlist(all.trues.t))*100

for(i in 1:10){
  print(sum(unlist(all.trues[[i]]))/length(unlist(all.trues[[i]]))*100)
}
for(i in 1:10){
  print(sum(unlist(all.trues.t[[i]]))/length(unlist(all.trues.t[[i]]))*100)
}

#plot(ecdf(tmpD$Dust_deposition))
#plot(ecdf(tmpP$Dust_deposition), col = "blue", add = T)

#Driest and wettest period

#In Blue Lake
duration = 100
span = 100

PMDI$PMDI.100  <- as.vector(rollapply(PMDI$PMDI, width = duration, FUN = mean, by = 1, align = "center", na.rm = T, fill = NA))

start.index <- c(which.min(PMDI$PMDI.100),which.max(PMDI$PMDI.100))

index.m <- seq(start.index[1]-(span/2),start.index[1]+(span/2),1)
index.p <- seq(start.index[2]-(span/2),start.index[2]+(span/2),1)

years.mega  <- PMDI$Age[index.m]
years.pluvial  <- PMDI$Age[index.p]

age.index.mega <- list()
dust.mega  <- list()
age.index.pluvial <- list()
dust.pluvial  <- list()

for(i in 1:max.ens){
  print(i)
  age.index.mega[[i]] <- which(sort(round(Blue_age_abs[complete.cases(Blue_age_abs),i+1])) %in% years.mega)
  c  <- age.index.mega[[i]]
  c  <- c[!is.na(c)]
  dust.mega[[i]] <- Blue_DMAR_sorted[c,i+1]
  
  age.index.pluvial[[i]] = which(sort(round(Blue_age_abs[complete.cases(Blue_age_abs),i+1])) %in% years.pluvial)
  d  <- age.index.pluvial[[i]]
  d  <- d[!is.na(d)]
  dust.pluvial[[i]] <- Blue_DMAR_sorted[d,i+1]
  
}

Blue.depo.mega.100 <- data.frame(Dust_deposition = c(unlist(dust.mega), unlist(dust.pluvial))*1000, Dataset = c(rep("Driest century", length(unlist(dust.mega))), rep("Wettest century", length(unlist(dust.pluvial)))))

Blue.depo.mega.100.plot <- ggplot(Blue.depo.mega.100, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  scale_x_continuous(limits = c(0,20))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Blue.depo.mega.100.plot

ggsave(Blue.depo.mega.100.plot, filename = paste0(results,"Blue_lake_",duration,"_",Sys.Date(),".pdf"))

Blue.depo.mega.100.all <- data.frame(Dust_deposition = c(unlist(dust.mega)*1000, unlist(dust.pluvial)*1000), Ensemble = as.factor(c(rep(1:1000,rapply(dust.mega, length, how = "unlist")[1:1000]),rep(1:1000,rapply(dust.pluvial, length, how = "unlist")[1:1000]))), Dataset = c(rep("Driest century",length(unlist(dust.mega))), rep("Wettest century", length(unlist(dust.pluvial)))))

Blue.depo.mega.100.all$Dust_deposition[Blue.depo.mega.100.all$Dust_deposition < 0]  <- 0

unpaired.100.plot <- na.omit(Blue.depo.mega.100.all) %>% 
  filter(Ensemble %in% as.character(sample(seq(1,1000,1),8))) %>%
  dabest(Dataset, Dust_deposition,
         idx = c("Driest century", "Wettest century"), paired = F)

unpaired.100.plot.viz <- plot(unpaired.100.plot, color.column = Ensemble)

ggsave(unpaired.100.plot.viz, filename = paste0(results,"Blue_lake_unpaired_",duration,"_",Sys.Date(),".pdf"))

#plot(ecdf(tmpD$Dust_deposition))
#plot(ecdf(tmpP$Dust_deposition), col = "blue", add = T)

trues <- vector()
trues.t  <- vector()

for(i in 1:1000){
  tmpD  <- Blue.depo.mega.100.all %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Driest century")
  tmpP  <- Blue.depo.mega.100.all %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Wettest century")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  t  <- t.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "greater")
  
  if(k$p.value <= 0.05){
    trues[i] <- T 
  }else{
    trues[i] <- F
  }
  
  if(t$p.value <= 0.05){
    trues.t[i] <- T
  }else{
    trues.t[i] <- F
  }
}
sum(trues)/length(trues)*100
sum(trues.t)/length(trues.t)*100

cent.Blue.stats <- na.omit(Blue.depo.mega.100.all) %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition),
            grp.IQR = IQR(Dust_deposition),
            grp.median = median(Dust_deposition))
cent.Blue.stats

#In the LME without uncertainty

cent.SW <- list()
ann.data.100 <- list()

for(i in 1:length(all.ens)){
  for(j in 2:length(names(ann.data[[1]]))){
    
    cent.SW[[j-1]] <- as.vector(rollapply(ann.data[[i]][[j]], width = duration, FUN = mean, by = 1, align = "center", na.rm = T, fill = NA))
    
  }
  names(cent.SW) <- names(ann.data[[1]])[-1]
  ann.data.100[[i]] <- cent.SW
  ann.data.100[[i]]$Years <- seq(850,2005,1)
}
names(ann.data.100) <- all.ens

start.index <- list()
for(i in 1:length(all.ens)){
  start.index[[i]] <- c(which.min(ann.data.100[[i]]$SOILWATER_10CM),which.max(ann.data.100[[i]]$SOILWATER_10CM))
}

index.m <- list()
index.p <- list()
for(i in 1:length(all.ens)){
  a = start.index[[i]][1]
  b = start.index[[i]][2]
  index.m[[i]] <- seq(a-(span/2),a+(span/2),1)
  index.p[[i]] <- seq(b-(span/2),b+(span/2),1)
}

dep_during_cent_mega <- list()
dep_during_cent_pluvial <- list()

for(i in 1:length(all.ens)){
  dep_during_cent_mega[[i]] <- LME[[i+1]][index.m[[i]]]*31536000000
  dep_during_cent_pluvial[[i]] <- LME[[i+1]][index.p[[i]]]*31536000000
}

trues <- vector()

for(i in 1:10){
  tmpD  <- dep_during_cent_mega[[i]]
  tmpP  <- dep_during_cent_pluvial[[i]]
  k  <- ks.test(tmpD,tmpP,alternative = "less")
  print(k$p.value)
  print(k$statistic)
  if(k$p.value <= 0.05){
    trues[i] <- T 
  }else{
    trues[i] <- F
  }
}
sum(trues, na.rm = T)/(length(trues))*100

#plot(ecdf(tmpD))
#plot(ecdf(tmpP), col = "blue", add = T)

dep_during_cent_mega.list <- list()
dep_during_cent_pluvial.list <- list()
for(i in 1:length(all.ens)){
  dep_during_cent_mega.list[[i]]  <- unlist(dep_during_cent_mega[[i]])
  dep_during_cent_pluvial.list[[i]]  <- unlist(dep_during_cent_pluvial[[i]])
}
dep_during_cent_mega <- unlist(dep_during_cent_mega)
dep_during_cent_pluvial <- unlist(dep_during_cent_pluvial)

rep.times.dep_during_mega <- vector()
rep.times.dep_during_pluvial <- vector()

for(i in 1: length(all.ens.names)){
  rep.times.dep_during_mega[i] <- length(dep_during_cent_mega.list[[i]])
  rep.times.dep_during_pluvial[i] <- length(dep_during_cent_pluvial.list[[i]])
}

dep.cent.nounc.mega <- data.frame(Dust_deposition = c(unlist(dep_during_cent_mega.list), unlist(dep_during_cent_pluvial.list)), Ensemble = c(rep("2", rep.times.dep_during_mega[1]), rep("3", rep.times.dep_during_mega[2]), rep("4", rep.times.dep_during_mega[3]), rep("5", rep.times.dep_during_mega[4]), rep("6", rep.times.dep_during_mega[5]), rep("7", rep.times.dep_during_mega[6]), rep("8", rep.times.dep_during_mega[7]), rep("9",rep.times.dep_during_mega[8]), rep("10", rep.times.dep_during_mega[9]), rep("11", rep.times.dep_during_mega[10]),rep("2", rep.times.dep_during_pluvial[1]), rep("3", rep.times.dep_during_pluvial[2]), rep("4", rep.times.dep_during_pluvial[3]), rep("5", rep.times.dep_during_pluvial[4]), rep("6", rep.times.dep_during_pluvial[5]), rep("7", rep.times.dep_during_pluvial[6]), rep("8", rep.times.dep_during_pluvial[7]), rep("9",rep.times.dep_during_pluvial[8]), rep("10", rep.times.dep_during_pluvial[9]), rep("11", rep.times.dep_during_pluvial[10])), Dataset = c(rep("Driest century",length(dep_during_cent_mega)), rep("Wettest century", length(dep_during_cent_pluvial))))

Depo.mega.cent.nounc.plot <- ggplot(dep.cent.nounc.mega, aes(Dust_deposition, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,20))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.65),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cent.nouncertainty.plot  <- ggplot(dep.cent.nounc.mega, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
cent.nouncertainty.plot

ggsave(cent.nouncertainty.plot, filename = paste0(results,"LME_nounc_",duration,"_",Sys.Date(),".pdf"))

cent.df.dep.nounc <- data.frame(Dust_deposition = c(dep_during_cent_mega,dep_during_cent_pluvial), Dataset = c(rep("Driest century",length(dep_during_cent_mega)), rep("Wettest century", length(dep_during_cent_pluvial))))

Nounc.mega.cent.plot <- ggplot(cent.df.dep.nounc, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  xlab("")+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Nounc.mega.cent.plot

ggsave(Nounc.mega.cent.plot, filename = paste0(results,"LME_nounc_ens",duration,"_",Sys.Date(),".pdf"))

cent.nounc.stats <- cent.df.dep.nounc %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition),
            grp.IQR = IQR(Dust_deposition),
            grp.median = median(Dust_deposition))
cent.nounc.stats


#In the LME with uncertainty

years.mega <- list()
years.pluvial  <- list()
for(i in 1:length(all.ens)){
  years.mega[[i]]  <- ann.data[[1]]$Years[index.m[[i]]]
  years.pluvial[[i]]  <- ann.data[[1]]$Years[index.p[[i]]]
}

dep_during_mega_unc <- list()
values <- list()

for(i in 1:length(all.ens)){
  for(j in 1:max.ens){
    
  values[[j]] <- LMEbins[[i]]$Dust[which(Blue_age_abs_com_sorted[,j] %in% years.mega[[i]])]
    
  }
  dep_during_mega_unc[[i]] <- values 
}

dep_during_pluvial_unc <- list()
values <- list()

for(i in 1:length(all.ens)){
  for(j in 1:100){
    values[[j]] <- LMEbins[[i]]$Dust[which(Blue_age_abs_com_sorted[,j] %in% years.pluvial[[i]])]
    
  }
  dep_during_pluvial_unc[[i]] <- values
}

#plot(ecdf(tmpD$Dust_deposition))
#plot(ecdf(tmpP$Dust_deposition), col = "blue", add = T)

dep_during_mega_unc.list <- list()
for(i in 1:length(all.ens)){
  dep_during_mega_unc.list[[i]]  <- unlist(dep_during_mega_unc[[i]])
}
dep_during_mega_unc <- unlist(dep_during_mega_unc)

dep_during_pluvial_unc.list <- list()
for(i in 1:length(all.ens)){
  dep_during_pluvial_unc.list[[i]]  <- unlist(dep_during_pluvial_unc[[i]])
}
dep_during_pluvial_unc <- unlist(dep_during_pluvial_unc)

rep.times.dep_during_mega <- vector()
rep.times.dep_during_pluvial <- vector()

for(i in 1: length(all.ens.names)){
  rep.times.dep_during_mega[i] <- length(dep_during_mega_unc.list[[i]])
  rep.times.dep_during_pluvial[i] <- length(dep_during_pluvial_unc.list[[i]])
}

dep.cent.unc.mega <- data.frame(Dust_deposition = c(unlist(dep_during_mega_unc.list)*31536000000, unlist(dep_during_pluvial_unc.list)*31536000000), Ensemble = c(rep("2", rep.times.dep_during_mega[1]), rep("3", rep.times.dep_during_mega[2]), rep("4", rep.times.dep_during_mega[3]), rep("5", rep.times.dep_during_mega[4]), rep("6", rep.times.dep_during_mega[5]), rep("7", rep.times.dep_during_mega[6]), rep("8", rep.times.dep_during_mega[7]), rep("9",rep.times.dep_during_mega[8]), rep("10", rep.times.dep_during_mega[9]), rep("11", rep.times.dep_during_mega[10]),rep("2", rep.times.dep_during_pluvial[1]), rep("3", rep.times.dep_during_pluvial[2]), rep("4", rep.times.dep_during_pluvial[3]), rep("5", rep.times.dep_during_pluvial[4]), rep("6", rep.times.dep_during_pluvial[5]), rep("7", rep.times.dep_during_pluvial[6]), rep("8", rep.times.dep_during_pluvial[7]), rep("9",rep.times.dep_during_pluvial[8]), rep("10", rep.times.dep_during_pluvial[9]), rep("11", rep.times.dep_during_pluvial[10])), Dataset = c(rep("Driest century",length(dep_during_mega_unc)), rep("Wettest century", length(dep_during_pluvial_unc))))

Depo.mega.cent.unc.plot <- ggplot(dep.cent.unc.mega, aes(Dust_deposition, linetype = Ensemble))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("Annual deposition (g/m2/yr)")+
  theme_bw()+
  ylab("")+
  scale_x_continuous(limits = c(0,20))+
  theme(legend.title = element_blank(), legend.position = c(0.8,0.65),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cent.uncertainty.plot  <- ggplot(dep.cent.unc.mega, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  facet_wrap(Ensemble ~.)+
  coord_cartesian(xlim = c(0,40))+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  xlab("")+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
cent.uncertainty.plot

ggsave(cent.uncertainty.plot, filename = paste0(results,"LME_unc_",duration,"_",Sys.Date(),".pdf"))

cent.df.dep.unc <- data.frame(Dust_deposition = c(dep_during_mega_unc*31536000000,dep_during_pluvial_unc*31536000000), Dataset = c(rep("Driest century",length(dep_during_mega_unc)), rep("Wettest century", length(dep_during_pluvial_unc))))

Unc.mega.cent.plot <- ggplot(cent.df.dep.unc, aes(x = Dust_deposition))+
  geom_density(aes(fill = Dataset, y = ..scaled..), alpha = 0.4, adjust = 2.5)+
  scale_fill_manual(values = c("#EFC000FF", "#0073C2FF"))+
  theme_bw()+
  scale_x_continuous(limits = c(0,40))+
  xlab("")+
  theme(legend.title = element_blank(), legend.position = c(0.75,0.75),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Unc.mega.cent.plot

ggsave(Unc.mega.cent.plot, filename = paste0(results,"LME_unc_ens",duration,"_",Sys.Date(),".pdf"))

trues <- vector()

for(i in c(2:4,6:11)){
  tmpD  <- dep.cent.unc.mega %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Driest century")
  tmpP  <- dep.cent.unc.mega %>% filter(Ensemble %in% as.character(i)) %>%
    filter(Dataset == "Wettest century")
  k  <- ks.test(tmpD$Dust_deposition,tmpP$Dust_deposition,alternative = "less")
  print(k$p.value)
  if(k$p.value <= 0.05){
    trues[i-1] <- T 
  }else{
    trues[i-1] <- F
  }
}
sum(trues, na.rm = T)/(length(trues)-1)*100

save.image("/scratch/sha59/LME/all_ensemble_part4_v2.RData")

cent.unc.stats <- dep.cent.unc.mega %>% 
  group_by(Dataset) %>%
  summarise(grp.mean = mean(Dust_deposition),
            grp.IQR = IQR(Dust_deposition),
            grp.median = median(Dust_deposition))
cent.unc.stats
