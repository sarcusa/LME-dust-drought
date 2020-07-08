########################################################################################
# Project: Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison
# Prepared by: S. Arcusa
# Version: 1
# Date: 07-10-2019
# Part 0: Creating figure 5
#######################################################################################

library(readxl)
library(ncdf4)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
library(tidyr)
library(lipdR)
library(geoChronR)
library(tidyquant)

# Functions
duration = 35
mav <- function(x,n=duration){stats::filter(x,rep(1/n,n), sides=2)}

# Preparing data

# Tree ring data
PMDI  <- as.data.frame(read_xlsx("/home/sha59/FourCornersPMDI.xlsx"))
PMDI.df <- data.frame(Time = PMDI$'Time (CE)', PMDI = PMDI$'Four Courners PMDI (35°N to 37°N, -107°E to -112°E')
PMDI.35 <- mav(PMDI.df$PMDI)

# LME data for ensembles 002 and 003

setwd("/scratch/sha59/LME/NetCDF") #Change to directory of surfdata file
ncin  <- nc_open("surfdata_1.9x2.5_simyr1350_c131018.nc")

LON <- ncvar_get(ncin,"LONGXY")[,1]
LAT <- ncvar_get(ncin,"LATIXY")[1,]
Time <- seq(as.Date("850-01-01"), by = "month", length.out = 13872)
Years <- as.numeric(format(as.Date(Time), format = "%Y"))
Months <- as.numeric(format(as.Date(Time), format = "%m"))
nsteps <- 13872

xmin  <- which.min(abs(240-LON))
xmax  <- which.min(abs(251-LON))
ymin  <- which.min(abs(31-LAT))
ymax  <- which.min(abs(41-LAT))

setwd("/projects/pd_lab/sha59/LME/run002")
ncin <- nc_open("SOILWATER.002.085001-200512.nc")
soilwater002 <- ncvar_get(ncin, "SOILWATER_10CM")
SW_002 <- soilwater002[xmin:xmax,ymin:ymax,]
SW_av_002 <- vector()
for(m in 1:nsteps){
  SW_av_002[m] <- mean(SW_002[,,m], na.rm = T)
}

setwd("/projects/pd_lab/sha59/LME/run003")
ncin <- nc_open("SOILWATER.003.085001-200512.nc")
soilwater003 <- ncvar_get(ncin, "SOILWATER_10CM")
SW_003 <- soilwater003[xmin:xmax,ymin:ymax,]
SW_av_003 <- vector()
for(m in 1:nsteps){
  SW_av_003[m] <- mean(SW_003[,,m], na.rm = T)
}

remove(ncin)

### Analysis

### LME

SW_monthly <- data.frame(r002 = SW_av_002, r003 = SW_av_003, Years = Years)
SW_annual  <- aggregate(. ~ Years, FUN= mean, SW_monthly)
SW_35years  <- as.data.frame(sapply(SW_annual[,2:3], function(x) mav(x)))
SW_35years$Years = SW_annual$Years
LME <- data.frame(Year = c(SW_35years$Years, SW_35years$Years), 
                  Soil_moisture = c(SW_35years$r002, SW_35years$r003), 
                  Ensemble = c(rep("002",1156),rep("003",1156)), 
                  Period = c(rep(NA,1156), rep(NA,1156)))

soil.stats <- matrix(data = NA, nrow = 4, ncol = 2)
for(i in 1:2){
  soil.stats[1,i] <- mean(SW_35years[[i]],na.rm = T)
  soil.stats[2,i] <- sd(SW_35years[[i]], na.rm = T)/2
  soil.stats[3,i] <- mean(SW_annual[[i+1]],na.rm = T)
  soil.stats[4,i] <- sd(SW_annual[[i+1]], na.rm = T)/2
}

index.mega <- list()
for(i in 1:2){
  index.mega[[i]] <- which(SW_35years[[i]] <= soil.stats[1,i] - soil.stats[2,i])
}

cumsum.mega <- vector()
soilw_during_mega <- list()
mega.list <- list()
non.mega.list <- list()

for(i in 1:2){
  cumsum.mega <- cumsum(c(1,abs(index.mega[[i]][-length(index.mega[[i]])] - index.mega[[i]][-1]) >1))
  result <- table(cumsum.mega)
  group <- as.integer(which(result >=35)) 
  want.index <- which(cumsum.mega %in% group)
  inlist.index <- index.mega[[i]][want.index]
  mega.list[[i]] <- inlist.index
  non.mega.list[[i]] <- seq(1,1156,1)[-mega.list[[i]]]
  
  soilw_during_mega[[i]] <- SW_35years[[i]][inlist.index]
}

LME$Period[mega.list[[1]]] <- "MD" 
LME$Period[-mega.list[[1]]] <- "Other"
LME$Period[mega.list[[2]]+1156] <- "MD" 
LME$Period[-mega.list[[2]]+1156] <- "Other"

LME.long <- gather(LME, key, value, -Year, -Period, -Ensemble) 

#### PALEO MODEL
# LME drought

LME.sub <- filter(LME.long, Ensemble == "002")

LME_plot  <- ggplot(LME.sub, aes(x = Year, y = value))+
  geom_line(col = "darkgrey")+
  geom_path(aes(color = factor(Period),group = 1))+
  scale_color_manual(values = c( "MD" = "red", "Other" = "darkgrey"))+
  #facet_grid(Ensemble ~., scales = "free")+
  scale_x_continuous(expand=c(0,0))+
  geom_hline(yintercept = soil.stats[1,1], linetype = "solid", color = "black")+
  geom_hline(yintercept = soil.stats[1,1]- soil.stats[2,1],linetype = "dashed", color = "black")+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.9,0.9),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlab("Year (A.D.)")+
  labs(y = expression("Soil moisture (kg m"^{-2}*")"))

LME_plot


# LME dust
setwd("C:/Users/sha59/OneDrive - Northern Arizona University/PhD thesis/Cornell/Data")

LME_dust <- read.csv(file = "LME_ann_dep_sum.csv")
LME_dust$Year <- seq(850,2005,1)

LME_dust_plot <- ggplot(LME_dust, aes(x = Year, y = DSTDEP_2))+
  geom_line(col = "darkgrey")+
  geom_ma(n = 35, size = 1, linetype = 1, col = "black")+
  scale_x_continuous(expand=c(0,0))+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = -0.2, unit = "cm"))+
  xlab("Year (A.D.)")+
  labs(y = expression("Dust deposition (g m"^{-2}*"yr"^{-1}*")"))
LME_dust_plot

# Paleo model

modelmodel <- cowplot::plot_grid(LME_plot,LME_dust_plot, rel_heights=c(1.25, 1.25), ncol = 1, align = "v", labels = c("(c)", "(d)"))
print(modelmodel)

ggsave(laketree, filename = "model_model_fig.pdf", height = 7, width = 5, dpi = 300)

##### PALEO DATA
# Blue Lake

BL <- readLipd("C:/Users/sha59/Dropbox/Blue Lake/BlueLake.Routson.2018.lpd")
B <- runBacon(BL,labIDVar = NULL, age14CVar = "age14c",
              age14CuncertaintyVar = "age14cuncertainty", 
              ageVar = "age", ageUncertaintyVar = "ageuncertainty",
              depthVar = "depth", reservoirAge14CVar = NULL, 
              reservoirAge14CUncertaintyVar = NULL,
              rejectedAgesVar = NULL, baconThick = 3)

start = 2011
end = 0

core <-  mapAgeEnsembleToPaleoData(B)

max.ens <- 1000

ae = convertBP2AD(selectData(L = core, varName = "ageensemble"))
DBD = selectData(L = core, varName = "dbd")
DBDm <- matrix(DBD$values, nrow = length(DBD$values), ncol = max.ens)
depth = selectData(L = core, varName = "depth")

ageDiff = diff(ae$values)
depthDiff = diff(depth$values)
replicated.depth.matrix = matrix(data=depthDiff,nrow = length(depthDiff),ncol = ncol(ageDiff))
srEns = replicated.depth.matrix/ageDiff
empty.values = as.matrix(t(rep(NaN, 1000)))
srEns = rbind(srEns, empty.values)

MAR = DBDm * srEns

df <- selectData(L = core, varName = "dust")
DMAR <- -df$values * MAR * 10000

med <- vector()
for(i in 1:length(DMAR[,1])){
  med[i] <- median(DMAR[i,], na.rm = T)
}

bestfit <- rowMeans(DMAR, na.rm = T)
bstfit.df <- data.frame(BestFit = bestfit, Time = rowMeans(ae$values, na.rm = T), Median = med)

Blue.plot = plotTimeseriesEnsRibbons(X = ae, Y = DMAR, 
                                     x.bin = seq(end,start, by= 1), 
                                     probs = c(.1,.25,.5,.75,.9), 
                                     colorLow = "azure", colorHigh = "dark grey")+
  scale_x_continuous(expand = c(0,0), breaks=seq(0,2000,500))+
  labs(title = "", y = expression("Dust deposition (g m"^{-2}*"yr"^{-1}*")"))+
  xlab("Year (A.D.)")+
  coord_cartesian(ylim = c(0, 120), xlim = c(end,start))+
  scale_y_continuous(expand=c(0,0))+
  theme( axis.title.y.right = element_text(angle = 90),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = -0.2, unit = "cm"))+
  #geom_line(data = bstfit.df, aes(x = Time, y = BestFit), col = "grey32")+
  #geom_ma(data = bstfit.df, aes(x = Time, y = BestFit), n = 30,colour = "black", linetype = 1, size = 1,ma_fun = SMA)+
  geom_line(data = bstfit.df, aes(x = Time, y = Median), col = "black")
Blue.plot

# Trees

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

Trees  <- data.frame(Year = PMDI$'Time (CE)', PMDI = PMDI.35, Period = rep(NA,2018))
Trees$Period[inlist.index] <- "MD"
Trees$Period[-inlist.index]  <- "Other"

Trees.long <- gather(Trees, key, value, -Year, -Period)

Trees.plot <- ggplot(Trees.long, aes(x = Year, y = value))+
  geom_line(col = "darkgrey")+
  scale_x_continuous(expand=c(0,0))+
  geom_path(aes(color = factor(Period),group = 1))+
  scale_color_manual(values = c( "MD" = "red", "Other" = "darkgrey"))+
  theme_bw()+
  xlab("Year (A.D.)")+
  ylab("PMDI")+
  geom_hline(yintercept = PMDI.stats[1], linetype = "solid", color = "black")+
  geom_hline(yintercept = PMDI.stats[1]- PMDI.stats[2],linetype = "dashed", color = "black")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.8,0.2))
Trees.plot

# Paleo data

laketree <- cowplot::plot_grid( Trees.plot,Blue.plot, rel_heights=c(1.00, 1.25), ncol = 1, align = "v", labels = c("(a)", "(b)"))
print(laketree)

ggsave(laketree, filename = "paleo_data_fig.pdf", height = 7, width = 5, dpi = 300)


##### Paleo data/ model comparison

compa_vertical <- cowplot::plot_grid(modelmodel, laketree, cols = 1)
compa_horiz <- cowplot::plot_grid(modelmodel, laketree, cols = 2, align = "h")

ggsave(compa_vertical, filename = "fig_5.pdf", height = 7, width = 5)

