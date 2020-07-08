########################################################################################
# Project: Dust-Drought Nexus in the Southwestern United States: A Proxy-Model Comparison
# Prepared by: S. Arcusa
# Version: 1
# Date: 07-10-2019
# Part 1b: Script to create ensemble map for results of multiple linear regression
#######################################################################################

R2.all <- array(c(R2.002, R2.003, R2.004, R2.005, R2.006, R2.007, R2.008, R2.009, R2.010, R2.011, R2.012, R2.013), dim = c(144,96,12))
coef_SOILWATER.all<- array(c(coef_SOILWATER.002, coef_SOILWATER.003, coef_SOILWATER.004, coef_SOILWATER.005, coef_SOILWATER.006, coef_SOILWATER.007, coef_SOILWATER.008, coef_SOILWATER.009, coef_SOILWATER.010, coef_SOILWATER.011, coef_SOILWATER.012, coef_SOILWATER.013), dim = c(144,96,12))
coef_U10cubed.all<- array(c(coef_U10cubed.002, coef_U10cubed.003, coef_U10cubed.004, coef_U10cubed.005, coef_U10cubed.006, coef_U10cubed.007, coef_U10cubed.008, coef_U10cubed.009, coef_U10cubed.010, coef_U10cubed.011, coef_U10cubed.012, coef_U10cubed.013), dim = c(144,96,12))
coef_f_m.all<- array(c(coef_f_m.002, coef_f_m.003, coef_f_m.004, coef_f_m.005, coef_f_m.006, coef_f_m.007, coef_f_m.008, coef_f_m.009, coef_f_m.010, coef_f_m.011, coef_f_m.012, coef_f_m.013), dim = c(144,96,12))
lmg_SOILWATER.all<- array(c(lmg_SOILWATER.002, lmg_SOILWATER.003, lmg_SOILWATER.004, lmg_SOILWATER.005, lmg_SOILWATER.006, lmg_SOILWATER.007, lmg_SOILWATER.008, lmg_SOILWATER.009, lmg_SOILWATER.010, lmg_SOILWATER.011, lmg_SOILWATER.012, lmg_SOILWATER.013), dim = c(144,96,12))
lmg_U10cubed.all<- array(c(lmg_U10cubed.002, lmg_U10cubed.003, lmg_U10cubed.004, lmg_U10cubed.005, lmg_U10cubed.006, lmg_U10cubed.007, lmg_U10cubed.008, lmg_U10cubed.009, lmg_U10cubed.010, lmg_U10cubed.011, lmg_U10cubed.012, lmg_U10cubed.013), dim = c(144,96,12))
lmg_f_m.all<- array(c(lmg_f_m.002, lmg_f_m.003, lmg_f_m.004, lmg_f_m.005, lmg_f_m.006, lmg_f_m.007, lmg_f_m.008, lmg_f_m.009, lmg_f_m.010, lmg_f_m.011, lmg_f_m.012, lmg_f_m.013), dim = c(144,96,12))

R2.all[R2.all == "NaN"] <- NA
coef_SOILWATER.all[coef_SOILWATER.all == "NaN"] <- NA
coef_U10cubed.all[coef_U10cubed.all == "NaN"] <- NA
coef_f_m.all[coef_f_m.all == "NaN"] <- NA
lmg_SOILWATER.all[lmg_SOILWATER.all == "NaN"] <- NA
lmg_U10cubed.all[lmg_U10cubed.all == "NaN"] <- NA
lmg_f_m.all[lmg_f_m.all == "NaN"] <- NA

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
    
  }
}

# Plotting
library(raster)
library(rgdal)
library(maps)
library(rasterVis)
library(viridis)
library(ggplot2)

Mybrks <- seq(0,1,by=0.1)
nb <- length(Mybrks)-1

XFrm <- -180
XTo <- 180
YFrm <- -90
YTo <- 90
XStp <- 60
YStp <- 30

R2.map <- raster(t(R2), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                 crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
R2.map <- flip(R2.map, direction = "y")
R2.map <- rotate(R2.map)
R2.map <- trim(R2.map)

SOILWATER.map <- raster(t(lmg_SOILWATER), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                        crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
SOILWATER.map <- flip(SOILWATER.map, direction = "y")
SOILWATER.map <- rotate(SOILWATER.map)
plot(SOILWATER.map)

U10cubed.map <- raster(t(lmg_U10cubed), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                       crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
U10cubed.map <- flip(U10cubed.map, direction = "y")
U10cubed.map <- rotate(U10cubed.map)

f_m.map <- raster(t(lmg_f_m), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                  crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
f_m.map <- flip(f_m.map, direction = "y")
f_m.map <- rotate(f_m.map)

R2sd.map <- raster(t(R2.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                 crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
R2sd.map <- flip(R2sd.map, direction = "y")
R2sd.map <- rotate(R2sd.map)

SOILWATERsd.map <- raster(t(lmg_SOILWATER.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                        crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
SOILWATERsd.map <- flip(SOILWATERsd.map, direction = "y")
SOILWATERsd.map <- rotate(SOILWATERsd.map)
plot(SOILWATERsd.map)

U10cubedsd.map <- raster(t(lmg_U10cubed.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                       crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
U10cubedsd.map <- flip(U10cubedsd.map, direction = "y")
U10cubedsd.map <- rotate(U10cubedsd.map)

f_msd.map <- raster(t(lmg_f_m.sd), xmn = min(lon), xmx = max(lon), ymx = max(lat), ymn = min(lat),
                  crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
f_msd.map <- flip(f_msd.map, direction = "y")
f_msd.map <- rotate(f_msd.map)

pdf("Ens_R2.pdf")
#png("Ens_R2.png", width = 350, height = 350)
plot(R2.map,col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_SOILWATER.pdf")
#png("Ens_SOILWATER.png", width = 350, height = 350)
plot(SOILWATER.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_U10cubed.pdf")
#png("Ens_U10cubed.png", width = 350, height = 350)
plot(U10cubed.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_f_m.pdf")
#png("Ens_f_m.png", width = 350, height = 350)
plot(f_m.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

pdf("Ens_corr_map.pdf")
#png("Ens_corr_map.png", width = 350, height = 350)
s <- stack(SOILWATER.map, f_m.map, U10cubed.map)
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

#### SD maps

#pdf("Ens_R2_sd.pdf")
png("Ens_R2_sd.png", width = 350, height = 350)
plot(R2sd.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

#pdf("Ens_SOILWATER_sd.pdf")
png("Ens_SOILWATER_sd.png", width = 350, height = 350)
plot(SOILWATERsd.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

#pdf("Ens_U10cubed_sd.pdf")
png("Ens_U10cubed_sd.png", width = 350, height = 350)
plot(U10cubedsd.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

#pdf("Ens_f_m_sd.pdf")
png("Ens_f_m_sd.png", width = 350, height = 350)
plot(f_msd.map, col = rev(viridis(nb)), axes=F,box=F, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
axis(1,tick=T,pos=YFrm, las =1, at=seq(XFrm,XTo,XStp))
axis(2,tick=T,pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=YTo, las =1, at = seq(XFrm,XTo,XStp))
axis(4,tick =T, labels= F, pos=XTo, at=seq(YFrm,YTo,YStp))
map("world", add=T, fill=F)
legend("right", x = XTo, y = YTo)
dev.off()

# Analysis for the Southwest
xmin  <- which.min(abs(240-lon))
xmax  <- which.min(abs(251-lon))
ymin  <- which.min(abs(31-lat))
ymax  <- which.min(abs(41-lat))

SW_soilwater <- lmg_SOILWATER[xmin:xmax,ymin:ymax]
mean(SW_soilwater, na.rm = T)
sd(SW_soilwater, na.rm = T)

SW_f_m <- lmg_f_m[xmin:xmax,ymin:ymax]
mean(SW_f_m, na.rm = T)
sd(SW_f_m, na.rm = T)

SW_U10cubed <- lmg_U10cubed[xmin:xmax,ymin:ymax]
mean(SW_U10cubed, na.rm = T)
sd(SW_U10cubed, na.rm = T)

SW_soilwater.sd <- sd(SW_soilwater.all, na.rm = T)
SW_soilwater.all <- lmg_SOILWATER.all[xmin:xmax,ymin:ymax,]

Gobi_Tak <- extent(90,110,30,40)
Gobi_Tak_map <- raster::crop(SOILWATER.map, Gobi_Tak)
mean(Gobi_Tak_map@data@values, na.rm = T)

high <- (length(which(R2 >= 0.5))/(length(which(is.na(R2) == F))))*100

# Plot the Southwest
library(maps)
SW.region = extent(-120,-109,31,41) #used to be 240,251 (use for the MEME)
SW.region = extent(240,251,31,41)
SW.SOILWATER <- raster::crop(SOILWATER.map, SW.region)
SW.R2 <- raster::crop(R2.map, SW.region)
SW.fm <- raster::crop(f_m.map, SW.region)
SW.U10 <- raster::crop(U10cubed.map, SW.region)
usa <- map("state", add =T)
#usa$x  <- 360 + usa$x #use for MEME

pdf("SW_LME_soil_moisture.pdf")
#png("SW_LME_soil_moisture.png", width = 350, height = 350)
plot(SW.SOILWATER, col = rev(viridis(8)), breaks = seq(from=0,to=0.14,by=0.02), lab.breaks = seq(from=0,to=0.14,by=0.02), zlim=c(0,0.14))
map(usa, add = T)
map("world", add = T, fill = F) #use world 2 for MEME
dev.off()

pdf("SW_LME_fm.pdf")
#png("SW_LME_fm.png", width = 350, height = 350)
plot(SW.fm, col = rev(viridis(10)), breaks = seq(from=0.7,to=1,by=0.05), lab.breaks = seq(from=0.7,to=1,by=0.05), zlim=c(0.7,1))
map(usa, add = T)
map("world", add = T, fill = F) #use world 2 for MEME
dev.off()

pdf("SW_LME_U10cubed.pdf")
png("SW_LME_U10cubed.png", width = 350, height = 350)
plot(SW.U10, col = rev(viridis(10)), breaks = seq(from=0,to=0.2,by=0.02), lab.breaks = seq(from=0,to=0.2,by=0.02), zlim=c(0,0.2))
map(usa, add = T)
map("world", add = T, fill = F) #use world 2 for MEME
dev.off()

pdf("SW_LME_R2.pdf")
#png("SW_LME_R2.png", width = 350, height = 350)
#plot(SW.SOILWATER, main = "LME Soil Moisture")
plot(SW.R2,col = rev(viridis(8)), axes=F,box=T, breaks = seq(from=0.65,to=1,by=0.05), lab.breaks = seq(from=0.65,to=1,by=0.05), zlim=c(0.65,1))
map(usa, add = T)
map("world", add = T, fill = F) #use world 2 for MEME
axis(1,tick=T,pos=31.88, las =1, at=seq(-118,-110,2))
axis(2,tick=T,pos=-119.90, las =1, at=seq(32,40,2))
dev.off()

png("SW_ens_corr_map.png", width = 350, height = 350)
#png("Ens_corr_map.png", width = 350, height = 350)
a <- stack(SW.SOILWATER, SW.fm, SW.U10)
cor.map <- plotRGB(a, stretch = 'hist', axes = T)
axis(1,tick=T,pos=31, las =1, at=seq(-120,-108,2))
#axis(2,tick=F, labels = F, pos=XFrm, las =1, at=seq(YFrm,YTo,YStp))
axis(3, tick =T, labels= F, pos=42, las =1, at = seq(-120,-109,2))
axis(4,tick =T, labels= T, pos=-109, at=seq(32,41,2))
map("world", add=T, fill=F)
map(usa, add = T)
dev.off()

# LME with 0-1 scale

par(mfrow=c(1,3))
pdf("SW_LME_all_var_0-1.pdf")
#png("SW_LME_all_var_0-1.png", width = 350, height = 350)

plot(SW.SOILWATER, col = rev(viridis(nb)), axes=T,box=T, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
title("Soil moisture")
map(usa, add = T)
map("world", add = T, fill = F) #use world 2 for MEME

plot(SW.fm, col = rev(viridis(nb)), axes=T,box=T, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
title("Bare ground exposure")
map(usa, add = T)
map("world", add = T, fill = F) #use world 2 for MEME

plot(SW.U10, col = rev(viridis(nb)), axes=T,box=T, breaks = Mybrks, lab.breaks = Mybrks, zlim=c(0,1))
title("Wind speed cubed")
map(usa, add = T)
map("world", add = T, fill = F) #use world 2 for MEME
dev.off()
