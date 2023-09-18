
rm(list = ls())
# Clear work space

setwd("//nmbu.no/my/home/Desktop/DATA")
#load data 

Bats45<-read.csv("45%.csv",sep=";", dec=",", header=TRUE)

## Add column for species ID ##
mystacinus<- Bats45$BatID %in% c("Amelia","Marie", "Turid", "Reeda", "Daisy", "Line", "Stine", "Ethel", "Ragnhild", "Ida", "Louise", "Aricia")

Bats45$Species<-ifelse(mystacinus == TRUE, "MYSTACINUS", "BRANDTII")

# Set new column with Species as factor instead of character 

Bats45$fSpecies<- as.factor(Bats45$Species)
str(Bats45)
names(Bats45)
#[1] "Date"     "Time"     "BatID"    "Easting"  "Northing" "BD"       "Dir.obs"  "Gain"     "Species"  "fSpecies"
summary(Bats45$BatID)
#Amelia   Aricia   Astrid    Dagny    Daisy    Ethel      Ida     Kaja     Line   Louise    Maren    Marie     Nora   Phoebe Ragnhild    Reeda 
#16       90       81       83      134      117      119       91       55       64      141       31       49       62       28       70 
#Sofia   Steffi    Stine     Thea    Turid 
#39       39       84      101        8 

##The number of onsite plots per individual bat 

table(Bats45$Dir.obs, Bats45$BatID)
#Amelia Aricia Astrid Dagny Daisy Ethel Ida Kaja Line Louise Maren Marie Nora Phoebe Ragnhild Reeda Sofia Steffi Stine Thea Turid
#No       11     49     59    40    98    79 106   52   37     29   110    25   41     56       28    51    26     34    58   43     8
#Yes       5     41     22    43    36    38  13   39   18     35    31     6    8      6        0    19    13      5    26   58     0
#Yes = onsite plots

table(Bats45$fSpecies, Bats45$BatID)


M.mystacinus <- subset(Bats45, fSpecies == "MYSTACINUS")
M.brandtii <- subset(Bats45, fSpecies == "BRANDTII")

par(mfrow=c(2,1))
plot(table(M.mystacinus$BatID), main = "Number of plots Mmys", xaxt = "n")
plot(table(M.brandtii$BatID), main = "Number of plots Mbra")

#Fix date/time

names(Bats45)[1:2]<-c("date.cap","time.cap")

Bats45$date.cap<-as.character(Bats45$date.cap)

Bats45$time.cap<-as.character(Bats45$time.cap)

#---note: looks like many records are missing date and time data (this is rectified further down in the script)

Bats45$datetime<-paste(Bats45$date.cap,Bats45$time.cap)

Bats45$datetime.posix<-as.POSIXct(strptime(Bats45$datetime,"%d.%m.%y %H:%M"),tz="GMT")

summary(Bats45$fSpecies)
summary(Bats45$Dir.obs=="Yes"|Bats45$BD=="Yes")
summary(Bats45$Dir.obs=="Yes")
summary(Bats45$BD=="Yes")


#Load packages

library(sp)
library(sf)
library(rgeos)
library(adehabitatMA)
library(digest)
library(mime)
library(jsonlite)
library(crs)
library(dplyr)
library(adehabitatHR)
library(Rcpp)
library(raster)
library(rgdal)
library(munsell)
library(tibble)
library(pillar)
library(ggplot2)
library(ResourceSelection)
library(mapview)
library(amt)
library(sjPlot)
library(Rgb)
library(nloptr)
library(colorspace)
library(rlang)
library(broom) 
library(lattice)
library(MASS)


## Create spatial data frame ##

coords <- coordinates(Bats45)<-Bats45[,c("Easting","Northing")]
#assigning coordinates 

proj4string(Bats45)<-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#Define the projection 

Bats45<-Bats45[!is.na(Bats45$Easting) | !is.na(Bats45$Northing),]
#Remove NA values 

#Run again to get per species datasets without NA's 
Mystacinus <- subset(Bats45, fSpecies == "MYSTACINUS")
Brandtii <- subset(Bats45, fSpecies == "BRANDTII")



###########################################################################################
######################      Import Raster Data (LiDAR) ####################################
###########################################################################################


#The variables are explained here:
#https://github.com/Jean-Romain/lidR/wiki/stdmetrics

##---------- 250 m2 
stack.r250<- stack("stdmetrics_z.grd")
object250<- stack.r250
names(object250) 


###############################################################################################
###   Used versus available habitat and Mmus versus Mbra - exploratory analyses             ###
###############################################################################################


bat.locs <- Bats45 

table(bat.locs$BatID, bat.locs$fSpecies)

#          BRANDTII MYSTACINUS
#Amelia          0         16
#Aricia          0         90
#Astrid         81          0
#Dagny          83          0
#Daisy           0        134
#Ethel           0        117
#Ida             0        119
#Kaja           91          0
#Line            0         55
#Louise          0         64
#Maren         141          0
#Marie           0         31
#Nora           49          0
#Phoebe         62          0
#Ragnhild        0         28
#Reeda           0         70
#Sofia          39          0
#Steffi         39          0
#Stine           0         84
#Thea          101          0
#Turid           0          8


xy.spatial <- SpatialPoints(bat.locs[,c("Easting","Northing")]) 
myobject <- crop(object250, extent(xy.spatial)) #Cropping the lidar data to only the area where we have bat observations

#Get the lidar variables from 'myobject'
zmax <- myobject$zmax
zmean   <- myobject$zmean  
zsd <- myobject$zsd
zskew    <- myobject$zskew  
zkurt    <- myobject$zkurt    
zentropy   <- myobject$zentropy
pzabovezmean <- myobject$pzabovezmean 
pzabove0.5   <- myobject$pzabove0.5   
zq10 <- myobject$zq10    #Drop z5, z15, z25, z35, z45, z55, z65, z75, z85 (do not need all "levels" - the are very similar)
zq20 <- myobject$zq20     
zq30 <- myobject$zq30   
zq40 <- myobject$zq40   
zq50 <- myobject$zq50   
zq60 <- myobject$zq60   
zq70 <- myobject$zq70   
zq80 <- myobject$zq80   
zq90 <- myobject$zq90   
#Dropped all the cum variables too

####################################

#Exploration of lidar variables - which show strong signals;
#i) either clear difference between used and available, and/or ii) clear difference between species

#Raster analysis - comparing Mystacinus and Brandtii - used and available for each lidar variable
#Population/species level (used an available for all individuals in each species pooled)
#Drop z5, z15, z25, z35, z45, z55, z65, z75; use 10-intervals only (likely high correlation among variables)


#----------------zmax

#Mystacinus
#xy.spatial.M <- SpatialPoints(Mystacinus[,c("Easting","Northing")])
#myobject.M <- crop(object250, extent(xy.spatial.M)) 
#myzmax.M <- crop(myobject.M$zmax, extent(xy.spatial.M))
#Mmys.mcp <- mcp(xy.spatial.M, 95) #I used 95% MCP here
#random.points.M <- spsample(Mmys.mcp, 1000, "random")
#Mystacinus$zmax <- extract(myzmax.M, xy.spatial.M)
#random.zmax.M <- extract(myzmax.M, random.points.M)

#Brandtii
#xy.spatial.B <- SpatialPoints(Brandtii[,c("Easting","Northing")])
#myobject.B <- crop(object250, extent(xy.spatial.B)) 
#myzmax.B <- crop(myobject.B$zmax, extent(xy.spatial.B))
#Mbra.mcp <- mcp(xy.spatial.B, 95) #I used 95% MCP here
#random.points.B <- spsample(Mbra.mcp, 1000, "random")
#Brandtii$zmax <- extract(myzmax.B, xy.spatial.B)
#random.zmax.B <- extract(myzmax.B, random.points.B)

# overlapping histograms?

#Compare distribution of
# 1. zmax values extracted from the whole area (square) defined by 
# the eastern/westermost "Easting" values and the northern and southernmost "Northing" values (one measure of availability)
# 2. zmax values in random points within the "pooled" MCP for all the Mystacinus bats (another measure of availability)
# 3. zmax values for actual observations of Mystacinus bats

#----------------zmax, max height
#windows()
#par(mfrow=c(1,2))
#hist(myzmax.M, col="grey", breaks = seq(0,40,1), freq=FALSE, bor="darkgrey", ylim = c(0,0.3), main = 'M. mystacinus', xlab ="zmax 250 m2") 
#hist(random.zmax.M, col=rgb(0,1,0,.5), breaks = seq(0,40,1), freq=FALSE, bor="darkgrey", ylim = c(0,0.3), add=TRUE, main = 'n')
#hist(Mystacinus$zmax, col=rgb(1,0,0,.5), breaks = seq(0,40,1), freq=FALSE, add=TRUE, bor = "red")

#hist(myzmax.B, col="grey", breaks = seq(0,40,1), freq=FALSE, bor="darkgrey", ylim = c(0,0.3), main = 'M. brandtii', xlab ="zmax 250 m2")
#hist(random.zmax.B, col=rgb(0,1,0,.5), breaks = seq(0,40,1), freq=FALSE, bor="darkgrey", ylim = c(0,0.3), add=TRUE, main = 'n')
#hist(Brandtii$zmax, col=rgb(1,0,0,.5), breaks = seq(0,40,1), freq=FALSE, add=TRUE, bor = "red")



    ##############################################################################################################
#### Figures with moth important LiDAR variables and pooled 95% MCPs for April's master  #####################
##############################################################################################################

#windows()
#zmax 

#plot(zmax)
#lines(Mbra.mcp, col="darkblue", lwd=2)
#lines(Mmys.mcp, col="darkred", lwd=2)
#mtext("zmax 250m2")



################################################################################################################
###------------- species and individual-level analyses -------------------------------------------##############
################################################################################################################



###############################################################
###-------------  M. mystacinus bats  ---------------------####
###############################################################

#------Amelia

Amelia.locs <- subset(bat.locs, BatID == "Amelia")

#Determine available space (simplistically: MCP)
Amelia.mcp <- mcp(SpatialPoints(Amelia.locs), 95)
Amelia.mcp.xy <- Amelia.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Amelia.mcp.xy <- fortify(Amelia.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Amelia.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Amelia.mcp, col="darkred", lwd=2)
mtext("Mmys Amelia zmax res250")

#Create data frame of covariates

names(Amelia.locs)

Amelia.sp<-SpatialPoints(Amelia.locs)


used.Amelia.df <- data.frame(
  used = rep(1, length(Amelia.sp)),
  Species = rep("Mys", length(Amelia.sp)),
  BatID = rep("Amelia", length(Amelia.sp)),
  zmax = extract(zmax, Amelia.sp),
  zmean = extract(zmean, Amelia.sp),
  zsd = extract(zsd, Amelia.sp),
  zentropy = extract(zentropy, Amelia.sp),
  pzabovezmean = extract(pzabovezmean, Amelia.sp),
  pzabove0.5 = extract(pzabove0.5, Amelia.sp),
  zq10 = extract(zq10, Amelia.sp),
  zq20 = extract(zq20, Amelia.sp),
  zq30 = extract(zq30, Amelia.sp),
  zq40 = extract(zq40, Amelia.sp),
  zq50 = extract(zq50, Amelia.sp),
  zq60 = extract(zq60, Amelia.sp),
  zq70 = extract(zq70, Amelia.sp),
  zq80 = extract(zq80, Amelia.sp),
  zq90 = extract(zq90, Amelia.sp))

head(used.Amelia.df)
summary(used.Amelia.df)
dim(used.Amelia.df)



#######################################
#Obtain some random points within mcp

random.points.Amelia <- spsample(Amelia.mcp, 16, "random")

#Calculate the covariates for the null set

avail.Amelia.df <- data.frame(
  used = rep(0, length(random.points.Amelia)),
  Species = rep("Mys", length(Amelia.sp)),
  BatID = rep("Amelia", length(Amelia.sp)),
  zmax = extract(zmax, random.points.Amelia),
  zmean = extract(zmean, random.points.Amelia),
  zsd = extract(zsd, random.points.Amelia),
  zentropy = extract(zentropy, random.points.Amelia),
  pzabovezmean = extract(pzabovezmean, random.points.Amelia),
  pzabove0.5 = extract(pzabove0.5, random.points.Amelia),
  zq10 = extract(zq10, random.points.Amelia),
  zq20 = extract(zq20, random.points.Amelia),
  zq30 = extract(zq30, random.points.Amelia),
  zq40 = extract(zq40, random.points.Amelia),
  zq50 = extract(zq50, random.points.Amelia),
  zq60 = extract(zq60, random.points.Amelia),
  zq70 = extract(zq70, random.points.Amelia),
  zq80 = extract(zq80, random.points.Amelia),
  zq90 = extract(zq90, random.points.Amelia))

head(avail.Amelia.df)
summary(avail.Amelia.df)
dim(avail.Amelia.df)



#------Aricia

Aricia.locs <- subset(bat.locs, BatID == "Aricia")

#Determine available space (simplistically: MCP)
Aricia.mcp <- mcp(SpatialPoints(Aricia.locs), 95)
Aricia.mcp.xy <- Aricia.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Aricia.mcp.xy <- fortify(Aricia.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Aricia.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Aricia.mcp, col="darkred", lwd=2)
mtext("Mmys Aricia zmax res250")

#Create data frame of covariates

names(Aricia.locs)

Aricia.sp<-SpatialPoints(Aricia.locs)


used.Aricia.df <- data.frame(
  used = rep(1, length(Aricia.sp)),
  Species = rep("Mys", length(Aricia.sp)),
  BatID = rep("Aricia", length(Aricia.sp)),
  zmax = extract(zmax, Aricia.sp),
  zmean = extract(zmean, Aricia.sp),
  zsd = extract(zsd, Aricia.sp),
  zentropy = extract(zentropy, Aricia.sp),
  pzabovezmean = extract(pzabovezmean, Aricia.sp),
  pzabove0.5 = extract(pzabove0.5, Aricia.sp),
  zq10 = extract(zq10, Aricia.sp),
  zq20 = extract(zq20, Aricia.sp),
  zq30 = extract(zq30, Aricia.sp),
  zq40 = extract(zq40, Aricia.sp),
  zq50 = extract(zq50, Aricia.sp),
  zq60 = extract(zq60, Aricia.sp),
  zq70 = extract(zq70, Aricia.sp),
  zq80 = extract(zq80, Aricia.sp),
  zq90 = extract(zq90, Aricia.sp))

head(used.Aricia.df)
summary(used.Aricia.df)
dim(used.Aricia.df)



#######################################
#Obtain some random points within mcp

random.points.Aricia <- spsample(Aricia.mcp, 90, "random")

#Calculate the covariates for the null set

avail.Aricia.df <- data.frame(
  used = rep(0, length(random.points.Aricia)),
  Species = rep("Mys", length(Aricia.sp)),
  BatID = rep("Aricia", length(Aricia.sp)),
  zmax = extract(zmax, random.points.Aricia),
  zmean = extract(zmean, random.points.Aricia),
  zsd = extract(zsd, random.points.Aricia),
  zentropy = extract(zentropy, random.points.Aricia),
  pzabovezmean = extract(pzabovezmean, random.points.Aricia),
  pzabove0.5 = extract(pzabove0.5, random.points.Aricia),
  zq10 = extract(zq10, random.points.Aricia),
  zq20 = extract(zq20, random.points.Aricia),
  zq30 = extract(zq30, random.points.Aricia),
  zq40 = extract(zq40, random.points.Aricia),
  zq50 = extract(zq50, random.points.Aricia),
  zq60 = extract(zq60, random.points.Aricia),
  zq70 = extract(zq70, random.points.Aricia),
  zq80 = extract(zq80, random.points.Aricia),
  zq90 = extract(zq90, random.points.Aricia))


head(avail.Aricia.df)
summary(avail.Aricia.df)
dim(avail.Aricia.df)



#------Daisy

Daisy.locs <- subset(bat.locs, BatID == "Daisy")

#Determine available space (simplistically: MCP)
Daisy.mcp <- mcp(SpatialPoints(Daisy.locs), 95)
Daisy.mcp.xy <- Daisy.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Daisy.mcp.xy <- fortify(Daisy.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Daisy.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Daisy.mcp, col="darkred", lwd=2)
mtext("Mmys Daisy zmax res250")

#Create data frame of covariates

names(Daisy.locs)

Daisy.sp<-SpatialPoints(Daisy.locs)


used.Daisy.df <- data.frame(
  used = rep(1, length(Daisy.sp)),
  Species = rep("Mys", length(Daisy.sp)),
  BatID = rep("Daisy", length(Daisy.sp)),
  zmax = extract(zmax, Daisy.sp),
  zmean = extract(zmean, Daisy.sp),
  zsd = extract(zsd, Daisy.sp),
  zentropy = extract(zentropy, Daisy.sp),
  pzabovezmean = extract(pzabovezmean, Daisy.sp),
  pzabove0.5 = extract(pzabove0.5, Daisy.sp),
  zq10 = extract(zq10, Daisy.sp),
  zq20 = extract(zq20, Daisy.sp),
  zq30 = extract(zq30, Daisy.sp),
  zq40 = extract(zq40, Daisy.sp),
  zq50 = extract(zq50, Daisy.sp),
  zq60 = extract(zq60, Daisy.sp),
  zq70 = extract(zq70, Daisy.sp),
  zq80 = extract(zq80, Daisy.sp),
  zq90 = extract(zq90, Daisy.sp))


head(used.Daisy.df)
summary(used.Daisy.df)
dim(used.Daisy.df)



#######################################
#Obtain some random points within mcp

random.points.Daisy <- spsample(Daisy.mcp, 134, "random")

#Calculate the covariates for the null set

avail.Daisy.df <- data.frame(
  used = rep(0, length(random.points.Daisy)),
  Species = rep("Mys", length(Daisy.sp)),
  BatID = rep("Daisy", length(Daisy.sp)),
  zmax = extract(zmax, random.points.Daisy),
  zmean = extract(zmean, random.points.Daisy),
  zsd = extract(zsd, random.points.Daisy),
  zentropy = extract(zentropy, random.points.Daisy),
  pzabovezmean = extract(pzabovezmean, random.points.Daisy),
  pzabove0.5 = extract(pzabove0.5, random.points.Daisy),
  zq10 = extract(zq10, random.points.Daisy),
  zq20 = extract(zq20, random.points.Daisy),
  zq30 = extract(zq30, random.points.Daisy),
  zq40 = extract(zq40, random.points.Daisy),
  zq50 = extract(zq50, random.points.Daisy),
  zq60 = extract(zq60, random.points.Daisy),
  zq70 = extract(zq70, random.points.Daisy),
  zq80 = extract(zq80, random.points.Daisy),
  zq90 = extract(zq90, random.points.Daisy))


head(avail.Daisy.df)
summary(avail.Daisy.df)
dim(avail.Daisy.df)


#------Ethel

Ethel.locs <- subset(bat.locs, BatID == "Ethel")

#Determine available space (simplistically: MCP)
Ethel.mcp <- mcp(SpatialPoints(Ethel.locs), 95)
Ethel.mcp.xy <- Ethel.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Ethel.mcp.xy <- fortify(Ethel.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Ethel.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Ethel.mcp, col="darkred", lwd=2)
mtext("Mmys Ethel zmax res250")

#Create data frame of covariates

names(Ethel.locs)

Ethel.sp<-SpatialPoints(Ethel.locs)


used.Ethel.df <- data.frame(
  used = rep(1, length(Ethel.sp)),
  Species = rep("Mys", length(Ethel.sp)),
  BatID = rep("Ethel", length(Ethel.sp)),
  zmax = extract(zmax, Ethel.sp),
  zmean = extract(zmean, Ethel.sp),
  zsd = extract(zsd, Ethel.sp),
  zentropy = extract(zentropy, Ethel.sp),
  pzabovezmean = extract(pzabovezmean, Ethel.sp),
  pzabove0.5 = extract(pzabove0.5, Ethel.sp),
  zq10 = extract(zq10, Ethel.sp),
  zq20 = extract(zq20, Ethel.sp),
  zq30 = extract(zq30, Ethel.sp),
  zq40 = extract(zq40, Ethel.sp),
  zq50 = extract(zq50, Ethel.sp),
  zq60 = extract(zq60, Ethel.sp),
  zq70 = extract(zq70, Ethel.sp),
  zq80 = extract(zq80, Ethel.sp),
  zq90 = extract(zq90, Ethel.sp))

head(used.Ethel.df)
summary(used.Ethel.df)
dim(used.Ethel.df)



#######################################
#Obtain some random points within mcp

random.points.Ethel <- spsample(Ethel.mcp, 117, "random")

#Calculate the covariates for the null set

avail.Ethel.df <- data.frame(
  used = rep(0, length(random.points.Ethel)),
  Species = rep("Mys", length(Ethel.sp)),
  BatID = rep("Ethel", length(Ethel.sp)),
  zmax = extract(zmax, random.points.Ethel),
  zmean = extract(zmean, random.points.Ethel),
  zsd = extract(zsd, random.points.Ethel),
  zentropy = extract(zentropy, random.points.Ethel),
  pzabovezmean = extract(pzabovezmean, random.points.Ethel),
  pzabove0.5 = extract(pzabove0.5, random.points.Ethel),
  zq10 = extract(zq10, random.points.Ethel),
  zq20 = extract(zq20, random.points.Ethel),
  zq30 = extract(zq30, random.points.Ethel),
  zq40 = extract(zq40, random.points.Ethel),
  zq50 = extract(zq50, random.points.Ethel),
  zq60 = extract(zq60, random.points.Ethel),
  zq70 = extract(zq70, random.points.Ethel),
  zq80 = extract(zq80, random.points.Ethel),
  zq90 = extract(zq90, random.points.Ethel))


head(avail.Ethel.df)
summary(avail.Ethel.df)
dim(avail.Ethel.df)




#------Ida

Ida.locs <- subset(bat.locs, BatID == "Ida")

#Determine available space (simplistically: MCP)
Ida.mcp <- mcp(SpatialPoints(Ida.locs), 95)
Ida.mcp.xy <- Ida.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Ida.mcp.xy <- fortify(Ida.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Ida.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Ida.mcp, col="darkred", lwd=2)
mtext("Mmys Ida zmax res250")

#Create data frame of covariates

names(Ida.locs)

Ida.sp<-SpatialPoints(Ida.locs)


used.Ida.df <- data.frame(
  used = rep(1, length(Ida.sp)),
  Species = rep("Mys", length(Ida.sp)),
  BatID = rep("Ida", length(Ida.sp)),
  zmax = extract(zmax, Ida.sp),
  zmean = extract(zmean, Ida.sp),
  zsd = extract(zsd, Ida.sp),
  zentropy = extract(zentropy, Ida.sp),
  pzabovezmean = extract(pzabovezmean, Ida.sp),
  pzabove0.5 = extract(pzabove0.5, Ida.sp),
  zq10 = extract(zq10, Ida.sp),
  zq20 = extract(zq20, Ida.sp),
  zq30 = extract(zq30, Ida.sp),
  zq40 = extract(zq40, Ida.sp),
  zq50 = extract(zq50, Ida.sp),
  zq60 = extract(zq60, Ida.sp),
  zq70 = extract(zq70, Ida.sp),
  zq80 = extract(zq80, Ida.sp),
  zq90 = extract(zq90, Ida.sp))


head(used.Ida.df)
summary(used.Ida.df)
dim(used.Ida.df)



#######################################
#Obtain some random points within mcp

random.points.Ida <- spsample(Ida.mcp, 119, "random")

#Calculate the covariates for the null set

avail.Ida.df <- data.frame(
  used = rep(0, length(random.points.Ida)),
  Species = rep("Mys", length(Ida.sp)),
  BatID = rep("Ida", length(Ida.sp)),
  zmax = extract(zmax, random.points.Ida),
  zmean = extract(zmean, random.points.Ida),
  zsd = extract(zsd, random.points.Ida),
  zentropy = extract(zentropy, random.points.Ida),
  pzabovezmean = extract(pzabovezmean, random.points.Ida),
  pzabove0.5 = extract(pzabove0.5, random.points.Ida),
  zq10 = extract(zq10, random.points.Ida),
  zq20 = extract(zq20, random.points.Ida),
  zq30 = extract(zq30, random.points.Ida),
  zq40 = extract(zq40, random.points.Ida),
  zq50 = extract(zq50, random.points.Ida),
  zq60 = extract(zq60, random.points.Ida),
  zq70 = extract(zq70, random.points.Ida),
  zq80 = extract(zq80, random.points.Ida),
  zq90 = extract(zq90, random.points.Ida))


head(avail.Ida.df)
summary(avail.Ida.df)
dim(avail.Ida.df)



#------Line

Line.locs <- subset(bat.locs, BatID == "Line")

#Determine available space (simplistically: MCP)
Line.mcp <- mcp(SpatialPoints(Line.locs), 95)
Line.mcp.xy <- Line.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Line.mcp.xy <- fortify(Line.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Line.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Line.mcp, col="darkred", lwd=2)
mtext("Mmys Line zmax res250")

#Create data frame of covariates

names(Line.locs)

Line.sp<-SpatialPoints(Line.locs)


used.Line.df <- data.frame(
  used = rep(1, length(Line.sp)),
  Species = rep("Mys", length(Line.sp)),
  BatID = rep("Line", length(Line.sp)),
  zmax = extract(zmax, Line.sp),
  zmean = extract(zmean, Line.sp),
  zsd = extract(zsd, Line.sp),
  zentropy = extract(zentropy, Line.sp),
  pzabovezmean = extract(pzabovezmean, Line.sp),
  pzabove0.5 = extract(pzabove0.5, Line.sp),
  zq10 = extract(zq10, Line.sp),
  zq20 = extract(zq20, Line.sp),
  zq30 = extract(zq30, Line.sp),
  zq40 = extract(zq40, Line.sp),
  zq50 = extract(zq50, Line.sp),
  zq60 = extract(zq60, Line.sp),
  zq70 = extract(zq70, Line.sp),
  zq80 = extract(zq80, Line.sp),
  zq90 = extract(zq90, Line.sp))


head(used.Line.df)
summary(used.Line.df)
dim(used.Line.df)



#######################################
#Obtain some random points within mcp

random.points.Line <- spsample(Line.mcp, 55, "random")

#Calculate the covariates for the null set

avail.Line.df <- data.frame(
  used = rep(0, length(random.points.Line)),
  Species = rep("Mys", length(Line.sp)),
  BatID = rep("Line", length(Line.sp)),
  zmax = extract(zmax, random.points.Line),
  zmean = extract(zmean, random.points.Line),
  zsd = extract(zsd, random.points.Line),
  zentropy = extract(zentropy, random.points.Line),
  pzabovezmean = extract(pzabovezmean, random.points.Line),
  pzabove0.5 = extract(pzabove0.5, random.points.Line),
  zq10 = extract(zq10, random.points.Line),
  zq20 = extract(zq20, random.points.Line),
  zq30 = extract(zq30, random.points.Line),
  zq40 = extract(zq40, random.points.Line),
  zq50 = extract(zq50, random.points.Line),
  zq60 = extract(zq60, random.points.Line),
  zq70 = extract(zq70, random.points.Line),
  zq80 = extract(zq80, random.points.Line),
  zq90 = extract(zq90, random.points.Line))

head(avail.Line.df)
summary(avail.Line.df)
dim(avail.Line.df)




#------Louise

Louise.locs <- subset(bat.locs, BatID == "Louise")

#Determine available space (simplistically: MCP)
Louise.mcp <- mcp(SpatialPoints(Louise.locs), 95)
Louise.mcp.xy <- Louise.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Louise.mcp.xy <- fortify(Louise.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Louise.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Louise.mcp, col="darkred", lwd=2)
mtext("Mmys Louise zmax res250")

#Create data frame of covariates

names(Louise.locs)

Louise.sp<-SpatialPoints(Louise.locs)


used.Louise.df <- data.frame(
  used = rep(1, length(Louise.sp)),
  Species = rep("Mys", length(Louise.sp)),
  BatID = rep("Louise", length(Louise.sp)),
  zmax = extract(zmax, Louise.sp),
  zmean = extract(zmean, Louise.sp),
  zsd = extract(zsd, Louise.sp),
  zentropy = extract(zentropy, Louise.sp),
  pzabovezmean = extract(pzabovezmean, Louise.sp),
  pzabove0.5 = extract(pzabove0.5, Louise.sp),
  zq10 = extract(zq10, Louise.sp),
  zq20 = extract(zq20, Louise.sp),
  zq30 = extract(zq30, Louise.sp),
  zq40 = extract(zq40, Louise.sp),
  zq50 = extract(zq50, Louise.sp),
  zq60 = extract(zq60, Louise.sp),
  zq70 = extract(zq70, Louise.sp),
  zq80 = extract(zq80, Louise.sp),
  zq90 = extract(zq90, Louise.sp))


head(used.Louise.df)
summary(used.Louise.df)
dim(used.Louise.df)



#######################################
#Obtain some random points within mcp

random.points.Louise <- spsample(Louise.mcp, 64, "random")

#Calculate the covariates for the null set

avail.Louise.df <- data.frame(
  used = rep(0, length(random.points.Louise)),
  Species = rep("Mys", length(Louise.sp)),
  BatID = rep("Louise", length(Louise.sp)),
  zmax = extract(zmax, random.points.Louise),
  zmean = extract(zmean, random.points.Louise),
  zsd = extract(zsd, random.points.Louise),
  zentropy = extract(zentropy, random.points.Louise),
  pzabovezmean = extract(pzabovezmean, random.points.Louise),
  pzabove0.5 = extract(pzabove0.5, random.points.Louise),
  zq10 = extract(zq10, random.points.Louise),
  zq20 = extract(zq20, random.points.Louise),
  zq30 = extract(zq30, random.points.Louise),
  zq40 = extract(zq40, random.points.Louise),
  zq50 = extract(zq50, random.points.Louise),
  zq60 = extract(zq60, random.points.Louise),
  zq70 = extract(zq70, random.points.Louise),
  zq80 = extract(zq80, random.points.Louise),
  zq90 = extract(zq90, random.points.Louise))


head(avail.Louise.df)
summary(avail.Louise.df)
dim(avail.Louise.df)



#------Marie

Marie.locs <- subset(bat.locs, BatID == "Marie")

#Determine available space (simplistically: MCP)
Marie.mcp <- mcp(SpatialPoints(Marie.locs), 95)
Marie.mcp.xy <- Marie.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Marie.mcp.xy <- fortify(Marie.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Marie.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Marie.mcp, col="darkred", lwd=2)
mtext("Mmys Marie zmax res250")

#Create data frame of covariates

names(Marie.locs)

Marie.sp<-SpatialPoints(Marie.locs)


used.Marie.df <- data.frame(
  used = rep(1, length(Marie.sp)),
  Species = rep("Mys", length(Marie.sp)),
  BatID = rep("Marie", length(Marie.sp)),
  zmax = extract(zmax, Marie.sp),
  zmean = extract(zmean, Marie.sp),
  zsd = extract(zsd, Marie.sp),
  zentropy = extract(zentropy, Marie.sp),
  pzabovezmean = extract(pzabovezmean, Marie.sp),
  pzabove0.5 = extract(pzabove0.5, Marie.sp),
  zq10 = extract(zq10, Marie.sp),
  zq20 = extract(zq20, Marie.sp),
  zq30 = extract(zq30, Marie.sp),
  zq40 = extract(zq40, Marie.sp),
  zq50 = extract(zq50, Marie.sp),
  zq60 = extract(zq60, Marie.sp),
  zq70 = extract(zq70, Marie.sp),
  zq80 = extract(zq80, Marie.sp),
  zq90 = extract(zq90, Marie.sp))


head(used.Marie.df)
summary(used.Marie.df)
dim(used.Marie.df)



#######################################
#Obtain some random points within mcp

random.points.Marie <- spsample(Marie.mcp, 31, "random")

#Calculate the covariates for the null set

avail.Marie.df <- data.frame(
  used = rep(0, length(random.points.Marie)),
  Species = rep("Mys", length(Marie.sp)),
  BatID = rep("Marie", length(Marie.sp)),
  zmax = extract(zmax, random.points.Marie),
  zmean = extract(zmean, random.points.Marie),
  zsd = extract(zsd, random.points.Marie),
  zentropy = extract(zentropy, random.points.Marie),
  pzabovezmean = extract(pzabovezmean, random.points.Marie),
  pzabove0.5 = extract(pzabove0.5, random.points.Marie),
  zq10 = extract(zq10, random.points.Marie),
  zq20 = extract(zq20, random.points.Marie),
  zq30 = extract(zq30, random.points.Marie),
  zq40 = extract(zq40, random.points.Marie),
  zq50 = extract(zq50, random.points.Marie),
  zq60 = extract(zq60, random.points.Marie),
  zq70 = extract(zq70, random.points.Marie),
  zq80 = extract(zq80, random.points.Marie),
  zq90 = extract(zq90, random.points.Marie))


head(avail.Marie.df)
summary(avail.Marie.df)
dim(avail.Marie.df) 



#------Ragnhild

Ragnhild.locs <- subset(bat.locs, BatID == "Ragnhild")

#Determine available space (simplistically: MCP)
Ragnhild.mcp <- mcp(SpatialPoints(Ragnhild.locs), 95)
Ragnhild.mcp.xy <- Ragnhild.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Ragnhild.mcp.xy <- fortify(Ragnhild.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Ragnhild.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Ragnhild.mcp, col="darkred", lwd=2)
mtext("Mmys Ragnhild zmax res250")

#Create data frame of covariates

names(Ragnhild.locs)

Ragnhild.sp<-SpatialPoints(Ragnhild.locs)


used.Ragnhild.df <- data.frame(
  used = rep(1, length(Ragnhild.sp)),
  Species = rep("Mys", length(Ragnhild.sp)),
  BatID = rep("Ragnhild", length(Ragnhild.sp)),
  zmax = extract(zmax, Ragnhild.sp),
  zmean = extract(zmean, Ragnhild.sp),
  zsd = extract(zsd, Ragnhild.sp),
  zentropy = extract(zentropy, Ragnhild.sp),
  pzabovezmean = extract(pzabovezmean, Ragnhild.sp),
  pzabove0.5 = extract(pzabove0.5, Ragnhild.sp),
  zq10 = extract(zq10, Ragnhild.sp),
  zq20 = extract(zq20, Ragnhild.sp),
  zq30 = extract(zq30, Ragnhild.sp),
  zq40 = extract(zq40, Ragnhild.sp),
  zq50 = extract(zq50, Ragnhild.sp),
  zq60 = extract(zq60, Ragnhild.sp),
  zq70 = extract(zq70, Ragnhild.sp),
  zq80 = extract(zq80, Ragnhild.sp),
  zq90 = extract(zq90, Ragnhild.sp))


head(used.Ragnhild.df)
summary(used.Ragnhild.df)
dim(used.Ragnhild.df)



#######################################
#Obtain some random points within mcp

random.points.Ragnhild <- spsample(Ragnhild.mcp, 28, "random")

#Calculate the covariates for the null set

avail.Ragnhild.df <- data.frame(
  used = rep(0, length(random.points.Ragnhild)),
  Species = rep("Mys", length(Ragnhild.sp)),
  BatID = rep("Ragnhild", length(Ragnhild.sp)),
  zmax = extract(zmax, random.points.Ragnhild),
  zmean = extract(zmean, random.points.Ragnhild),
  zsd = extract(zsd, random.points.Ragnhild),
  zentropy = extract(zentropy, random.points.Ragnhild),
  pzabovezmean = extract(pzabovezmean, random.points.Ragnhild),
  pzabove0.5 = extract(pzabove0.5, random.points.Ragnhild),
  zq10 = extract(zq10, random.points.Ragnhild),
  zq20 = extract(zq20, random.points.Ragnhild),
  zq30 = extract(zq30, random.points.Ragnhild),
  zq40 = extract(zq40, random.points.Ragnhild),
  zq50 = extract(zq50, random.points.Ragnhild),
  zq60 = extract(zq60, random.points.Ragnhild),
  zq70 = extract(zq70, random.points.Ragnhild),
  zq80 = extract(zq80, random.points.Ragnhild),
  zq90 = extract(zq90, random.points.Ragnhild))


head(avail.Ragnhild.df)
summary(avail.Ragnhild.df)
dim(avail.Ragnhild.df)



#------Reeda

Reeda.locs <- subset(bat.locs, BatID == "Reeda")

#Determine available space (simplistically: MCP)
Reeda.mcp <- mcp(SpatialPoints(Reeda.locs), 95)
Reeda.mcp.xy <- Reeda.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Reeda.mcp.xy <- fortify(Reeda.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Reeda.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Reeda.mcp, col="darkred", lwd=2)
mtext("Mmys Reeda zmax res250")

#Create data frame of covariates

names(Reeda.locs)

Reeda.sp<-SpatialPoints(Reeda.locs)


used.Reeda.df <- data.frame(
  used = rep(1, length(Reeda.sp)),
  Species = rep("Mys", length(Reeda.sp)),
  BatID = rep("Reeda", length(Reeda.sp)),
  zmax = extract(zmax, Reeda.sp),
  zmean = extract(zmean, Reeda.sp),
  zsd = extract(zsd, Reeda.sp),
  zentropy = extract(zentropy, Reeda.sp),
  pzabovezmean = extract(pzabovezmean, Reeda.sp),
  pzabove0.5 = extract(pzabove0.5, Reeda.sp),
  zq10 = extract(zq10, Reeda.sp),
  zq20 = extract(zq20, Reeda.sp),
  zq30 = extract(zq30, Reeda.sp),
  zq40 = extract(zq40, Reeda.sp),
  zq50 = extract(zq50, Reeda.sp),
  zq60 = extract(zq60, Reeda.sp),
  zq70 = extract(zq70, Reeda.sp),
  zq80 = extract(zq80, Reeda.sp),
  zq90 = extract(zq90, Reeda.sp))


head(used.Reeda.df)
summary(used.Reeda.df)
dim(used.Reeda.df)



#######################################
#Obtain some random points within mcp

random.points.Reeda <- spsample(Reeda.mcp, 70, "random")

#Calculate the covariates for the null set

avail.Reeda.df <- data.frame(
  used = rep(0, length(random.points.Reeda)),
  Species = rep("Mys", length(Reeda.sp)),
  BatID = rep("Reeda", length(Reeda.sp)),
  zmax = extract(zmax, random.points.Reeda),
  zmean = extract(zmean, random.points.Reeda),
  zsd = extract(zsd, random.points.Reeda),
  zentropy = extract(zentropy, random.points.Reeda),
  pzabovezmean = extract(pzabovezmean, random.points.Reeda),
  pzabove0.5 = extract(pzabove0.5, random.points.Reeda),
  zq10 = extract(zq10, random.points.Reeda),
  zq20 = extract(zq20, random.points.Reeda),
  zq30 = extract(zq30, random.points.Reeda),
  zq40 = extract(zq40, random.points.Reeda),
  zq50 = extract(zq50, random.points.Reeda),
  zq60 = extract(zq60, random.points.Reeda),
  zq70 = extract(zq70, random.points.Reeda),
  zq80 = extract(zq80, random.points.Reeda),
  zq90 = extract(zq90, random.points.Reeda))


head(avail.Reeda.df)
summary(avail.Reeda.df)
dim(avail.Reeda.df)




#------Stine

Stine.locs <- subset(bat.locs, BatID == "Stine")

#Determine available space (simplistically: MCP)
Stine.mcp <- mcp(SpatialPoints(Stine.locs), 95)
Stine.mcp.xy <- Stine.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Stine.mcp.xy <- fortify(Stine.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Stine.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Stine.mcp, col="darkred", lwd=2)
mtext("Mmys Stine zmax res250")

#Create data frame of covariates

names(Stine.locs)

Stine.sp<-SpatialPoints(Stine.locs)


used.Stine.df <- data.frame(
  used = rep(1, length(Stine.sp)),
  Species = rep("Mys", length(Stine.sp)),
  BatID = rep("Stine", length(Stine.sp)),
  zmax = extract(zmax, Stine.sp),
  zmean = extract(zmean, Stine.sp),
  zsd = extract(zsd, Stine.sp),
  zentropy = extract(zentropy, Stine.sp),
  pzabovezmean = extract(pzabovezmean, Stine.sp),
  pzabove0.5 = extract(pzabove0.5, Stine.sp),
  zq10 = extract(zq10, Stine.sp),
  zq20 = extract(zq20, Stine.sp),
  zq30 = extract(zq30, Stine.sp),
  zq40 = extract(zq40, Stine.sp),
  zq50 = extract(zq50, Stine.sp),
  zq60 = extract(zq60, Stine.sp),
  zq70 = extract(zq70, Stine.sp),
  zq80 = extract(zq80, Stine.sp),
  zq90 = extract(zq90, Stine.sp))


head(used.Stine.df)
summary(used.Stine.df)
dim(used.Stine.df)



#######################################
#Obtain some random points within mcp

random.points.Stine <- spsample(Stine.mcp, 84, "random")

#Calculate the covariates for the null set

avail.Stine.df <- data.frame(
  used = rep(0, length(random.points.Stine)),
  Species = rep("Mys", length(Stine.sp)),
  BatID = rep("Stine", length(Stine.sp)),
  zmax = extract(zmax, random.points.Stine),
  zmean = extract(zmean, random.points.Stine),
  zsd = extract(zsd, random.points.Stine),
  zentropy = extract(zentropy, random.points.Stine),
  pzabovezmean = extract(pzabovezmean, random.points.Stine),
  pzabove0.5 = extract(pzabove0.5, random.points.Stine),
  zq10 = extract(zq10, random.points.Stine),
  zq20 = extract(zq20, random.points.Stine),
  zq30 = extract(zq30, random.points.Stine),
  zq40 = extract(zq40, random.points.Stine),
  zq50 = extract(zq50, random.points.Stine),
  zq60 = extract(zq60, random.points.Stine),
  zq70 = extract(zq70, random.points.Stine),
  zq80 = extract(zq80, random.points.Stine),
  zq90 = extract(zq90, random.points.Stine))


head(avail.Stine.df)
summary(avail.Stine.df)
dim(avail.Stine.df)



#------Turid

Turid.locs <- subset(bat.locs, BatID == "Turid")

#Determine available space (simplistically: MCP)
Turid.mcp <- mcp(SpatialPoints(Turid.locs), 95)
Turid.mcp.xy <- Turid.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Turid.mcp.xy <- fortify(Turid.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Turid.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Turid.mcp, col="darkred", lwd=2)
mtext("Mmys Turid zmax res250")

#Create data frame of covariates

names(Turid.locs)

Turid.sp<-SpatialPoints(Turid.locs)


used.Turid.df <- data.frame(
  used = rep(1, length(Turid.sp)),
  Species = rep("Mys", length(Turid.sp)),
  BatID = rep("Turid", length(Turid.sp)),
  zmax = extract(zmax, Turid.sp),
  zmean = extract(zmean, Turid.sp),
  zsd = extract(zsd, Turid.sp),
  zentropy = extract(zentropy, Turid.sp),
  pzabovezmean = extract(pzabovezmean, Turid.sp),
  pzabove0.5 = extract(pzabove0.5, Turid.sp),
  zq10 = extract(zq10, Turid.sp),
  zq20 = extract(zq20, Turid.sp),
  zq30 = extract(zq30, Turid.sp),
  zq40 = extract(zq40, Turid.sp),
  zq50 = extract(zq50, Turid.sp),
  zq60 = extract(zq60, Turid.sp),
  zq70 = extract(zq70, Turid.sp),
  zq80 = extract(zq80, Turid.sp),
  zq90 = extract(zq90, Turid.sp))


head(used.Turid.df)
summary(used.Turid.df)
dim(used.Turid.df)




#######################################
#Obtain some random points within mcp

random.points.Turid <- spsample(Turid.mcp, 8, "random")

#Calculate the covariates for the null set

avail.Turid.df <- data.frame(
  used = rep(0, length(random.points.Turid)),
  Species = rep("Mys", length(Turid.sp)),
  BatID = rep("Turid", length(Turid.sp)),
  zmax = extract(zmax, random.points.Turid),
  zmean = extract(zmean, random.points.Turid),
  zsd = extract(zsd, random.points.Turid),
  zentropy = extract(zentropy, random.points.Turid),
  pzabovezmean = extract(pzabovezmean, random.points.Turid),
  pzabove0.5 = extract(pzabove0.5, random.points.Turid),
  zq10 = extract(zq10, random.points.Turid),
  zq20 = extract(zq20, random.points.Turid),
  zq30 = extract(zq30, random.points.Turid),
  zq40 = extract(zq40, random.points.Turid),
  zq50 = extract(zq50, random.points.Turid),
  zq60 = extract(zq60, random.points.Turid),
  zq70 = extract(zq70, random.points.Turid),
  zq80 = extract(zq80, random.points.Turid),
  zq90 = extract(zq90, random.points.Turid))


head(avail.Turid.df)
summary(avail.Turid.df)
dim(avail.Turid.df)


#Combine Used data for all Mmys individuals
Mmys.data.used <- rbind(
  used.Amelia.df, 
  used.Aricia.df, 
  used.Daisy.df, 
  used.Ethel.df,
  used.Ida.df, 
  used.Line.df, 
  used.Louise.df, 
  used.Marie.df, 
  used.Ragnhild.df, 
  used.Reeda.df, 
  used.Stine.df, 
  used.Turid.df)
head(Mmys.data.used)
summary(Mmys.data.used)
dim(Mmys.data.used)
table(Mmys.data.used$BatID)

#Combine Avilable data for all Mmys individuals
Mmys.data.avail <- rbind(
  avail.Amelia.df, 
  avail.Aricia.df, 
  avail.Daisy.df, 
  avail.Ethel.df,
  avail.Ida.df, 
  avail.Line.df, 
  avail.Louise.df, 
  avail.Marie.df, 
  avail.Ragnhild.df, 
  avail.Reeda.df, 
  avail.Stine.df, 
  avail.Turid.df)
head(Mmys.data.avail)
summary(Mmys.data.avail)
dim(Mmys.data.avail)
table(Mmys.data.avail$BatID)

#Combine Used and Avilable data for all Mmys individuals
Mmys.data <- rbind(Mmys.data.used, Mmys.data.avail)
head(Mmys.data)
summary(Mmys.data)
dim(Mmys.data)
table(Mmys.data$BatID)

###############################################################
###-------------  M. bandtii bats  ------------------------####
###############################################################

#------Astrid

Astrid.locs <- subset(bat.locs, BatID == "Astrid")

#Determine available space (simplistically: MCP)
Astrid.mcp <- mcp(SpatialPoints(Astrid.locs), 95)
Astrid.mcp.xy <- Astrid.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Astrid.mcp.xy <- fortify(Astrid.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Astrid.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Astrid.mcp, col="blue", lwd=2)
mtext("Mbra Astrid zmax res250")

#Create data frame of covariates

names(Astrid.locs)

Astrid.sp<-SpatialPoints(Astrid.locs)


used.Astrid.df <- data.frame(
  used = rep(1, length(Astrid.sp)),
  Species = rep("Mbra", length(Astrid.sp)),
  BatID = rep("Astrid", length(Astrid.sp)),
  zmax = extract(zmax, Astrid.sp),
  zmean = extract(zmean, Astrid.sp),
  zsd = extract(zsd, Astrid.sp),
  zentropy = extract(zentropy, Astrid.sp),
  pzabovezmean = extract(pzabovezmean, Astrid.sp),
  pzabove0.5 = extract(pzabove0.5, Astrid.sp),
  zq10 = extract(zq10, Astrid.sp),
  zq20 = extract(zq20, Astrid.sp),
  zq30 = extract(zq30, Astrid.sp),
  zq40 = extract(zq40, Astrid.sp),
  zq50 = extract(zq50, Astrid.sp),
  zq60 = extract(zq60, Astrid.sp),
  zq70 = extract(zq70, Astrid.sp),
  zq80 = extract(zq80, Astrid.sp),
  zq90 = extract(zq90, Astrid.sp))


head(used.Astrid.df)
summary(used.Astrid.df)
dim(used.Astrid.df) 



#######################################
#Obtain some random points within mcp

random.points.Astrid <- spsample(Astrid.mcp, 81, "random")

#Calculate the covariates for the null set

avail.Astrid.df <- data.frame(
  used = rep(0, length(random.points.Astrid)),
  Species = rep("Mbra", length(Astrid.sp)),
  BatID = rep("Astrid", length(Astrid.sp)),
  zmax = extract(zmax, random.points.Astrid),
  zmean = extract(zmean, random.points.Astrid),
  zsd = extract(zsd, random.points.Astrid),
  zentropy = extract(zentropy, random.points.Astrid),
  pzabovezmean = extract(pzabovezmean, random.points.Astrid),
  pzabove0.5 = extract(pzabove0.5, random.points.Astrid),
  zq10 = extract(zq10, random.points.Astrid),
  zq20 = extract(zq20, random.points.Astrid),
  zq30 = extract(zq30, random.points.Astrid),
  zq40 = extract(zq40, random.points.Astrid),
  zq50 = extract(zq50, random.points.Astrid),
  zq60 = extract(zq60, random.points.Astrid),
  zq70 = extract(zq70, random.points.Astrid),
  zq80 = extract(zq80, random.points.Astrid),
  zq90 = extract(zq90, random.points.Astrid))


head(avail.Astrid.df)
summary(avail.Astrid.df)
dim(avail.Astrid.df)


#------Dagny

Dagny.locs <- subset(bat.locs, BatID == "Dagny")

#Determine available space (simplistically: MCP)
Dagny.mcp <- mcp(SpatialPoints(Dagny.locs), 95)
Dagny.mcp.xy <- Dagny.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Dagny.mcp.xy <- fortify(Dagny.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Dagny.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Dagny.mcp, col="blue", lwd=2)
mtext("Mbra Dagny zmax res250")

#Create data frame of covariates

names(Dagny.locs)

Dagny.sp<-SpatialPoints(Dagny.locs)


used.Dagny.df <- data.frame(
  used = rep(1, length(Dagny.sp)),
  Species = rep("Mbra", length(Dagny.sp)),
  BatID = rep("Dagny", length(Dagny.sp)),
  zmax = extract(zmax, Dagny.sp),
  zmean = extract(zmean, Dagny.sp),
  zsd = extract(zsd, Dagny.sp),
  zentropy = extract(zentropy, Dagny.sp),
  pzabovezmean = extract(pzabovezmean, Dagny.sp),
  pzabove0.5 = extract(pzabove0.5, Dagny.sp),
  zq10 = extract(zq10, Dagny.sp),
  zq20 = extract(zq20, Dagny.sp),
  zq30 = extract(zq30, Dagny.sp),
  zq40 = extract(zq40, Dagny.sp),
  zq50 = extract(zq50, Dagny.sp),
  zq60 = extract(zq60, Dagny.sp),
  zq70 = extract(zq70, Dagny.sp),
  zq80 = extract(zq80, Dagny.sp),
  zq90 = extract(zq90, Dagny.sp))


head(used.Dagny.df)
summary(used.Dagny.df)
dim(used.Dagny.df) 



#######################################
#Obtain some random points within mcp

random.points.Dagny <- spsample(Dagny.mcp, 83, "random")

#Calculate the covariates for the null set

avail.Dagny.df <- data.frame(
  used = rep(0, length(random.points.Dagny)),
  Species = rep("Mbra", length(Dagny.sp)),
  BatID = rep("Dagny", length(Dagny.sp)),
  zmax = extract(zmax, random.points.Dagny),
  zmean = extract(zmean, random.points.Dagny),
  zsd = extract(zsd, random.points.Dagny),
  zentropy = extract(zentropy, random.points.Dagny),
  pzabovezmean = extract(pzabovezmean, random.points.Dagny),
  pzabove0.5 = extract(pzabove0.5, random.points.Dagny),
  zq10 = extract(zq10, random.points.Dagny),
  zq20 = extract(zq20, random.points.Dagny),
  zq30 = extract(zq30, random.points.Dagny),
  zq40 = extract(zq40, random.points.Dagny),
  zq50 = extract(zq50, random.points.Dagny),
  zq60 = extract(zq60, random.points.Dagny),
  zq70 = extract(zq70, random.points.Dagny),
  zq80 = extract(zq80, random.points.Dagny),
  zq90 = extract(zq90, random.points.Dagny))


head(avail.Dagny.df)
summary(avail.Dagny.df)
dim(avail.Dagny.df)


#------Kaja

Kaja.locs <- subset(bat.locs, BatID == "Kaja")

#Determine available space (simplistically: MCP)
Kaja.mcp <- mcp(SpatialPoints(Kaja.locs), 95)
Kaja.mcp.xy <- Kaja.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Kaja.mcp.xy <- fortify(Kaja.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Kaja.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Kaja.mcp, col="blue", lwd=2)
mtext("Mbra Kaja zmax res250")

#Create data frame of covariates

names(Kaja.locs)

Kaja.sp<-SpatialPoints(Kaja.locs)


used.Kaja.df <- data.frame(
  used = rep(1, length(Kaja.sp)),
  Species = rep("Mbra", length(Kaja.sp)),
  BatID = rep("Kaja", length(Kaja.sp)),
  zmax = extract(zmax, Kaja.sp),
  zmean = extract(zmean, Kaja.sp),
  zsd = extract(zsd, Kaja.sp),
  zentropy = extract(zentropy, Kaja.sp),
  pzabovezmean = extract(pzabovezmean, Kaja.sp),
  pzabove0.5 = extract(pzabove0.5, Kaja.sp),
  zq10 = extract(zq10, Kaja.sp),
  zq20 = extract(zq20, Kaja.sp),
  zq30 = extract(zq30, Kaja.sp),
  zq40 = extract(zq40, Kaja.sp),
  zq50 = extract(zq50, Kaja.sp),
  zq60 = extract(zq60, Kaja.sp),
  zq70 = extract(zq70, Kaja.sp),
  zq80 = extract(zq80, Kaja.sp),
  zq90 = extract(zq90, Kaja.sp))


head(used.Kaja.df)
summary(used.Kaja.df)
dim(used.Kaja.df) 



#######################################
#Obtain some random points within mcp

random.points.Kaja <- spsample(Kaja.mcp, 91, "random")

#Calculate the covariates for the null set

avail.Kaja.df <- data.frame(
  used = rep(0, length(random.points.Kaja)),
  Species = rep("Mbra", length(Kaja.sp)),
  BatID = rep("Kaja", length(Kaja.sp)),
  zmax = extract(zmax, random.points.Kaja),
  zmean = extract(zmean, random.points.Kaja),
  zsd = extract(zsd, random.points.Kaja),
  zentropy = extract(zentropy, random.points.Kaja),
  pzabovezmean = extract(pzabovezmean, random.points.Kaja),
  pzabove0.5 = extract(pzabove0.5, random.points.Kaja),
  zq10 = extract(zq10, random.points.Kaja),
  zq20 = extract(zq20, random.points.Kaja),
  zq30 = extract(zq30, random.points.Kaja),
  zq40 = extract(zq40, random.points.Kaja),
  zq50 = extract(zq50, random.points.Kaja),
  zq60 = extract(zq60, random.points.Kaja),
  zq70 = extract(zq70, random.points.Kaja),
  zq80 = extract(zq80, random.points.Kaja),
  zq90 = extract(zq90, random.points.Kaja))

head(avail.Kaja.df)
summary(avail.Kaja.df)
dim(avail.Kaja.df)


#------Maren

Maren.locs <- subset(bat.locs, BatID == "Maren")

#Determine available space (simplistically: MCP)
Maren.mcp <- mcp(SpatialPoints(Maren.locs), 95)
Maren.mcp.xy <- Maren.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Maren.mcp.xy <- fortify(Maren.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Maren.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Maren.mcp, col="blue", lwd=2)
mtext("Mbra Maren zmax res250")

#Create data frame of covariates

names(Maren.locs)

Maren.sp<-SpatialPoints(Maren.locs)


used.Maren.df <- data.frame(
  used = rep(1, length(Maren.sp)),
  Species = rep("Mbra", length(Maren.sp)),
  BatID = rep("Maren", length(Maren.sp)),
  zmax = extract(zmax, Maren.sp),
  zmean = extract(zmean, Maren.sp),
  zsd = extract(zsd, Maren.sp),
  zentropy = extract(zentropy, Maren.sp),
  pzabovezmean = extract(pzabovezmean, Maren.sp),
  pzabove0.5 = extract(pzabove0.5, Maren.sp),
  zq10 = extract(zq10, Maren.sp),
  zq20 = extract(zq20, Maren.sp),
  zq30 = extract(zq30, Maren.sp),
  zq40 = extract(zq40, Maren.sp),
  zq50 = extract(zq50, Maren.sp),
  zq60 = extract(zq60, Maren.sp),
  zq70 = extract(zq70, Maren.sp),
  zq80 = extract(zq80, Maren.sp),
  zq90 = extract(zq90, Maren.sp))

head(used.Maren.df)
summary(used.Maren.df)
dim(used.Maren.df) 



#######################################
#Obtain some random points within mcp

random.points.Maren <- spsample(Maren.mcp, 141, "random")

#Calculate the covariates for the null set

avail.Maren.df <- data.frame(
  used = rep(0, length(random.points.Maren)),
  Species = rep("Mbra", length(Maren.sp)),
  BatID = rep("Maren", length(Maren.sp)),
  zmax = extract(zmax, random.points.Maren),
  zmean = extract(zmean, random.points.Maren),
  zsd = extract(zsd, random.points.Maren),
  zentropy = extract(zentropy, random.points.Maren),
  pzabovezmean = extract(pzabovezmean, random.points.Maren),
  pzabove0.5 = extract(pzabove0.5, random.points.Maren),
  zq10 = extract(zq10, random.points.Maren),
  zq20 = extract(zq20, random.points.Maren),
  zq30 = extract(zq30, random.points.Maren),
  zq40 = extract(zq40, random.points.Maren),
  zq50 = extract(zq50, random.points.Maren),
  zq60 = extract(zq60, random.points.Maren),
  zq70 = extract(zq70, random.points.Maren),
  zq80 = extract(zq80, random.points.Maren),
  zq90 = extract(zq90, random.points.Maren))


head(avail.Maren.df)
summary(avail.Maren.df)
dim(avail.Maren.df)


#------Nora

Nora.locs <- subset(bat.locs, BatID == "Nora")

#Determine available space (simplistically: MCP)
Nora.mcp <- mcp(SpatialPoints(Nora.locs), 95)
Nora.mcp.xy <- Nora.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Nora.mcp.xy <- fortify(Nora.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Nora.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Nora.mcp, col="blue", lwd=2)
mtext("Mbra Nora zmax res250")

#Create data frame of covariates

names(Nora.locs)

Nora.sp<-SpatialPoints(Nora.locs)


used.Nora.df <- data.frame(
  used = rep(1, length(Nora.sp)),
  Species = rep("Mbra", length(Nora.sp)),
  BatID = rep("Nora", length(Nora.sp)),
  zmax = extract(zmax, Nora.sp),
  zmean = extract(zmean, Nora.sp),
  zsd = extract(zsd, Nora.sp),
  zentropy = extract(zentropy, Nora.sp),
  pzabovezmean = extract(pzabovezmean, Nora.sp),
  pzabove0.5 = extract(pzabove0.5, Nora.sp),
  zq10 = extract(zq10, Nora.sp),
  zq20 = extract(zq20, Nora.sp),
  zq30 = extract(zq30, Nora.sp),
  zq40 = extract(zq40, Nora.sp),
  zq50 = extract(zq50, Nora.sp),
  zq60 = extract(zq60, Nora.sp),
  zq70 = extract(zq70, Nora.sp),
  zq80 = extract(zq80, Nora.sp),
  zq90 = extract(zq90, Nora.sp))

head(used.Nora.df)
summary(used.Nora.df)
dim(used.Nora.df) 



#######################################
#Obtain some random points within mcp

random.points.Nora <- spsample(Nora.mcp, 49, "random")

#Calculate the covariates for the null set

avail.Nora.df <- data.frame(
  used = rep(0, length(random.points.Nora)),
  Species = rep("Mbra", length(Nora.sp)),
  BatID = rep("Nora", length(Nora.sp)),
  zmax = extract(zmax, random.points.Nora),
  zmean = extract(zmean, random.points.Nora),
  zsd = extract(zsd, random.points.Nora),
  zentropy = extract(zentropy, random.points.Nora),
  pzabovezmean = extract(pzabovezmean, random.points.Nora),
  pzabove0.5 = extract(pzabove0.5, random.points.Nora),
  zq10 = extract(zq10, random.points.Nora),
  zq20 = extract(zq20, random.points.Nora),
  zq30 = extract(zq30, random.points.Nora),
  zq40 = extract(zq40, random.points.Nora),
  zq50 = extract(zq50, random.points.Nora),
  zq60 = extract(zq60, random.points.Nora),
  zq70 = extract(zq70, random.points.Nora),
  zq80 = extract(zq80, random.points.Nora),
  zq90 = extract(zq90, random.points.Nora))

head(avail.Nora.df)
summary(avail.Nora.df)
dim(avail.Nora.df)


#------Phoebe

Phoebe.locs <- subset(bat.locs, BatID == "Phoebe")

#Determine available space (simplistically: MCP)
Phoebe.mcp <- mcp(SpatialPoints(Phoebe.locs), 95)
Phoebe.mcp.xy <- Phoebe.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Phoebe.mcp.xy <- fortify(Phoebe.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Phoebe.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Phoebe.mcp, col="blue", lwd=2)
mtext("Mbra Phoebe zmax res250")

#Create data frame of covariates

names(Phoebe.locs)

Phoebe.sp<-SpatialPoints(Phoebe.locs)


used.Phoebe.df <- data.frame(
  used = rep(1, length(Phoebe.sp)),
  Species = rep("Mbra", length(Phoebe.sp)),
  BatID = rep("Phoebe", length(Phoebe.sp)),
  zmax = extract(zmax, Phoebe.sp),
  zmean = extract(zmean, Phoebe.sp),
  zsd = extract(zsd, Phoebe.sp),
  zentropy = extract(zentropy, Phoebe.sp),
  pzabovezmean = extract(pzabovezmean, Phoebe.sp),
  pzabove0.5 = extract(pzabove0.5, Phoebe.sp),
  zq10 = extract(zq10, Phoebe.sp),
  zq20 = extract(zq20, Phoebe.sp),
  zq30 = extract(zq30, Phoebe.sp),
  zq40 = extract(zq40, Phoebe.sp),
  zq50 = extract(zq50, Phoebe.sp),
  zq60 = extract(zq60, Phoebe.sp),
  zq70 = extract(zq70, Phoebe.sp),
  zq80 = extract(zq80, Phoebe.sp),
  zq90 = extract(zq90, Phoebe.sp))


head(used.Phoebe.df)
summary(used.Phoebe.df)
dim(used.Phoebe.df) 



#######################################
#Obtain some random points within mcp

random.points.Phoebe <- spsample(Phoebe.mcp, 62, "random")

#Calculate the covariates for the null set

avail.Phoebe.df <- data.frame(
  used = rep(0, length(random.points.Phoebe)),
  Species = rep("Mbra", length(Phoebe.sp)),
  BatID = rep("Phoebe", length(Phoebe.sp)),
  zmax = extract(zmax, random.points.Phoebe),
  zmean = extract(zmean, random.points.Phoebe),
  zsd = extract(zsd, random.points.Phoebe),
  zentropy = extract(zentropy, random.points.Phoebe),
  pzabovezmean = extract(pzabovezmean, random.points.Phoebe),
  pzabove0.5 = extract(pzabove0.5, random.points.Phoebe),
  zq10 = extract(zq10, random.points.Phoebe),
  zq20 = extract(zq20, random.points.Phoebe),
  zq30 = extract(zq30, random.points.Phoebe),
  zq40 = extract(zq40, random.points.Phoebe),
  zq50 = extract(zq50, random.points.Phoebe),
  zq60 = extract(zq60, random.points.Phoebe),
  zq70 = extract(zq70, random.points.Phoebe),
  zq80 = extract(zq80, random.points.Phoebe),
  zq90 = extract(zq90, random.points.Phoebe))

head(avail.Phoebe.df)
summary(avail.Phoebe.df)
dim(avail.Phoebe.df)


#------Sofia

Sofia.locs <- subset(bat.locs, BatID == "Sofia")

#Determine available space (simplistically: MCP)
Sofia.mcp <- mcp(SpatialPoints(Sofia.locs), 95)
Sofia.mcp.xy <- Sofia.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Sofia.mcp.xy <- fortify(Sofia.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Sofia.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Sofia.mcp, col="blue", lwd=2)
mtext("Mbra Sofia zmax res250")

#Create data frame of covariates

names(Sofia.locs)

Sofia.sp<-SpatialPoints(Sofia.locs)


used.Sofia.df <- data.frame(
  used = rep(1, length(Sofia.sp)),
  Species = rep("Mbra", length(Sofia.sp)),
  BatID = rep("Sofia", length(Sofia.sp)),
  zmax = extract(zmax, Sofia.sp),
  zmean = extract(zmean, Sofia.sp),
  zsd = extract(zsd, Sofia.sp),
  zentropy = extract(zentropy, Sofia.sp),
  pzabovezmean = extract(pzabovezmean, Sofia.sp),
  pzabove0.5 = extract(pzabove0.5, Sofia.sp),
  zq10 = extract(zq10, Sofia.sp),
  zq20 = extract(zq20, Sofia.sp),
  zq30 = extract(zq30, Sofia.sp),
  zq40 = extract(zq40, Sofia.sp),
  zq50 = extract(zq50, Sofia.sp),
  zq60 = extract(zq60, Sofia.sp),
  zq70 = extract(zq70, Sofia.sp),
  zq80 = extract(zq80, Sofia.sp),
  zq90 = extract(zq90, Sofia.sp))


head(used.Sofia.df)
summary(used.Sofia.df)
dim(used.Sofia.df) 



#######################################
#Obtain some random points within mcp

random.points.Sofia <- spsample(Sofia.mcp, 39, "random")

#Calculate the covariates for the null set

avail.Sofia.df <- data.frame(
  used = rep(0, length(random.points.Sofia)),
  Species = rep("Mbra", length(Sofia.sp)),
  BatID = rep("Sofia", length(Sofia.sp)),
  zmax = extract(zmax, random.points.Sofia),
  zmean = extract(zmean, random.points.Sofia),
  zsd = extract(zsd, random.points.Sofia),
  zentropy = extract(zentropy, random.points.Sofia),
  pzabovezmean = extract(pzabovezmean, random.points.Sofia),
  pzabove0.5 = extract(pzabove0.5, random.points.Sofia),
  zq10 = extract(zq10, random.points.Sofia),
  zq20 = extract(zq20, random.points.Sofia),
  zq30 = extract(zq30, random.points.Sofia),
  zq40 = extract(zq40, random.points.Sofia),
  zq50 = extract(zq50, random.points.Sofia),
  zq60 = extract(zq60, random.points.Sofia),
  zq70 = extract(zq70, random.points.Sofia),
  zq80 = extract(zq80, random.points.Sofia),
  zq90 = extract(zq90, random.points.Sofia))



head(avail.Sofia.df)
summary(avail.Sofia.df)
dim(avail.Sofia.df)


#------Steffi

Steffi.locs <- subset(bat.locs, BatID == "Steffi")

#Determine available space (simplistically: MCP)
Steffi.mcp <- mcp(SpatialPoints(Steffi.locs), 95)
Steffi.mcp.xy <- Steffi.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Steffi.mcp.xy <- fortify(Steffi.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Steffi.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Steffi.mcp, col="blue", lwd=2)
mtext("Mbra Steffi zmax res250")

#Create data frame of covariates

names(Steffi.locs)

Steffi.sp<-SpatialPoints(Steffi.locs)


used.Steffi.df <- data.frame(
  used = rep(1, length(Steffi.sp)),
  Species = rep("Mbra", length(Steffi.sp)),
  BatID = rep("Steffi", length(Steffi.sp)),
  zmax = extract(zmax, Steffi.sp),
  zmean = extract(zmean, Steffi.sp),
  zsd = extract(zsd, Steffi.sp),
  zentropy = extract(zentropy, Steffi.sp),
  pzabovezmean = extract(pzabovezmean, Steffi.sp),
  pzabove0.5 = extract(pzabove0.5, Steffi.sp),
  zq10 = extract(zq10, Steffi.sp),
  zq20 = extract(zq20, Steffi.sp),
  zq30 = extract(zq30, Steffi.sp),
  zq40 = extract(zq40, Steffi.sp),
  zq50 = extract(zq50, Steffi.sp),
  zq60 = extract(zq60, Steffi.sp),
  zq70 = extract(zq70, Steffi.sp),
  zq80 = extract(zq80, Steffi.sp),
  zq90 = extract(zq90, Steffi.sp))

head(used.Steffi.df)
summary(used.Steffi.df)
dim(used.Steffi.df) 



#######################################
#Obtain some random points within mcp

random.points.Steffi <- spsample(Steffi.mcp, 39, "random")

#Calculate the covariates for the null set

avail.Steffi.df <- data.frame(
  used = rep(0, length(random.points.Steffi)),
  Species = rep("Mbra", length(Steffi.sp)),
  BatID = rep("Steffi", length(Steffi.sp)),
  zmax = extract(zmax, random.points.Steffi),
  zmean = extract(zmean, random.points.Steffi),
  zsd = extract(zsd, random.points.Steffi),
  zentropy = extract(zentropy, random.points.Steffi),
  pzabovezmean = extract(pzabovezmean, random.points.Steffi),
  pzabove0.5 = extract(pzabove0.5, random.points.Steffi),
  zq10 = extract(zq10, random.points.Steffi),
  zq20 = extract(zq20, random.points.Steffi),
  zq30 = extract(zq30, random.points.Steffi),
  zq40 = extract(zq40, random.points.Steffi),
  zq50 = extract(zq50, random.points.Steffi),
  zq60 = extract(zq60, random.points.Steffi),
  zq70 = extract(zq70, random.points.Steffi),
  zq80 = extract(zq80, random.points.Steffi),
  zq90 = extract(zq90, random.points.Steffi))

head(avail.Steffi.df)
summary(avail.Steffi.df)
dim(avail.Steffi.df)


#------Thea

Thea.locs <- subset(bat.locs, BatID == "Thea")

#Determine available space (simplistically: MCP)
Thea.mcp <- mcp(SpatialPoints(Thea.locs), 95)
Thea.mcp.xy <- Thea.mcp@polygons[[1]]@Polygons[[1]]@coords

# FORTIYFY
require(ggplot2)
Thea.mcp.xy <- fortify(Thea.mcp)

par(mfrow=c(1,1))
plot(zmax)
points(Thea.locs, asp=1, pch=19, col=rgb(0,0,.4,.5), cex=0.5)
lines(Thea.mcp, col="blue", lwd=2)
mtext("Mbra Thea zmax res250")

#Create data frame of covariates

names(Thea.locs)

Thea.sp<-SpatialPoints(Thea.locs)


used.Thea.df <- data.frame(
  used = rep(1, length(Thea.sp)),
  Species = rep("Mbra", length(Thea.sp)),
  BatID = rep("Thea", length(Thea.sp)),
  zmax = extract(zmax, Thea.sp),
  zmean = extract(zmean, Thea.sp),
  zsd = extract(zsd, Thea.sp),
  zentropy = extract(zentropy, Thea.sp),
  pzabovezmean = extract(pzabovezmean, Thea.sp),
  pzabove0.5 = extract(pzabove0.5, Thea.sp),
  zq10 = extract(zq10, Thea.sp),
  zq20 = extract(zq20, Thea.sp),
  zq30 = extract(zq30, Thea.sp),
  zq40 = extract(zq40, Thea.sp),
  zq50 = extract(zq50, Thea.sp),
  zq60 = extract(zq60, Thea.sp),
  zq70 = extract(zq70, Thea.sp),
  zq80 = extract(zq80, Thea.sp),
  zq90 = extract(zq90, Thea.sp))


head(used.Thea.df)
summary(used.Thea.df)
dim(used.Thea.df) 



#######################################
#Obtain some random points within mcp

random.points.Thea <- spsample(Thea.mcp, 101, "random")

#Calculate the covariates for the null set

avail.Thea.df <- data.frame(
  used = rep(0, length(random.points.Thea)),
  Species = rep("Mbra", length(Thea.sp)),
  BatID = rep("Thea", length(Thea.sp)),
  zmax = extract(zmax, random.points.Thea),
  zmean = extract(zmean, random.points.Thea),
  zsd = extract(zsd, random.points.Thea),
  zentropy = extract(zentropy, random.points.Thea),
  pzabovezmean = extract(pzabovezmean, random.points.Thea),
  pzabove0.5 = extract(pzabove0.5, random.points.Thea),
  zq10 = extract(zq10, random.points.Thea),
  zq20 = extract(zq20, random.points.Thea),
  zq30 = extract(zq30, random.points.Thea),
  zq40 = extract(zq40, random.points.Thea),
  zq50 = extract(zq50, random.points.Thea),
  zq60 = extract(zq60, random.points.Thea),
  zq70 = extract(zq70, random.points.Thea),
  zq80 = extract(zq80, random.points.Thea),
  zq90 = extract(zq90, random.points.Thea))



head(avail.Thea.df)
summary(avail.Thea.df)
dim(avail.Thea.df)

###################
###################

#--------------zmax 250 plots pooled MCPs for each species
#windows()
#plot(zmax)
#lines(Amelia.mcp, col="darkred", lwd=2)
#lines(Aricia.mcp, col="darkred", lwd=2)
#lines(Daisy.mcp, col="darkred", lwd=2)
#lines(Ethel.mcp, col="darkred", lwd=2)
#lines(Ida.mcp, col="darkred", lwd=2)
#lines(Line.mcp, col="darkred", lwd=2)
#lines(Louise.mcp, col="darkred", lwd=2)
#lines(Marie.mcp, col="darkred", lwd=2)
#lines(Ragnhild.mcp, col="darkred", lwd=2)
#lines(Reeda.mcp, col="darkred", lwd=2)
#lines(Stine.mcp, col="darkred", lwd=2)
#lines(Turid.mcp, col="darkred", lwd=2)
#mtext("95% MCP M. mystacinus zmax250")

#windows()
#plot(zmax)
#lines(Astrid.mcp, col="blue", lwd=2)
#lines(Dagny.mcp, col="blue", lwd=2)
#lines(Kaja.mcp, col="blue", lwd=2)
#lines(Maren.mcp, col="blue", lwd=2)
#lines(Nora.mcp, col="blue", lwd=2)
#lines(Phoebe.mcp, col="blue", lwd=2)
#lines(Sofia.mcp, col="blue", lwd=2)
#lines(Steffi.mcp, col="blue", lwd=2)
#lines(Thea.mcp, col="blue", lwd=2)
#mtext("95% MCP M. brandtii zmax250")


#################################
#################################

#Combine used data for all Mbra individuals
Mbra.data.used <- rbind(
  used.Astrid.df,
  used.Dagny.df,
  used.Kaja.df,
  used.Maren.df,
  used.Nora.df,
  used.Phoebe.df,
  used.Sofia.df,
  used.Steffi.df,
  used.Thea.df)
head(Mbra.data.used)
summary(Mbra.data.used)
dim(Mbra.data.used)
table(Mbra.data.used$BatID)

#Combine available data for all Mbra individuals
Mbra.data.avail <- rbind(
  avail.Astrid.df, 
  avail.Dagny.df, 
  avail.Kaja.df, 
  avail.Maren.df, 
  avail.Nora.df, 
  avail.Phoebe.df, 
  avail.Sofia.df, 
  avail.Steffi.df, 
  avail.Thea.df)
head(Mbra.data.avail)
summary(Mbra.data.avail)
dim(Mbra.data.avail)
table(Mbra.data.avail$BatID)


#Combine Used and Avilable data for all Mbra individuals
Mbra.data <- rbind(Mbra.data.used, Mbra.data.avail)
head(Mbra.data)
summary(Mbra.data)
dim(Mbra.data)
table(Mbra.data$BatID)

#####################################################
#####################################################


#Combine Used and Avilable data for all Mmys and Mbra individuals
Bats.data <- rbind(Mmys.data, Mbra.data)
head(Bats.data)
summary(Bats.data) #There are a few NA's in the lidar data - not sure why and how it affects results - look into this later
dim(Bats.data)
#[1] 3004   18

Bats.data.used <- subset(Bats.data, used == "1")
Bats.data.avail <- subset(Bats.data, used == "0")

dim(Bats.data.used)
#[1] 1502   18
dim(Bats.data.avail)
#[1] 1502   18


#I am going to export this as a csv file to look at the NAs
#write.table(Bats.data, file = "Bats_data_incl_NA_250.csv", sep = ",", quote = FALSE, append = FALSE, na ="NA")

subBats.data = subset(Bats.data, select = -c(zentropy) ) #657 observations had NA for zentropy
dim(subBats.data)
#[1] 3004   17

data_250 <- na.omit(subBats.data)
summary(data_250)
dim(data_250)
#[1] 2998   17    #5 used observastions (Mmys Ida (1), Mbra Astrid(3), Mbra Stef(3) removed because NA for all lidar variables)

data_250.used <- subset(data_250, used == "1")
data_250.avail <- subset(data_250, used == "0")

dim(data_250.used)
#[1] 1496   17
dim(data_250.avail)
#[1] 1502   17

#I am going to export this as a csv file so that you can run RSFs later without going through all the previous code
#write.table(data_250, file = "data_95MCP250_excl_zentropy.csv", sep = ",", quote = FALSE, append = FALSE, na ="NA")


#Moved the RSF analyses to separate script

