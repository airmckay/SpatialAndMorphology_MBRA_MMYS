rm(list = ls())
# Clear work space

## Load, inspect and format observations ## *****************************************************************

setwd("//nmbu.no/my/home/Desktop/DATA")
#load data

Bats45<-read.csv("45%.csv",sep=";", dec=",", header=TRUE)
## THIS IS A SEMICOLON / COMMA SEPERATED CSV - DO NOT USE PERIODS TO DEFINE SEP OR COMMAS TO DEFINE DEC!!! 

## Add column for species ID ##
mystacinus<- Bats45$BatID %in% c("Amelia","Marie", "Turid", "Reeda", "Daisy", "Line", "Stine", "Ethel", "Ragnhild", "Ida", "Louise", "Aricia")

Bats45$Species<-ifelse(mystacinus == TRUE, "MYSTACINUS", "BRANDTII")

# Set new column with Species as factor isntead of character 

Bats45$fSpecies<- as.factor(Bats45$Species)
str(Bats45)
summary(Bats45$BatID)

##Attempting to quantify the number of onsite plots per individual bat 
##This did not work very well - I ended up just altering my CSV and making it useless

#onsite<-Bats45 %>% count(Bats45$BatID, Bats45$BD=="Yes"|Bats45$Dir.obs=="Yes", sort=TRUE)
#summary(onsite)

  
#Fix date/time

names(Bats45)[1:2]<-c("date.cap","time.cap")

Bats45$date.cap<-as.character(Bats45$date.cap)

Bats45$time.cap<-as.character(Bats45$time.cap)

Bats45$datetime<-paste(Bats45$date.cap,Bats45$time.cap)

Bats45$datetime.posix<-as.POSIXct(strptime(Bats45$datetime,"%d.%m.%y %H:%M"),tz="GMT")

#Important summaries
summary(Bats45$fSpecies)
summary(Bats45$Dir.obs=="Yes"|Bats45$BD=="Yes")
summary(Bats45$Dir.obs=="Yes")
summary(Bats45$BD=="Yes")

#onsite<-filter(Bats45, BD=="Yes"|Dir.obs=="Yes")
#summary(onsite$BatID)
#summary(onsite$fSpecies)

#Load packages

library(sp)
library(sf)
library(rgeos)
library(adehabitatMA)
library(crs)
library(dplyr)
library(adehabitatHR)
library(Rcpp)
library(raster)
library(rgdal)
library(ggplot2)
library(ResourceSelection)
library(mapview)
library(amt)
library(sjPlot)
#library(broom) does not cooperate

## Create spatial data frame ## *****************************************************************

UTMs <- coordinates(Bats45)<-Bats45[,c("Easting","Northing")]
#assigning coordinates 

proj4string(Bats45)<-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#Define the projection 

Bats45<-Bats45[!is.na(Bats45$Easting) | !is.na(Bats45$Northing),]
#Remove NA values 

## Adding relevant shape file(s) ## *****************************************************************

system("unzip Nittedal.zip")

require(rgdal)
shape<- readOGR(dsn=".", layer= "32_0233AR5_ArealressursGrense_KURVE")
plot(shape, col='khaki1')

## Import Raster Data ##  *****************************************************************

system("unzip Metrics250")
r250<-raster("stdmetrics_z.grd")

proj4string(r250)<-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

plot(r250)
str(r250)

r16<-raster("stdmetrics_z_16.grd")
proj4string(r16)<-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
plot(r16)
str(r16)


#Cannot see attributes... 
r16@data@attributes[1]


## Or...
r <- raster()
r[] <- 1:ncell(r)
r <- writeRaster(r,"Metrics250.grd", overwrite=TRUE) %>% projectRaster(crs="+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
str(r)
crs(r)
plot(r, col = terrain.colors(3), main="try this")

r

b<-brick(r250)
str(b)
nlayers(b)
plot(b)

## What is likely happening here is that all of the variable are
## stacked on top of each other and this makes it difficult to 
## visualize anything meaningful. Need to isloate variables and plot
## them separately. 
library(ggplot2)
library(Rgb)
library(munsell)
library(colorspace)
library(rlang)

rsst <- raster('Metrics250')
plot(rsst)

dat_grid <- data.frame(xyFromCell(r16, 1:ncell(r16)),
                       vals = r16[])

ggplot(dat_grid, aes(x = x, y = y, fill = vals)) +
  geom_tile() 

#See no difference between b and b1.... 
#I am a bit worried about proceeding in case I have somehow 
#imported these files incorrectly... but here I go... 

## Plotting points ## *****************************************************************
## Distinguishing between species as well as individuals

#Plotting two species together 
Species.col<- c("blue","yellow")
names(Species.col)<-levels(Bats45$fSpecies)
plot(Bats45, col = Species.col[Bats45$fSpecies],pch=20, cex=0.5)

#MBRA blue, MMYS yellow

# Plotting all individuals with different colors 
BatID.col<- c("red","gold","blue","tomato","sienna1","indianred","orangered","lightseagreen","orangered3","lightcoral","cornflowerblue","salmon4","cyan","brown1","brown4","chocolate1","skyblue","steelblue4","darkgoldenrod1","azure","brown3")
plot(Bats45, col=BatID.col[Bats45$BatID],pch=20, cex=0.5)
#MBRA blue/green, MMYS red/yellow

#Plotting species in parts
#plot(Bats45[Bats45$Species=="MYSTACINUS",], col="yellow")
#plot(Bats45[Bats45$Species=="BRANDTII",], col= "blue", add=TRUE)
#summary(Bats45$Species=="MYSTACINUS")
#686 MBRA, 816 MMYS 


## MCPS ##  *****************************************************************

# MCP of both species groups 
#temp <- Bats45
#temp@data <- data.frame(id = Bats45@data$fSpecies)#--
#Bats45.mcp <- mcp(temp, percent = 95)
#plot(Bats45, col = Species.col[Bats45$fSpecies],pch=20, cex=0.5)
#plot(Bats45.mcp, border = adjustcolor(Species.col), col = adjustcolor(Species.col, alpha = 0.3), add = TRUE)

# MCP of all bats individually 
temp <- Bats45
temp@data <- data.frame(id = Bats45@data$BatID)#---
Bats.mcp <- mcp(temp, percent = 95)
plot(Bats45, col = BatID.col[Bats45$BatID], pch=20, cex=0.5)
plot(Bats.mcp, border = adjustcolor(BatID.col), col = adjustcolor(BatID.col, alpha = 0.3), add = TRUE)

# MCP for each individual bat - Nora and Thea examples
#temp@data <- data.frame(id = Bats45@data$BatID=="Nora")#--
#Nora.mcp <- mcp(temp, percent = 95)
#Nora.col<-("black")
#names(Nora.col)<-levels(Bats45$BatID=="Nora")
#plot(Bats45, col = Nora.col[Bats45$BatID=="Nora"], pch=20, cex=0.5)
#plot(Nora.mcp, border = adjustcolor(Nora.col), col = adjustcolor(Nora.col, alpha = 0.3), add = TRUE)
## MCP loads but not the points ## help here **** 

#temp@data <- data.frame(id = Bats45@data$BatID=="Thea")
#Thea.mcp <- mcp(temp, percent = 95)
#Thea.col<-("orange")
#names(Thea.col)<-levels(Bats45$BatID=="Thea")
#plot(Bats45, col = Thea.col[Bats45$BatID=="Thea"],pch=20, cex=0.5)
#plot(Thea.mcp, border = adjustcolor(Thea.col), col = adjustcolor(Thea.col, alpha = 0.3), add = TRUE)

## Results of MCPs 
#Is it necessary to export the MCPs as shape files or otherwise export their results?

area=Bats.mcp$area


## RSF ##  *****************************************************************

## The basic set up for a RSF should look like:
## 1. Preparing used vs. available points
## 2. Prearing linear measures (not sure if that would be Nitelva or roads or both in our case)
## 3. Then prepare raster covariates (canopy height, vegetation density, etc.)
## 4. THEN build the actual RSF with logistic regressions 

## 1
str(Bats45)
summary(Bats.mcp)
## The "Used" space has already been defined through mcps
## Need to create available points next. 

## From: https://www.danaseidel.com/MovEco-R-Workshop/Materials/Day5/SelectionFunctions/
ID <- unique(Bats45$BatID)
availables <- list()
for(i in 1:length(ID)){
  st_sample(st_as_sf(Bats.mcp)[i,], 10*nrow(filter(as(Bats45,"Spatial")[1], ID == ID[i]))) %>%
    
    st_sf(geometry = .) %>%
    mutate(ID = ID[i], 
           timestamp = NA, 
           Used = 0) -> availables[[i]] 
}

availables %>% 
  do.call(rbind,.) %>%
  rbind(Bats45, .) -> Bats_all

## Issues with Bats45 not being subsettable


## Code from: https://www.r-bloggers.com/home-range-estimation-mcp/
## This runs but it does not help with the next step
#ID<-(Bats45$BatID)
#xy<-data.frame(Bats45$Easting, Bats45$Northing, "id"=ID)
#centroid <- apply(xy[, 1:2], 2, mean)
#d <- sqrt(((xy[, 1] - centroid[1])^2) + ((xy[, 2] - centroid[2])^2))
#indx <- 1:length(d)
#pct <- indx[d <= quantile(d, .8)]
#mcp.pts <- xy[pct, ]

#https://ecosystems.psu.edu/research/labs/walter-lab/manual/chapter-8-resource-selection/8-1-minimum-convex-polygon
##Another way of making MCPs that may be better for RSF, may not be necessary
#xysp<-SpatialPoints(coords)
#proj4string(xysp)<-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#sppt<-data.frame(xysp)
#idsp<-data.frame(Bats45[3])
#merge<-data.frame(idsp)
#coordinates(merge)<-sppt
#plot(merge)
#str(merge)

#cp <- mcp(merge[,1], percent=100)
## The home-range size
#as.data.frame(cp)
## Plot the home ranges
#plot(cp)
#plot(cp[2,])#only plot deer D8
## ... And the relocations (Fig. 8.2)
#plot(merge, col=as.data.frame(merge)[,1], add=TRUE)

#Next; https://terpconnect.umd.edu/~egurarie/teaching/SpatialModelling_AKTWS2018/6_RSF_SSF.html
## This seems to work the best so far but still has hang ups. 

xy.obs<-Bats45[sample(1:nrow(Bats45), 1502), c("Easting","Northing", "BatID")]
summary(xy.obs)
xy.random<-spsample(Bats.mcp, 200, "regular")

## I cannot sample any more than 3 random points and I do not know why
## Does not change if I change the mcp to 95 or 100% 
## OR if I increase the number of observations to include in xy.obs
##Is it because I have multiple polygons?

plot(xy.random, asp = 1, col = "darkblue", pch = 19, cex = 0.5)
points(xy.obs, pch = 19, col = "orange", cex = 0.5)
summary(xy.obs)

## 
library(amt)
rp<-random_points(Bats.mcp, n=200, type="regular")
plot(rp)



## BONUS SCRIPT ## *****************************************************************


## Plotting buffers 
#coords<-Bats45[,c("Easting","Northing")]
#Bats.spdf<-SpatialPointsDataFrame(coords=coords, data=Bats45, proj4string = CRS(crs))
#Bat_buffer<-st_buffer(temp, dist=22,nQuadSegs = 1, endCapStyle = "ROUND",joinStyle = "ROUND",mitreLimit = 1)
#plot(Bat_buffer)


#nQuadSegs needs to be >/= 22, and then it insists on plotting the first 10 attributs 
#nQuadSegs=integer; number of segments per quadrant (fourth of a circle), for all or perfeature
#Need to redesign Thea_object 



