rm(list = ls())
# Clear work space

setwd("//nmbu.no/my/home/Desktop/DATA")
#load data

dataset<-read.csv("batdata (1).csv", sep=";", dec=",", header=TRUE)

#---QUICK INITIAL LOOK AT THE 'dataset' DATA OBJECT
str(dataset)
head(dataset)
dim(dataset)
names(dataset)
summary(dataset$RFA) 


source(file = "HighstatLibV4.R")
######################################################################
#Load all packages
library(lattice)
library(MASS)
library(ggplot2)

######################################################################
#Data exploration


histogram(dataset$Weight, breaks=10) #just inspecing the distribution of the weight data
histogram(dataset$RFA, breaks=5) #just inspecing the distribution of the RFA data
histogram(dataset$RFA, breaks=10) #just inspecing the distribution of the RFA data


#Outliers?              
dotchart(dataset$Weight, 
         main = "Weight",
         xlab = "Values of variable",
         ylab = "Order of the data")       #No extreme outliers

dotchart(dataset$RFA, 
         main = "RFA",
         xlab = "Values of variable",
         ylab = "Order of the data")       #No extreme outliers

#Break it down into species

xyplot(RFA ~ Weight|Species, data = dataset, 
       layout = c(2,1), ylab = "RFA",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       panel = function(x,y) {
         panel.points(x,y, pch = 16, cex = 0.5, col =1)
         if (length(x) > 5){
           tmp <- lm(y~x)
           panel.abline(tmp, col = 1, lwd = 3)}
       })

#Break it down into species and sex
xyplot(RFA ~ Weight|Sex*Species, data = dataset, 
       layout = c(2,2), ylab = "RFA",
       strip = function(bg = 'white', ...) 
         strip.default(bg = 'white', ...),
       panel = function(x,y) {
         panel.points(x,y, pch = 16, cex = 0.5, col =1)
         if (length(x) > 5){
           tmp <- lm(y~x)
           panel.abline(tmp, col = 1, lwd = 3)}
       })

#####################################################################
## Aim model the relationhip between Weight and RFA, and check if it depend on species and sex



m1 <- lm(RFA ~ Sex*Species, data = dataset)
summary(m1)
anova(m1)

#Can you simplify the model?

drop1(m1, test="Chi")

#No, because the Species*Sexinteraction is significant

bwplot(RFA ~  Sex | Species, 
       strip = strip.custom(bg = 'white'),
       cex = .5,
       data = dataset,
       ylab= "Right forearm length (mm)",
       main= "Right Forearm Length (mm) Distribution")

# Violin Plot alternative
v<-ggplot(dataset, aes(x=Sex, y=RFA, fill=Species)) + 
        geom_violin(trim=FALSE) + 
        labs(x = "Sex", y = "Right Forearm Length (mm)")+
        geom_boxplot(width=0.1,position=position_dodge(0.9))+ 
        theme(legend.position="none")
v

m2 <- lm(RFA ~ Species, data = dataset)
summary(m2)
anova(m2)



## ______________________________________________________________ ##

# Average forearm length of each species, together and by sex with confidence intervals 
library(dplyr)
library(scales)
library(ggplot2)



summary(dataset$Species)
str(dataset)
head(dataset)
summary(dataset$Sex)


BF <-filter(B, Sex == "Female")
BM<-filter(B, Sex == "Male")
summary(BF$RFA)
summary(BM$RFA)

summary(M$RFA, Sex == "Male")
summary(M$RFA, Sex == "Female")
sd(M$RFA)
summary(B$RFA)
sd(B$RFA)

summary(B$Sex)
summary(M$Sex)


