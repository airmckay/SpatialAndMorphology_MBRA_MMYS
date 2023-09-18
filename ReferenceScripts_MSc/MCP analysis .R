rm(list = ls()) 
# Clear work space

setwd("~/Desktop/R - Thesis /DATA")

library(ggplot2) 
library(dplyr)

mcp.all<-read.csv("mcp_total.csv", header=TRUE)
head(mcp.all)

hr<-area$mcp.all~lm(fSpecies+obs,data=mcp.all)
summary(hr)
# Length   Class    Mode 
# 3 formula    call 

lm1<-lm(area~fSpecies+obs,data=mcp.all)
summary(lm1)

#Estimate Std. Error t value Pr(>|t|)   
#(Intercept)         86.87963   24.89248   3.490  0.00261 **
#  fSpeciesMYSTACINUS -61.91432   19.64918  -3.151  0.00553 **
#  obs                 -0.01378    0.26294  -0.052  0.95877   

drop1(lm1, test="Chi")

#area ~ fSpecies + obs
#Df Sum of Sq   RSS    AIC Pr(>Chi)   
#<none>                35308 161.97            
#fSpecies  1   19475.9 54784 169.20 0.002387 **
#  obs       1       5.4 35314 159.98 0.954848   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Species is significant but # of fixes is not 

lm2<-lm(area~fSpecies*obs,data=mcp.all)
summary(lm2)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)             145.7897    35.4146   4.117  0.00072 ***
#  fSpeciesMYSTACINUS     -144.5000    42.1269  -3.430  0.00319 ** 
#  obs                      -0.7867     0.4298  -1.830  0.08481 .  
#fSpeciesMYSTACINUS:obs    1.1210     0.5177   2.166  0.04485 *  

drop1(lm2, test="Chi")

#Model:
#  area ~ fSpecies * obs
#Df Sum of Sq   RSS    AIC Pr(>Chi)  
#<none>                    27674 158.86           
#fSpecies:obs  1    7634.5 35308 161.97   0.0237 *

#weakly significant interaction between species and observations

anova(lm1)

#Analysis of Variance Table

#Response: area
#             Df Sum Sq Mean Sq F value   Pr(>F)   
#  fSpecies   1  19642 19642.4 10.0136 0.005366 **
#  obs        1      5     5.4  0.0027 0.958770   
#  Residuals 18  35308  1961.6                    

anova(lm2)

#Response: area
#                Df  Sum Sq Mean Sq F value   Pr(>F)   
#  fSpecies      1 19642.4 19642.4 12.0664 0.002905 **
#  obs           1     5.4     5.4  0.0033 0.954783   
#  fSpecies:obs  1  7634.5  7634.5  4.6899 0.044850 * 
#  Residuals    17 27673.7  1627.9 


boxplot(area ~ fSpecies, strip = strip.custom(bg = 'white'),
        cex = .5,
        data = mcp.all,
        ylab = " ",
        xlab = "Species", 
        names = c("M. brandtii", "M. mystacinus"),
        col = c("indianred1", "turquoise3")
        ) 



plot(area~ fSpecies+obs, data=mcp.all)
library(ggplot2)
v<-ggplot(mcp.all, aes(x=fSpecies, y=area, fill=fSpecies)) + 
  geom_violin(trim=FALSE) + 
  labs(x = " ", y = " ")+
  theme(axis.text.x=element_blank()) +
  geom_boxplot(width=0.1) +
  theme(legend.position="none")
v

summary(mcp.all$area, fSpecies=="MYSTACINUS")
summary(mcp.all$area, fSpecies=="BRANDTII")

m<-filter(mcp.all, fSpecies=="MYSTACINUS")
b<-filter(mcp.all, fSpecies=="BRANDTII")
summary(m$area)
summary(b$area)
head(b$area)
summary(mcp.all$id)
str(mcp.all)
summary(b$area)
summary(m$area)

#plotting the predicted response 
library(Lahman)
if (require("Lahman")) {
  # Find 3 bats with larges home range
  tbl_df(mcp.all) %>%
    group_by(id) %>%
    tally(area) %>%
    top_n(3) 
}

#1 Astrid  189.
#2 Nora    129.
#3 Phoebe  154.

lm2<-lm(area~fSpecies-obs,data=mcp.all)
require(ggplot2)

t <- ggplot(mcp.all, aes(x=obs, y=area, col=fSpecies)) + geom_point() + ylab(make_unit_label("area", area)) + geom_smooth(method=lm, se=FALSE, fullrange=TRUE) 


t 
library(ggforc)
ggplot(mcp.all) + geom_point(aes(x = as.numeric(obs), y = as.numeric(area), col=fSpecies)) +  xlab(make_unit_label("obs", obs)) +ylab(make_unit_label("area", area))
rlang::last_error()

citation(package="MuMIn")
 
help(mcp)
