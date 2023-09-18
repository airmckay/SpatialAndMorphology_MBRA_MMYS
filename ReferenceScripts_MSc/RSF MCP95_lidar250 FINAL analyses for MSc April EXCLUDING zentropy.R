#################################################################################################################################
##################      RSF 95% individual MCPs LiDAR resolution 250   ##########################################################
#################################################################################################################################

rm(list = ls())
# Clear work space
 
setwd("~/Desktop/R - Thesis /DATA")

library(ggplot2)
library(ResourceSelection)
library(MASS)
library(lme4)

data_2<-read.csv("data_95MCP250_excl_zentropy.csv",sep=",", dec=".", header=TRUE)
head(data_2)
 
#Make new variable zcv
data_2$zcv<-data_2$zsd/data_2$zmean
head(data_2)
dim(data_2)
#[1] 2998   18
summary(data_2) 

data_2$Species=as.factor(data_2$Species)
data_2$BatID=as.factor(data_2$BatID)
data_2[,4:18]=scale(data_2[,4:18],scale=TRUE)#standardize lidar data to mean of zero
str(data_2)
summary(data_2)

#Check correlations between pairs of lidar variables
lidar.data <- data_2[, 4:length(data_2)]
round(cor(lidar.data), 2)

#              zmax zmean   zsd pzabovezmean pzabove0.5  zq10  zq20  zq30  zq40  zq50  zq60  zq70  zq80  zq90   zcv
#zmax          1.00  0.73  0.89         0.07       0.72  0.01  0.04  0.17  0.37  0.49  0.58  0.66  0.72  0.80 -0.06
#zmean         0.73  1.00  0.91         0.65       0.91  0.08  0.16  0.40  0.72  0.86  0.94  0.97  0.96  0.93 -0.46
#zsd           0.89  0.91  1.00         0.41       0.86  0.00  0.04  0.19  0.48  0.64  0.76  0.87  0.92  0.97 -0.40
#pzabovezmean  0.07  0.65  0.41         1.00       0.60  0.05  0.13  0.32  0.54  0.65  0.69  0.68  0.64  0.54 -0.69
#pzabove0.5    0.72  0.91  0.86         0.60       1.00  0.09  0.20  0.39  0.60  0.73  0.82  0.87  0.90  0.89 -0.50
#zq10          0.01  0.08  0.00         0.05       0.09  1.00  0.59  0.21  0.10  0.07  0.06  0.05  0.04  0.03 -0.03
#zq20          0.04  0.16  0.04         0.13       0.20  0.59  1.00  0.46  0.22  0.16  0.14  0.11  0.10  0.08 -0.07
#zq30          0.17  0.40  0.19         0.32       0.39  0.21  0.46  1.00  0.56  0.43  0.36  0.31  0.28  0.25 -0.16
#zq40          0.37  0.72  0.48         0.54       0.60  0.10  0.22  0.56  1.00  0.84  0.73  0.64  0.57  0.51 -0.26
#zq50          0.49  0.86  0.64         0.65       0.73  0.07  0.16  0.43  0.84  1.00  0.91  0.81  0.74  0.66 -0.34
#zq60          0.58  0.94  0.76         0.69       0.82  0.06  0.14  0.36  0.73  0.91  1.00  0.93  0.86  0.78 -0.40
#zq70          0.66  0.97  0.87         0.68       0.87  0.05  0.11  0.31  0.64  0.81  0.93  1.00  0.96  0.88 -0.45
#zq80          0.72  0.96  0.92         0.64       0.90  0.04  0.10  0.28  0.57  0.74  0.86  0.96  1.00  0.95 -0.48
#zq90          0.80  0.93  0.97         0.54       0.89  0.03  0.08  0.25  0.51  0.66  0.78  0.88  0.95  1.00 -0.49
#zcv          -0.06 -0.46 -0.40        -0.69      -0.50 -0.03 -0.07 -0.16 -0.26 -0.34 -0.40 -0.45 -0.48 -0.49  1.00

#The lidar variables are highly correlated!
#Test each of them separately, in combination with Species and the lidar variable ? Species interaction

Mmys <- subset(data_2, Species == "Mys")
Mbra <- subset(data_2, Species == "Mbra")


#----------------Exploring lidar variables influencing probability of use in M_mys

#We will fit mixed models with BatID as random effect to account for among-individuals variation
fit1 = glmer(used ~ 1|BatID, data=Mmys, family=binomial(link="logit"),
             nAGQ = 0)#Intercept model
fit2 = glmer(used ~ zmax + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#zmax model - not taking into account which species
fit3 = glmer(used ~ zmean + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit4 = glmer(used ~ zsd + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit5 = glmer(used ~ pzabovezmean + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit6 = glmer(used ~ pzabove0.5 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit7 = glmer(used ~ zq10 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit8 = glmer(used ~ zq20 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit9 = glmer(used ~ zq30 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit10 = glmer(used ~ zq40 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit11 = glmer(used ~ zq50 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit12 = glmer(used ~ zq60 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit13 = glmer(used ~ zq70 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit14 = glmer(used ~ zq80 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit15 = glmer(used ~ zq90 + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)#
fit16 = glmer(used ~ zcv + (1|BatID), data=Mmys,
              family=binomial(link="logit"),nAGQ = 0)
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5)
summary(fit6)
summary(fit7)
summary(fit8)
summary(fit9)
summary(fit10)
summary(fit11)
summary(fit12)
summary(fit13)
summary(fit14)
summary(fit15)
summary(fit16)


#                                         Fixed effects:
#                                  Estimate  Std. Error    z value Pr(>|z|)  
#summary(fit1)       (Intercept)  -0.001226    0.049523    -0.025     0.98
#summary(fit2)       zmax          0.56845     0.05678     10.012   <2e-16 ***     
#summary(fit3)       zmean         0.28960     0.04925      5.881 4.08e-09 ***    
#summary(fit4)       zsd           0.37371     0.05079      7.358 1.87e-13 ***       
#summary(fit5)       pzabovezmean -0.104977    0.050002    -2.099     0.0358 *
#summary(fit6)       pzabove0.5    0.23349     0.05022      4.649 3.33e-06 ***
#summary(fit7)       zq10          1.08648    28.07220      0.039     0.969
#summary(fit8)       zq20          0.067619    0.082128     0.823     0.410      
#summary(fit9)       zq30          0.062207    0.048318     1.287     0.198      
#summary(fit10)      zq40          0.11303     0.04620      2.447     0.0144 *      
#summary(fit11)      zq50          0.20342     0.04687      4.340 1.43e-05 ***      
#summary(fit12)      zq60          0.26371     0.04744      5.559 2.71e-08 ***     
#summary(fit13)      zq70          0.26479     0.04819      5.495 3.92e-08 ***      
#summary(fit14)      zq80          0.27160     0.04910      5.531 3.18e-08 ***      
#summary(fit15)      zq90          0.29352     0.04998      5.872 4.29e-09 ***      
#summary(fit16)      zcv           0.36400     0.06341      5.740 9.46e-09 ***       
                                    



#----------------Exploring lidar variables influencing probability of use in M_bra

#We will fit mixed models with BatID as random effect to account for among-individuals variation
fit1 = glmer(used ~ 1|BatID, data=Mbra, family=binomial(link="logit"),
             nAGQ = 0)#Intercept model
fit2 = glmer(used ~ zmax + (1|BatID), data=Mbra,
             family=binomial(link="logit"),nAGQ = 0)#zmax model - not taking into account which species
fit3 = glmer(used ~ zmean + (1|BatID), data=Mbra,
             family=binomial(link="logit"),nAGQ = 0)#
fit4 = glmer(used ~ zsd + (1|BatID), data=Mbra,
             family=binomial(link="logit"),nAGQ = 0)#
fit5 = glmer(used ~ pzabovezmean + (1|BatID), data=Mbra,
             family=binomial(link="logit"),nAGQ = 0)#
fit6 = glmer(used ~ pzabove0.5 + (1|BatID), data=Mbra,
             family=binomial(link="logit"),nAGQ = 0)#
fit7 = glmer(used ~ zq10 + (1|BatID), data=Mbra,
             family=binomial(link="logit"),nAGQ = 0)#
fit8 = glmer(used ~ zq20 + (1|BatID), data=Mbra,
             family=binomial(link="logit"),nAGQ = 0)#
fit9 = glmer(used ~ zq30 + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#
fit10 = glmer(used ~ zq40 + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#
fit11 = glmer(used ~ zq50 + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#
fit12 = glmer(used ~ zq60 + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#
fit13 = glmer(used ~ zq70 + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#
fit14 = glmer(used ~ zq80 + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#
fit15 = glmer(used ~ zq90 + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#
fit16 = glmer(used ~ zcv + (1|BatID), data=Mbra,
              family=binomial(link="logit"),nAGQ = 0)#

summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5)
summary(fit6)
summary(fit7)
summary(fit8)
summary(fit9)
summary(fit10)
summary(fit11)
summary(fit12)
summary(fit13)
summary(fit14)
summary(fit15)
summary(fit16)

#                                         Fixed effects:
#                                 Estimate Std. Error z value Pr(>|z|)  
#summary(fit1)      (Intercept)   -0.007315   0.054094  -0.135    0.892 
#summary(fit2)       zmax          0.59866    0.05753   10.406   <2e-16 ***    
#summary(fit3)       zmean         0.33770    0.06338    5.328 9.91e-08 ***   
#summary(fit4)       zsd           0.49433    0.06143    8.047 8.46e-16 ***      
#summary(fit5)       pzabovezmean -0.42778    0.05665   -7.551 4.31e-14 ***
#summary(fit6)       pzabove0.5    0.33931    0.05875    5.775 7.68e-09 ***
#summary(fit7)       zq10          2.55180    1.68871    1.511    0.131      
#summary(fit8)       zq20          0.139719   0.069896   1.999    0.0456 *     
#summary(fit9)       zq30          0.091476   0.060772   1.505    0.132     
#summary(fit10)      zq40          0.033818   0.062137   0.544    0.586     
#summary(fit11)      zq50          0.081969   0.062729   1.307    0.191     
#summary(fit12)      zq60          0.18005    0.06366    2.828    0.00468 **    
#summary(fit13)      zq70          0.31328    0.06315    4.961 7.01e-07 ***     
#summary(fit14)      zq80          0.36014    0.06187    5.821 5.85e-09 ***     
#summary(fit15)      zq90          0.40007    0.06042    6.622 3.55e-11 ***     
#summary(fit16)      zcv           0.28229    0.05259    5.367 7.99e-08 ***      




#------------Construct models for both species simultanously
#------------use significant level p<0.01
#------------avoiding collinearity (r<0.6)


#For M_mys the following variables are significant at p<0.001 level:
#                                  Estimate  Std. Error    z value Pr(>|z|)  
#summary(fit2)       zmax          0.56845     0.05678     10.012   <2e-16 ***     
#summary(fit3)       zmean         0.28960     0.04925      5.881 4.08e-09 ***    
#summary(fit4)       zsd           0.37371     0.05079      7.358 1.87e-13 ***       
#summary(fit6)       pzabove0.5    0.23349     0.05022      4.649 3.33e-06 ***
#summary(fit11)      zq50          0.20342     0.04687      4.340 1.43e-05 ***      
#summary(fit12)      zq60          0.26371     0.04744      5.559 2.71e-08 ***     
#summary(fit13)      zq70          0.26479     0.04819      5.495 3.92e-08 ***      
#summary(fit14)      zq80          0.27160     0.04910      5.531 3.18e-08 ***      
#summary(fit15)      zq90          0.29352     0.04998      5.872 4.29e-09 ***      
#summary(fit16)      zcv           0.36400     0.06341      5.740 9.46e-09 *** 

#For M_bra the following variables are significant at p<0.001 level:
#                                 Estimate Std. Error z value Pr(>|z|)  
#summary(fit2)       zmax          0.59866    0.05753   10.406   <2e-16 ***    
#summary(fit3)       zmean         0.33770    0.06338    5.328 9.91e-08 ***   
#summary(fit4)       zsd           0.49433    0.06143    8.047 8.46e-16 ***      
#summary(fit5)       pzabovezmean -0.42778    0.05665   -7.551 4.31e-14 ***
#summary(fit6)       pzabove0.5    0.33931    0.05875    5.775 7.68e-09 ***
#summary(fit13)      zq70          0.31328    0.06315    4.961 7.01e-07 ***     
#summary(fit14)      zq80          0.36014    0.06187    5.821 5.85e-09 ***     
#summary(fit15)      zq90          0.40007    0.06042    6.622 3.55e-11 ***     
#summary(fit16)      zcv           0.28229    0.05259    5.367 7.99e-08 ***     

#Start with LiDAR variable with strongest "signal" on the response (highest z-value)
#Take correlations between pairs of variables into account (r > 0.6 too strongy correlated to be included in same model)
#We end up with a final model with Species and the following non-correlated (r < 0.6) LiDAR variables:
#zmax, zcv and zq60
#Include all SpeciesxLiDAR variable interaction because that is how we test for difference between the two species

#This is the full (most complex) model:
Mfull <- glmer(used ~ Species + zmax + zq60 + zcv + Species:zmax +  Species:zq60 + Species:zcv +  (1|Species/BatID), data=data_2, family=binomial(link="logit"),nAGQ = 0)
drop1(Mfull, test = "Chi") #Likelihood ratio test
#Single term deletions
#Model:
#  used ~ Species + zmax + zq60 + zcv + Species:zmax + Species:zq60 + 
#  Species:zcv + (1 | Species/BatID)
#
#Df    AIC     LRT   Pr(Chi)    
#<none>          3855.4                      
#Species:zmax  1 3856.3  2.8419 0.0918369 .    #Drop this term; not significant (p>0.05)
#Species:zq60  1 3864.7 11.2983 0.0007758 ***
#Species:zcv   1 3870.6 17.1746  3.41e-05 ***

#Fit reduced model without Species:zmax interaction
M2<- glmer(used ~ Species + zmax + zq60 + zcv + Species:zq60 + Species:zcv +  (1|Species/BatID), data=data_2, family=binomial(link="logit"),nAGQ = 0)
drop1(M2, test = "Chi") #Likelihood ratio test
#Single term deletions - NOTE! Only terms that are candidates for dropping are listed
#Main effects Species, zq60 and zcv are not cadidates for dropping because they are included in interactions
#
#Model:
#  used ~ Species + zmax + zq60 + zcv + Species:zq60 + Species:zcv + 
#  (1 | Species/BatID)
#Df    AIC     LRT   Pr(Chi)    
#<none>          3856.3                      
#zmax          1 3986.3 132.009 < 2.2e-16 ***
#Species:zq60  1 3862.8   8.521  0.003511 ** 
#Species:zcv   1 3869.5  15.247 9.434e-05 *** 
#
#All terms are significant and thus the model cannot be further reduced, so this is your final model

summary(M2)
#Random effects:
#  Groups        Name        Variance  Std.Dev.
#BatID:Species (Intercept) 0.0007192 0.02682 
#Species       (Intercept) 0.0000000 0.00000 
#Number of obs: 2998, groups:  BatID:Species, 21; Species, 2
#
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)      0.05584    0.05873   0.951 0.341726    
#SpeciesMys      -0.10827    0.08118  -1.334 0.182290    
#zmax             0.56227    0.05021  11.198  < 2e-16 ***
#zq60            -0.06224    0.07798  -0.798 0.424781    
#zcv              0.24791    0.05863   4.228 2.35e-05 ***
#SpeciesMys:zq60  0.26531    0.08975   2.956 0.003117 ** 
#SpeciesMys:zcv   0.38513    0.10030   3.840 0.000123 *** 

#This is your final model (NOTE: this is a GLMM model WITH batID nested within Species as random effect. 
#Present these results in a table in your thesis

############################################################
#Sketch model fit for model  to understand what the model is telling us

##refitting model with GLM to scetch results
fit <- glm(used ~ Species + zmax + zq60 + zcv + Species:zq60 + Species:zcv  , family=binomial(), data=data_2)
summary(fit) #gives approximately the same parameter estimates as the GLMM, thus the GLMM-based and GLM-based figure will look the same
#Note! No need to mention in your thesis that you fitted GLMs to make the figures 

#Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)      0.05540    0.05793   0.956 0.338929    
#SpeciesMys      -0.10791    0.08009  -1.347 0.177862    
#zmax             0.56158    0.05017  11.194  < 2e-16 ***
#zq60            -0.06163    0.07795  -0.791 0.429170    
#zcv              0.24744    0.05850   4.230 2.34e-05 ***
#SpeciesMys:zq60  0.26348    0.08968   2.938 0.003304 ** 
#SpeciesMys:zcv   0.38521    0.10018   3.845 0.000121 ***



#-----------------------------Figure zmax

##getting predicted use as a function of zmax, at average zq60 and zcv
range.zmax <- range(data_2$zmax)
range.zmax
#[1] -1.450767  2.243304

plotting_dfm <- expand.grid(zmax = seq(from=-1.42, to = 2.24, by=0.001),
                            Species    = c("Mbra","Mys"),
                            zq60 = mean(data_2$zq60),
                            zcv = mean(data_2$zcv))
plotting_dfm$preds <- plogis( predict(fit , newdata=plotting_dfm))

 
##plotting the predicted response on the two covariates
windows()
pl <- ggplot(plotting_dfm, aes(x=zmax, y =preds, color=as.factor(Species)))
pl + 
  geom_point( ) +
  ggtitle("Predicted Use by zmax and Species") + 
  ggplot2::ylab("Predicted Use")


#-----------------------------Figure zq60

##getting predicted use as a function of zq60, at average zmax and zcv
range.zq60 <- range(data_2$zq60)
range.zq60
#[1] -0.5654901  3.7842825

plotting_dfm <- expand.grid(zq60 = seq(from=-0.56, to = 3.78, by=0.001),
                            Species    = c("Mbra","Mys"),
                            zmax = mean(data_2$zmax),
                            zcv = mean(data_2$zcv))
plotting_dfm$preds <- plogis( predict(fit , newdata=plotting_dfm))


##plotting the predicted response on the two covariates
windows()
pl <- ggplot(plotting_dfm, aes(x=zq60, y =preds, color=as.factor(Species)))
pl + 
  geom_point( ) +
  ggtitle("Predicted Use by zq60 and Species") + 
  ggplot2::ylab("Predicted Use")



#------------------------------Figure zcv

##getting predicted use as a function of zcv, at average zmax and zq60
range.zcv <- range(data_2$zcv)
range.zcv
#[1] -1.039593  6.961540

plotting_dfm <- expand.grid(zcv = seq(from=-1.04, to = 5.32, by=0.001),
                            Species    = c("Mbra","Mys"),
                            zmax = mean(data_2$zmax),
                            zq60 = mean(data_2$zq60))
plotting_dfm$preds <- plogis( predict(fit , newdata=plotting_dfm))


##plotting the predicted response on the two covariates
windows()
pl <- ggplot(plotting_dfm, aes(x=zcv, y =preds, color=as.factor(Species)))
pl + 
  geom_point( ) +
  ggtitle("Predicted Use by zcv and Species") + 
  ggplot2::ylab("Predicted Use")

