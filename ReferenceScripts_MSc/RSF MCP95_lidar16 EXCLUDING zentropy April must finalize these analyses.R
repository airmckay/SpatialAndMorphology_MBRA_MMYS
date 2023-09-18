#################################################################################################################################
##################      RSF 95% individual MCPs LiDAR resolution 16   ##########################################################
#################################################################################################################################

rm(list = ls())
# Clear work space

setwd("~/Desktop/R - Thesis /DATA")

library(ggplot2)
library(ResourceSelection)
library(MASS)
library(lme4)

data_2<-read.csv("data_95MCP16_final.csv",sep=",", dec=".", header=TRUE)
head(data_2)
 
#Make new variable zcv (coefficient of variation in height) 
data_2$zcv<-data_2$zsd/data_2$zmean
head(data_2)
dim(data_2)
#[1] 2970   18
summary(data_2) 

data_2 <- na.omit(data_2)
dim(data_2)
#[1] 2970   18  #There were no more NAs

data_2$Species=as.factor(data_2$Species)
data_2$BatID=as.factor(data_2$BatID)
data_2[,4:18]=scale(data_2[,4:18],scale=TRUE)#standardize lidar data to mean of zero
str(data_2)
summary(data_2)

#Check correlations between pairs of lidar variables
lidar.data <- data_2[, 4:length(data_2)]
round(cor(lidar.data), 2)

#              zmax zmean   zsd pzabovezmean pzabove0.5  zq10  zq20  zq30  zq40  zq50  zq60  zq70  zq80  zq90   zcv
#zmax          1.00  0.86  0.96         0.36       0.80  0.25  0.38  0.51  0.64  0.74  0.80  0.85  0.89  0.94 -0.18
#zmean         0.86  1.00  0.84         0.65       0.90  0.46  0.64  0.78  0.89  0.94  0.96  0.97  0.96  0.93 -0.46
#zsd           0.96  0.84  1.00         0.42       0.76  0.09  0.23  0.39  0.58  0.72  0.81  0.87  0.92  0.96 -0.28
#pzabovezmean  0.36  0.65  0.42         1.00       0.61  0.25  0.41  0.55  0.65  0.68  0.69  0.67  0.63  0.55 -0.81
#pzabove0.5    0.80  0.90  0.76         0.61       1.00  0.43  0.58  0.69  0.77  0.82  0.85  0.86  0.87  0.86 -0.47
#zq10          0.25  0.46  0.09         0.25       0.43  1.00  0.75  0.58  0.47  0.40  0.37  0.33  0.31  0.28 -0.23
#zq20          0.38  0.64  0.23         0.41       0.58  0.75  1.00  0.83  0.68  0.59  0.54  0.50  0.47  0.43 -0.30
#zq30          0.51  0.78  0.39         0.55       0.69  0.58  0.83  1.00  0.86  0.76  0.70  0.65  0.61  0.57 -0.36
#zq40          0.64  0.89  0.58         0.65       0.77  0.47  0.68  0.86  1.00  0.91  0.85  0.80  0.76  0.71 -0.41
#zq50          0.74  0.94  0.72         0.68       0.82  0.40  0.59  0.76  0.91  1.00  0.95  0.91  0.86  0.81 -0.44
#zq60          0.80  0.96  0.81         0.69       0.85  0.37  0.54  0.70  0.85  0.95  1.00  0.97  0.93  0.88 -0.46
#zq70          0.85  0.97  0.87         0.67       0.86  0.33  0.50  0.65  0.80  0.91  0.97  1.00  0.97  0.93 -0.46
#zq80          0.89  0.96  0.92         0.63       0.87  0.31  0.47  0.61  0.76  0.86  0.93  0.97  1.00  0.97 -0.46
#zq90          0.94  0.93  0.96         0.55       0.86  0.28  0.43  0.57  0.71  0.81  0.88  0.93  0.97  1.00 -0.42
#zcv          -0.18 -0.46 -0.28        -0.81      -0.47 -0.23 -0.30 -0.36 -0.41 -0.44 -0.46 -0.46 -0.46 -0.42  1.00

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

# Use a hypothesis test to determine which variables have strong influence on M.mys 
# Pass = p < 0.001

#                                    Fixed effects:
#                         Estimate Std. Error z value Pr(>|z|)  
summary(fit1)#(Intercept)0.007389   0.049630   0.149    0.882
summary(fit2)#zmax         0.37634    0.05142   7.319  2.5e-13 ***
summary(fit3)#zmean        0.26117    0.04992   5.232 1.68e-07 *** 
summary(fit4)#zsd          0.33848    0.05056   6.694 2.17e-11 ***
summary(fit5)#pzabovezm   -0.11812    0.04833  -2.444   0.0145 *
summary(fit6)#pzabove0.5   0.20394    0.04965   4.107    4e-05 ***
summary(fit7)#zq10         0.06303    0.04927   1.279    0.201
summary(fit8)#zq20         0.097599   0.049418   1.975   0.0483 *
summary(fit9)#zq30         0.1275532  0.0478417   2.666  0.00767 **
summary(fit10)#zq40        0.161213   0.048158   3.348 0.000815 ***
summary(fit11)#zq50        0.20971    0.04834   4.339 1.43e-05 ***
summary(fit12)#zq60        0.23878    0.04858   4.915 8.86e-07 ***
summary(fit13)#zq70        0.24015    0.04884   4.917 8.79e-07 ***
summary(fit14)#zq80        0.26150    0.04921   5.313 1.08e-07 ***
summary(fit15)#zq90        0.29616    0.04999   5.924 3.14e-09 ***
summary(fit16)#zcv         0.18268    0.05124   3.565 0.000364 ***

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

# Use a hypothesis test to determine which variables have strong influence on M.mys 
# Pass = p < 0.001

#                                   Fixed effects:
#                         Estimate Std. Error z value Pr(>|z|)  
summary(fit1)#(Intercept) -0.02972    0.05452  -0.545    0.586
summary(fit2)#zmax         0.29726    0.05737   5.182  2.2e-07 ***
summary(fit3)#zmean        0.220823   0.058899   3.749 0.000177 ***
summary(fit4)#zsd          0.31289    0.05936   5.271 1.35e-07 *** 
summary(fit5)#pzabovezm   -0.25370    0.05798  -4.376 1.21e-05 ***
summary(fit6)#pzabove0.5   0.195872   0.057376   3.414 0.000641 ***
summary(fit7)#zq10         0.03147    0.05658   0.556    0.578
summary(fit8)#zq20        -0.002611   0.056265  -0.046    0.963
summary(fit9)#zq30         0.06071    0.05933   1.023    0.306
summary(fit10)#zq40        0.16376    0.05993   2.733  0.00628 **
summary(fit11)#zq50        0.197435   0.060407   3.268  0.00108 **
summary(fit12)#zq60        0.203734   0.060078   3.391 0.000696 ***
summary(fit13)#zq70        0.225719   0.059554   3.790 0.000151 ***
summary(fit14)#zq80        0.26209    0.05952   4.403 1.07e-05 ***
summary(fit15)#zq90        0.27635    0.05853   4.722 2.34e-06 ***
summary(fit16)#zcv         0.18216    0.05622   3.240  0.00119 **

#------------Construct models for both species simultanously
#------------use significant level p<0.01
#------------avoiding collinearity (r<0.6)

#Passed Hypothesis test for M. mys: 

#zmax, zmean, zsd, pzabove0.5, zq40, zq50, zq60, zq70, zq80, zq90, zcv

#Passed Hypothesis test for M. bra:

#zmax, zmean, zsd, pzabovezmean, pzabove0.5, zq60, zq70, zq80, zq90

#Excluding variables that are highly correlated, we are left with:

#zmax, zcv

#This is the full (most complex) model:
Mfull <- glmer(used ~ Species + zmax + zcv + Species:zmax + Species:zcv +  (1|Species/BatID), data=data_2, family=binomial(link="logit"),nAGQ = 0)
drop1(Mfull, test = "Chi") #Likelihood ratio test
#Single term deletions
#Df    AIC     LRT Pr(Chi)
#<none>          4003.5                
#Species:zmax  1 4003.6 2.09621  0.1477
#Species:zcv   1 4002.1 0.57034  0.4501

#Neither of the species / variable interactions were significant, drop both. Cannot be dropped further

M2 <- glmer(used ~ Species + zmax + zcv + (1|Species/BatID), data=data_2, family=binomial(link="logit"),nAGQ = 0)
summary(M2)

#Present output in a table in results with 250 resolution 
#Generalized linear mixed model fit by maximum likelihood (Adaptive
 #                                                         Gauss-Hermite Quadrature, nAGQ = 0) [glmerMod]
#Family: binomial  ( logit )
#Formula: used ~ Species + zmax + zcv + (1 | Species/BatID)
#Data: data_2

#AIC      BIC   logLik deviance df.resid 
#4001.8   4037.8  -1994.9   3989.8     2964 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.3069 -0.8934 -0.7584  0.9455  1.4779 

#Random effects:
#  Groups        Name        Variance  Std.Dev. 
#BatID:Species (Intercept) 0.000e+00 0.000e+00
#Species       (Intercept) 1.246e-18 1.116e-09
#Number of obs: 2970, groups:  BatID:Species, 21; Species, 2

# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.04401    0.05635   0.781    0.435    
# SpeciesMys  -0.09197    0.07667  -1.200    0.230    
# zmax         0.39282    0.03930   9.995  < 2e-16 ***
#   zcv          0.26375    0.04037   6.533 6.44e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) SpcsMy zmax  
# SpeciesMys -0.746              
# zmax        0.146 -0.185       
# zcv         0.034 -0.021  0.215
##refitting model with GLM to scetch results
fit <- glm(used ~ Species + zmax + zcv, family=binomial(), data=data_2)
summary(fit) #gives approximately the same parameter estimates as the GLMM,

#------------------------------Figure zcv

##getting predicted use as a function of zcv, at average zmax 
range.zcv <- range(data_2$zcv)
range.zcv
#[1] -1.407101  5.535054

plotting_dfm <- expand.grid(zcv = seq(from=-1.41, to = 5.54, by=0.001),
                            Species    = c("Mbra","Mys"),
                            zmax = mean(data_2$zmax))
plotting_dfm$preds <- plogis( predict(fit , newdata=plotting_dfm))


##plotting the predicted response on the two covariates
windows()
pl <- ggplot(plotting_dfm, aes(x=zcv, y =preds, color=as.factor(Species)))
pl + 
  geom_point( ) +
  ggtitle("Predicted Use by zcv and Species") + 
  ggplot2::ylab("Predicted Use")

#------------------------------Figure zmax

##getting predicted use as a function of zmax, at average zcv
range.zmax <- range(data_2$zmax)
range.zmax
#[1] -0.9218276  2.6942289

plotting_dfm <- expand.grid(zmax = seq(from=-0.922, to = 2.69, by=0.001),
                            Species    = c("Mbra","Mys"),
                            zcv = mean(data_2$zcv))
plotting_dfm$preds <- plogis( predict(fit , newdata=plotting_dfm))


##plotting the predicted response on the two covariates
windows()
pl <- ggplot(plotting_dfm, aes(x=zmax, y =preds, color=as.factor(Species)))
pl + 
  geom_point( ) +
  ggtitle("Predicted Use by zmax and Species") + 
  ggplot2::ylab("Predicted Use")


