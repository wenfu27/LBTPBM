# Lyme paper
# A two-part Bayesian spatio-temporal modelling
# Author Wen Fu
# 25 May 2020

library(sp)
library(gstat)
library(dplyr)
library(purrr)
library(maps) 
library(RColorBrewer)
# Spatial plotting
library(rgdal)
library(raster)
library(rNOMADS)
library(RSQLite)
library(dbplyr)
#library(chron)
library(lattice)
library(ncdf4)
#library(plot3D)
library(ecmwfr)
library(keyring)
library(maptools)
library(ggplot2)
library(plyr)
#install.packages("INLA", dep = TRUE, repos = "https://inla.r-inla-download.org/R/stable")
library("INLA")
#install.packages('foreach')
#install.packages(c('spdep','spatstat','latticeExtra','JointModel'))
library("spdep")
library("rgdal")
library("RColorBrewer")
library("spatstat")
library("sp")
library("maptools")
library("latticeExtra")
library("gridExtra")
library("gstat")
library("ggplot2")
library("MASS")
library("JointModel")

# Input data 
setwd('C:\\Users\\Wen\\')

# IDarea1 = ID for each pixel 
# long= longitude, lat= latitude
# inc= LB seasonal incidence estimated by kriging interpolation
# deer= index of deer presence 
# ndvi= vegetation indice (Normalized difference vegetation index)
# maxST= averaged maximum soil temperature at season Q
# maxrh= averaged maximum relative humidity at season Q
# maxst_bf= averaged maximum soil temperature at season Q-1
# sd= saturation deficit, calculated by maxST and maxrh
# ndd14= number of dry days (14-day increments per unit)
# citique = proportion of tick bite reports
# rodent = species richeness of rodents
# season = Winter Jan to Mar, Spring Apr to Jun, Summer Jul to Sept, Autumn Oct to Dec
# IDtime1 = 1, 2, 3 to 16 quarters during 2016-19 


############### classify the covariates 

## LB dataset 2016-19 4 years
Dat<- read.csv2('Dat1619.csv')
Dat$deer_cat <- cut(Dat$deer, c(0, 0.4, 0.6, 0.8, 1), include.lowest = TRUE)
Dat$ndvi_cat<- cut(Dat$ndvi, c(-1, 0.59, 1), include.lowest = TRUE)
Dat$sd_rd<- round(Dat$sd)
Dat$sd_cat<- cut(Dat$sd_rd, c(1, 3, 5, 9, 23), include.lowest = TRUE)
Dat$maxst_bf_rd<- round(Dat$maxst_bf)
Dat$maxstbf_cat<- cut(Dat$maxst_bf_rd, c(-1, 10, 15, 22, 36), include.lowest = T)
Dat$cit_rd<-round(Dat$citique,digit=2)
quantile(Dat$cit_rd, probs=seq(0, 1, 0.33))
Dat$cit_rd<-round(Dat$citique,digit=2)
Dat$cit_cat<- cut(Dat$cit_rd, c(0, 0.05, 0.25, 4), include.lowest = T)
Dat$rod <- round(Dat$rodent5,digits = 2)
Dat$rod_cat<- cut(Dat$rod, c(1, 3.27, 5), include.lowest = TRUE)

################## 2020 prediction
###  LB dataset 2020   
Dat20<- read.csv2('Dat20.csv')
Dat20$deer_cat <- cut(Dat20$deer, c(0, 0.4, 0.6, 0.8, 1), include.lowest = TRUE)
Dat20$ndvi_cat<- cut(Dat20$ndvi, c(-1, 0.59, 1), include.lowest = TRUE)
Dat20$sd_rd<- round(Dat20$sd)
Dat20$sd_cat<- cut(Dat20$sd_rd, c(1, 3, 5, 9, 23), include.lowest = TRUE)
Dat20$maxst_bf_rd<- round(Dat20$maxst_bf)
Dat20$maxstbf_cat<- cut(Dat20$maxst_bf_rd, c(-1, 10, 15, 22, 36), include.lowest = T)
Dat20$cit_rd<-round(Dat20$citique,digit=2)
Dat20$cit_cat<- cut(Dat20$cit_rd, c(0, 0.05, 0.25, 4), include.lowest = T)
Dat20$rod <- round(Dat20$rodent5,digits = 2)
Dat20$rod_cat<- cut(Dat20$rod, c(1, 3.27, 5), include.lowest = TRUE)


# Data 2016-2019
# create dataset with positive values only for the continuous part
Datlog <- Dat[Dat$inc>0,]
nb <- length(Dat$inc) # length of binary part  
nc <- length(Datlog$inc) # length of continuous part 
np= 1573 # number of grid cells 
Dat$b <- ifelse(Dat$inc==0,0,1) # zero value indicator (binary part outcome)
Datlog$c <- Datlog$inc # positive values only (continuous part outcome)

# add outcome NA of 2020 for the two parts + link column for prediction link
nd=6292 # predict row for 2020
yy <- matrix(NA, ncol = 3, nrow = nb+nc+nd+nd)
yy[1:(nb+nd),1] <- c(Dat$b,rep(NA, nd)) # binary outcome
yy[nb+nd+(1:(nc+nd)),2] <- c(Datlog$inc,rep(NA, nd)) # continuous outcome
yy[1:(nb+nd),3] <- 1 # logistic link
yy[nb+nd+(1:(nc+nd)),3] <- 2 # gamma link 

yb = yy[,1]
yc = yy[,2] 
link=yy[,3]

#### Add Prediction part in 2016-19 dataset 
# binary outcome for 2020
Dat2<- Dat20
Dat2$b<- NA
Dat<- rbind.data.frame(Dat, Dat2)
# continuous outcome for 2020
Dat22<- Dat20
Dat22$c<- NA
Datlog<- rbind.data.frame(Datlog, Dat22)

# set up unique identifiers for the random-effects
Datlog$idl <- as.integer(Datlog$IDarea1)
Datlog$idl2 <- as.integer(Datlog$IDarea1)
Datlog$IDTime <- as.integer(Datlog$IDtime1)

# fixed effects
linear.covariate <- data.frame(
  InteB = c(rep(1,nb+nd), rep(0,nc+nd)), # intercept (binary part)
  InteC = c(rep(0,nb+nd), rep(1,nc+nd)), # intercept (continuous part)
  ndvi_b = c(as.character(Dat$ndvi_cat), rep(0,nc+nd)),
  rodcatb = c(as.character(Dat$rod_cat), rep(0,nc+nd)),
  deer_c = c(rep(0,nb+nd), as.character(Datlog$deer_cat)),
  citcat_c = c(rep(0,nb+nd), as.character(Datlog$cit_cat)),
  sd_c = c(rep(0,nb+nd), as.character(Datlog$sd_cat)),
  nddcut_c = c(rep(0,nb+nd), Datlog$ndd14),
  maxstbf_c = c(rep(0,nb+nd), as.character(Datlog$maxstbf_cat)))

linear.covariate$ndvi_b<- as.factor(as.character(linear.covariate$ndvi_b))
linear.covariate$rodcatb<- as.factor(as.character(linear.covariate$rodcatb))
linear.covariate$deer_c<- as.factor(as.character(linear.covariate$deer_c))
linear.covariate$sd_c<- as.factor(as.character(linear.covariate$sd_c))
linear.covariate$maxstbf_c<- as.factor(as.character(linear.covariate$maxstbf_c))
linear.covariate$citcat_c<- as.factor(as.character(linear.covariate$citcat_c))
linear.covariate$ndvi_b<- relevel(linear.covariate$ndvi_b, ref='[-1,0.59]')
linear.covariate$rodcatb <- relevel(linear.covariate$rodcatb, ref='[1,3.27]')
linear.covariate$citcat_c<- relevel(linear.covariate$citcat_c, ref='[0,0.05]')
linear.covariate$maxstbf_c <- relevel(linear.covariate$maxstbf_c, ref='[-1,10]')
linear.covariate$deer_c <- relevel(linear.covariate$deer_c, ref='[0,0.4]')
linear.covariate$sd_c <- relevel(linear.covariate$sd_c, ref='[1,3]')

# random-effects
random.covariate<-list(IDl=c(rep(NA,nb+nd),Datlog$IDarea1), # random intercept (continuous)
                       IDb=c(as.integer(Dat$IDarea1),rep(NA,nc+nd)), # random intercept (binary)
                       IDTl=c(rep(NA,nb+nd),Datlog$IDtime1), # random intercept (continuous)
                       IDTb=c(as.integer(Dat$IDtime1),rep(NA,nc+nd))) # random intercept (binary)
jointdf = data.frame(linear.covariate, random.covariate, yb, yc,link)
joint.data <- as.list(inla.rbind.data.frames(jointdf))
Yjoint = cbind(joint.data$yb, joint.data$yc) # outcomes
joint.data$Y <- Yjoint

###
joint.data$ndvi_b_06_1 <- ifelse(joint.data$ndvi_b=="(0.59,1]", 1, 0)
joint.data$rodcatb_3_5 <- ifelse(joint.data$rodcatb=="(3.27,5]", 1, 0)
joint.data$deer_c_4_6 <- ifelse(joint.data$deer_c=="(0.4,0.6]", 1, 0)
joint.data$deer_c_6_8 <- ifelse(joint.data$deer_c=="(0.6,0.8]", 1, 0)
joint.data$deer_c_8_10 <- ifelse(joint.data$deer_c=="(0.8,1]", 1, 0)
joint.data$citcat_c_05_25 <- ifelse(joint.data$citcat_c=="(0.05,0.25]", 1, 0)
joint.data$citcat_c_25_40 <- ifelse(joint.data$citcat_c=="(0.25,4]", 1, 0)
joint.data$sd_c_3_5 <- ifelse(joint.data$sd_c=="(3,5]", 1, 0)
joint.data$sd_c_5_9 <- ifelse(joint.data$sd_c=="(5,9]", 1, 0)
joint.data$sd_c_9_23 <- ifelse(joint.data$sd_c=="(9,23]", 1, 0)
joint.data$maxstbf_c_10_15 <- ifelse(joint.data$maxstbf_c=="(10,15]", 1, 0)
joint.data$maxstbf_c_15_22 <- ifelse(joint.data$maxstbf_c=="(15,22]", 1, 0)
joint.data$maxstbf_c_22_36 <- ifelse(joint.data$maxstbf_c=="(22,36]", 1, 0)

# Define the neighbour structure 
map<- readOGR('map/template_fr.shp')  
nb <- poly2nb(map)
head(nb)
tail(nb) # 1573 cell
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

# INLA formula

formulaJ <- Y ~  -1 + InteB + InteC + ndvi_b_06_1 + rodcatb_3_5 + deer_c_4_6 + deer_c_6_8 + deer_c_8_10 +
  maxstbf_c_10_15 + maxstbf_c_15_22 + maxstbf_c_22_36 + sd_c_3_5 + sd_c_5_9 + sd_c_9_23 +
  citcat_c_05_25 + citcat_c_25_40+
  f(IDb, model = "bym", graph = g, adjust.for.con.comp = TRUE, constr = TRUE) +
  f(IDl,copy="IDb") +
  f(IDTb,  model="seasonal",  season.length=4) +
  f(IDTl, copy="IDTb")

#Fit model with INLA ()
TPinla <- inla(formulaJ, family=c("binomial", "gamma"), 
               data=joint.data,
               Ntrials=c(rep(1,length(Dat$inc)),rep(NA,length(Datlog$inc))),
               control.predictor=list(compute=TRUE, link=link ),
               control.family=list(list(control.link = list(model = "logit")),
                                   list(control.link = list(model = "log"))),
               control.inla=list(strategy="adaptive"),
               control.fixed =list(expand.factor.strategy='inla', remove.names="(Intercept)"), 
               control.compute=list(waic=T,cpo = T, dic = T), 
               verbose=TRUE)

summary(TPinla)
options("scipen"=100, "digits"=2)
exp(TPinla$summary.fixed[, 3:5])

#PIT
First_positive_value <- length(na.omit(joint.data$Y[,1]))
Last_positive_value <- First_positive_value + length(na.omit(joint.data$Y[,2]))
hist(TPinla$cpo$pit[First_positive_value:Last_positive_value], main="", xlab="PIT (positive values of the two-part model)")
