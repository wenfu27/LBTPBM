# This the model code for paper 
# 'Spatial and seasonal determinants of Lyme borreliosis incidence in France, 2016 to 2021'
# Upload date: 14 Feb 2023

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

# 1-  Import datasets
# Data from the period 2016-19 used for parameter inference
# Data from the year 2020 and 2021 used separately for model validation

setwd('C:\\Users\\Wen\\')

# IDarea1 = ID for each pixel 
# long= longitude, lat= latitude
# inc= LB seasonal incidence estimated by kriging interpolation
# deer= index of deer presence (0-100%, fixed in time)
# ndvi= vegetation indice (Normalized difference vegetation index)
# maxST= averaged maximum soil temperature at quarter Q
# maxrh= averaged maximum relative humidity at quarter Q
# sd_mm= averaged mean saturation deficit at quarter Q
# ndd14= number of dry days (14-day increments per unit) at quarter Q
# citique = proportion of tick bite reports adjusted by risk-population weight at quarter Q
# rodent5 = species richeness of rodents (fixed in time)
# season = Winter Jan to Mar, Spring Apr to Jun, Summer Jul to Sept, Autumn Oct to Dec
# IDtime1 = 1, 2, 3 to 24 quarters during 2016-21

# import dataset
###  LB my dataset
Dat1619<- read.csv2('Dat1619.csv')
Baz20<- read.csv2('Baz2020.csv')
Baz21<- read.csv2('Baz2021.csv')

### Catergorise the covariates (except ndd14 as continuous variable), by biological relevance or by quadratic methods 
Dat <- Dat1619
Dat$deer_cat <- cut(Dat$deer, c(0, 0.4, 0.6, 0.8, 1), include.lowest = TRUE)
Dat$ndvi_cat<- cut(Dat$ndvi, c(-1,0.59, 1), include.lowest = TRUE)
Dat$sd_cat<- cut(Dat$sd_mm, c(0, 1.5, 3, 5, 13), include.lowest = TRUE)
Dat$maxst_cat<- cut(Dat$maxST, c(-2, 7, 15, 22, 37), include.lowest = T)
Dat$cit_cat<- cut(Dat$citique, c(0, 0.05, 0.25, 4), include.lowest = T)
Dat$rod_cat<- cut(Dat$rodent5, c(1, 3.7, 5), include.lowest = TRUE)

## 2020  + 2021 prediction 
Dat20<- rbind(Dat2020, Dat2021)
Dat20$deer_cat <- cut(Dat20$deer, c(0, 0.4, 0.6, 0.8, 1), include.lowest = TRUE)
Dat20$ndvi_cat<- cut(Dat20$ndvi, c(-1, 0.59, 1), include.lowest = TRUE)
Dat20$sd_cat<- cut(Dat20$sd_mm, c(0, 1.5, 3, 5, 13), include.lowest = TRUE)
Dat20$maxst_cat<- cut(Dat20$maxST, c(-2, 7, 15, 22, 37), include.lowest = T)
Dat20$cit_cat<- cut(Dat20$citique, c(0, 0.05, 0.25, 4), include.lowest = T)
Dat20$rod_cat<- cut(Dat20$rodent5, c(1, 3.7, 5), include.lowest = TRUE)

# import the grid map
# define the spatial structure 
map <- readRDS('map.rsd')
nb <- poly2nb(map)
tail(nb) # 1573 cell
nb2INLA("map.adj", nb)
FR.adj <- paste(getwd(),"/map.adj",sep="")
g <- inla.read.graph(filename = "map.adj")

# Data 2016-2019
# create dataset with positive values only for the continuous part
Datlog <- Dat[Dat$inc>0,]
nb <- length(Dat$inc) # length of binary part  25168 units of analysis 
nc <- length(Datlog$inc) # length of continuous part 20138 units of anlaysis 
np= 1573 # number of pixels
Dat$b <- ifelse(Dat$inc==0,0,1) # zero value indicator (binary part outcome)
Datlog$c <- Datlog$inc # positive values only (continuous part outcome)

# add outcome NA of 2020+2021 for the two parts + a columne of link for prediction link
nd=6292*2 # predict row for 2020+2021
yy <- matrix(NA, ncol = 3, nrow = nb+nc+nd+nd)
yy[1:(nb+nd),1] <- c(Dat$b,rep(NA, nd)) # binary outcome
yy[nb+nd+(1:(nc+nd)),2] <- c(Datlog$inc,rep(NA, nd)) # continuous outcome
yy[1:(nb+nd),3] <- 1 # logistic link
yy[nb+nd+(1:(nc+nd)),3] <- 2 # gamma link 

yb = yy[,1]
yc = yy[,2] 
link=yy[,3]

# add binary outcome for 2020+2021
Dat2<- Dat20
Dat2$b <-NA
Dat<- rbind.data.frame(Dat, Dat2)

# add continuous outcome for 2020+2021
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
  season_b =c(as.character(Dat$season),rep(0,nc+nd)),
  deer_c = c(rep(0,nb+nd), as.character(Datlog$deer_cat)),
  citcat_c = c(rep(0,nb+nd), as.character(Datlog$cit_cat)),
  sd_c = c(rep(0,nb+nd), as.character(Datlog$sd_cat)),
  maxstbf_c = c(rep(0,nb+nd), as.character(Datlog$maxstbf_cat)),
  maxst_c = c(rep(0,nb+nd), as.character(Datlog$maxst_cat)),
  season_c =c(rep(0,nb+nd), as.character(Datlog$season)),
  ndd_c = c(rep(0,nb+nd), Datlog$ndd14))

linear.covariate$ndvi_b<- as.factor(as.character(linear.covariate$ndvi_b))
linear.covariate$rodcatb<- as.factor(as.character(linear.covariate$rodcatb))
linear.covariate$season_b<- as.factor(as.character(linear.covariate$season_b))
linear.covariate$deer_c<- as.factor(as.character(linear.covariate$deer_c))
linear.covariate$sd_c<- as.factor(as.character(linear.covariate$sd_c))
linear.covariate$maxstbf_c<- as.factor(as.character(linear.covariate$maxstbf_c))
linear.covariate$citcat_c<- as.factor(as.character(linear.covariate$citcat_c))
linear.covariate$season_c<- as.factor(as.character(linear.covariate$season_c))
linear.covariate$ndvi_b<- relevel(linear.covariate$ndvi_b, ref='[-1,0.59]')
linear.covariate$rodcatb <- relevel(linear.covariate$rodcatb, ref='[1,3.7]')
linear.covariate$season_b <- relevel(linear.covariate$season_b, ref='winter')
linear.covariate$citcat_c<- relevel(linear.covariate$citcat_c, ref='[0,0.05]')
linear.covariate$maxstbf_c <- relevel(linear.covariate$maxstbf_c, ref='[-2,7]')
linear.covariate$maxst_c <- relevel(factor(linear.covariate$maxst_c), ref='[-2,7]')
linear.covariate$deer_c <- relevel(linear.covariate$deer_c, ref='[0,0.4]')
linear.covariate$sd_c <- relevel(linear.covariate$sd_c, ref='[0,1.5]')


# random-effects correlated spatial random effects and seasonal random effects 
random.covariate<-list(IDbl=c(as.integer(Dat$IDarea1),Datlog$IDarea1),
                       IDgp=c(rep(1, nb+nd), rep(2, nc+nd)), # groups for correlation
                       IDTl=c(rep(NA,nb+nd),Datlog$IDtime1), 
                       IDTb=c(as.integer(Dat$IDtime1),rep(NA,nc+nd)))

jointdf = data.frame(linear.covariate, random.covariate, yb, yc)
joint.data <- as.list(inla.rbind.data.frames(jointdf))
Yjoint = cbind(joint.data$yb, joint.data$yc) # outcomes
joint.data$Y <- Yjoint


# factorise the categorised covariates
joint.data$ndvi_b_06_1 <- ifelse(joint.data$ndvi_b=="(0.59,1]", 1, 0)
joint.data$rodcatb_3_5 <- ifelse(joint.data$rodcatb=="(3.7,5]", 1, 0)
joint.data$season_b_spring <- ifelse(joint.data$season_b=="spring", 1, 0)
joint.data$season_b_summer <- ifelse(joint.data$season_b=="summer", 1, 0)
joint.data$season_b_fall <- ifelse(joint.data$season_b=="fall", 1, 0)
joint.data$deer_c_4_6 <- ifelse(joint.data$deer_c=="(0.4,0.6]", 1, 0)
joint.data$deer_c_6_8 <- ifelse(joint.data$deer_c=="(0.6,0.8]", 1, 0)
joint.data$deer_c_8_10 <- ifelse(joint.data$deer_c=="(0.8,1]", 1, 0)
joint.data$citcat_c_05_25 <- ifelse(joint.data$citcat_c=="(0.05,0.25]", 1, 0)
joint.data$citcat_c_25_40 <- ifelse(joint.data$citcat_c=="(0.25,4]", 1, 0)
joint.data$sd_c_1_3 <- ifelse(joint.data$sd_c=="(1.5,3]", 1, 0)
joint.data$sd_c_3_5 <- ifelse(joint.data$sd_c=="(3,5]", 1, 0)
joint.data$sd_c_5 <- ifelse(joint.data$sd_c=="(5,13]", 1, 0)
joint.data$maxstbf_c_7_15 <- ifelse(joint.data$maxstbf_c=="(7,15]", 1, 0)
joint.data$maxstbf_c_15_22 <- ifelse(joint.data$maxstbf_c=="(15,22]", 1, 0)
joint.data$maxstbf_c_22_37 <- ifelse(joint.data$maxstbf_c=="(22,37]", 1, 0)
joint.data$maxst_c_7_15 <- ifelse(joint.data$maxst_c=="(7,15]", 1, 0)
joint.data$maxst_c_15_22 <- ifelse(joint.data$maxst_c=="(15,22]", 1, 0)
joint.data$maxst_c_22_37 <- ifelse(joint.data$maxstbf_c=="(22,37]", 1, 0)


inla.setOption(inla.mode="experimental") # reduce the run out of memory 

## Model formula

formulaJ <- Y ~   -1 + InteB + InteC + ndvi_b_06_1 + rodcatb_3_5 + deer_c_4_6 + deer_c_6_8 + deer_c_8_10 +
                  maxst_c_7_15 + maxst_c_15_22 + maxst_c_22_37 + sd_c_1_3 + sd_c_3_5 + sd_c_5 + citcat_c_05_25 + citcat_c_25_40 + 
                  f(IDbl, model = "bym", graph = g, group=IDgp, adjust.for.con.comp = TRUE, constr = TRUE) +
                  f(IDTb, model="seasonal", season.length=4) +
                  f(IDTl, model="seasonal", season.length=4)
   
  
##Fit model with INLA ()

TPinla <- inla(formulaJ, family=c("binomial", "gamma"), 
               data=joint.data,
               Ntrials=c(rep(1,length(Dat$inc)),rep(NA,length(Datlog$inc))),
               control.predictor=list(compute=TRUE, link=link ),
               control.family=list(list(control.link = list(model = "logit")),
                                   list(control.link = list(model = "log"))),
               # control.inla=list(int.strategy="eb"), # reduce the thread 
               control.inla=list(strategy="adaptive"),
               control.fixed =list(expand.factor.strategy='inla', remove.names="(Intercept)"), 
               # control.compute=list(config= TRUE), 
               control.compute=list(waic=T,cpo = T, dic = T),
               verbose=TRUE)


summary(TPinla)
options("scipen"=100, "digits"=2)
exp(TPinla$summary.fixed[, 3:5]) 
plot(TPinla)

# plot the random effect 2016-2019
plot(TPinla$summary.random$IDbl$sd, ylab='spatial residual (BYM model)')
plot(TPinla$summary.random$IDTb$mean, ylab='seasonal variation (logistic part)')
plot(TPinla$summary.random$IDTb$sd, ylab='seasonal variation (logistic part)')
plot(TPinla$summary.random$IDTl$mean, ylab='mean of seasonal random effects')
plot(TPinla$summary.random$IDTl$sd, ylab='standard deviation of seasonal random effects')

# plot the prediction 
par(mfrow=c(2,1))
plot(TPinla$summary.fitted.values$mean[37753:70474], ylab = 'incidence')
plot(TPinla$summary.fitted.values$sd[37753:70474], ylab = 'standard deviation')

## Model validation using histogram of PIT values for 2020-2021
ns <- 1000
# 
snb <- inla.posterior.sample(ns, TPinla)
snb[[1]]$hyperpar
# observation 2020-2021
inctt20<- read.csv2('inc2020.csv')
inctt21<- read.csv2('inc2021.csv')
obs20 <- inctt20$kginc
obs20_ga <- which(obs20 >0)
obs21 <- inctt21$var1.pred
obs21_ga <- which(obs21 >0)
nb <- 25168 # length of binary part  25168 units of analysis 
nc <- 20138 # length of continuous part 20138 units of anlaysis 

# add outcome NA of 2020+2021 for the two parts + link col for prediction link
nd=6292*2 # predict row for 2020+2021
#
snb_l <- sapply(snb, function(x){
  rgamma(10632, scale = exp(x$latent[(nb+nd+nc+1):(nb+nc+nd+nd)][c(obs20_ga, obs21_ga)])/x$hyperpar[1],
         shape= x$hyperpar[1])
        } )
pit_ga <- rowMeans(snb_l<c(obs20[obs20_ga],obs21[obs21_ga]))
length(pit_ga)
pit_ga_20 <- pit_ga[1: 5332]
pit_ga_21 <- pit_ga[5333: 10632]
hist(pit_ga_20, xlab= 'PIT (positive value of the two-part model), 2020', main ='', ylim= c(0, 800))
hist(pit_ga_21, xlab= 'PIT (positive value of the two-part model), 2021', main = '', ylim = c(0, 800))
