
#Load packages (if packages need to be installed first, do: install.packages("packagename")
require(psych)    
require(nlme)     	
require(lme4)     
require(haven)     
require(lmerTest)
require(plyr)
require(dplyr)

#########################################################
#R-code for Factor analysis on complete data (Section 2)#
#########################################################

#Load data
data.factor <- read_spss("EB82 Multilevel Analysis Full Data.sav")
data.factor <- as.data.frame(data.factor)

items       <- data.factor[,c(5:9)]
covmat      <- cov(items, use="complete.obs")

#Factor analysis extracting 1 factor using Principle Axis Factoring
fit0        <- fa(covmat, nfactors=1, n.obs=1495, covar=TRUE, fm="pa")
fit0

#Remove objects
rm(list = ls())


#################################################
# R-code for analyses performed in Section 3.1.1# 
#################################################

#Models Section 3.1

#Load data
data     <- read_spss("EB82 Multilevel Analysis Section 3.1.1.sav")

#Delete patients with missing data for analysis analogous to BWS ANOVA

data2 <- ddply(data, "id",
      function(df)if(any(is.na(df[, 3]))) NULL else df)

fit1  <- lmer(SexualFunction~StudyPeriod*Treatment +  
             (1|id), REML=FALSE ,data=data2)
summary(fit1)

#67% prediction interval random intercept placebo and T+S patients
lo.p <- fixef(fit1)[1] - sqrt(13.47)
hi.p <- fixef(fit1)[1] + sqrt(13.47)


lo.t <- (fixef(fit1)[1] + fixef(fit1)[3]) - sqrt(13.47)
hi.t <- (fixef(fit1)[1] + fixef(fit1)[3]) + sqrt(13.47)


#ML on aggregated scores (use all observed patients, including the ones with missings)
fit2  <- lmer(SexualFunction~StudyPeriod*Treatment +  
              (1 + StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)

#Model estimation runs into errors as the variance-covariance matrix of the random effects
#is unidentifiable

#Fit the same model, but then constrain the covariance between the random intercept and 
#random slope to zero

fit2a  <- lmer(SexualFunction~StudyPeriod*Treatment +  
              (0 + StudyPeriod|id) + (1|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit2a)

#Remove objects
rm(list = ls())


#################################################
# R-code for analyses performed in Section 3.1.2# 
#################################################

#Read data
data     <- read_spss("EB82 Multilevel Analysis Section 3.1.2 - 3.2 - 3.3.sav")

#ML model on individual events. (1 + StudyPeriod|id) is the random part of the model-> random intercept + random slope for StudyPeriod

fit1  <- lmer(SexualFunction~ StudyPeriod*Treatment+
              (1 + StudyPeriod|id), data=data, REML=FALSE, na.action=na.omit)   
summary(fit1)

#67% predictive interval placebo and T+S patients
lo.p <- fixef(fit1)[2] - sqrt(18.01)
hi.p <- fixef(fit1)[2] + sqrt(18.01)

lo.t <- (fixef(fit1)[2]+fixef(fit1)[4]) - sqrt(18.01)
hi.t <- (fixef(fit1)[2]+fixef(fit1)[4]) + sqrt(18.01)

#Model for age (grand-mean centered), BMI (grand-mean centered) and menopausal status

#Recode Menopausal status to 0-1
data$MenoStat2 <- ifelse(data$MenoStat==1, 0, 1)

#Grand-mean center Age and BMI
data$Age_c <- data$Age - mean(data$Age) 
data$BMI_c <- data$BMI - mean(data$BMI)

fit2  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                Age_c + BMI_c + MenoStat2 +
                (1 +StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit2)

#Deviance difference test
anova(fit1, fit2)

###############################################
# R-code for analyses performed in Section 3.2# 
###############################################

#Linear effect (Include EventCount)
fit3  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                EventCount +
                (1 +StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit3)

#Random slope EventCount significant?
fit3a  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                 EventCount +
                 (1 + StudyPeriod + EventCount |id),
               data=data, REML=FALSE, na.action=na.omit)
summary(fit3a)

#Deviance difference test
anova(fit3, fit3a)

#Quadratic effect (EventCount + I(EventCount^2))
fit4  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                EventCount + I(EventCount^2) + 
                (1 +StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit4)

#Percentage of periods that exceed the 8 events
counts <- data %>% group_by(id,StudyPeriod) %>% summarise(Freq=n())
counts <- as.data.frame(counts)
maxc   <- subset(counts, counts[,3]>8)
perc   <- nrow(maxc)/nrow(counts)

#Random quadratic effect
fit4a  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                 EventCount + I(EventCount^2) + 
                 (1 +StudyPeriod+(EventCount+ I(EventCount^2))|id), 
               data=data, REML=FALSE, na.action=na.omit)
summary(fit4a)
#Model does not converge

#Interactions with study period
fit5  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                (EventCount + I(EventCount^2))*StudyPeriod + 
                (1 +StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit5)

#Deviance difference test model fit5 with model fit4
anova(fit4, fit5)


###############################################
# R-code for analyses performed in Section 3.3# 
###############################################

#Function to calculate person means
#Input:
  #data = data set
  #group_by = the variable that indicates the person or group (id)
  #X = variable of which the means needs to be calculated

person.mean <- function(data, group_by, X){
  id.means  <- ddply(data, group_by, function(x, ind){mean(x[,ind], na.rm=TRUE)}, X)
  id.col    <- which(colnames(data)==group_by)
  n.rows    <- as.vector(table(data[,id.col]))
  y.id.mean <- t(as.vector(0))
  
  for (i in 1:nrow(id.means)){
    v         <- rep(id.means[i,2], n.rows[i])
    y.id.mean <- append(y.id.mean, v)
  }
  
  y.id.mean   <- as.matrix(y.id.mean[-1],,1)
  return(y.id.mean)
}

#Patient means study period
data$StudyPeriod_id.means <- person.mean(data, "id", "StudyPeriod")

#Center the variable around zero
data$StudyPeriod_id.meansC <- data$StudyPeriod_id.means - mean(data$StudyPeriod_id.means)

#Average proportion of ATP events
mean(data$StudyPeriod_id.means)

#Run a model with group means of study period x treatment interaction
fit1a  <- lmer(SexualFunction~ StudyPeriod*Treatment + EventCount + I(EventCount^2) +
                StudyPeriod_id.meansC*Treatment
              +(1 + StudyPeriod|id), data=data, REML=FALSE, na.action=na.omit)   
summary(fit1a)

#Patient means of event count
data$EventCount_id.means <- person.mean(data, "id", "EventCount")

#Center the variable around zero
data$EventCount_id.meansC <- data$EventCount_id.means - mean(data$EventCount_id.means)

#Average of EventCount
mean(data$EventCount_id.means)

#Run a model with group means EventCount

fit4b  <- lmer(SexualFunction~ StudyPeriod*Treatment + EventCount + I(EventCount^2) +
                EventCount_id.meansC + I(EventCount_id.meansC^2) + 
              +(1 + StudyPeriod|id), data=data, REML=FALSE, na.action=na.omit)   
summary(fit4b)

#Remove objects
rm(list = ls())








