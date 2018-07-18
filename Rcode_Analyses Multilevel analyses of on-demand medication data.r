
#Load packages
require(nlme)
require(lme4)
require(haven)
require(lmerTest)
require(plyr)
require(dplyr)
require(DataCombine)

#########################################################
#R-code for Factor analysis on complete data (Section 2)#
#########################################################

#Load data
data.factor <- read_spss("EB82 Multilevel Analysis Section 4.sav")

items       <- data.factor[,c(6:10)]
items       <- as.data.frame(items)
covmat      <- cov(items, use="complete.obs")

fit0        <- factanal(covmat=covmat, factors=1, n.obs=1495)
fit0

#Remove objects
rm(list = ls())


###############################################
# R-code for analyses performed in Section 3.1# 
###############################################

#Models Section 3.1

#Load data
data     <- read_spss("EB82 Multilevel Analysis Section 3.1.sav")

#Delete patients with missing data for analysis analogous to BWS ANOVA

data2 <- ddply(data, "id",
      function(df)if(any(is.na(df[, 3]))) NULL else df)

fit1  <- lme(SexualFunction~StudyPeriod*Treatment, 
             random = ~1|id, data=data2, method="ML")
summary(fit1)

#Variance components
VarCorr(fit1)

#Deviance
AIC(fit1) - (6*2)

#67% prediction interval random intercept placebo and T+S patients
lo.p <- fixef(fit1)[1] - sqrt(getVarCov(fit1)[1])
hi.p <- fixef(fit1)[1] + sqrt(getVarCov(fit1)[1])

lo.t <- (fixef(fit1)[1] + fixef(fit1)[3]) - sqrt(getVarCov(fit1)[1])
hi.t <- (fixef(fit1)[1] + fixef(fit1)[3]) + sqrt(getVarCov(fit1)[1])


#ML on aggregated scores (use all observed patients, including the ones with missings)
fit2  <- lme(SexualFunction~StudyPeriod*Treatment, 
              random = ~StudyPeriod|id, 
              data=data, na.action=na.omit, method="ML")
summary(fit2)

#Variance components
VarCorr(fit2)

#Deviance
AIC(fit2) - (8*2)

#Remove objects
rm(list = ls())


###############################################
# R-code for analyses performed in Section 3.2# 
###############################################

#Read data
data     <- read_spss("EB82 Multilevel Analysis Section 3.2 and 3.3.sav")

#ML model on individual events.  

fit1  <- lmer(SexualFunction~ StudyPeriod*Treatment
              +(1 + StudyPeriod|id), data=data, REML=FALSE, na.action=na.omit)
summary(fit1)

#67% predictive interval placebo and T+S patients
lo.p <- fixef(fit1)[2] - sqrt(18.01)
hi.p <- fixef(fit1)[2] + sqrt(18.01)

lo.t <- (fixef(fit1)[2]+fixef(fit1)[4]) - sqrt(18.01)
hi.t <- (fixef(fit1)[2]+fixef(fit1)[4]) + sqrt(18.01)

#Model for age (grand-mean centered), BMI (grand-mean centered) and menopausal status

#Grand-mean center Age and BMI
data$Age_c <- data$Age - mean(data$Age) 
data$BMI_c <- data$BMI - mean(data$BMI)

fit2  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                Age_c + BMI_c + MenoStat +
                (1 +StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit2)

###############################################
# R-code for analyses performed in Section 3.3# 
###############################################

#Linear effect
fit3  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                EventCount +
                (1 +StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit3)

#Random slope EventCount significant?
fit3a  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                 EventCount +
                 (1 +StudyPeriod|id) + 
                 (0 + EventCount|id),
               data=data, REML=FALSE, na.action=na.omit)
summary(fit3a)

#Likelihood ratio test
anova(fit3, fit3a)

#Quadratic effect
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

#Likelihood ratio test between model with and without quadratic effect
anova(fit3, fit4)

#Random quadratic effect
fit4a  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                 EventCount + I(EventCount^2) + 
                 (1 +StudyPeriod+(EventCount+ I(EventCount^2))|id), 
               data=data, REML=FALSE, na.action=na.omit)
summary(fit4a)
#No convergence

#Interactions with study period
fit5  <- lmer(SexualFunction~ StudyPeriod*Treatment + 
                (EventCount + I(EventCount^2))*StudyPeriod + 
                (1 +StudyPeriod|id), 
              data=data, REML=FALSE, na.action=na.omit)
summary(fit5)

#Likelihood ratio test
anova(fit4, fit5)

#Remove objects
rm(list = ls())


###############################################
# R-code for analyses performed in Section 4.1# 
###############################################

data      <- read_spss("EB82 Multilevel Analysis Section 4.sav")
data      <- as.data.frame(data)

#Function to group-mean center variables.
group.mean <- function(data, group_by, X){
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
#Select only baseline and ATP data
data1 <- subset(data, (StudyPeriod==0|StudyPeriod==2))
data1$StudyPeriod <- ifelse(data1$StudyPeriod==2, 1, 0)

#Assign each first day period a missing value
for (i in 1:nrow(data1)){
  if(data1$EventCount[i] == 0 & data1$EventCount[i]==0){
    data1$DaysInt[i] <- NA
  }
}

#Lag the SexualFunction predictor
data1$id_reg <- paste0(data1$id, data1$StudyPeriod)

data1$gr.mean.Y <- group.mean(data1, "id", "SexualFunction")
data1$GMC.Y     <- data1$SexualFunction - data1$gr.mean.Y

data2 <- slide(data1, Var="GMC.Y", GroupVar="id_reg", NewVar="lagged.gmc.y", slideBy=-1)
data1 <- as.data.frame(data2)

data1$gr.mean.Y_c <- data1$gr.mean.Y - mean(data1$gr.mean.Y)

#Model
fit1<- lmer(DaysInt~ StudyPeriod*Treatment +lagged.gmc.y +
              gr.mean.Y_c +
              (1+StudyPeriod+lagged.gmc.y|id),
            data=data1, REML = FALSE, na.action = na.omit)
summary(fit1)

#Random slope lagged predictor significant?
fit1a<- lmer(DaysInt~ StudyPeriod*Treatment 
             +lagged.gmc.y + gr.mean.Y_c + (1+StudyPeriod|id),
             data=data1, REML = FALSE, na.action = na.omit)

fit1b<- lmer(DaysInt~ StudyPeriod*Treatment 
             +lagged.gmc.y + gr.mean.Y_c + (1+StudyPeriod|id)+
               (0+lagged.gmc.y|id) ,
             data=data1, REML = FALSE, na.action = na.omit)

#Likelihood ratio test
anova(fit1a, fit1b)

#Model with medication intake. 
data2 <- slide(data1, Var="SSEQ2", GroupVar="id_reg", NewVar="lagged.SSEQ2", slideBy=-1)
data1 <- as.data.frame(data2)

fit2 <- lmer(DaysInt~ StudyPeriod*Treatment +lagged.gmc.y*lagged.SSEQ2 +
               gr.mean.Y_c +
               (1+StudyPeriod+lagged.gmc.y|id),
             data=data1, REML = FALSE, na.action = na.omit)
summary(fit2)


###############################################
# R-code for analyses performed in Section 4.2# 
###############################################

#Intercept-only model for Empirical Bayes estimates of the means
fit0 <- lmer(SexualFunction~1 + (1|id), 
             data=data1, REML=FALSE, na.action=na.omit)

#Extract Empirical Bayes Estimates
data1$e.bi[!is.na(data1$SexualFunction)] <- resid(fit0)

#Lag empirical Bayes patient-centered predictor
data2 <- slide(data1, Var="e.bi", GroupVar="id_reg", NewVar="lev1pred.BE", slideBy=-1)
data1 <- as.data.frame(data2)

fit3 <- lmer(SexualFunction~1 +StudyPeriod*Treatment+
               lev1pred.BE + (1+StudyPeriod+lev1pred.BE|id), 
             data=data1, REML=FALSE, na.action=na.omit)
summary(fit3)

#Random Slope significant? 
fit3a <- lmer(SexualFunction~1 +StudyPeriod*Treatment+
                lev1pred.BE + (1+StudyPeriod|id), 
              data=data1, REML=FALSE, na.action=na.omit)

fit3b <- lmer(SexualFunction~1 +StudyPeriod*Treatment+
                lev1pred.BE + (1+StudyPeriod|id) + (0+lev1pred.BE|id), 
              data=data1, REML=FALSE, na.action=na.omit)

#Likelihood ratio test
anova(fit3a, fit3b)

#67% prediction interval
lo <- fixef(fit3)[4]-sqrt(0.02145)
hi <- fixef(fit3)[4]+sqrt(0.02145)

#Include additional variables.
#Interaction with treatment group
fit4 <- lmer(SexualFunction~1 +StudyPeriod*Treatment+
               lev1pred.BE*Treatment + (1+StudyPeriod+lev1pred.BE|id), 
             data=data1, REML=FALSE, na.action=na.omit)
summary(fit4)

#Interaction lagged predictor Days interval, medication intake, Study Period
#Lag 1st level variable medication intake: 
data2 <- slide(data1, Var="SSEQ2", GroupVar="id_reg", NewVar="lagged.SSEQ2", slideBy=-1)
data1 <- as.data.frame(data2)

fit5 <- lmer(SexualFunction~1 +StudyPeriod*Treatment
             +lev1pred.BE*DaysInt +
               lev1pred.BE*lagged.SSEQ2+
               lev1pred.BE*StudyPeriod+
               (1+StudyPeriod+lev1pred.BE|id), data=data1, REML=FALSE, na.action=na.omit)
summary(fit5)

#Remove objects
rm(list = ls())






