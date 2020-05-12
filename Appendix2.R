######################################
#        R CODE FOR ANALYSING        #      
#       PROBITSIM DATA               #
#           2018-09-12               #
######################################


#install.packages("sm")
#install.packages("lme4")
#install.packages("Matrix")
#install.packages("Mass")
#install.packages("Hmisc")
#install.packages("Matching")
#install.packages("MASS")

#detach(probitsim)
rm(list=ls()) 
library(foreign)
library(Hmisc)
library(sm)
library(lme4)
library(Matching)
library(MASS) 
library(boot)
library(plyr)
library(cobalt)
library(tableone)
library(survey)
library(WeightIt)
library(ggplot2)


#setwd("C:/Users/ingwa901/Dropbox/Stratos-causality/simulations/CURRENT VERSION SIMULATION + DATA ANALYSIS/R-code for analyses")

probitsim<-as.data.frame(read.dta("PROBITsim2018_v12.dta", convert.factors = FALSE))

#Some transformations

probitsim$cage<-probitsim$age-mean(probitsim$age)
probitsim$cage2<-probitsim$cage^2
probitsim$cwgt0<-probitsim$wgt0-mean(probitsim$wgt0)
probitsim$cwgt02<-probitsim$cwgt0^2
attach(probitsim)

# 2. DATA DESCRIPTION with Hmisc-package
describe(probitsim)
summary(probitsim)

#TABLE 1

tr.a1<-subset(probitsim,probitsim$a1==1)
co.a1<-subset(probitsim,probitsim$a1==0)

#
table(a2,a1)
prop.table(table(tr.a1$a2))

table(tr.a1$a3,tr.a1$a2)
prop.table(table(tr.a1$a3[tr.a1$a2==1]))
prop.table(table(tr.a1$a3[tr.a1$a2==0]))

prop.table(table(smoke))

# TABLE 2
#OVERALL
mean(sex)
mean(wgt0)
sd(wgt0)
mean(wgt3)
sd(wgt3)
mean(age)
sd(age)
#STRATA
#TREATMENT ARM

table(tr.a1$sex)
prop.table(table(tr.a1$sex))

mean(tr.a1$wgt0)
sd(tr.a1$wgt0)

mean(tr.a1$wgt3)
sd(tr.a1$wgt3)

mean(tr.a1$age)
sd(tr.a1$age)

table(tr.a1$smoke)
prop.table(table(tr.a1$smoke))

table(factor(tr.a1$educ))
prop.table(table(factor(tr.a1$educ))

table(tr.a1$allergy)
prop.table(table(tr.a1$allergy))


#CONTROL ARM

table(factor(co.a1$sex))
prop.table(table(factor(co.a1$sex)))

mean(co.a1$wgt3)
sd(co.a1$wgt3)

table(co.a1$smoke)
prop.table(table(co.a1$smoke))

mean(co.a1$age)
sd(co.a1$age)

table(factor(co.a1$educ))
prop.table(table(factor(co.a1$educ)))

table(factor(co.a1$allergy))
prop.table(table(factor(co.a1$allergy)))


#boys
boys<-subset(probitsim,probitsim$sex==1)
mean(boys$wgt0)
#girls
girls<-subset(probitsim,probitsim$sex==0)
mean(girls$wgt0)
#low education
low.educ<-subset(probitsim,probitsim$educ==1)
mean(low.educ$wgt0)
#smoking
smoking<-subset(probitsim,probitsim$smoke==1)
mean(smoking$wgt0)
#control arm
co<-subset(probitsim,probitsim$a1==0)
#intervention arm
int<-subset(probitsim,probitsim$a1==1)

#EFFECT OF A1
#
mod.ate.a1<-lm(wgt3~a1,data=probitsim)
summary(mod.ate.a1)
ate.a1<-coef(mod.ate.a1)[2]
se.a1<-9.074

ate.a1-1.96*se.a1
ate.a1+1.96*se.a1



# 3. NUC ESTIMATION APPROACHES, ALL BOOTSTRAP FUNCTIONS IN THE END OF THE SCRIPT 

# TO OBTAIN ESTIMATES IN TABLE 4 TREATMENT EFFECTS OF A2
#GÖR OM TABLE 2....
# First fit a propenstiy score model 

#ps.mod.a2<-glm(a2~cage+factor(educ)+factor(smoke),family=binomial,data=probitsim)

ps.mod.a2<-glm(a2~factor(location)+factor(educ)+cage+cage2+factor(smoke)
+factor(allergy),family=binomial,data=probitsim)

summary(ps.mod.a2)

probitsim$ps.a2<-fitted.values(ps.mod.a2)
 
attach(probitsim)

# Balance checks for the ps-model 
weight<-a2/ps.a2+(1-a2)/(1-ps.a2)

mean(a2*weight*age)
mean((1-a2)*weight*age)
mean(age)

mean(a2*weight*smoke)
mean((1-a2)*weight*smoke)
mean(smoke)

mean(a2*weight*allergy)
mean((1-a2)*weight*allergy)
mean(allergy)


#          


# Overlap

sm.density.compare(probitsim$ps.a2, probitsim$a2, xlab="Propensity score")
title(main="Comparing PS distributions of treated and controls, A2")

# Crude regression

mod.crude.a2<-lm(wgt3~a2)
ate.crude.a2<-coef(mod.crude.a2)[2]
se.crude.a2<-summary(mod.crude.a2)[[4]][4]

ate.crude.a2
se.crude.a2

#95%confidence interval
ate.crude.a2-1.96*se.crude.a2
ate.crude.a2+1.96*se.crude.a2

# Regression adjustment (without interactions) 

adj.a2<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location))
summary(adj.a2)
ate.adj.a2<-coef(adj.a2)[2]
se.adj.a2<-summary(adj.a2)[[4]][2,2]
ate.adj.a2
se.adj.a2

#95%confidence interval
ate.adj.a2-1.96*se.adj.a2
ate.adj.a2+1.96*se.adj.a2

# Regression adjustment (with interaction), separate regressions for the a2 groups
treated.a2<-subset(probitsim,probitsim$a2==1)
controls.a2<-subset(probitsim,probitsim$a2==0)

control.mod.a2<-lm(wgt3~factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location),data=controls.a2)
probitsim$mY0 <- predict(control.mod.a2,probitsim,type="response")

treated.mod.a2<-lm(wgt3~factor(educ)+smoke+factor(allergy)+cage+cage2+factor(location),data=treated.a2)
probitsim$mY1 <- predict(treated.mod.a2,probitsim,type="response")

attach(probitsim)

ate.reg<-mean(mY1-mY0)
ate.reg

# Regression with interactions ATT

mY0.treated <- predict(control.mod.a2,newdata=treated.a2,type="response")

treated.mod.a2<-lm(wgt3~factor(educ)+smoke+factor(allergy)+cage+cage2+factor(location),data=treated.a2)
mY1.treated <- predict(treated.mod.a2,type="response")

att.reg<-mean(mY1.treated)-mean(mY0.treated)


# Regression with PS

reg.PS<-lm(wgt3~a2+ps.a2,data=probitsim)
reg.PS

ate.reg.PS<-coef(reg.PS)[2]
se.reg.PS<-summary(reg.PS)[[4]][2,2]
ate.reg.PS
se.reg.PS


#95%confidence interval, ignoring estimation of PS, see bootstrap SE below
ate.reg.PS-1.96*se.reg.PS
ate.reg.PS+1.96*se.reg.PS

# PS stratification

#ATE

probitsim$strata<-cut(probitsim$ps.a2,quantile(probitsim$ps,probs = 
 seq(0, 1, 1/6)))
fits <- lmList(wgt3~ a2 | strata, data=probitsim)

strata.effect<-mean(c(coef(fits)[2])$a2)
strata.effect

#ATT

probitsim$strata.att<-cut(probitsim$ps.a2,quantile(treated.a2$ps,probs = 
 seq(0, 1, 1/6)))
fits.att <- lmList(wgt3~ a2 | strata.att, data=probitsim)

strata.effect<-mean(c(coef(fits.att)[2])$a2)
strata.effect

# PS matching (1 match)

match.ate<-Match(Y=probitsim$wgt3,Tr=probitsim$a2,X=probitsim$ps.a2,M=1,estimand="ATE",caliper=0.1)
match.ate$est 
match.ate$se

love.plot(match.ate, formula= a2~factor(location)+factor(educ)+cage+cage2+factor(smoke)+
factor(allergy),family=binomial,data=probitsim,order="alphabetical",grid=TRUE)+
ggtitle("Matching balance")

# In the paper 

match.att<-Match(Y=probitsim$wgt3,Tr=probitsim$a2,X=probitsim$ps.a2,M=1,estimand="ATT",caliper=0.1)
match.att$est 
match.att$se



#95%confidence interval
match.ate$est-1.96*match.ate$se
match.ate$est+1.96*match.ate$se
# balance table

# PS matching (3 matches)

match.ate<-Match(Y=probitsim$wgt3,Tr=probitsim$a2,X=probitsim$ps.a2,M=3,estimand="ATE",caliper=0.1)
match.ate$est 
match.ate$se


match.ate<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATE",caliper=0.1)
match.ate$est 
match.ate$se

match.att<-Match(Y=probitsim$wgt3,Tr=probitsim$a2,X=probitsim$ps.a2,M=3,estimand="ATT",caliper=0.1)
match.att$est 
match.att$se



#95%confidence interval

match.ate$est-1.96*match.ate$se
match.ate$est+1.96*match.ate$se

# PS IPW 

ipw.ate<-1/sum(a2/ps.a2)*sum(a2*wgt3/ps.a2)-1/sum((1-a2)/(1-ps.a2))*sum((1-a2)*wgt3/(1-ps.a2))
ipw.ate

w.out1 <- weightit(a2 ~ factor(location)+factor(educ)+cage+cage2+factor(smoke)
+factor(allergy),family=binomial,data=probitsim,estimand = "ATE", method = "ps")

set.cobalt.options(binary = "std")
love.plot(w.out1)


ipw.att<-mean(treated.a2$wgt3)-sum((1-a2)*wgt3*(ps.a2/(1-ps.a2)))/sum((1-a2)*(ps.a2/(1-ps.a2)))
ipw.att


ipw.att<-mean(treated.a2$wgt3)-1/nrow(treated.a2)*sum((1-a2)*wgt3*(ps.a2)/(1-ps.a2)+(a2-ps.a2)*mY0*(ps.a2)/(1-ps.a2))
ipw.att



# PS DR IPW 

dr.ate<-(1/nrow(probitsim))*sum((a2*wgt3-(a2-ps.a2)*mY1)/ps.a2-((1-a2)*wgt3+(a2-ps.a2)*mY0)/(1-ps.a2))

dr.att<-mean(treated.a2$wgt3)-1/nrow(treated.a2)*sum((1-a2)*wgt3*(ps.a2)/(1-ps.a2)+(a2-ps.a2)*mY0*(ps.a2)/(1-ps.a2))
dr.att


#############################################################################
# TO OBTAIN ESTIMATES IN TABLE 3, TREATMENT EFFECTS OF A3 for A1=0 AND A1=1 #
#############################################################################
#probitsim<-probitsim[a1==0,]
probitsim<-probitsim[a1==1,]
attach(probitsim)
nrow(probitsim)
# 8377 individuals A1=0
# 8667 individuals A1=1
# first fit a propenstiy score model and evaluate overlap assumption


#ps.mod.a3<-glm(a3~factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+
factor(sex)+cwgt0+cwgt02,family=binomial,data=probitsim)

ps.mod.a3<-glm(a3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+
factor(sex)+cwgt0+cwgt02,family=binomial,data=probitsim)

#note that a2 is excluded for the stratum A1=0.

summary(ps.mod.a3)

probitsim$ps.a3<-fitted.values(ps.mod.a3)
attach(probitsim)

# balance checks for the ps-model

weight<-a3/ps.a3+(1-a3)/(1-ps.a3)

mean(a3*weight*age)
mean((1-a3)*weight*age)
mean(age)

mean(a3*weight*smoke)
mean((1-a3)*weight*smoke)
mean(smoke)

mean(a3*weight*allergy)
mean((1-a3)*weight*allergy)
mean(allergy)

mean(a3*weight*cesarean)
mean((1-a3)*weight*cesarean)
mean(cesarean)

mean(a3*weight*sex)
mean((1-a3)*weight*sex)
mean(sex)

mean(a3*weight*wgt0)
mean((1-a3)*weight*wgt0)
mean(wgt0)

vars1 <- c("a2","educ","smoke","allergy","cage","cage2","location","cesarean","sex","cwgt0","cwgt02")
probitsim$weight<-(a3/ps.a3)+(1-a3)/(1-ps.a3)
attach(probitsim)

#UNADJUSTED
IPW <- svydesign(ids=~0, data=probitsim, weights=weight)
IPW.1 <- CreateTableOne(vars = vars1, strata = "a3",
                          data = probitsim, test = FALSE)
print(IPW.1, smd = TRUE)

tabIPW <- svyCreateTableOne(vars = vars1,strata = "a3",
                             data = IPW.1, test = FALSE)

obj<-round(cbind(ExtractSmd(IPW.1),ExtractSmd(tabIPW)),3)
obj2<-sort(obj[,1], decreasing = TRUE)

plot(obj2,obj[,1],ylim=c(0,0.3),pch=1)
lines(obj2,obj[,1])



# Overlap

sm.density.compare(probitsim$ps.a3, probitsim$a3, xlab="Propensity score")
title(main="Comparing PS distributions of treated and controls, A3 for A1=0")

# Crude regression

crude.a3<-lm(wgt3~a3)
summary(crude.a3)

ate.crude.a3<-coef(crude.a3)[2]
se.crude.a3<-summary(crude.a3)[[4]][2,2]

ate.crude.a3
se.crude.a3

#95% confidence interval

ate.crude.a3-1.96*se.crude.a3
ate.crude.a3+1.96*se.crude.a3

# Regression adjustment (simple) 
#ATE and ATT
adj.a3<-lm(wgt3~a3+a2+factor(educ)+smoke+allergy+cage+cage2+factor(location)+cesarean+sex+cwgt0+cwgt02)
summary(adj.a3)

#remove a2 for strata a1=0

#95% confidence interval

ate.adj.a3<-coef(adj.a3)[2]
se.adj.a3<-summary(adj.a3)[[4]][2,2]
ate.adj.a3
se.adj.a3

#95% confidence interval

ate.adj.a3-1.96*se.adj.a3
ate.adj.a3+1.96*se.adj.a3

# Regression adjustment (with interaction)
# separate regressions for the a3 groups
#ATE

treated.a3<-subset(probitsim,probitsim$a3==1)
controls.a3<-subset(probitsim,probitsim$a3==0)

control.mod.a3<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+factor(sex)+cwgt0+cwgt02,data=controls.a3)
#remove a2 for strata a1=0
probitsim$mY0 <- predict(control.mod.a3,probitsim,type="response")

treated.mod.a3<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+factor(sex)+cwgt0+cwgt02,data=treated.a3)
#remove a2 for strata a1=0
probitsim$mY1 <- predict(treated.mod.a3,probitsim,type="response")

attach(probitsim)

mean(mY1)-mean(mY0)

#ATT

mY0.treated <- predict(control.mod.a3,treated.a3,type="response")
mY1.treated <- predict(treated.mod.a3,treated.a3,type="response")

mean(mY1.treated)-mean(mY0.treated)

# PS stratification
#ATE

probitsim$strata<-cut(probitsim$ps.a3,quantile(probitsim$ps.a3,probs = 
 seq(0, 1, 1/6)))
#checking up on na

summary(probitsim$strata)
number<-which(is.na(probitsim$strata))

#probitsim$strata[number]<-"(0.0909,0.335]" 
probitsim$strata[number,]<-"(0.0794,0.413]" 

summary(probitsim$strata)
fits <- lmList(wgt3~ a3 | strata, data=probitsim)

strata.effect<-mean(c(coef(fits)[2])$a3)
strata.effect

#ATT

probitsim$strata.treated<-cut(probitsim$ps.a3,quantile(treated.a3$ps.a3,probs = 
 seq(0, 1, 1/6)))
fits.att <- lmList(wgt3~ a3 | strata.treated, data=probitsim)

strata.effect.att<-mean(c(coef(fits.att)[2])$a3)
strata.effect.att


# Regression with PS

reg.PS<-lm(wgt3~a3+ps.a3,data=probitsim)
reg.PS

ate.reg.PS<-coef(reg.PS)[2]
se.reg.PS<-summary(reg.PS)[[4]][2,2]
ate.reg.PS
se.reg.PS


#95%confidence interval, ignoring estimation of PS, see bootstrap SE below
ate.reg.PS-1.96*se.reg.PS
ate.reg.PS+1.96*se.reg.PS

# PS matching (1 match)

match.ate<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATE",caliper=0.1)
match.ate$est 
match.ate$se


#BALANCE CHECK 
#some transformations
probitsim$educ<-as.factor(probitsim$educ)
probitsim$smoke<-as.factor(probitsim$smoke)
probitsim$allergy<-as.factor(probitsim$allergy)
probitsim$location<-as.factor(probitsim$location)
probitsim$cesarean<-as.factor(probitsim$cesarean)
probitsim$sex<-as.factor(probitsim$sex)


#in the paper
love.plot(match.ate, formula=a3~educ+smoke+allergy+cage+cage2+location+cesarean+
sex+cwgt0+cwgt02, data=probitsim,order="alphabetical",grid=TRUE)+ggtitle("Matching balance, A1=0")

#in the paper

love.plot(match.ate, formula=a3~a2+educ+smoke+allergy+cage+cage2+location+cesarean+
sex+cwgt0+cwgt02, data=probitsim,order="alphabetical",grid=TRUE)+ggtitle("Matching balance, A1=1")

#or alternatively by dropping caliper

match.ate.raw<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATE")
match.ate.raw$est 
match.ate.raw$se

match.att<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATT",caliper=0.1)
match.att$est 
match.att$se


# PS matching (3 matches)

match.ate<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=3,estimand="ATE", caliper=0.1)
match.ate$est 
match.ate$se


match.att<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=3,estimand="ATT",caliper=0.1)
match.att$est 
match.att$se


# PS IPW

ipw.ate<-(1/nrow(probitsim))*(sum(a3*wgt3/ps.a3-(1-a3)*wgt3/(1-ps.a3)))
ipw.att<-(1/sum(a3))*(sum(a3*wgt3-(1-a3)*wgt3*ps.a3/(1-ps.a3)))

# PS DR IPW 

dr.ate<-(1/nrow(probitsim))*(sum((a3*wgt3)/ps.a3-(a3-ps.a3)*mY1)/ps.a3-(1-a3)*wgt3/(1-ps.a3)+(a3-ps.a3)*mY0/(1-ps.a3)))
dr.att<-(1/nrow(probitsim))*(sum((a3*wgt3)-(1-a3)*wgt3*ps.a3/(1-ps.a3)+(a3-ps.a3)*mY0*ps.a3/(1-ps.a3)))


#################################################################################
# TO OBTAIN ESTIMATES IN TABLE 5, TREATMENT EFFECTS OF A3 for A1=0 AND EDUC=LOW #
#################################################################################

probitsim<-subset(probitsim,probitsim$educ==1)
attach(probitsim)

sm.density.compare(probitsim$ps.a3, probitsim$a3, xlab="Propensity score")
title(main="Comparing PS distributions of treated and controls for low educated mothers, A3")

MatchBalance
# PS matching (1 match)

match.ate<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATE")
match.ate$est 
match.ate$se

#or alternatively by dropping caliper

match.ate.raw<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATE")
match.ate.raw$est 
match.ate.raw$se

match.att<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATT",caliper=0.1)
match.att$est 
match.att$se


# PS matching (3 matches)

match.ate<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=3,estimand="ATE", caliper=0.1)
match.ate$est 
match.ate$se

summary(match.ate)


match.att<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=3,estimand="ATT",caliper=0.1)
match.att$est 
match.att$se

#Checking balance

balance<-MatchBalance(a3~a2+educ+smoke+allergy+cage+cage2+location+cesarean+sex+cwgt0+cwgt02,match.out=match.ate,data=probitsim)
balance
#recover matched data

matched <- rbind(probitsim[match.ate$index.treated,], probitsim[match.ate$index.control,])
nrow(matched)

#without ties

match.ate.noties<-Match(Y=probitsim$wgt3,Tr=probitsim$a3,X=probitsim$ps.a3,M=1,estimand="ATE",ties=FALSE)
match.ate.noties$est
match.ate.noties$se

matched <- rbind(probitsim[match.ate$index.treated,], probitsim[match.ate$index.control,])
nrow(matched)

matched.noties <- data.frame(match.ate.noties$index.treated,match.ate.noties$index.control)
matched.noties


#########################
#BOOTSTRAP FUNCTIONS A2 #
#########################

#regression adjustment with interactions
fc.reg <- function(datat,i){
d2 <- datat[i,]
	
treated.a2<-subset(d2,d2$a2==1)
controls.a2<-subset(d2,d2$a2==0)

control.mod.a2<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location),data=controls.a2)
mY0 <- predict(control.mod.a2,d2,type="response")

treated.mod.a2<-lm(wgt3~a2+factor(educ)+factor(smoke)+allergy+cage+cage2+factor(location),data=treated.a2)
mY1 <- predict(treated.mod.a2,d2,type="response")

ate.reg<-mean(mY1-mY0)

	return(ate.reg)
}

set.seed(626)

boot.reg<- boot(probitsim, fc.reg, R=1000)
boot.reg

boot.ci(boot.out = boot.reg, type = "norm")

#regression adjustment with interactions
fc.reg.att <- function(datat,i){
d2 <- datat[i,]
	
treated.a2<-subset(d2,d2$a2==1)
controls.a2<-subset(d2,d2$a2==0)

control.mod.a2<-lm(wgt3~factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location),data=controls.a2)
mY0.treated <- predict(control.mod.a2,treated.a2,type="response")

treated.mod.a2<-lm(wgt3~factor(educ)+factor(smoke)+allergy+cage+cage2+factor(location),data=treated.a2)
mY1.treated <- predict(treated.mod.a2,type="response")

att.reg<-mean(mY1.treated-mY0.treated)

	return(att.reg)
}

set.seed(626)

boot.reg<- boot(probitsim, fc.reg.att, R=1000)
boot.reg

boot.ci(boot.out = boot.reg, type = "norm")



#PS regression

fc.psreg <- function(datat,i){
d2 <- datat[i,]

ps.mod.a2<-glm(a2~factor(location)+factor(educ)+cage+cage2+factor(smoke)
+factor(allergy),family=binomial,data=d2)

boot.ps.a2<-fitted.values(ps.mod.a2)
reg.PS<-lm(wgt3~a2+boot.ps.a2,data=d2)
reg.PS

ate.psreg<-coef(reg.PS)[2]

	return(ate.psreg)
}

set.seed(6791)

boot.psreg<- boot(probitsim, fc.psreg, R=1000)
boot.psreg

boot.ci(boot.out = boot.reg, type = "norm")


# PS-stratification six strata, ATE

fc.strat <- function(datat,i){
d2 <- datat[i,]

d2$strata<-cut(d2$ps.a2,quantile(d2$ps.a2,probs = 
 seq(0, 1, 1/6)))
fits <- lmList(wgt3~ a2 | strata, data=d2)

strata.effect<-mean(c(coef(fits)[2])$a2)
return(strata.effect)
}

set.seed(321)

boot.psstrat<- boot(probitsim, fc.strat, R=1000)
boot.psstrat
boot.ci(boot.out = boot.psstrat, type = "norm")

# PS-stratification six strata, ATT

fc.strat.att <- function(datat,i){
d2 <- datat[i,]
treated<-subset(d2,d2$a2==1)

d2$strata<-cut(d2$ps.a2,quantile(treated$ps.a2,probs = 
 seq(0, 1, 1/6)))
fits <- lmList(wgt3~ a2 | strata, data=d2)

strata.effect<-mean(c(coef(fits)[2])$a2)
return(strata.effect)
}

set.seed(321)

boot.psstrat<- boot(probitsim, fc.strat.att, R=1000)
boot.psstrat
boot.ci(boot.out = boot.psstrat, type = "norm")



# IPW
#ATE
fc.ipw <- function(datat,i){
d2 <- datat[i,]

ipw.ate<-1/sum(d2$a2/d2$ps.a2)*sum(d2$a2*d2$wgt3/d2$ps.a2)-1/sum((1-d2$a2)/(1-d2$ps.a2))*sum((1-d2$a2)*d2$wgt3/(1-d2$ps.a2))
return(ipw.ate)
}

set.seed(867)

boot.ipw<- boot(probitsim, fc.ipw, R=1000)
boot.ipw
boot.ci(boot.out = boot.ipw, type = "norm")

#ATT
fc.ipw.att<- function(datat,i){
d2 <- datat[i,]

treated.a2<-subset(d2,d2$a2==1)
ipw.att<-mean(treated.a2$wgt3)-sum((1-d2$a2)*d2$wgt3*(d2$ps.a2/(1-d2$ps.a2)))/sum((1-d2$a2)*(d2$ps.a2/(1-d2$ps.a2)))

return(ipw.att)
}

set.seed(867)

boot.ipw<- boot(probitsim, fc.ipw.att, R=1000)
boot.ipw
boot.ci(boot.out = boot.ipw, type = "norm")

# IPW DR
#ATE
fc.dr <- function(datat,i){
d2 <- datat[i,]

dr.ate<-(1/nrow(d2))*sum((d2$a2*d2$wgt3-(d2$a2-d2$ps.a2)*d2$mY1)/d2$ps.a2-((1-d2$a2)*d2$wgt3+(d2$a2-d2$ps.a2)*d2$mY0)/(1-d2$ps.a2))


return(dr.ate)
}

set.seed(534)

boot.dr<- boot(probitsim, fc.dr, R=1000)
boot.dr
boot.ci(boot.out = boot.dr, type = "norm")

#ATT

fc.dr.att <- function(datat,i){
d2 <- datat[i,]
treated.a2<-subset(d2,d2$a2==1)
dr.att<-mean(treated.a2$wgt3)-1/nrow(treated.a2)*sum((1-d2$a2)*d2$wgt3*(d2$ps.a2)/(1-d2$ps.a2)+(d2$a2-d2$ps.a2)*d2$mY0*(d2$ps.a2)/(1-d2$ps.a2))
dr.att

return(dr.att)
}

set.seed(534)

boot.dr<- boot(probitsim, fc.dr.att, R=1000)
boot.dr
boot.ci(boot.out = boot.dr, type = "norm")


###########################################
#BOOTSTRAP FUNCTIONS A3 for A1=0 and A1=1 #
###########################################

#regression adjustment with interactions
#checking data should be 8377 individuals for A1=0 and 8667 for A1=1

nrow(probitsim)

#ATE

fc.reg.a3 <- function(datat,i){
d2 <- datat[i,]
	
treated.a3<-subset(d2,d2$a3==1)
controls.a3<-subset(d2,d2$a3==0)

control.mod.a3<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+factor(sex)+cwgt0+cwgt02,data=controls.a3)

mY0 <- predict(control.mod.a3,d2,type="response")

treated.mod.a3<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+factor(sex)+cwgt0+cwgt02,data=treated.a3)

mY1 <- predict(treated.mod.a3,d2,type="response")

ate.reg<-mean(mY1-mY0)

return(ate.reg)
}

set.seed(754)

boot.reg<- boot(probitsim, fc.reg.a3, R=1000)
boot.reg

boot.ci(boot.out = boot.reg, type = "norm")


#ATT

fc.regatt.a3 <- function(datat,i){
d2 <- datat[i,]
	
treated.a3<-subset(d2,d2$a3==1)
controls.a3<-subset(d2,d2$a3==0)

control.mod.a3<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+factor(sex)+cwgt0+cwgt02,data=controls.a3)
mY0.treated <- predict(control.mod.a3,treated.a3,type="response")

treated.mod.a3<-lm(wgt3~a2+factor(educ)+factor(smoke)+factor(allergy)+cage+cage2+factor(location)+factor(cesarean)+factor(sex)+cwgt0+cwgt02,data=treated.a3)
mY1.treated <- predict(treated.mod.a3,type="response")

att.reg<-mean(mY1.treated-mY0.treated)

return(att.reg)
}

set.seed(9064)

boot.reg<- boot(probitsim, fc.regatt.a3, R=1000)
boot.reg

boot.ci(boot.out = boot.reg, type = "norm")


#PS regression

fc.psreg <- function(datat,i){
d2 <- datat[i,]
	
reg.PS<-lm(wgt3~a3+ps.a3,data=d2)
reg.PS

ate.psreg<-coef(reg.PS)[2]

	return(ate.psreg)
}

set.seed(61)

boot.psreg<- boot(probitsim, fc.psreg, R=1000)
boot.psreg

boot.ci(boot.out = boot.psreg, type = "norm")

# PS-stratification six strata

fc.strat <- function(datat,i){
d2 <- datat[i,]

d2$strata<-cut(d2$ps.a3,quantile(d2$ps.a3,probs = 
 seq(0, 1, 1/6)))
fits <- lmList(wgt3~ a3 | strata, data=d2)

strata.effect<-mean(c(coef(fits)[2])$a3)
return(strata.effect)
}

set.seed(79283)

boot.psstrat<- boot(probitsim, fc.strat, R=1000)
boot.psstrat
boot.ci(boot.out = boot.psstrat, type = "norm")

# PS-stratification six strata, ATT

fc.strat.att <- function(datat,i){
d2 <- datat[i,]
treated<-subset(d2,d2$a3==1)

d2$strata<-cut(d2$ps.a3,quantile(treated$ps.a3,probs = 
 seq(0, 1, 1/6)))
fits <- lmList(wgt3~ a3 | strata, data=d2)

strata.effect<-mean(c(coef(fits)[2])$a3)
return(strata.effect)
}

set.seed(6059)

boot.psstrat<- boot(probitsim, fc.strat.att, R=1000)
boot.psstrat
boot.ci(boot.out = boot.psstrat, type = "norm")


# IPW
#ATE
fc.ipw <- function(datat,i){
d2 <- datat[i,]

ipw.ate<-1/sum(d2$a3/d2$ps.a3)*sum(d2$a3*d2$wgt3/d2$ps.a3)-1/sum((1-d2$a3)/(1-d2$ps.a3))*sum((1-d2$a3)*d2$wgt3/(1-d2$ps.a3))
return(ipw.ate)
}

set.seed(867)

boot.ipw<- boot(probitsim, fc.ipw, R=1000)
boot.ipw
boot.ci(boot.out = boot.ipw, type = "norm")

#ATT

fc.ipw.att<- function(datat,i){
d2 <- datat[i,]

treated.a2<-subset(d2,d2$a3==1)
ipw.att<-mean(treated.a3$wgt3)-sum((1-d2$a3)*d2$wgt3*(d2$ps.a3/(1-d2$ps.a3)))/sum((1-d2$a3)*(d2$ps.a3/(1-d2$ps.a3)))

return(ipw.att)
}

set.seed(867)

boot.ipw<- boot(probitsim, fc.ipw.att, R=1000)
boot.ipw
boot.ci(boot.out = boot.ipw, type = "norm")

# IPW DR
#ATE
fc.dr <- function(datat,i){
d2 <- datat[i,]

dr.ate<-(1/nrow(d2))*sum((d2$a3*d2$wgt3-(d2$a3-d2$ps.a3)*d2$mY1)/d2$ps.a3-((1-d2$a3)*d2$wgt3+(d2$a3-d2$ps.a3)*d2$mY0)/(1-d2$ps.a3))

return(dr.ate)
}

set.seed(534)

boot.dr<- boot(probitsim, fc.dr, R=1000)
boot.dr
boot.ci(boot.out = boot.dr, type = "norm")


#ATT

fc.dr.att <- function(datat,i){
d2 <- datat[i,]
treated.a3<-subset(d2,d2$a3==1)
dr.att<-mean(treated.a3$wgt3)-1/nrow(treated.a3)*sum((1-d2$a3)*d2$wgt3*(d2$ps.a3)/(1-d2$ps.a3)+(d2$a3-d2$ps.a3)*d2$mY0*(d2$ps.a3)/(1-d2$ps.a3))
dr.att

return(dr.att)
}

set.seed(544)

boot.dr<- boot(probitsim, fc.dr.att, R=1000)
boot.dr
boot.ci(boot.out = boot.dr, type = "norm")



###########################################
#BOOTSTRAP FUNCTIONS A3 for A1=0 and A1=1 #
###########################################


# IPW

#ATE
fc.ipw <- function(datat,i){
d2 <- datat[i,]

ipw.ate<-1/sum(d2$a3/d2$ps.a3)*sum(d2$a3*d2$wgt3/d2$ps.a3)-1/sum((1-d2$a3)/(1-d2$ps.a3))*sum((1-d2$a3)*d2$wgt3/(1-d2$ps.a3))
return(ipw.ate)
}

set.seed(867)

boot.ipw<- boot(probitsim, fc.ipw, R=1000)
boot.ipw
boot.ci(boot.out = boot.ipw, type = "norm")












