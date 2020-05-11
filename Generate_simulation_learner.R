setwd("C://Users//slecessie//Dropbox//Stratos-causality//simulations//CURRENT VERSION SIMULATION + DATA ANALYSIS//data")
#########################################################################
#			    				                                                      #
#  This code produces the dataset "PROBITsim_april2018.csv  ,           #
#   which is the simulated dataset inspired by the  probit trial        #
#                                                                       #
########################################################################

##########################
#                        #
#DATA GENERATING PROCESS #
#                        #
##########################

# We will define a set of covariates and four different treatments A1-A4

rm(list=ls()) 
inv.logit <- function(x) { exp(x)/(1+exp(x)) }
set.seed(1)

n <- 17044

#################
#               #         
#COVARIATES     #
#               #  
#################

#Four locations are generated
Location <- sample(c(1,2,3,4),n,replace=T,prob=c(0.33,0.16,0.26,0.25))
Location <- as.factor(Location)

# Age lognormal 
Age <- round(exp(rnorm(n,3.17,0.19)))
# Any Age < 13 set to 13
Age[Age < 13] <- 13

# Sex: boys coded with 1
Sex <- rbinom(n,1,0.52)

# Educ a factor with three levels, lower = less educated, depends on Location,
# The education is lower in location 1 and 3, rural areas
Educ <- sample(c(1,2,3),n,replace=T,prob=c(0.31,0.54,0.15)) 
Educ[Location == 2] <- sample(c(1,2,3),sum(Location == 2),replace=T,prob=c(0.44,0.44,0.12))
Educ[Location == 3] <- sample(c(1,2,3),sum(Location == 3),replace=T,prob=c(0.33,0.51,0.16))
Educ[Location == 4] <- sample(c(1,2,3),sum(Location == 4),replace=T,prob=c(0.41,0.48,0.11))
Educ <- as.factor(Educ)

# Smoking, Allergy, Cesarean modelled on Educ

# Smoke split: 0.60, 0.80, 0.95 non-smokers
Smoke <- rbinom(n,1,as.numeric(Educ == 1)*0.40 + as.numeric(Educ == 2)*0.25 + as.numeric(Educ == 3)*0.10)
# Allergy split: 0.97, 0.95, 0.93 non-allergy
Allergy <- rbinom(n,1,as.numeric(Educ == 1)*0.03 + as.numeric(Educ == 2)*0.05 + as.numeric(Educ == 3)*0.07)
# Cesarean: 0.90, 0.88, 0.84 not
Cesarean <- rbinom(n,1,as.numeric(Educ == 1)*0.1 + as.numeric(Educ == 2)*0.12 + as.numeric(Educ == 3)*0.16)

# birthweight dependent on Sex of child, education Smoking , standard deviation depends on Sex.
Wgt0 <- round(rnorm(n,2950 + 140*Sex + 80*(Educ==2) + 160*(Educ==3) -200*Smoke, 390 + 30*Sex),0)


##########################
#               	       #         
# TREATMENTS  A1 AND A2  #
#                        #  
##########################

#There are four definitions of treatment, A1,A2,A3,A4 (Intervention,intervention actually received, starting to breastfeed and breastfeeding for three months)) according to the four definitions in the paper


A1 <- rbinom(n,1,0.5)
A2.prob<-inv.logit(-1.9+ 0.1*Age + 0.5*as.numeric(Educ == 2)+1.0*as.numeric(Educ == 3)-1*Smoke) # Probability to participate in programme depending on Educ, smoke and age
summary( A2.prob)
mean(A2.prob[Educ==1])
mean(A2.prob[Educ==2])
mean(A2.prob[Educ==3])
mean(A2.prob[Smoke==1])
mean(A2.prob[Smoke==0])

A2potential<-rbinom(n,1,A2.prob)  # Would intervention been taken up if it was offered?
A2.observed<-ifelse(A1==1,A2potential,0)


rm(A2.prob)
#######################################################################
#                                                                     #
## PRINCIPAL STRATA AND A3=STARTING TO BREASTFEED                     #
#######################################################################

#coefficients for principal strata probabilities

new.tmtmod <- c(-2.5, #constant
                1.5, #difference between always takers and always/compliers
                0.25, # educ =2
                0.5, # educ =3
                0.1, #age
                0.008,#sex
                -0.5, #smoke extra jan 2018
                0.0006) #birthweigt

# principle strata for BF

X.always <- cbind(1,0,Educ==2,Educ==3,Age,Sex,Smoke,(Wgt0-3000)) #defining always takers
# probability to be an always taker
pA3.always <- inv.logit(X.always%*%new.tmtmod)
X.always.comp <- cbind(1,1,Educ==2,Educ==3,Age,Sex,Smoke, (Wgt0-3000))  #defining always takers and starting when encouraged
# probability to be an always BF taker or complier if A1 =1  
pA3.always.comp <- inv.logit(X.always.comp%*%new.tmtmod)
summary(pA3.always)
summary(pA3.always.comp)

# draw princ.strata.variable which is a random number between 0 and 1. 
# If princ.strata.variable is smaller than pA3.always, women is always BF taker
# We assume that only those with A2potential = 1 can be compliers, because those are the ones who actually go to 
# the programme and can be influenced.
# So if A2potential = 0 and princ.strata.variable > pA3.always, the women is never taker
# If A2potential = 1 and  pA3.always <= princ.strata.variable  <= pA3.1always.comp, the women is a complier 
# If A2potential = 1 and  princ.strata.variable  > pA3.1always.comp, the women is a a never taker. 

princ.strata.variable <- as.numeric(runif(n))

# Define compliers (1), never-takers(0) and always-takers(2);  we assume monotonicity
princ.strat <- 2 * (princ.strata.variable < pA3.always ) + 
               1 * (A2potential == 1& princ.strata.variable >= pA3.always  & princ.strata.variable < pA3.always.comp) + 0

# define potential A3 for A1=1 and A1 =0
A3pot.A1.1 <- 1*((princ.strat==2) | (princ.strat==1))
A3pot.A1.0 <- 1*( princ.strat==2)

# Now define what the potential A3 would be if A1=1 and A2=1
# Then we do not have to look whether A2potential = 1 
A3pot.A1.1.A2.1 <- 1*( princ.strata.variable < pA3.always.comp) 

# define observed A3 depending on principle strata and observed value of A1
A3.observed <- ifelse ( A1==1, A3pot.A1.1, A3pot.A1.0)

# list of the different potential treatments per women
cbind(Educ, A1, A2potential, pA3.always, pA3.always.comp,princ.strat,A3pot.A1.1,A3pot.A1.0, A3pot.A1.1.A2.1)[1:30,]


rm(new.tmtmod,X.always, X.always.comp)
########################################################################################################



#############################
#                           #
#                           #
# DURATION OF BREASTFEEDING #
#                           #
#############################


# We define the duration of breastfeeding, given that you start.
# It depends on Educ, Wgt0, Allergy, Age, Sex, Infc and Hosp


parameters.duration <- c(1.5, 
                         0.001, #centered birthweight
                         1, # education medium
                         1.5, #education high
                           1, #A2
                         0.5,#allergy
                         0.05, #age
                         -1,  #cesarean
                         -.5) #sex
      
X.dur.A2.1 <- cbind(1,(Wgt0-3000),Educ==2,Educ==3,1,Allergy,Age,Cesarean,Sex) #covariates with A2 set to 1
X.dur.A2.0 <- cbind(1,(Wgt0-3000),Educ==2,Educ==3,0,Allergy,Age,Cesarean,Sex) #covariates with A2 set to 0

# duration if A2 =1 and A3=1 (women starts breastfeeding)
durpot.A2.1.A3.1 <- pmin(rpois(n,pmax(X.dur.A2.1%*%parameters.duration,0)),3)
durpot.A1.0.A3.1 <- pmin(rpois(n,pmax(X.dur.A2.0%*%parameters.duration,0)),3)

round(prop.table(table(durpot.A2.1.A3.1)),2)
round(prop.table(table(durpot.A1.0.A3.1)),2)

# potential duration  if A2=A2potential and A3=1
durpot.A2pot.A3.1 <- ifelse(A2potential==1, durpot.A2.1.A3.1, durpot.A1.0.A3.1)  

# potential duration  if A1=1
durpot.A1.1 <- ifelse(A3pot.A1.1==0, 0, durpot.A2pot.A3.1)  

# potential duration  if A1=0
durpot.A1.0 <- ifelse(A3pot.A1.0==0, 0, durpot.A1.0.A3.1)  

# potential duration  if A1=1 and A2=1 ,  everyone follows programme
durpot.A1.1.A2.1 <- ifelse(A3pot.A1.1.A2.1==0, 0, durpot.A2.1.A3.1)  

# actual duration observed
dur.observed <-ifelse(A1==1, durpot.A1.1, durpot.A1.0)

# This table is giving the different potental durations for different potential treatments
cbind(A1, A2potential, A3pot.A1.1, A3pot.A1.0,A3pot.A1.1.A2.1, durpot.A2.1.A3.1, durpot.A1.0.A3.1, dur.observed)[1:30,]

# Observed A4: the mother started breastfeeding and continued for the full 3 months
A4.observed <-ifelse(dur.observed >=3,1,0)

rm(parameters.duration,X.dur.A2.1,X.dur.A2.0)

################
#              #
#  OUTCOMES    #
#              #
################

# Defining the potential outcome models under the four different durations of breastfeeding

#  define random subject effect
subject <- rnorm(n,0,40)

new.wtmod <- c(5800,
               1,# Birthweigth
               16,# Location==2
               -20,# Location==3
               -15, # Location==4
               10, # Education==2
               20,  # Education==3
               -50, # Smoke
               -25, # Allergy
               -10, # Age
               -40, # Cesarean
               500, # Sex
               100, # Duration (per month)
               50, # Duration and medium education interaction
               100, # Duration and lower education interaction
               -0.02,# Birthweight and duration interaction, lower BW more effect of BF
               50 ) # Smoking and duration interaction

weight.3mnd <- function(duration)
{X.wt3 <- cbind(1,I(Wgt0-3000),Location==2,Location==3,Location==4,Educ==2,Educ==3,
                Smoke,Allergy,Age,Cesarean,Sex,duration, duration*(Educ==2), duration*(Educ==1),
                duration*(Wgt0-3000),duration*Smoke) 
round(X.wt3%*%new.wtmod + subject + rnorm(n,0,1),0) 
}

# potential outcome if BF is started and duration is fixed.
Wgt3pot.dur0 <-  weight.3mnd(0)
Wgt3pot.dur1 <-  weight.3mnd(1)
Wgt3pot.dur2 <-  weight.3mnd(2)
Wgt3pot.dur3 <-  weight.3mnd(3)


mean(Wgt3pot.dur3[Educ==1]) - mean(Wgt3pot.dur0[Educ==1])
mean(Wgt3pot.dur3[Educ==2]) - mean(Wgt3pot.dur0[Educ==2])
mean(Wgt3pot.dur3[Educ==3]) - mean(Wgt3pot.dur0[Educ==3])
mean(Wgt3pot.dur3[Smoke==1]) - mean(Wgt3pot.dur0[Smoke==1])
mean(Wgt3pot.dur3[Smoke==0]) - mean(Wgt3pot.dur0[Smoke==0])
mean(Wgt3pot.dur3[Wgt0<2000]) - mean(Wgt3pot.dur0[Wgt0<2000])
mean(Wgt3pot.dur3[Wgt0<3000]) - mean(Wgt3pot.dur0[Wgt0<3000])
mean(Wgt3pot.dur3[Wgt0>4000]) - mean(Wgt3pot.dur0[Wgt0>4000])
mean(Wgt3pot.dur3[Wgt0>4500]) - mean(Wgt3pot.dur0[Wgt0>4500])



# potential outcome if programme is followed and BF is started.
# Use duration we would have observed if A2 =1
Wgt3pot.A2.1.A3.1 <- (durpot.A2.1.A3.1==0)*Wgt3pot.dur0 + (durpot.A2.1.A3.1==1)*Wgt3pot.dur1 +
                       (durpot.A2.1.A3.1==2)* Wgt3pot.dur2 + (durpot.A2.1.A3.1==3)* Wgt3pot.dur3

# potential outcome if programme is not followed and BF is started.
# Use duration we would have observed if A2 =0
Wgt3pot.A1.0.A3.1 <- (durpot.A1.0.A3.1==0)*Wgt3pot.dur0 + (durpot.A1.0.A3.1==1)*Wgt3pot.dur1 +
                       (durpot.A1.0.A3.1==2)* Wgt3pot.dur2 + (durpot.A1.0.A3.1==3)* Wgt3pot.dur3

# potential outcome if intervention if offered and BF is started.
# Use duration we would have observed if A2 =A2 potential
Wgt3pot.A1.1.A3.1 <- ifelse((A2potential==1), Wgt3pot.A2.1.A3.1,Wgt3pot.A1.0.A3.1)
  
# potential outcome if everyone would have followed programme but not necessarily started BF
# we use thatA1=1, A2=1 and that if A3 potential = 0, then weight = Wgt3pot.dur0, else  Wgt3pot.A2.1.A3.1
Wgt3pot.A2.1 <- ifelse((A3pot.A1.1.A2.1==0), Wgt3pot.dur0, Wgt3pot.A2.1.A3.1)

# potential outcome if no-one would have followed programme (ie under randomisation for control)
# we use A3 potential under A1=0 and duration under A2=0
Wgt3pot.A1.0 <- ifelse((A3pot.A1.0==0), Wgt3pot.dur0, Wgt3pot.A1.0.A3.1)

# potential outcome under randomisation for intervention
# we use weight under A2=potential
Wgt3pot.A1.1 <- ifelse((A2potential==1), Wgt3pot.A2.1, Wgt3pot.A1.0)


# Observed weight
Wgt3.observed <-ifelse( A1==0, Wgt3pot.A1.0, Wgt3pot.A1.1)

rm(subject,new.wtmod,weight.3mnd)

###############################
# WRITE DATA TO FILE         #
#                             #
##############################



PROBIT.sim <- data.frame(A1,A2.observed,Location,Age,Educ,Smoke, Allergy, Cesarean, 
                Wgt0, Sex, A3.observed, dur.observed, A4.observed, Wgt3.observed,
                A2potential,  A3pot.A1.0, A3pot.A1.1, A3pot.A1.1.A2.1, 
                durpot.A1.0, durpot.A1.1, durpot.A1.1.A2.1, durpot.A1.0.A3.1, durpot.A2.1.A3.1, durpot.A2pot.A3.1,
                Wgt3pot.A1.0, Wgt3pot.A1.1,  Wgt3pot.A2.1, Wgt3pot.A1.0.A3.1, Wgt3pot.A1.1.A3.1 , Wgt3pot.A2.1.A3.1, 
                Wgt3pot.dur0,  Wgt3pot.dur1, Wgt3pot.dur2, Wgt3pot.dur3)          



write.csv(PROBIT.sim,"PROBITsim_april2018.csv",row.names = FALSE)

