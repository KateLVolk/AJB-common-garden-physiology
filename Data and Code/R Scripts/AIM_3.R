##GEUTRIL PHYSIO MANUSCRIP:::::AIM 3
#triat heritabilities, sem and evolvability

setwd("D:/GEUTRI_PHYSIO/Submission Documents/Data and Code/Data")

#load data
stoma<-read.csv("stomata_units.csv")
stoma$Block.ID <- as.factor(stoma$Block.ID)
#Visualize Normality
qqnorm(stoma$GCLab);qqline(stoma$GCLab)
qqnorm(stoma$GCLad);qqline(stoma$GCLad)
qqnorm(stoma$SDab);qqline(stoma$SDab)
qqnorm(stoma$Sdad);qqline(stoma$Sdad)
qqnorm(stoma$SAIab);qqline(stoma$SAIab)
qqnorm(stoma$SAIad);qqline(stoma$SAIad)

#load libraries
library(lme4); library(lmerTest); library(effects);library(permute); library(knitr)
library(tidyverse); library(rptR);library(emmeans);library(ggplot2); library(tidyverse)
library(data.table);library(afex)

#-------------------------------------------------------------------------------
#-----------------------------ABAXIAL GCL---------------------------------------
#-------------------------------------------------------------------------------
stoma = subset(stoma,!is.na(GCLab))
##SUBSET and do this 
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 

PRA <- subset(stoma, Region == "Prairie")
GLA<- subset(stoma, Region =="GL Alvar")
MBA <- subset(stoma, Region =="M Alvar")
stoma = subset(stoma,!is.na(GCLab))

##Prairie Model
pra.mod <- lmer_alt(GCLab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)
anova(pra.mod)
print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family h2 == VA/VP 
prair.fam.h2<-(2.5*vars[3,4]) / (vars[3,4]+vars[1,4]+vars[5,4]) ##added both Vg:e and residual to Vg here
prair.fam.h2 #0.2164769

#calculate standard error from standard deviation (sd / sqrt(n))
prah2sem<-(2.5*vars[3,5])/sqrt(nrow(PRA))
prah2sem#0.1464105

#CVa 
pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$GCLab))
pracva #0.03632229

##Manitoba Model ##something wrong with nesting I think..model will not run
mba.mod <- lmer_alt(GCLab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = MBA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                     check.nobs.vs.rankZ = "ignore",
                                                     check.nobs.vs.nRE="ignore"))
summary(mba.mod) ## returns 0 for pop:block vcov
vars<-as.data.frame(VarCorr(mba.mod))

#Manitoba Family h2

mani.fam.h2<-(2.5*vars[3,4])/(vars[3,4]+vars[1,4]+vars[5,4])
mani.fam.h2 #0.3905943

#calculate standard error from standard deviation (sd / sqrt(n))
mbah2sem<-(2.5*vars[3,5])/sqrt(nrow(MBA))
mbah2sem#0.1951921

#CVa 
mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$GCLab))
mbacva #0.04246351



##Great Lake Model ##something wrong with nesting I think..model will not run unless I add lmerControl
gla.mod<-lmer_alt(GCLab~
                    (1| Population) +
                    (1| Family.Unique)+
                    (1|Population:Block.ID)+
                    (1|Family.Unique:Block.ID), 
                  data = GLA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                   check.nobs.vs.rankZ = "ignore",
                                                   check.nobs.vs.nRE="ignore"))
summary(gla.mod)
vars<-as.data.frame(VarCorr(gla.mod))

##Great Lake h2
great.fam.h2<-(2.5*vars[3,4]) /(vars[3,4]+vars[1,4]+vars[5,4])
great.fam.h2 #0.4269977


#calculate standard error from standard deviation (sd / sqrt(n))
glah2sem<-(2.5*vars[3,5])/sqrt(nrow(GLA))
glah2sem#0.1044503

#CVa 
glacva<-((sqrt(2.5*vars[3,4]))/mean(GLA$GCLab))
glacva #0.04488784

#-------------------------------------------------------------------------------
#-----------------------------ADAXIAL GCL---------------------------------------
#-------------------------------------------------------------------------------
stoma<-read.csv("stomata_units.csv")
stoma$Block.ID <- as.factor(stoma$Block.ID)
stoma = subset(stoma,!is.na(GCLad))
##SUBSET and do this 
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 

PRA <- subset(stoma, Region == "Prairie")
GLA<- subset(stoma, Region =="GL Alvar")
MBA <- subset(stoma, Region =="M Alvar")


##Prairie Model
pra.mod <- lmer_alt(GCLad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)
anova(pra.mod)
print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family h2 == VA/VP 
prair.fam.h2<-(2.5*vars[3,4]) / (vars[3,4]+vars[1,4]+vars[5,4]) ##added both Vg:e and residual to Vg here
prair.fam.h2 #0.207254

#calculate standard error from standard deviation (sd / sqrt(n))
prah2sem<-(2.5*vars[3,5])/sqrt(nrow(PRA))
prah2sem#0.1464105

#CVa 
pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$GCLad))
pracva #0.03632229




##Manitoba Model ##something wrong with nesting I think..model will not run
mba.mod <- lmer_alt(GCLad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = MBA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                     check.nobs.vs.rankZ = "ignore",
                                                     check.nobs.vs.nRE="ignore"))
summary(mba.mod) ## returns 0 for pop:block vcov
vars<-as.data.frame(VarCorr(mba.mod))

#Manitoba Family h2

mani.fam.h2<-(2.5*vars[3,4])/(vars[3,4]+vars[1,4]+vars[5,4])
mani.fam.h2 #0

#calculate standard error from standard deviation (sd / sqrt(n))
mbah2sem<-(2.5*vars[3,5])/sqrt(nrow(MBA))
mbah2sem#0

#CVa 
mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$GCLad))
mbacva #0



##Great Lake Model ##something wrong with nesting I think..model will not run unless I add lmerControl
gla.mod<-lmer_alt(GCLad~
                    (1| Population) +
                    (1| Family.Unique)+
                    (1|Population:Block.ID)+
                    (1|Family.Unique:Block.ID), 
                  data = GLA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                   check.nobs.vs.rankZ = "ignore",
                                                   check.nobs.vs.nRE="ignore"))
summary(gla.mod)
vars<-as.data.frame(VarCorr(gla.mod))

##Great Lake h2
great.fam.h2<-(2.5*vars[3,4]) /(vars[3,4]+vars[1,4]+vars[5,4])
great.fam.h2 #0.3911679


#calculate standard error from standard deviation (sd / sqrt(n))
glah2sem<-(2.5*vars[3,5])/sqrt(nrow(GLA))
glah2sem# 0.09704961

#CVa 
glacva<-((sqrt(2.5*vars[3,4]))/mean(GLA$GCLad))
glacva #0.04228569
#-------------------------------------------------------------------------------
#-----------------------------ABAXIAL SD---------------------------------------
#-------------------------------------------------------------------------------
stoma<-read.csv("stomata_units.csv")
stoma$Block.ID <- as.factor(stoma$Block.ID)
stoma = subset(stoma,!is.na(SDab))
##SUBSET and do this 
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 

PRA <- subset(stoma, Region == "Prairie")
GLA<- subset(stoma, Region =="GL Alvar")
MBA <- subset(stoma, Region =="M Alvar")


##Prairie Model
pra.mod <- lmer_alt(SDab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)
anova(pra.mod)
print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family h2 == VA/VP 
prair.fam.h2<-(2.5*vars[3,4]) / (vars[3,4]+vars[1,4]+vars[5,4]) ##added both Vg:e and residual to Vg here
prair.fam.h2 #0.133251

#calculate standard error from standard deviation (sd / sqrt(n))
prah2sem<-(2.5*vars[3,5])/sqrt(nrow(PRA))
prah2sem#2.925845

#CVa 
pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SDab))
pracva #0.09931627


##Manitoba Model ##something wrong with nesting I think..model will not run
mba.mod <- lmer_alt(SDab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = MBA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                     check.nobs.vs.rankZ = "ignore",
                                                     check.nobs.vs.nRE="ignore"))
summary(mba.mod) ## returns 0 for pop:block vcov
vars<-as.data.frame(VarCorr(mba.mod))

#Manitoba Family h2

mani.fam.h2<-(2.5*vars[3,4])/(vars[3,4]+vars[1,4]+vars[5,4])
mani.fam.h2 #0.1320467

#calculate standard error from standard deviation (sd / sqrt(n))
mbah2sem<-(2.5*vars[3,5])/sqrt(nrow(MBA))
mbah2sem#4.192403

#CVa 
mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SDab))
mbacva #0.09924224


##Great Lake Model ##something wrong with nesting I think..model will not run unless I add lmerControl
gla.mod<-lmer_alt(SDab~
                    (1| Population) +
                    (1| Family.Unique)+
                    (1|Population:Block.ID)+
                    (1|Family.Unique:Block.ID), 
                  data = GLA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                   check.nobs.vs.rankZ = "ignore",
                                                   check.nobs.vs.nRE="ignore"))
summary(gla.mod)
vars<-as.data.frame(VarCorr(gla.mod))

##Great Lake h2
great.fam.h2<-(2.5*vars[3,4]) /(vars[3,4]+vars[1,4]+vars[5,4])
great.fam.h2 #0.2692473


#calculate standard error from standard deviation (sd / sqrt(n))
glah2sem<-(2.5*vars[3,5])/sqrt(nrow(GLA))
glah2sem# 2.785598

#CVa 
glacva<-((sqrt(2.5*vars[3,4]))/mean(GLA$SDab))
glacva #0.1250742
#-------------------------------------------------------------------------------
#-----------------------------ADAXIAL SD---------------------------------------
#-------------------------------------------------------------------------------
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SDad))
##SUBSET and do this 
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 

PRA <- subset(stoma, Region == "Prairie")
GLA<- subset(stoma, Region =="GL Alvar")
MBA <- subset(stoma, Region =="M Alvar")


##Prairie Model
pra.mod <- lmer_alt(SDad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)
anova(pra.mod)
print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family h2 == VA/VP 
prair.fam.h2<-(2.5*vars[3,4]) / (vars[3,4]+vars[1,4]+vars[5,4]) ##added both Vg:e and residual to Vg here
prair.fam.h2 #0.3651705

#calculate standard error from standard deviation (sd / sqrt(n))
prah2sem<-(2.5*vars[3,5])/sqrt(nrow(PRA))
prah2sem#3.491505

#CVa 
pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SDad))
pracva #0.1732927


##Manitoba Model ##something wrong with nesting I think..model will not run
mba.mod <- lmer_alt(SDad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = MBA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                     check.nobs.vs.rankZ = "ignore",
                                                     check.nobs.vs.nRE="ignore"))
summary(mba.mod) ## returns 0 for pop:block vcov
vars<-as.data.frame(VarCorr(mba.mod))

#Manitoba Family h2

mani.fam.h2<-(2.5*vars[3,4])/(vars[3,4]+vars[1,4]+vars[5,4])
mani.fam.h2 #0.271087

#calculate standard error from standard deviation (sd / sqrt(n))
mbah2sem<-(2.5*vars[3,5])/sqrt(nrow(MBA))
mbah2sem#4.353598

#CVa 
mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SDad))
mbacva #0.1297914


##Great Lake Model ##something wrong with nesting I think..model will not run unless I add lmerControl
gla.mod<-lmer_alt(SDad~
                    (1| Population) +
                    (1| Family.Unique)+
                    (1|Population:Block.ID)+
                    (1|Family.Unique:Block.ID), 
                  data = GLA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                   check.nobs.vs.rankZ = "ignore",
                                                   check.nobs.vs.nRE="ignore"))
summary(gla.mod)
vars<-as.data.frame(VarCorr(gla.mod))

##Great Lake h2
great.fam.h2<-(2.5*vars[3,4]) /(vars[3,4]+vars[1,4]+vars[5,4])
great.fam.h2 #0.4453984


#calculate standard error from standard deviation (sd / sqrt(n))
glah2sem<-(2.5*vars[3,5])/sqrt(nrow(GLA))
glah2sem# 2.558508

#CVa 
glacva<-((sqrt(2.5*vars[3,4]))/mean(GLA$SDad))
glacva #0.1681944

#-------------------------------------------------------------------------------
#-----------------------------ABAXIAL SAI---------------------------------------
#-------------------------------------------------------------------------------
stoma<-read.csv("stomata_units.csv")

stoma = subset(stoma,!is.na(SAIab))
##SUBSET and do this 
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 

PRA <- subset(stoma, Region == "Prairie")
GLA<- subset(stoma, Region =="GL Alvar")
MBA <- subset(stoma, Region =="M Alvar")


##Prairie Model
pra.mod <- lmer_alt(SAIab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)
anova(pra.mod)
print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family h2 == VA/VP 
prair.fam.h2<-(2.5*vars[3,4]) / (vars[3,4]+vars[1,4]+vars[5,4]) ##added both Vg:e and residual to Vg here
prair.fam.h2 #0.1674781

#calculate standard error from standard deviation (sd / sqrt(n))
prah2sem<-(2.5*vars[3,5])/sqrt(nrow(PRA))
prah2sem#0.06192677

#CVa 
pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SAIab))
pracva #0.1005247


##Manitoba Model ##something wrong with nesting I think..model will not run
mba.mod <- lmer_alt(SAIab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = MBA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                     check.nobs.vs.rankZ = "ignore",
                                                     check.nobs.vs.nRE="ignore"))
summary(mba.mod) ## returns 0 for pop:block vcov
vars<-as.data.frame(VarCorr(mba.mod))

#Manitoba Family h2

mani.fam.h2<-(2.5*vars[3,4])/(vars[3,4]+vars[1,4]+vars[5,4])
mani.fam.h2 #0

#calculate standard error from standard deviation (sd / sqrt(n))
mbah2sem<-(2.5*vars[3,5])/sqrt(nrow(MBA))
mbah2sem#0

#CVa 
mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SAIab))
mbacva #0



##Great Lake Model ##something wrong with nesting I think..model will not run unless I add lmerControl
gla.mod<-lmer_alt(SAIab~
                    (1| Population) +
                    (1| Family.Unique)+
                    (1|Population:Block.ID)+
                    (1|Family.Unique:Block.ID), 
                  data = GLA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                   check.nobs.vs.rankZ = "ignore",
                                                   check.nobs.vs.nRE="ignore"))
summary(gla.mod)
vars<-as.data.frame(VarCorr(gla.mod))

##Great Lake h2
great.fam.h2<-(2.5*vars[3,4]) /(vars[3,4]+vars[1,4]+vars[5,4])
great.fam.h2 #0.1702273


#calculate standard error from standard deviation (sd / sqrt(n))
glah2sem<-(2.5*vars[3,5])/sqrt(nrow(GLA))
glah2sem# 0.04300333

#CVa 
glacva<-((sqrt(2.5*vars[3,4]))/mean(GLA$SAIab))
glacva #0.09497059

#-------------------------------------------------------------------------------
#-----------------------------ADAXIAL SAI---------------------------------------
#-------------------------------------------------------------------------------
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SAIad))
##SUBSET and do this 
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 

PRA <- subset(stoma, Region == "Prairie")
GLA<- subset(stoma, Region =="GL Alvar")
MBA <- subset(stoma, Region =="M Alvar")


##Prairie Model
pra.mod <- lmer_alt(SAIad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)
anova(pra.mod)
print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family h2 == VA/VP 
prair.fam.h2<-(2.5*vars[3,4]) / (vars[3,4]+vars[1,4]+vars[5,4]) ##added both Vg:e and residual to Vg here
prair.fam.h2 #0.4423379

#calculate standard error from standard deviation (sd / sqrt(n))
prah2sem<-(2.5*vars[3,5])/sqrt(nrow(PRA))
prah2sem#0.1036254

#CVa 
pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SAIad))
pracva #0.1841476


##Manitoba Model ##something wrong with nesting I think..model will not run
mba.mod <- lmer_alt(SAIad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = MBA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                     check.nobs.vs.rankZ = "ignore",
                                                     check.nobs.vs.nRE="ignore"))
summary(mba.mod) ## returns 0 for pop:block vcov
vars<-as.data.frame(VarCorr(mba.mod))

#Manitoba Family h2

mani.fam.h2<-(2.5*vars[3,4])/(vars[3,4]+vars[1,4]+vars[5,4])
mani.fam.h2 #0.5246773

#calculate standard error from standard deviation (sd / sqrt(n))
mbah2sem<-(2.5*vars[3,5])/sqrt(nrow(MBA))
mbah2sem#0.1517817

#CVa 
mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SAIad))
mbacva #0.1735809


##Great Lake Model ##something wrong with nesting I think..model will not run unless I add lmerControl
gla.mod<-lmer_alt(SAIad~
                    (1| Population) +
                    (1| Family.Unique)+
                    (1|Population:Block.ID)+
                    (1|Family.Unique:Block.ID), 
                  data = GLA,  control=lmerControl(check.nobs.vs.nlev = "ignore",
                                                   check.nobs.vs.rankZ = "ignore",
                                                   check.nobs.vs.nRE="ignore"))
summary(gla.mod)
vars<-as.data.frame(VarCorr(gla.mod))

##Great Lake h2
great.fam.h2<-(2.5*vars[3,4]) /(vars[3,4]+vars[1,4]+vars[5,4])
great.fam.h2 #0.4095687


#calculate standard error from standard deviation (sd / sqrt(n))
glah2sem<-(2.5*vars[3,5])/sqrt(nrow(GLA))
glah2sem# 0.06020287

#CVa 
glacva<-((sqrt(2.5*vars[3,4]))/mean(GLA$SAIad))
glacva #0.1495959
