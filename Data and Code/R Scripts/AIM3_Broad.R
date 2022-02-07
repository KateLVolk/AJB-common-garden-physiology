##GEUTRIL PHYSIO MANUSCRIP:::::AIM 3.1
#BROAD SENSE HERITABILITY 
## VG/VP

setwd("D:/GEUTRI_PHYSIO/Submission Documents/Data and Code/Data")

#load data
stoma<-read.csv("stomata_units.csv")
stoma$Block.ID <- as.factor(stoma$Block.ID)

#load libraries
library(lme4); library(lmerTest); library(effects);library(permute); library(knitr)
library(tidyverse); library(rptR);library(emmeans);library(ggplot2); library(tidyverse)
library(data.table);library(afex)


## WUE?
wue<-read.csv("Refined_WUE.csv")
wue$Region<-as.factor(wue$Region)
wue$Population<-as.factor(wue$Population)
wue$Family.Unique<-as.factor(wue$Family.Unique)
wue$BLOCK<-as.factor(wue$BLOCK)
PRA <- subset(wue, Region == "Prairie")
PRA = subset(PRA,!is.na(dC13))
GLA<- subset(wue, Region =="GL Alvar")
GLA = subset(GLA,!is.na(dC13))
MBA <- subset(wue, Region =="MB_Alvar")
MBA = subset(MBA,!is.na(dC13))

pra.mod <- lmer_alt(dC13~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:BLOCK)+
                      (1|Family.Unique:BLOCK), 
                    data = PRA)
summary(pra.mod)

print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family H2 == VG/VP 
prair.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
prair.fam.H2 #0.9932033

#calculate standard error from standard deviation (sd / sqrt(n))
praH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
praH2sem#0.2262908
 ##########ISSUE 
mba.mod <- lmer_alt(dC13~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:BLOCK)+
                      (1|Family.Unique:BLOCK), 
                    data = MBA)
summary(mba.mod)

print(VarCorr(mba.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(mba.mod))

#H2 == VG/VP 
mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 #0.9913402

#calculate standard error from standard deviation (sd / sqrt(n))
mbaH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
mbaH2sem#0.1708571

gla.mod <- lmer_alt(dC13~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:BLOCK)+
                      (1|Family.Unique:BLOCK), 
                    data = GLA)
summary(gla.mod)

print(VarCorr(gla.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(gla.mod))

H2 == VG/VP 
mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 #0.9913402

#calculate standard error from standard deviation (sd / sqrt(n))
mbaH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
mbaH2sem#0.1708571
#-------------------------------------------------------------------------------
#-----------------------------ABAXIAL GCL---------------------------------------
#-------------------------------------------------------------------------------
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 

PRA <- subset(stoma, Region == "Prairie")
PRA = subset(PRA,!is.na(GCLab))
GLA<- subset(stoma, Region =="GL Alvar")
GLA = subset(GLA,!is.na(GCLab))
MBA <- subset(stoma, Region =="M Alvar")
MBA = subset(MBA,!is.na(GCLab))


##Prairie Model
pra.mod <- lmer_alt(GCLab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)

print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family H2 == VG/VP 
prair.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
prair.fam.H2 #0.08659076

#calculate standard error from standard deviation (sd / sqrt(n))
praH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
praH2sem#0.03788302

#CVa 
#pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$GCLab))
#pracva #0.03632229

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

#Manitoba Family H2

mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 #0.8049687

#calculate standard error from standard deviation (sd / sqrt(n))
ManH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
ManH2sem
#CVa 
#mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$GCLab))
#mbacva #0.04246351



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

##Great Lake H2

greatl.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
greatl.fam.H2 #0.8536896

#calculate standard error from standard deviation (sd / sqrt(n))
GreatH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
GreatH2sem

#-------------------------------------------------------------------------------
#-----------------------------ADAXIAL GCL---------------------------------------
#-------------------------------------------------------------------------------
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 
stoma<-read.csv("stomata_units.csv")
PRA <- subset(stoma, Region == "Prairie")
PRA = subset(PRA,!is.na(GCLad))
GLA<- subset(stoma, Region =="GL Alvar")
GLA = subset(GLA,!is.na(GCLad))
MBA <- subset(stoma, Region =="M Alvar")
MBA = subset(MBA,!is.na(GCLad))

##Prairie Model
pra.mod <- lmer_alt(GCLad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)

print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family H2 == VG/VP 
prair.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
prair.fam.H2 #0.08659076

#calculate standard error from standard deviation (sd / sqrt(n))
praH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
praH2sem#

#CVa 
#pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$GCLad))
#pracva #0.03632229

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

#Manitoba Family H2

mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 #0.8015926

#calculate standard error from standard deviation (sd / sqrt(n))
ManH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
ManH2sem#
#CVa 
#mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$GCLad))
#mbacva #0.04246351

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

##Great Lake H2

greatl.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
greatl.fam.H2 #0.856333

#calculate standard error from standard deviation (sd / sqrt(n))
GreatH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
GreatH2sem#


#-------------------------------------------------------------------------------
#-----------------------------ABAXIAL SD---------------------------------------
#-------------------------------------------------------------------------------
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 
stoma<-read.csv("stomata_units.csv")
PRA <- subset(stoma, Region == "Prairie")
PRA = subset(PRA,!is.na(SDab))
GLA<- subset(stoma, Region =="GL Alvar")
GLA = subset(GLA,!is.na(SDab))
MBA <- subset(stoma, Region =="M Alvar")
MBA = subset(MBA,!is.na(SDab))

##Prairie Model
pra.mod <- lmer_alt(SDab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)

print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family H2 == VG/VP 
prair.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
prair.fam.H2 #0.3201125

#calculate standard error from standard deviation (sd / sqrt(n))
praH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
praH2sem#0.03788302

#CVa 
#pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SDab))
#pracva #0.03632229

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

#Manitoba Family H2

mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 

#calculate standard error from standard deviation (sd / sqrt(n))
ManH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
ManH2sem
#CVa 
#mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SDab))
#mbacva #0.04246351

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

##Great Lake H2

greatl.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
greatl.fam.H2

#calculate standard error from standard deviation (sd / sqrt(n))
GreatH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
GreatH2sem

#-------------------------------------------------------------------------------
#-----------------------------ADAXIAL SD---------------------------------------
#-------------------------------------------------------------------------------
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 
stoma<-read.csv("stomata_units.csv")
PRA <- subset(stoma, Region == "Prairie")
PRA = subset(PRA,!is.na(SDad))
GLA<- subset(stoma, Region =="GL Alvar")
GLA = subset(GLA,!is.na(SDad))
MBA <- subset(stoma, Region =="M Alvar")
MBA = subset(MBA,!is.na(SDad))

##Prairie Model
pra.mod <- lmer_alt(SDad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)

print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family H2 == VG/VP 
prair.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
prair.fam.H2 #0.3453555

#calculate standard error from standard deviation (sd / sqrt(n))
praH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
praH2sem#4.169132

#CVa 
#pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SDad))
#pracva #0.03632229

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

#Manitoba Family H2

mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 #0.9000177

#calculate standard error from standard deviation (sd / sqrt(n))
ManH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
ManH2sem #5.808519
#CVa 
#mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SDad))
#mbacva #0.04246351

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

##Great Lake H2

greatl.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
greatl.fam.H2 #0.845774

#calculate standard error from standard deviation (sd / sqrt(n))
GreatH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
GreatH2sem#7.069845

#-------------------------------------------------------------------------------
#-----------------------------ABAXIAL SAI---------------------------------------
#-------------------------------------------------------------------------------
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 
stoma<-read.csv("stomata_units.csv")
PRA <- subset(stoma, Region == "Prairie")
PRA = subset(PRA,!is.na(SAIab))
GLA<- subset(stoma, Region =="GL Alvar")
GLA = subset(GLA,!is.na(SAIab))
MBA <- subset(stoma, Region =="M Alvar")
MBA = subset(MBA,!is.na(SAIab))

##Prairie Model
pra.mod <- lmer_alt(SAIab~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)

print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family H2 == VG/VP 
prair.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
prair.fam.H2 #0.3432657

#calculate standard error from standard deviation (sd / sqrt(n))
praH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
praH2sem#0.08684502

#CVa 
#pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SAIab))
#pracva #0.03632229

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

#Manitoba Family H2

mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 #0.803068

#calculate standard error from standard deviation (sd / sqrt(n))
ManH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
ManH2sem #0.1171684
#CVa 
#mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SAIab))
#mbacva #0.04246351

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

##Great Lake H2

greatl.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
greatl.fam.H2 #0.8577337

#calculate standard error from standard deviation (sd / sqrt(n))
GreatH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
GreatH2sem#0.1777078

#-------------------------------------------------------------------------------
#-----------------------------ADAXIAL SAI---------------------------------------
#-------------------------------------------------------------------------------
#(1}Population )+ (1| Population: Block)+ (same for family +1| Family) or add residula variance to each term if Family: block 
stoma<-read.csv("stomata_units.csv")
PRA <- subset(stoma, Region == "Prairie")
PRA = subset(PRA,!is.na(SAIad))
GLA<- subset(stoma, Region =="GL Alvar")
GLA = subset(GLA,!is.na(SAIad))
MBA <- subset(stoma, Region =="M Alvar")
MBA = subset(MBA,!is.na(SAIad))

##Prairie Model
pra.mod <- lmer_alt(SAIad~
                      (1| Population) +
                      (1| Family.Unique)+
                      (1|Population:Block.ID)+
                      (1|Family.Unique:Block.ID), 
                    data = PRA)
summary(pra.mod)

print(VarCorr(pra.mod), comp=c("Variance", "Std.Dev."))
vars<-as.data.frame(VarCorr(pra.mod))

#Prairie Family H2 == VG/VP 
prair.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
prair.fam.H2 #0.3330234

#calculate standard error from standard deviation (sd / sqrt(n))
praH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
praH2sem#0.107087

#CVa 
#pracva<-((sqrt(2.5*vars[3,4]))/mean(PRA$SAIad))
#pracva #0.03632229

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

#Manitoba Family H2

mani.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
mani.fam.H2 #0.8423006

#calculate standard error from standard deviation (sd / sqrt(n))
ManH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
ManH2sem #0.1483975
#CVa 
#mbacva<-((sqrt(2.5*vars[3,4]))/mean(MBA$SAIad))
#mbacva #0.04246351

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

##Great Lake H2

greatl.fam.H2<-(vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]) / 
  (vars[1,4]+vars[2,4]+vars[3,4]+vars[4,4]+vars[5,4]) 
greatl.fam.H2 #0.85135

#calculate standard error from standard deviation (sd / sqrt(n))
GreatH2sem<-(vars[1,5]+vars[2,5]+vars[3,5]+vars[4,5])/sqrt(nrow(PRA))
GreatH2sem#0.1638455
