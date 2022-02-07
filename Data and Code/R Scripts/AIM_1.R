#AIM 1: does prairie smoke exhibit genetic differentiation in physiological traits?

setwd("D:/GEUTRI_PHYSIO/spreadsheets/finals")

#load data
stoma<-read.csv("stomata_units.csv")
wue<- read.csv("Refined_WUE.csv")
##Trait means for WUE
PRA <- subset(wue, Region == "Prairie")
GLA<- subset(wue, Region =="GL Alvar")
MBA <- subset(wue, Region =="MB_Alvar")

#Visualize Normality
qqnorm(stoma$GCLab);qqline(stoma$GCLab)
qqnorm(stoma$GCLad);qqline(stoma$GCLad)
qqnorm(stoma$SDab);qqline(stoma$SDab)
qqnorm(stoma$Sdad);qqline(stoma$Sdad)
qqnorm(stoma$SAIab);qqline(stoma$SAIab)
qqnorm(stoma$SAIad);qqline(stoma$SAIad)
qqnorm(wue$dC13);qqline(wue$dC13)
#load libraries
library(tidyverse)

##set up models, test ANVOA, and post-hoc tests

# GCL ab
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(GCLab))
#model
gclab.mod<-lm((GCLab~Region), data=stoma)
summary(gclab.mod)
#anova
anova(gclab.mod)
#post-hoc 
TukeyHSD(aov(gclab.mod))

# GCL ad
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(GCLad))
#model
gclad.mod<-lm((GCLad~Region), data=stoma)
summary(gclad.mod)
#anova
anova(gclad.mod)
#post-hoc 
TukeyHSD(aov(gclad.mod))

# SDab
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SDab))
#model
SDab.mod<-lm((SDab~Region), data=stoma)
summary(SDab.mod)
#anova
anova(SDab.mod)
#post-hoc 
TukeyHSD(aov(SDab.mod))

# SDad
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(Sdad))
#model
Sdad.mod<-lm((Sdad~Region), data=stoma)
summary(Sdad.mod)
#anova
anova(Sdad.mod)
#post-hoc 
TukeyHSD(aov(Sdad.mod))

# SAIab
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SAIab))
#model
SAIab.mod<-lm((SAIab~Region), data=stoma)
summary(SAIab.mod)
#anova
anova(SAIab.mod)
#post-hoc 
TukeyHSD(aov(SAIab.mod))

# SAIad
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SAIad))
#model
SAIad.mod<-lm((SAIad~Region), data=stoma)
summary(SAIad.mod)
#anova
anova(SAIad.mod)
#post-hoc 
TukeyHSD(aov(SAIad.mod))

##WUE
# SAIab
wue = subset(wue,!is.na(dC13))
#model
dC13.mod<-lm((dC13~Region), data=wue)
summary(dC13.mod)
#anova
anova(dC13.mod)
#post-hoc 
TukeyHSD(aov(dC13.mod))


## Means and Standard Errors for traits
#GCLab
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(GCLab))
GCLab1<- stoma %>%
  group_by(Region)%>%
  summarise(
    regmean=mean(GCLab), 
    regsme=sd(GCLab)/sqrt(length(GCLab)))
GCLab1

#GCLad
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(GCLad))
GCLad1<- stoma %>%
  group_by(Region)%>%
  summarise(
    regmean=mean(GCLad), 
    regsme=sd(GCLad)/sqrt(length(GCLad)))
GCLad1

#SDab
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SDab))
SDab1<- stoma %>%
  group_by(Region)%>%
  summarise(
    regmean=mean(SDab), 
    regsme=sd(SDab)/sqrt(length(SDab)))
SDab1

#SDad
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SDad))
SDad1<- stoma %>%
  group_by(Region)%>%
  summarise(
    regmean=mean(SDad), 
    regsme=sd(SDad)/sqrt(length(SDad)))
SDad1

#SAIab
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SAIab))
SAIab1<- stoma %>%
  group_by(Region)%>%
  summarise(
    regmean=mean(SAIab), 
    regsme=sd(SAIab)/sqrt(length(SAIab)))
SAIab1

#SAIad
stoma<-read.csv("stomata_units.csv")
stoma = subset(stoma,!is.na(SAIad))
SAIad1<- stoma %>%
  group_by(Region)%>%
  summarise(
    regmean=mean(SAIad), 
    regsme=sd(SAIad)/sqrt(length(SAIad)))
SAIad1

#WUE
wue<- read.csv("Refined_WUE.csv")
wue=subset(wue, !is.na(dC13))
wue1<- wue %>%
  group_by(Region)%>%
  summarise(
    regmean=mean(dC13), 
    regsme=sd(dC13)/sqrt(length(dC13)))
wue1
