#AIM 2: Climate PC regression models and figures

setwd("D:/GEUTRI_PHYSIO/Submission Documents/Data and Code/Data")
#load data
clim <- read.csv("climatedata.csv")
pops <- read.csv("geutri_pop_means_sem_final.csv")
#merge data 
df <- merge(pops, clim, by="Population")
names(df)
#grab climate variables for pca
climmy<-df[,22:ncol(df)]
names(climmy)

#load libraries
library(ggbiplot)
library(effects)

#remove columns with NAs (MAR)
clim2<-climmy[,-25]

#do the PCA
env<-prcomp(clim2,scale.=T)
summary(env) #first 2 PCs explain 77%, with the third explainign 95% 

#plotting the PCA
ggbiplot(env, 
         ellipse=T, 
         groups=df$Region.x, 
         var.scale = 2)+
  theme_classic(16)+
  labs(color="Regions")+
  xlab('PC1 (47.0%)')+
  ylab('PC2 (29.5%)')+
  scale_color_manual(values=c("black", "gray43", "grey63"))+
  theme(legend.position = "top") ####NEED TO CHANGE LOADINGS COLOR ?



#extract loadings of the PCA to determine greatest climatic variable impacting differences in pops
#the highest absolute value represents the biggest driver (abs function)

loadings1<-as.data.frame(env$rotation[, 1:3])
str(loadings)
write.csv(loadings, "loadings.csv")
loadings<-read.csv("loadings.csv")
library(tidyverse)
library(knitr)
loadings %>%
  arrange(desc(abs(PC3)) 
  ) %>%
  knitr::kable()
#PC1: FFP, DD_18, MAT
#PC2: CMD, AHM, PAS
#PC3: MWMT, TD, MSP

#view variance explained by each PC
env_pca_var = env$sdev^2/sum(env$sdev^2)
barplot(env_pca_var, ylab = "Variance Explained", 
        xlab = "Principal Component", cex.axis = 1.5, cex.lab = 1.4)

#take data from pca and merge to data
env_pca_values = as.data.frame(env$x, stringsAsFactors = F)
df1<-cbind(df, env_pca_values)                       
write.csv(df1, "popmeans_climate.csv")

#regression analyses 
#LOAD DATA
pc<-read.csv("popmeans_climate.csv")

##GCLab----------------------------------------------
pc = subset(pc,!is.na(GCLab))
#first mod
abgc1<-lm(GCLab~PC1 + PC2, data=pc)
summary(abgc1)# r2=.413

#second mod
abgc2<-lm(GCLab~PC1, data=pc)
summary(abgc2)#r2=0.001

#third mod
abgc3<-lm(GCLab~PC2, data=pc)
summary(abgc3)#r2=.39

#GCLab AIC Values
AIC(abgc1)#45.25474
AIC(abgc2)#53.9836
AIC(abgc3)#45.09221

anova(abgc1, abgc3)

a<-ggplot(pc, aes(x=pc$PC2, y=pc$GCLab))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Abaxial guard cell length (µm)\n")+
  xlab("\nPC 2")+
  theme(legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2 = 0.39**'",x=5,y=31,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.95, .3),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

a

##GCLad----------------------------------------------
pc<-read.csv("popmeans_climate.csv")
pc = subset(pc,!is.na(GCLad))
#first mod
adgc1<-lm(GCLad~PC1 + PC2, data=pc)
summary(adgc1)#0.3524

#second mod
adgc2<-lm(GCLad~PC1, data=pc)
summary(adgc2)#-0.05146 

#third mod
adgc3<-lm(GCLad~PC2, data=pc)
summary(adgc3)#0.3838

#GCLad AIC Values
AIC(adgc1)#57.7528
AIC(adgc2)#65.63787
AIC(adgc3)#56.02054

anova(adgc1, adgc3)

b<-ggplot(pc, aes(x=pc$PC2, y=pc$GCLad))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Adaxial guard cell length (µm)\n")+
  xlab("\nPC 2")+
  theme(legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2 = 0.39**'",x=5,y=31,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.95, .3),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

b

##SDab-----------------------------------------------
pc<-read.csv("popmeans_climate.csv")
pc = subset(pc,!is.na(SDab))
#first mod
absd1<-lm(SDab~PC1 + PC2, data=pc)
summary(absd1)# r2=0.4024 

#second mod
absd2<-lm(SDab~PC1, data=pc)
summary(absd2)#r2=-0.02851 

#third mod
absd3<-lm(SDab~PC2, data=pc)
summary(absd3)#r2=0.4021

#SDab AIC Values
AIC(absd1)#181.4311
AIC(absd2)#190.3672
AIC(absd3)#180.6046


c<-ggplot(pc, aes(x=pc$PC2, y=pc$SDab))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Abaxial stomatal density (mm²)\n")+
  xlab("\nPC 2")+
  theme(legend.position = "top", legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2 = 0.40**'",x=5,y=320,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.25, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

c

##SDad-----------------------------------------------
pc<-read.csv("popmeans_climate.csv")
pc = subset(pc,!is.na(SDad))
#first mod
adsd1<-lm(SDad~PC1 + PC2, data=pc)
summary(adsd1)# r2=0.3454

#second mod
adsd2<-lm(SDad~PC1, data=pc)
summary(adsd2)#r2= -0.0421 

#third mod
adsd3<-lm(SDad~PC2, data=pc)
summary(adsd3)#r2=0.3685 

#SDad AIC Values
AIC(adsd1)#169.9535
AIC(adsd2)#177.4842
AIC(adsd3)#168.4673

d<-ggplot(pc, aes(x=pc$PC2, y=pc$SDad))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Adaxial stomatal density (mm²)\n")+
  xlab("\nPC 2")+
  theme(legend.position = "top", legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2 = 0.37**'",x=5,y=220,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.25, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))

d

##SAIab-----------------------------------------------
pc<-read.csv("popmeans_climate.csv")
pc = subset(pc,!is.na(SAIab))
#first mod
absai1<-lm(SAIab~PC1 + PC2, data=pc)
summary(absai1)# r2=0.4124 

#second mod
absai2<-lm(SAIab~PC1, data=pc)
summary(absai2)#0.006068  

#third mod
absai3<-lm(SAIab~PC2, data=pc)
summary(absai3)#r2=0.3755

#SAIab AIC Values
AIC(absai1)#36.52273
AIC(absai2)#45.14682
AIC(absai3)#36.78321

e<-ggplot(pc, aes(x=pc$PC2, y=pc$SAIab))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Abaxial stomatal area index\n")+
  xlab("\nPC 2")+
  theme(legend.position = "top", legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2 = 0.41**'",x=5,y=7.4,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.25, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))
e

##SAIad-----------------------------------------------
pc<-read.csv("popmeans_climate.csv")
pc = subset(pc,!is.na(SAIad))
#first mod
adsai1<-lm(SAIad~PC1 + PC2, data=pc)
summary(adsai1)# 0.3057   

#second mod
adsai2<-lm(SAIad~PC1, data=pc)
summary(adsai2)#-0.04452    

#third mod
adsai3<-lm(SAIad~PC2, data=pc)
summary(adsai3)#r2=0.3334  

#SAIad AIC Values
AIC(adsai1)#33.21256
AIC(adsai2)#39.72449
AIC(adsai3)#31.63938

ff<-ggplot(pc, aes(x=pc$PC2, y=pc$SAIad))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Adaxial stomatal area index\n")+
  xlab("\nPC 2")+
  theme(legend.position = "top", legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2 = 0.33**'",x=5,y=5.5,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.25, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))
ff

##WUE-----------------------------------------------
pc<-read.csv("popmeans_climate.csv")
pc = subset(pc,!is.na(WUE))
#first mod
wue1<-lm(WUE~PC1 + PC2, data=pc)
summary(wue1)# 0.3707  

#second mod
wue2<-lm(WUE~PC1, data=pc)
summary(wue2)# -0.02201  

#third mod
wue3<-lm(WUE~PC2, data=pc)
summary(wue3)#r2=0.3656  

#WUE AIC Values
AIC(wue1)#21.06286
AIC(wue2)#28.95255
AIC(wue3)#20.36933

g<-ggplot(pc, aes(x=pc$PC2, y=pc$WUE))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Carbon isotope composition\n")+
  xlab("\nPC 2")+
  theme(legend.position = "top", legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2 = 0.37**'",x=5.5,y=-28.5,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.25, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))
g

g+annotate("text",label="High WUE",x=-2,y=-28.5,parse=TRUE,size=4,colour="black")#+
  annotate("text",label="Low WUE",x=-2,y=-29,parse=TRUE,size=4,colour="black")
  

## Figures ##

a<-ggplot(pc, aes(x=pc$PC2, y=pc$GCLab))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  theme_classic(14)+
  geom_point(aes(shape=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  ylab(bquote('Abaxial guard cell length ('*mu~'m'*')'))+
  xlab("")+
  ylim(25,30.5)+
  #xlab("\nPC 2")+
  annotate("text",label="'R^2==0.39**'",x=3,y=30,parse=TRUE,size=4,colour="black")+
  annotate("text",label="A",x=-2.5,y=30,size=6,fontface="bold")+
  theme(legend.title = element_blank())
b<-ggplot(pc, aes(x=pc$PC2, y=pc$GCLad))+
  geom_smooth(method='lm',color="black", alpha=0.15)+
  theme_classic(14)+
  geom_point(aes(shape=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  ylab(bquote('Adaxial guard cell length ('*mu~'m'*')'))+
  xlab("")+
  ylim(25,31)+
  #xlab("\nPC 2")+
  annotate("text",label="'R^2==0.38**'",x=3,y=30,parse=TRUE,size=4,colour="black")+
  annotate("text",label="B",x=-2.5,y=30,size=6,fontface="bold")+
  theme(legend.title = element_blank())
c<-ggplot(pc, aes(x=pc$PC2, y=pc$SDab))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  theme_classic(14)+  
  geom_point(aes(shape=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  ylab(bquote('Abaxial stomatal density ('*mm^2*')'))+
  xlab("")+
  ylim(100,300)+
  #xlab("\nPC 2")+
  annotate("text",label="'R^2==0.40**'",x=3,y=290,parse=TRUE,size=4,colour="black")+
  annotate("text",label="C",x=-2.5,y=290,size=6,fontface="bold")+
  theme(legend.title = element_blank())

d<-ggplot(pc, aes(x=pc$PC2, y=pc$SDad))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  theme_classic(14)+
  geom_point(aes(shape=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  ylab(bquote('Adaxial stomatal density ('*mm^2*')'))+
  xlab("")+
  ylim(100,300)+
  #xlab("\nPC 2")+
  annotate("text",label="'R^2==0.37**'",x=3,y=290,parse=TRUE,size=4,colour="black")+
  annotate("text",label="D",x=-2.5,y=290,size=6,fontface="bold")+
  theme(legend.title = element_blank())

e<-ggplot(pc, aes(x=pc$PC2, y=pc$SAIab))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  theme_classic(14)+
  geom_point(aes(shape=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  ylab(bquote('Abaxial stomatal area index ('*mm^2*')'))+
  xlab("\nPC 2")+
  ylim(2.5,7.5)+  
  annotate("text",label="'R^2==0.38**'",x=3,y=7,parse=TRUE,size=4,colour="black")+
  annotate("text",label="E",x=-2.5,y=7,size=6,fontface="bold")+
  theme(legend.title = element_blank())

f<-ggplot(pc, aes(x=pc$PC2, y=pc$SAIad))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  theme_classic(14)+
  geom_point(aes(shape=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  ylab(bquote('Adaxial stomatal area index ('*mm^2*')'))+
  xlab("\nPC 2")+
  ylim(2.5,7.5)+
  annotate("text",label="'R^2==0.33**'",x=3,y=7,parse=TRUE,size=4,colour="black")+
  annotate("text",label="F",x=-2.5,y=7,size=6,fontface="bold")+
  theme(legend.title = element_blank())

ggarrange(a,b,c,d,e,f,nrow=3,ncol=2, 
          align="hv", common.legend = TRUE, 
          legend="right")

g<-ggplot(pc, aes(x=pc$PC2, y=pc$WUE))+
  geom_smooth(method='lm', color="black", alpha=0.15)+
  geom_point(aes(shape=pc$Region, size=pc$Region))+
  scale_shape_manual(values=c(0,1,2))+
  theme_classic(14) +
  ylab("Carbon isotope composition\n")+
  xlab("\nPC 2")+
  theme(legend.position = "top", legend.title=element_blank())+
  scale_size_manual(values=c(4,4,4))+
  annotate("text",label="'R^2==0.37**'",x=3.9,y=-28,parse=TRUE,size=4,colour="black")+
  theme(
    legend.position = c(.25, .08),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))
g+
  annotate("text",label="High WUE",x=-2.3,y=-28,size=4,fontface="bold")+
  annotate("text",label="Low WUE",x=-2.3,y=-31,size=4,fontface="bold")
