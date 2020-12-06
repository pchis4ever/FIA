library(lme4)
library(ggplot2)
library(car)
library(MASS)
library(survival)
library(lmtest)
library(coxme)
library(visreg)
library(lavaan)
library(data.table)
library(piecewiseSEM)
library(mgcv)
library(dplyr)
library(binr)
library(sjstats)

setwd("~/R/R files")


tree <- read.csv("SD_TREE.csv")
tree <- tree[c("CN",	"PLT_CN",	"PREV_TRE_CN",	"INVYR",	"STATECD",	"COUNTYCD",	"PLOT",	"SUBP",	"TREE",	
               "CONDID",	"AZIMUTH",	"DIST",	"PREVCOND",	"STATUSCD",	"SPCD",	"SPGRPCD",	"DIA", 
               "DAMAGE_AGENT_CD1",	"DAMAGE_AGENT_CD2",	"DAMAGE_AGENT_CD3", "MIST_CL_CD", "TOTAGE", "HT",
               "ACTUALHT", "AGENTCD", "CR",	"CCLCD", "MORTYR")]
tree <- subset(tree, COUNTYCD == 33 | COUNTYCD == 81 | COUNTYCD == 93 | COUNTYCD == 103 | COUNTYCD == 47)
tree$TreeID <- tree$TREE + tree$PLOT*1000 + tree$SUBP*100 + tree$COUNTYCD*100000000
tree$PlotID <- tree$PLOT + tree$COUNTYCD*100000000
table(unique(tree[c("SPCD", "TreeID")])$SPCD)


cond <- read.csv("SD_COND.csv")
cond <- cond[c("GSSTK", "PLT_CN",	"SICOND", "SUBPPROP_UNADJ", "MICRPROP_UNADJ", "CARBON_STANDING_DEAD", "PHYSCLCD", "INVYR",	"STATECD", "COUNTYCD",	"PLOT",	"CONDID",	"COND_STATUS_CD",	
               "RESERVCD",	"OWNCD",	"OWNGRPCD", "FORTYPCD",	"STDAGE",	"STDSZCD", "DSTRBCD1",	"DSTRBYR1",	
               "DSTRBCD2",	"DSTRBYR2",	"DSTRBCD3",	"DSTRBYR3",	"TRTCD1",	"TRTYR1", "TRTCD2",	"TRTYR2",	"TRTCD3",	
               "TRTYR3", "BALIVE", "HABTYPCD1",	"HABTYPCD1_PUB_CD", "LIVE_CANOPY_CVR_PCT", "NBR_LIVE_STEMS", "LIVE_MISSING_CANOPY_CVR_PCT")]

plot <- read.csv("SD_PLOT.csv")
plot <- plot[c("PLOT", "STATECD", "COUNTYCD", "LAT", "LON", "MEASYEAR", "INVYR", "CN", "ELEV")]
plot <- subset(plot, (COUNTYCD == 33 | COUNTYCD == 81 | COUNTYCD == 93 | COUNTYCD == 103 | COUNTYCD == 47))


subp <- read.csv("SD_SUBPLOT.csv")
subp <- subp[c("INVYR",	"STATECD",	"COUNTYCD",	"PLOT",	"SUBP", "SLOPE", "ASPECT")]


#adding plot, subplot, and condition information
tree <- merge(plot, tree, by=c("STATECD", "COUNTYCD", "PLOT", "INVYR"))
tree <- merge(subp, tree, by=c("STATECD", "COUNTYCD", "PLOT", "SUBP", "INVYR"))
tree <- merge(cond, tree, by=c("STATECD", "COUNTYCD", "PLOT", "INVYR", "CONDID", "PLT_CN"))

allpanels <- tree

write.csv(tree, "allpanels.csv")

#denote if an entry is the first or last entry for a tree by putting a "1" in the appropriate field for "First" or "Last". Do in excel and save as .csv
firstlast <- read.csv("firstlast.csv")

first <- subset(firstlast, First == 1 & STATUSCD == 1 & Last == 0)
last <- subset(firstlast, Last == 1 & First == 0)
firstlast <- merge(first, last, by=c("STATECD", "COUNTYCD", "PLOT", "SUBP", "TREE", "TreeID", "PlotID"))


#getting only entries were dead at final measurement. Excluding entries that were harvested at final remeasurement, or killed by something other than insects.
firstlast <- subset(firstlast, STATUSCD.x == 1 & (STATUSCD.y == 1 | STATUSCD.y == 2))

ponderosa <-subset(firstlast, SPCD.x == 122)


#only non-private plots
#ponderosa <- subset(ponderosa, OWNGRPCD.x != 40)

#getting only trees alive or killed by insects
ponderosa <- subset(ponderosa, STATUSCD.y ==1 | (STATUSCD.y == 2 & AGENTCD.y ==10))


#export plot coordinates
names(ponderosa)[names(ponderosa) == "LAT.x"] <- "Latitude"
names(ponderosa)[names(ponderosa) == "LON.x"] <- "Longitude"
UniquePlotList <- unique(ponderosa[c("PlotID","Latitude","Longitude")])
write.csv(UniquePlotList, "UniquePlots.csv")


#adding climate data
prism <- read.csv("PRISM_4K_Noninterpolated.csv")
#prism <- prism[-c(1,2,3,4,5,6,7,8,9,10), ]
colnames(prism) <- c("PlotID","Latitude","Longitude","ElevationPRISM","Date","Precip","Temp")
prism <- subset(prism, Date == "Annual")
prism <- prism[c('PlotID', 'Precip', 'Temp', 'ElevationPRISM')]
#climate variables are exported by PRISM as factors for some reason, so resetting to numeric
prism$Precip <- (as.numeric(as.character(prism$Precip)))
prism$Temp <- (as.numeric(as.character(prism$Temp)))
prism$ElevationPRISM <- (as.numeric(as.character(prism$ElevationPRISM)))
#visualizing relationships between climate and elevation, based on PRISM provided elevations
qplot(x= ElevationPRISM, y = Precip, data = prism)
cor.test(prism$ElevationPRISM, prism$Precip, method="pearson")
qplot(x= ElevationPRISM, y = Temp, data = prism)
cor.test(prism$ElevationPRISM, prism$Temp, method="pearson")

#merging with main dataframe
ponderosa <- merge(ponderosa, prism, by ='PlotID')


#coding mortality
ponderosa$InsectDeath <- ifelse(ponderosa$STATUSCD.y == 2, 1, 0)



#coding forest type as factor
ponderosa$FORTYPCD.x <- as.factor(ponderosa$FORTYPCD.x)

#coding physio class as a factor
ponderosa$PHYSCLCD.x <- as.numeric(ponderosa$PHYSCLCD.x)

ponderosa$Mesic <- ifelse(ponderosa$PHYSCLCD.y<30 & ponderosa$PHYSCLCD.y >20, 1, 0)
ponderosa$Xeric <- ifelse(ponderosa$PHYSCLCD.y< 20, 1, 0)
ponderosa$Hydric <- ifelse(ponderosa$PHYSCLCD.y > 30, 1, 0)

#specifying ponderosa forest vs non-ponderosa
ponderosa$PonderosaForest <- ifelse(ponderosa$FORTYPCD.x == '221', 1, 0)

#linearizing slope by converting to degrees
ponderosa$Slope <- atan(ponderosa$SLOPE.x/100) * (180/pi)

#getting basal area of ponderosa per plot
pondo <- subset(allpanels, STATUSCD == 1 & SPCD == 122 & DIA != 'NA')

#basal area of individual tree, in sq ft
pondo$BasalArea <- ((0.5*(pondo$DIA/12))^2 * pi)
pondo$SubplotBasalArea <- ifelse (pondo$DIA > 4.9, pondo$BasalArea, 0)
pondo$MicroplotBasalArea <- ifelse (pondo$DIA < 5.0, pondo$BasalArea, 0)

#sum of all basal areas, per condition
groups <- group_by(pondo, COUNTYCD, PLOT, INVYR, CONDID, SUBPPROP_UNADJ, MICRPROP_UNADJ)
pondoBA <- summarise(groups, sum(SubplotBasalArea), sum(MicroplotBasalArea))

colnames(pondoBA) <- c('COUNTYCD', 'PLOT', 'INVYR.x', 'CONDID.x', 'SUBPPROP_UNADJ', 'MICRPROP_UNADJ', 'SubplotBasalArea', 'MicroplotBasalArea')

#getting proportion of subplots/microplots occupied for each condition 
pondoBA$subplotacres <- ((pi * (24^2)) / 43560) *4
pondoBA$subplotacresadj <- pondoBA$subplotacres * pondoBA$SUBPPROP_UNADJ

pondoBA$microplotacres <- ((pi * (6.8^2)) / 43560) *4
pondoBA$microplotacresadj <- pondoBA$microplotacres * pondoBA$MICRPROP_UNADJ

pondoBA$PondoBA <- pondoBA$SubplotBasalArea/pondoBA$subplotacresadj + pondoBA$MicroplotBasalArea/pondoBA$microplotacresadj
summary(pondoBA$PondoBA)
pondoBA <- pondoBA[c('COUNTYCD', 'PLOT', 'INVYR.x', 'CONDID.x', 'PondoBA')]

#merging
ponderosa <- merge(pondoBA, ponderosa, by = c('COUNTYCD', 'PLOT', 'INVYR.x', 'CONDID.x'), all.y=T)



#getting basal area of non-ponderosa per plot
notpondo <- subset(allpanels, STATUSCD == 1 & SPCD != 122 & DIA != 'NA')

#basal area of individual tree, in sq ft
notpondo$BasalArea <- ((0.5*(notpondo$DIA/12))^2 * pi)
notpondo$SubplotBasalArea <- ifelse (notpondo$DIA > 4.9, notpondo$BasalArea, 0)
notpondo$MicroplotBasalArea <- ifelse (notpondo$DIA < 5.0, notpondo$BasalArea, 0)

#sum of all basal areas, per condition
groups <- group_by(notpondo, COUNTYCD, PLOT, INVYR, CONDID, SUBPPROP_UNADJ, MICRPROP_UNADJ)
notpondoBA <- summarise(groups, sum(SubplotBasalArea), sum(MicroplotBasalArea))

colnames(notpondoBA) <- c('COUNTYCD', 'PLOT', 'INVYR.x', 'CONDID.x', 'SUBPPROP_UNADJ', 'MICRPROP_UNADJ', 'SubplotBasalArea', 'MicroplotBasalArea')

#getting proportion of subplots/microplots occupied for each condition 
notpondoBA$subplotacres <- ((pi * (24^2)) / 43560) *4
notpondoBA$subplotacresadj <- notpondoBA$subplotacres * notpondoBA$SUBPPROP_UNADJ

notpondoBA$microplotacres <- ((pi * (6.8^2)) / 43560) *4
notpondoBA$microplotacresadj <- notpondoBA$microplotacres * notpondoBA$MICRPROP_UNADJ

notpondoBA$NotPondoBA <- notpondoBA$SubplotBasalArea/notpondoBA$subplotacresadj + notpondoBA$MicroplotBasalArea/notpondoBA$microplotacresadj
summary(notpondoBA$NotPondoBA)
notpondoBA <- notpondoBA[c('COUNTYCD', 'PLOT', 'INVYR.x', 'CONDID.x', 'NotPondoBA')]

#merging
ponderosa <- merge(notpondoBA, ponderosa, by = c('COUNTYCD', 'PLOT', 'INVYR.x', 'CONDID.x'), all.y=T)

#setting non-ponderosa basal area that was NA to 0
ponderosa$NotPondoBA[is.na(ponderosa$NotPondoBA)] <- 0

#total BA
ponderosa$BA <- ponderosa$NotPondoBA + ponderosa$PondoBA

#proportion non pondo
ponderosa$PropNotPondo <- ponderosa$NotPondoBA / ponderosa$BA

#proportion pondo
ponderosa$PropPondo <- 1 - ponderosa$PropNotPondo


#replacing 0 aspect (from plots with slope <5%) with neutral aspect of 90
ponderosa$AspectAdjusted <- ifelse(ponderosa$ASPECT.x == 0, 90, ponderosa$ASPECT.x)

#converting aspect to southness
ponderosa$Southness <- ifelse(ponderosa$AspectAdjusted > 180, 360-ponderosa$AspectAdjusted, ponderosa$AspectAdjusted)
table(ponderosa$Southness)

#heat load index
ponderosa$SlopeRadians <- ponderosa$Slope*pi/180
ponderosa$LatitudeRadians <- ponderosa$LAT.y*pi/180
ponderosa$AspectRadians <- ponderosa$ASPECT.x*pi/180
ponderosa$FoldedAspectNorth <- pi-abs(ponderosa$AspectRadians-pi)
ponderosa$FoldedAspect <- abs(pi - abs(ponderosa$AspectRadians-(5*pi/4)))
ponderosa$FoldedAspectDegrees <- ponderosa$FoldedAspect*180/pi
ponderosa$HLI <- 0.339 + 0.808*cos(ponderosa$LatitudeRadians)*cos(ponderosa$SlopeRadians) - 0.196*sin(ponderosa$LatitudeRadians)*sin(ponderosa$SlopeRadians) - 0.482*cos(ponderosa$FoldedAspect)*sin(ponderosa$SlopeRadians)
ponderosa$HLI2 <- 0.339 + 0.808*cos(ponderosa$LatitudeRadians)*cos(ponderosa$SlopeRadians) - 0.196*sin(ponderosa$LatitudeRadians)*sin(ponderosa$SlopeRadians) - 0.482*cos(ponderosa$FoldedAspectNorth)*sin(ponderosa$SlopeRadians)

#Folded aspect around a SW-NE axis, flat plots adjusted to neutral aspect of 90
ponderosa$FoldedAspectDegreesNESW_adjusted <- ifelse(ponderosa$ASPECT.x > 0,  ponderosa$FoldedAspectDegrees, 90)

#converting to metric
ponderosa$Dia <- ponderosa$DIA.x * 2.54
ponderosa$LiveBasalArea <- ponderosa$BALIVE.x * 0.229568
ponderosa$Elevation <- ponderosa$ELEV.x*0.3048
ponderosa$AvgPrecip <- ponderosa$Precip
ponderosa$AvgTemp <- ponderosa$Temp
ponderosa$PondoBasal <- ponderosa$PondoBA * 0.229568
ponderosa$NotPondoBasal <- ponderosa$NotPondoBA * 0.229568



#renaming clunky variable headers
setnames(ponderosa, old = c('LIVE_CANOPY_CVR_PCT.x', 'CCLCD.x', 'STDAGE.x', 'SICOND.x', 'CR.x', 'CARBON_STANDING_DEAD.x'), new = c('CanopyCover', 'CrownClass', 'StandAge', 'SiteIndex', 'CrownRatio', 'StandingDeadCarbon'))

#Coding disturbances- fire, weather, disease, animal, general
ponderosa$Fire <- ifelse((ponderosa$DSTRBCD1.x >29 & ponderosa$DSTRBCD1.x <33)|(ponderosa$DSTRBCD2.x >29 & ponderosa$DSTRBCD2.x <33)|(ponderosa$DSTRBCD3.x >29 & ponderosa$DSTRBCD3.x <33)|(ponderosa$DSTRBCD1.y >29 & ponderosa$DSTRBCD1.y <33)|(ponderosa$DSTRBCD2.y >29 & ponderosa$DSTRBCD2.y <33)|(ponderosa$DSTRBCD3.y >29 & ponderosa$DSTRBCD3.y <33), 1, 0)
table(ponderosa$Fire)
ponderosa$Weather <- ifelse((ponderosa$DSTRBCD1.x >49 & ponderosa$DSTRBCD1.x <55)|(ponderosa$DSTRBCD2.x >49 & ponderosa$DSTRBCD2.x <55)|(ponderosa$DSTRBCD3.x >49 & ponderosa$DSTRBCD3.x <55)|(ponderosa$DSTRBCD1.y >49 & ponderosa$DSTRBCD1.y <55)|(ponderosa$DSTRBCD2.y >49 & ponderosa$DSTRBCD2.y <55)|(ponderosa$DSTRBCD3.y >49 & ponderosa$DSTRBCD3.y <55), 1, 0)
table(ponderosa$Weather)
ponderosa$Disease <- ifelse((ponderosa$DSTRBCD1.x >19 & ponderosa$DSTRBCD1.x<23)|(ponderosa$DSTRBCD2.x >19 & ponderosa$DSTRBCD2.x<23)|(ponderosa$DSTRBCD3.x >19 & ponderosa$DSTRBCD3.x<23)|(ponderosa$DSTRBCD1.y >19 & ponderosa$DSTRBCD1.y<23)|(ponderosa$DSTRBCD2.y >19 & ponderosa$DSTRBCD2.y<23)|(ponderosa$DSTRBCD3.y >19 & ponderosa$DSTRBCD3.y<23), 1, 0)
table(ponderosa$Disease)
ponderosa$Animal <- ifelse((ponderosa$DSTRBCD1.x >39 & ponderosa$DSTRBCD1.x<47)|(ponderosa$DSTRBCD2.x >39 & ponderosa$DSTRBCD2.x<47)|(ponderosa$DSTRBCD3.x >39 & ponderosa$DSTRBCD3.x<47) | (ponderosa$DSTRBCD1.y >39 & ponderosa$DSTRBCD1.y<47)|(ponderosa$DSTRBCD2.y >39 & ponderosa$DSTRBCD2.y<47)|(ponderosa$DSTRBCD3.y >39 & ponderosa$DSTRBCD3.y<47), 1, 0)
table(ponderosa$Animal)
ponderosa$Disturbance <- ifelse(ponderosa$DSTRBCD1.x>19 | ponderosa$DSTRBCD1.y>19 | ponderosa$DSTRBCD2.y>19 | ponderosa$DSTRBCD3.y>19, 1, 0)
table(ponderosa$Disturbance)
ponderosa$Cutting <- ifelse(ponderosa$TRTCD1.x==10, 1, 0)

#total time tree observed in the study interval
pondallpanels <- subset(allpanels, SPCD == 122)
plotmeasyr1 <- unique(pondallpanels[,c('PlotID', 'MEASYEAR')])
#write.csv(plotmeasyr1, 'plotsbyyear.csv')
#denote first and last entry in excel here
plotsbyyear<-read.csv('plotsbyyear.csv')
plotsbyyearlast <- subset(plotsbyyear, LAST == 1)
plotsbyyearlast <- plotsbyyearlast[c('PlotID', 'MEASYEAR')] 

colnames(plotsbyyearlast) <- c("PlotID", "LastYearMeasured")

ponderosa <- merge(plotsbyyearlast, ponderosa, by = 'PlotID')

ponderosa$Time <- ponderosa$LastYearMeasured-ponderosa$MEASYEAR.x
ponderosa <- subset(ponderosa, STATUSCD.y ==1 | (STATUSCD.y == 2 & AGENTCD.y == 10))


#writing formatted data file
#write.csv(ponderosa, "ponderosalogisticWITHPRIVATE.csv")

#start here with formatted data file
ponderosa <- read.csv("ponderosalogisticWITHPRIVATE.csv")

ponderosa$ActualHeight <- ifelse(is.na(ponderosa$ACTUALHT.x), ponderosa$HT.x, ponderosa$ACTUALHT.x)
ponderosa$HDratio <- ponderosa$ActualHeight/ (ponderosa$DIA.x/12)

ponderosa$Broken <- ifelse(ponderosa$HT.x>ponderosa$ActualHeight, 1, 0)
sum(ponderosa$Broken)

#mean centering and scaling by sd
ponderosacenter<-ponderosa


center_scale <- function(x) {
  scale(x, scale = T)
}


ponderosacenter$Dia <- center_scale(ponderosa$Dia)
ponderosacenter$Elevation <- center_scale(ponderosa$Elevation)
ponderosacenter$LiveBasalArea <- center_scale(ponderosa$LiveBasalArea)
ponderosacenter$PondoBasal <- center_scale(ponderosa$PondoBasal)
ponderosacenter$NotPondoBasal <- center_scale(ponderosa$NotPondoBasal)
ponderosacenter$CrownClass <- center_scale(ponderosa$CrownClass)
ponderosacenter$CrownRatio <- center_scale(ponderosa$CrownRatio)
ponderosacenter$StandAge <- center_scale(ponderosa$StandAge)
ponderosacenter$SiteIndex <- center_scale(ponderosa$SiteIndex)
ponderosacenter$AvgPrecip <- center_scale(ponderosa$AvgPrecip)
ponderosacenter$AvgTemp <- center_scale(ponderosa$AvgTemp)
ponderosacenter$PropPondo <- center_scale(ponderosa$PropPondo)
ponderosacenter$FoldedAspectDegreesNESW_adjusted <- center_scale(ponderosa$FoldedAspectDegreesNESW_adjusted)
ponderosacenter$Slope <- center_scale(ponderosa$Slope)
ponderosacenter$Southness <- center_scale(ponderosa$Southness)
ponderosacenter$Time <- center_scale(ponderosa$Time)
ponderosacenter$HLI <- center_scale(ponderosa$HLI)
ponderosacenter$StandSize <- center_scale(ponderosa$STDSZCD.x)
ponderosacenter$HD_ratio <- center_scale(ponderosa$HDratio)
ponderosacenter$Broken <- (ponderosa$Broken)
ponderosacenter$Height <- center_scale(ponderosa$HT.x)



#logistic regression
survive <- glmer (InsectDeath  ~  Time + StandAge + SiteIndex + HD_ratio + FoldedAspectDegreesNESW_adjusted + Slope + CrownClass + HLI + PondoBasal + NotPondoBasal + Dia + CrownRatio + Elevation
                  + (1 | PlotID), data = ponderosacenter, family=binomial, control=glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=80000)))
summary(survive)
vif(survive)

survive <- glmer (InsectDeath  ~  Time + PondoBasal + NotPondoBasal + Dia + CrownRatio + Elevation
                  + (1 | PlotID), data = ponderosacenter, family=binomial, control=glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=80000)))
summary(survive)
vif(survive)

#add climate
survive <- glmer (InsectDeath  ~  AvgTemp + AvgPrecip + PondoBasal + NotPondoBasal + Dia + CrownRatio + Elevation
                  + (1 | PlotID) + Time, data = ponderosacenter, family=binomial, control=glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=80000)))
summary(survive)
vif(survive)

survive <- glmer (InsectDeath  ~  Dia + Time + AvgPrecip + PondoBasal + NotPondoBasal + CrownRatio + Elevation
                  + (1 | PlotID), data = ponderosacenter, family=binomial, control=glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=100000)))
summary(survive)
vif(survive)

visreg(survive, 'CrownRatio', xlab= 'Crown Ratio (%)', type = 'contrast', scale = 'linear', ylab = expression(paste(Delta,'log odds (mortality)')), band = F)
visreg(survive, 'AvgPrecip', scale = 'linear', type = 'contrast', xlab= 'Average Annual Precipitation (cm)', ylab = expression(paste(Delta,'log odds (mortality)')), band = F)
visreg(survive, 'Elevation', scale = 'linear', type = 'contrast', xlab= 'Elevation (m)', ylab = expression(paste(Delta,'log odds (mortality)')), band = F)
visreg(survive, 'PondoBasal', type = 'contrast', scale = 'linear', xlab= expression(paste("Ponderosa Basal Area (", cm^2/ha, ")")), ylab = expression(paste(Delta,'log odds (mortality)')), band = F)
visreg(survive, 'NotPondoBasal', type = 'contrast', scale = 'linear', xlab= expression(paste("Non-Ponderosa Basal Area (", cm^2/ha, ")")), ylab = expression(paste(Delta,'log odds (mortality)')), band = F)
visreg(survive, 'Dia', type = 'contrast', scale = 'linear', xlab= 'Diameter (cm)', ylab = expression(paste(Delta,'log odds (mortality)')), band = F)

sd(ponderosa$Dia)
mean(ponderosa$Dia)
sd(ponderosa$AvgPrecip)
mean(ponderosa$AvgPrecip)
sd(ponderosa$CrownRatio)
mean(ponderosa$CrownRatio)
sd(ponderosa$Elevation)
mean(ponderosa$Elevation)
sd(ponderosa$PondoBasal)
mean(ponderosa$PondoBasal)
sd(ponderosa$NotPondoBasal)
mean(ponderosa$NotPondoBasal)



#INTERACTION MODEL
survive <- glmer (InsectDeath  ~  Time + AvgPrecip* (StandAge + LiveBasalArea + Dia + CrownRatio + Elevation) 
                  + (1 | PlotID), data = ponderosacenter, family=binomial, control=glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=80000)))
summary(survive)
vif(survive)

survive <- glmer (InsectDeath  ~  Dia + PondoBasal + NotPondoBasal + CrownRatio + AvgPrecip * (Elevation) 
                  + (1 | PlotID) + Time, data = ponderosacenter, family=binomial, control=glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=80000)))
summary(survive)
vif(survive)





#VISUALIZATIONS
visreg(survive, 'Elevation', band = F, ylab = expression(paste(Delta,'log odds (mortality)')), by = 'AvgPrecip', breaks = c(-1,1), layout=c(2,1), type = 'contrast', xlab = "Elevation (m)")
visreg(survive, 'FoldedAspectDegreesNESW_adjusted', band = F, ylab = expression(paste(Delta,'log odds (mortality)')), by = 'AvgPrecip',breaks = 2, layout=c(3,1), type='contrast', xlab = expression(paste("Non-Ponderosa Basal Area (", m^2/ha, ")")))


visreg(survive, 'LiveBasalArea', type = 'contrast', scale = 'linear')
visreg(survive, 'CrownRatio', by = 'AvgPrecip', breaks = 4)

visreg(survive, 'AvgPrecip', by = 'LogDia', type = 'contrast', breaks = 3, scale= 'response')
visreg(survive, 'Elevation', by = 'AvgPrecip', breaks = 2, type = 'contrast', scale = 'response',band = F)
vif(survive)

#binning

ponderosabig <- ponderosa

ponderosabig$ElevBins <- cut(ponderosabig$ELEV.x, breaks=c(0,5000,7000), include.lowest=TRUE)
ponderosabig$PrecipBins <- cut(ponderosabig$AvgPrecip, breaks=c(0,1400,1600,2000), include.lowest=TRUE)

ponderosabig$ElevBins <- cut(ponderosabig$ELEV.x, breaks=c(0,5500,7000), include.lowest=TRUE)
ponderosabig$PrecipBins <- cut(ponderosabig$AvgPrecip, breaks=c(0,650,900), include.lowest=TRUE)

ponderosabig %>%
  group_by(ElevBins, PrecipBins) %>%
  summarize(PropKilled = mean(InsectDeath, na.rm = TRUE), N = sum(STATUSCD.x), Killed = sum(InsectDeath), BA = mean(LiveBasalArea), DBH = mean(DIA.x))


cuts <- bins(ponderosa$ELEV.x, target.bins = 3, minpts = 500)
cuts$breaks <- bins.getvals(cuts)
cuts$binct

cuts <- bins(ponderosa$AvgPrecip, target.bins = 2, minpts = 1500)
cuts$breaks <- bins.getvals(cuts)
cuts$binct

#confidence intervals
cc <- confint(survive,parm="beta_")
ctab <- cbind(est=fixef(survive),cc)
ctab <- as.data.frame(ctab)

ctab[2,1] <- ctab[2,1]/sd(ponderosa$PondoBasal)
ctab[2,2] <- ctab[2,2]/sd(ponderosa$PondoBasal)
ctab[2,3] <- ctab[2,3]/sd(ponderosa$PondoBasal)

ctab[3,1] <- ctab[3,1]/sd(ponderosa$NotPondoBasal)
ctab[3,2] <- ctab[3,2]/sd(ponderosa$NotPondoBasal)
ctab[3,3] <- ctab[3,3]/sd(ponderosa$NotPondoBasal)

ctab[4,1] <- ctab[4,1]/sd(ponderosa$Dia)
ctab[4,2] <- ctab[4,2]/sd(ponderosa$Dia)
ctab[4,3] <- ctab[4,3]/sd(ponderosa$Dia)

ctab[5,1] <- ctab[5,1]/sd(ponderosa$CrownRatio)
ctab[5,2] <- ctab[5,2]/sd(ponderosa$CrownRatio)
ctab[5,3] <- ctab[5,3]/sd(ponderosa$CrownRatio)

ctab[6,1] <- ctab[6,1]/sd(ponderosa$Elevation)
ctab[6,2] <- ctab[6,2]/sd(ponderosa$Elevation)
ctab[6,3] <- ctab[6,3]/sd(ponderosa$Elevation)

ctab[7,1] <- ctab[7,1]/sd(ponderosa$Time)
ctab[7,2] <- ctab[7,2]/sd(ponderosa$Time)
ctab[7,3] <- ctab[7,3]/sd(ponderosa$Time)










#GAMs
survive2 <- gam (InsectDeath  ~ s(PropPondo, sp = 0.8) + s(FoldedAspectDegreesNESW_adjusted, sp = 0.8) + s(SiteIndex, sp = 0.8) + s(LiveBasalArea, sp=0.8) + s(Dia, sp= 0.8) + s(CrownRatio, sp = 0.8) + s(AvgPrecip, sp=0.8) + s(Elevation, sp= 0.8)  
                    + s(PlotID, bs= "re") + s(MEASYEAR.x, bs = "re"), data = ponderosa, family="binomial", link = "logit")
summary(survive2)
plot(survive2, trans=plogis)
visreg(survive2, "Dia", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)', xlab = 'Dia', breaks=2)
visreg(survive2, "CrownRatio", by = "AvgPrecip", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)', xlab = 'Crown Ratio', breaks=2)
visreg(survive2, "SiteIndex", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)')
visreg(survive2, "FoldedAspectDegreesNESW_adjusted", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)')
visreg(survive2, "LiveBasalArea", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)')
visreg(survive2, "Elevation", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)')
visreg(survive2, "AvgPrecip", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)')
visreg(survive2, "PropPondo", scale = "linear", band = F, partial = T, ylab ='Log odds (insect death)')




survive1 <- gam (InsectDeath  ~ Time  + s(SiteIndex) + s(LiveBasalArea) + s(Dia) + s(CrownRatio) + s(AvgPrecip) + s(Elevation)  
                 + s(PlotID, bs= "re") + s(MEASYEAR.x, bs="re"), data = ponderosa, select = TRUE, method = "REML", family="binomial", link = "logit")
summary(survive1)
plot(survive1, trans=plogis)


survive2 <- gam (InsectDeath  ~ s(SiteIndex, sp = 0.6) + s(LiveBasalArea, sp = 0.6) + s(Dia, sp = 0.6) + s(CrownRatio, sp = 1) + s(AvgPrecip, sp = 0.6) + s(Elevation, sp = 0.6) 
                 + s(PlotID, MEASYEAR.x, Time, bs="re"), data = ponderosa, family="binomial", link = "logit")
summary(survive2)
plot(survive2, trans=plogis)

ponderosapublic <- subset(ponderosa, OWNGRPCD.x != 40)
GAMsurvive <- gam (InsectDeath  ~ Time + s(PropNotPondo, sp = 0.8) + s(HLI, sp = 0.8) + s(AspectAdjusted, sp=0.8) + s(LiveBasalArea, sp=0.8) + s(Dia, sp=0.8) + s(CrownRatio, sp=0.8) + s(Elevation, sp=0.8) + s(AvgPrecip, sp = 0.8)
                    + s(PlotID, bs = "re") , data = ponderosa, family="binomial", link="logit")
summary(GAMsurvive)

#interactions with precip
GAMsurvive2 <- gam (InsectDeath  ~ Time + ti(AvgPrecip, PropNotPondo)+ti(AvgPrecip, HLI)+ ti(AvgPrecip, Dia)+ ti(AvgPrecip, Elevation) +ti(AvgPrecip, CrownRatio)+ ti(AvgPrecip, LiveBasalArea) + s(PropNotPondo, sp = 1.2) + s(HLI, sp = 1.2) + s(LiveBasalArea, sp=1.2) + s(Dia, sp=1.2) + s(CrownRatio, sp=1.2) + s(Elevation, sp=1.2) + s(AvgPrecip, sp = 1.2)
                   + s(PlotID, bs = "re") , data = ponderosapublic, family="binomial", link="logit")
summary(GAMsurvive2)

visreg(GAMsurvive2, "CrownRatio", by = "AvgPrecip", scale = "linear", band = T, partial = T, ylab ='Log odds (insect death)', xlab = 'Crown Ratio', breaks=2)
visreg(GAMsurvive2, "HLI", by = "AvgPrecip", scale = "linear", band = T, partial = T, ylab ='Log odds (insect death)', xlab = 'Crown Ratio', breaks=2)
visreg(GAMsurvive2, "LiveBasalArea", by = "AvgPrecip", scale = "linear", band = T, partial = T, ylab ='Log odds (insect death)', xlab = 'Crown Ratio', breaks=2)
visreg(GAMsurvive2, "Elevation", by = "AvgPrecip", scale = "linear", band = T, partial = T, ylab ='Log odds (insect death)', xlab = 'Crown Ratio', breaks=2)
visreg(GAMsurvive2, "Dia", by = "AvgPrecip", scale = "linear", band = T, partial = T, ylab ='Log odds (insect death)', xlab = 'Crown Ratio', breaks=2)
visreg(GAMsurvive2, "PropNotPondo", by = "AvgPrecip", scale = "linear", band = T, partial = T, ylab ='Log odds (insect death)', xlab = 'Crown Ratio', breaks=2)
visreg(GAMsurvive2, "PropNotPondo")
visreg(GAMsurvive2, "CrownRatio")
visreg(GAMsurvive2, "HLI")
visreg(GAMsurvive2, "LiveBasalArea")
visreg(GAMsurvive2, "Dia")
visreg(GAMsurvive2, "Elevation")
visreg(GAMsurvive2, "AvgPrecip")




#interactions
survive <- gam (InsectDeath  ~ s(SiteIndex, sp = 0.8) + s(PropNotPondo, sp = 0.8) + s(HLI, sp = 0.8) + s(LiveBasalArea, sp=0.8) + s(Dia, sp=0.8) + s(CrownRatio, sp=0.8) + s(Elevation, sp=0.8) + s(AvgPrecip, sp = 0.8) + ti(SiteIndex, AvgPrecip) + ti(AvgPrecip, HLI) + ti (AvgPrecip, LiveBasalArea) + ti(AvgPrecip, Dia) + ti(AvgPrecip, Elevation) + ti(AvgPrecip, PropNotPondo) + ti(AvgPrecip, CrownRatio)
                + s(PlotID, bs = "re") + s(MEASYEAR.x, bs = "re"), data = ponderosa, family="binomial", link="logit", method = "REML")
summary(survive)
plot(survive, trans=plogis)



#adding different climate data
climate <- read.csv('Chisholm_SD_PRISM2.csv')
climate$CN <- (as.numeric(as.character(climate$CN)))

colnames(climate)[which(names(climate) == "prism_tmean")] <- "Temp"
colnames(climate)[which(names(climate) == "PRISM_ppt")] <- "Precip"
colnames(climate)[which(names(climate) == "ELEV")] <- "ELEV2"

plot2 <- merge(climate,plot, by = 'CN')
plot2 <- plot2[c("Temp", "Precip", "PLOT", "COUNTYCD")]
plot2$Precip <- (as.numeric(as.character(plot2$Precip)))
plot2$Temp <- (as.numeric(as.character(plot2$Temp)))
plot2$PlotID <- plot2$PLOT + plot2$COUNTYCD*100000000
plot3 <- plot2[!duplicated(plot2$PlotID), ]

ponderosa2 <- merge(plot3, ponderosa, by =c('PlotID', 'PLOT', 'COUNTYCD'), all= F)

ponderosa <- ponderosa2

ponderosa3 <-anti_join(ponderosa, plot3, by =c('PlotID', 'PLOT', 'COUNTYCD'), all = F)



ponderosa %>%
  group_by(OWN) %>%
  summarize(PropKilled = mean(InsectDeath, na.rm = TRUE), N = sum(STATUSCD.x), Killed = sum(InsectDeath), BA = mean(LiveBasalArea), DBH = mean(DIA.x))

