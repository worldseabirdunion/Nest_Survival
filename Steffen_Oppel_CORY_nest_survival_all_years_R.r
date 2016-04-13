#################################################################################
###  Analysis of Cory's Shearwater nest survival rates using RMARK   ############
#################################################################################

## Nest monitoring data from 2009 - 2011
## Analysis of weekly nest survival over 21 encounter occassions

## PLEASE CITE:
##Hervias, S., Henriques, A., Oliveira, N., Pipa, T., Cowen, H., Ramos, J., Nogales, M.,
##Geraldes, P., Silva, C., Ruiz de Ybanez, R., Oppel, S., 2013. Studying the effects of 
##multiple invasive mammals on Cory's shearwater nest survival. Biological Invasions 15: 143-155.





setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Corvo\\Science\\Analysis\\CORY_Nest_survival")
setwd("C:\\STEFFEN\\RSPB\\Corvo\\Science\\Analysis\\CORY_Nest_survival")
library(RMark)

## load data from summary table

CORY <- read.table("Summary_NestSurvival_allyears.csv", sep=',', header=T) # this imports the data from the query to an object in R
dim(CORY)
head(CORY)

## extract the names of the colonies so you can match them
cols<-levels(CORY$site)


## determine number of missing data for each year
missing_data<-data.frame(year=c(2009,2010,2011))
for (y in 2009:2011){
CORYy<-subset(CORY, year==y)
for (i in 4:31){
x<-CORYy[!is.na(CORYy[,i]),]
missing_data[missing_data$year==y,i-2]<-(dim(CORYy)[1])-(dim(x)[1])
}
}
names(missing_data)<-names(CORY[3:31])
missing_data[,c(1,23:26)]
fix(CORY)
head(CORY)


#############################################################################################################################
####################### DATA EXPLORATION ##############################################################################
#############################################################################################################################

### Test for correlation among predictor variables

#### CORRELATION OF ENVIRONMENTAL VARIABLES #######

cor.matrix<-as.matrix(CORY[,c(24:28,32,33)])
cor.matrix<-rcorr(cor.matrix, type="pearson")
P.val<-p.adjust(cor.matrix$P, method = 'bonferroni')

write.table(cor.matrix$r,"Correlation_matrix.csv", sep=",", row.names=F, quote=F)
write.table(P.val,"Correlation_matrix_p_values.csv", sep=",", row.names=F, quote=F)





#############################################################################################################################
####################### DATA EXPLORATION ##############################################################################
#############################################################################################################################

### the approach below takes a long time, once run you can load the results and plot them here:
setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\Corvo\\Science\\Analysis\\CORY_Nest_survival")
setwd("C:\\STEFFEN\\RSPB\\Corvo\\Science\\Analysis\\CORY_Nest_survival")
RFvarImp<-read.table('CORY_nest_success_variable_importance_RandomForest.csv', header=T, sep=',')
par(mar=c(8,4,1,1))
barplot(RFvarImp$relative_importance, names.arg=RFvarImp$variable, las=2, ylab="relative variable importance")



##########################################
#### THE CODE RF.MODELSEL.R IS FOUND IN THE GITHUB REPOSITORY: www.github.com/worldseabirdunion/Modeling.git
############################################################################################################

###### RANDOM FOREST MODEL SELECTION based on Murphy et al. (2010), Ecology 91:252
library(randomForest)
source('C:/STEFFEN/RSPB/Statistics/rf.modelSel.R')  ######### Found at www.github.com/worldseabirdunion/Modeling.git
source('A:/RSPB/Statistics/rf.modelSel.R')
rf.classif <- rf.modelSel(CORY[,c(3,8:15,22:33)], factor(CORY[,7]), imp.scale="mir")
rf.classif$IMPORTANCE

varimp<-data.frame()
for (i in 1:dim(CORY)[1]){
CORYboot<-CORY[-i,]
rf.regress <- rf.modelSel(CORYboot[,c(3,8:15,22:33)], CORYboot[,7], imp.scale="se")
imp<-rf.regress$IMPORTANCE
imp$RUN<-i
imp$variable<-row.names(imp)
varimp<-rbind(varimp,imp, deparse.level = 0)
}

VarSelect<-aggregate(importance~variable, data=varimp, FUN='mean')
VarSelect$sum<-aggregate(importance~variable, data=varimp, FUN='sum')[,2]
VarSelect$mean_importance<-VarSelect$sum/i
VarSelect<-VarSelect[order(VarSelect$mean_importance, decreasing=T),] 
VarSelect$relative_importance<-(VarSelect$mean_importance/VarSelect$mean_importance[1])*100
VarSelect

write.table(VarSelect,"CORY_nest_success_variable_importance_RandomForest.csv", sep=',', row.names=F)


### visualise the effects of rats on nest success
RF<-randomForest(fate~year+site+rat7+chamber+elevation+habitat+rat5+rat6+rat8+rat9+rat10+nest_height+nest_lenght+substrate+walls+proximate_village, data=CORY, ntree=1500, mtry=4)
par(mfrow=c(3,2))
partialPlot(RF, CORY, elevation)
partialPlot(RF, CORY, site)
partialPlot(RF, CORY, rat7)
partialPlot(RF, CORY, chamber)
partialPlot(RF, CORY, habitat)
partialPlot(RF, CORY, nest_height)




#############################################################################################################################
############## DATA PROCESSING TO MEET INPUT REQUIREMENTS OF RMARK  ##########################################################
#############################################################################################################################


CORY$Year<-CORY$year
CORY$year<-NULL
CORY$ID<-NULL
names(CORY)[2:6]<-c('FirstFound','LastPresent','LastChecked','Fate','AgeDay1')




# CREATING time-dependent covariates for "Incubation", and "Rats"

create.stage.var<-function(data,agevar.name,stagevar.name,time.intervals,cutoff)
{
nocc=length(time.intervals)
age.mat=matrix(data[,agevar.name],nrow=dim(data)[1],ncol=nocc-1)
age.mat=t(t(age.mat)+cumsum(c(0,time.intervals[1:(nocc-2)])))
stage.mat=t(apply(age.mat,1,function(x) as.numeric(x<cutoff)))
stage.mat=data.frame(stage.mat)
names(stage.mat)=paste(stagevar.name,1:(nocc-1),sep="")
return(stage.mat)
}
CORY$AgeDay1<-1-CORY$FirstFound
x<-create.stage.var(CORY,"AgeDay1","Incubation",rep(1,21),9)
CORY<-cbind(CORY,x)



CORY$Rat1<-CORY$rat5
CORY$Rat2<-CORY$rat6
CORY$Rat3<-CORY$rat6
CORY$Rat4<-CORY$rat6
CORY$Rat5<-CORY$rat6
CORY$Rat6<-CORY$rat7
CORY$Rat7<-CORY$rat7
CORY$Rat8<-CORY$rat7
CORY$Rat9<-CORY$rat7
CORY$Rat10<-CORY$rat8
CORY$Rat11<-CORY$rat8
CORY$Rat12<-CORY$rat8
CORY$Rat13<-CORY$rat8
CORY$Rat14<-CORY$rat9
CORY$Rat15<-CORY$rat9
CORY$Rat16<-CORY$rat9
CORY$Rat17<-CORY$rat9
CORY$Rat18<-CORY$rat10
CORY$Rat19<-CORY$rat10
CORY$Rat20<-CORY$rat10
CORY$Rat21<-CORY$rat10


CORY$Mouse1<-CORY$mouse5
CORY$Mouse2<-CORY$mouse6
CORY$Mouse3<-CORY$mouse6
CORY$Mouse4<-CORY$mouse6
CORY$Mouse5<-CORY$mouse6
CORY$Mouse6<-CORY$mouse7
CORY$Mouse7<-CORY$mouse7
CORY$Mouse8<-CORY$mouse7
CORY$Mouse9<-CORY$mouse7
CORY$Mouse10<-CORY$mouse8
CORY$Mouse11<-CORY$mouse8
CORY$Mouse12<-CORY$mouse8
CORY$Mouse13<-CORY$mouse8
CORY$Mouse14<-CORY$mouse9
CORY$Mouse15<-CORY$mouse9
CORY$Mouse16<-CORY$mouse9
CORY$Mouse17<-CORY$mouse9
CORY$Mouse18<-CORY$mouse10
CORY$Mouse19<-CORY$mouse10
CORY$Mouse20<-CORY$mouse10
CORY$Mouse21<-CORY$mouse10


### reformat variables
## to analyse groups you create an indicator variable for n-1 levels of the group
levels(CORY$site)
CORY$name_pico<-ifelse(CORY$site=="Cancela do Pico",1,0)
CORY$name_faja<-ifelse(CORY$site=="Faj? Rochosa",1,0)
CORY$name_velha<-ifelse(CORY$site=="Fonte Velha",1,0)
CORY$name_miradouro<-ifelse(CORY$site=="Miradouro do Portal",1,0)
CORY$name_acucar<-ifelse(CORY$site=="P?o de A??car",1,0)
CORY$site<-NULL

levels(CORY$habitat)
CORY$hab_arundo<-ifelse(CORY$habitat=="Arundo donax",1,0)
CORY$hab_coast<-ifelse(CORY$habitat=="coast with rocks",1,0)
CORY$hab_pasture<-ifelse(CORY$habitat=="grazing land",1,0)
CORY$hab_road<-ifelse(CORY$habitat=="road",1,0)
CORY$habitat<-NULL

CORY$year2009<-ifelse(CORY$Year==2009,1,0)
CORY$year2010<-ifelse(CORY$Year==2010,1,0)
CORY$Year<-NULL

CORY $ curved_entrance<-ifelse(CORY $ curved_entrance=="Yes",1,0)
CORY $ chamber<-ifelse(CORY $ chamber=="yes",1,0)
CORY $ walls<-ifelse(CORY $ walls=="yes",1,0)
CORY $ substrate<-ifelse(CORY $ substrate=="rock",1,0)


######## REMOVE THE COLUMNS NO LONGER NEEDED ##########
for (i in 1:12){
CORY[,8]<-NULL
}




### LOOK AT THE DATA YOU CREATED
str(CORY)
head(CORY) # This command means you can see the first six line
summary(CORY) # This command show you mean, median, min., max. for each covariable





#############################################################################################################################
####################### NEST SURVIVAL ANALYSIS ##############################################################################
#############################################################################################################################


########## FOR ALL YEARS COMBINED ###################################
########## ignoring cat measurements, because not measured in 2009 ##

CORY.process<-process.data(CORY, model="Nest", begin.time=1, nocc=21)
CORY.ddl<-make.design.data(CORY.process)



########## STEP 1: determine the time structure for further analysis ###################################

run.CORY<-function()
{
# 0. Constant weekly survival rate (DSR) - the null model
m0<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~1)))

# 1. Linear temporal trend
m1<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~Time)))

# 2. Quadratic temporal trend
m2<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~I(Time^2))))

# 3. Nest survival varies every week
m3<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time)))

# 4. Nest survival varies between incubation and chick rearing
m4<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~Incubation)))

# 5. Constant weekly survival rate (DSR) but different among years
m5<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~year2009+year2010)))

# 6. Linear temporal trend varying in each year
m6<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~Time*year2009+Time*year2010)))

# 7. Quadratic temporal trend varying in each year
m7<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~I(Time^2)*year2009+I(Time^2)*year2010)))

# 8. Nest survival varies every week and in each year
m8<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010)))

# 9. Nest survival varies between incubation and chick rearing, differently in each year
m9<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~Incubation*year2009+Incubation*year2010)))

return(collect.models() )
}

CORY.results<-run.CORY()  
CORY_AIC_table<-print.marklist(CORY.results)
NULL_DEVIANCE<-CORY_AIC_table$Deviance[3]


### overwhelming evidence for a time*year structure!


########## STEP 2: evaluating a suite of candidate models for all years ###################################
########## using the temporal structure identified in Step 1: Nest survival varies every week in every year

run.CORY<-function()
{
# 0. Nest survival is different every week - the null model
m0<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010)))  # this command is to run a very simple model just to check any error

# 1. Nest survival follows relationship with elevation and distance to village
m1_elevation<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~elevation+proximate_village+time*year2009+time*year2010)))

# 2. Nest survival varies among colonies
m2_colony<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_pico+name_faja+name_velha+name_miradouro+name_acucar)))

# 3. Nest survival varies among habitats
m3_habitat<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+hab_arundo+hab_coast+hab_pasture+hab_road)))

# 5. Nest survival varies with burrow characteristics that may limit access for predators
m5_burrow<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+nest_height+chamber)))

# 6. Nest survival varies with soil and vegetation characteristics
m6_substrate<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+substrate+vegetation_cover)))


############ the following 3 models were removed because we include a different analysis to test for the effect of burrow characteristics
# 7. Nest survival varies with colony, and a nest specific parameter that limits access for predators (chamber)
#m7_colony_chamber<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_pico+name_faja+name_velha+name_miradouro+name_acucar+chamber)))

# 8. Nest survival varies with colony, and a nest specific parameter that limits access for predators (walls)
#m8_colony_height<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_pico+name_faja+name_velha+name_miradouro+name_acucar+walls)))

# 9. Nest survival varies with colony, and a nest specific parameter that limits access for predators (curved_entrance)
#m9_colony_length<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_pico+name_faja+name_velha+name_miradouro+name_acucar+curved_entrance)))
############



# 10. Nest survival varies with rat abundance in each month
m10_rat<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~Rat+time*year2009+time*year2010)))

# 11. Nest survival varies with mouse abundance in each month
m11_mouse<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~Mouse+time*year2009+time*year2010)))

# 12. Nest survival varies by colony and Rat abundance in each month
m12_rat_colony<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+Rat+name_pico+name_faja+name_velha+name_miradouro+name_acucar)))

# 13. Nest survival varies by colony and Mouse abundance in each month
m13_mouse_colony<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+Mouse+name_pico+name_faja+name_velha+name_miradouro+name_acucar)))

# 15. Nest survival varies by July rat abundance
#m15_July_rat<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+Rat7)))

# 16. Nest survival varies by rat abundance and between Pico and other colonies
# Pico has highest rat number, because old farm with lots of fruit trees  and few cats -> 
#m16_pico_rat<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_pico + Rat)))

# 17. Nest survival varies by rat abundance and between Fonte Velha and other colonies
# has strangely high hatching rate, but low fledging success, so lot of chick predation, but no egg predation -> has medium cat activity
#m17_velha_rat<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_velha + Rat)))

return(collect.models() )
}

CORY.results<-run.CORY()  
CORY_AIC_table<-print.marklist(CORY.results)


### estimate the proportion of explained variation
m2_colony<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_pico+name_faja+name_velha+name_miradouro+name_acucar)))
m_saturated<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time*year2009+time*year2010+name_pico+name_faja+name_velha+name_miradouro+name_acucar+chamber+Rat+Mouse+substrate+vegetation_cover+nest_width+nest_height+nest_lenght+curved_entrance+hab_arundo+hab_coast+hab_pasture+hab_road)))
null<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~1)))
VarExpl<-(-(m2_colony$results$lnl) + null$results$lnl)/(-(m_saturated$results$lnl) + null$results$lnl)   ## R^2 according to http://www.phidot.org/forum/viewtopic.php?f=1&t=168&p=292&hilit=%25+variation+explained#p292
VarExpl




########## STEP 3: summarising the results and drawing a plot ###################################
########## requires re-calculation of models without the year interaction to reduce variability in parameter estimates


run.CORY<-function()
{
m2_colony<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time+name_pico+name_faja+name_velha+name_miradouro+name_acucar)))
m3_habitat<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time+hab_arundo+hab_coast+hab_pasture+hab_road)))
m12_rat_colony<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time+Rat+name_pico+name_faja+name_velha+name_miradouro+name_acucar)))
m13_mouse_colony<-mark(CORY,nocc=21,model="Nest",model.parameters=list(S=list(formula=~time+Mouse+name_pico+name_faja+name_velha+name_miradouro+name_acucar)))
return(collect.models() )
}

CORY.results<-run.CORY()  

## model average parameter estimates
CORY_averaged_parameter_estimates<-model.average(CORY.results, vcv=TRUE)
b<-CORY_averaged_parameter_estimates$estimates
b$lcl[c(1,19)]<-1   ## the weeks with no recorded losses are set to 1

## plot weekly nest survival
library(Hmisc)
pdf("Cory_weekly_nest_survival_all_years.pdf", width=10, height=6)
errbar(b$par.index, b$estimate, b$ucl, b$lcl, cap=0.015, xlab="week of the breeding season", ylab="weekly nest survival probability",add=FALSE, lty=1, type='b', ylim=c(0.7, 1.01),lwd=1, pch=16, cex.lab=1.5,axes=F)
#axis(1, at=c(1:20), labels=c("29 May","5 Jun", "12 Jun","19 Jun","26 Jun","3 Jul", "10 Jul","17 Jul", "24 Jul", "31 Jul", "7 Aug", "14 Aug","21 Aug","28 Aug", "4 Sep", "11 Sep", "18 Sep", "25 Sep", "2 Oct", "9 Oct"), tck=-0.02, cex.axis=1.3)
axis(1, at=c(1:20), labels=c(1:20), tck=-0.02, cex.axis=1.3)
axis(2, at=c(0.7,0.8,0.9,1.0), tck=-0.02, labels=T, cex.axis=1.3, las=1)
#title(main="Cory's Shearwater Nest Survival on Corvo 2009 - 2011", sub = NULL, xlab = NULL, ylab = NULL, line = NA, outer = FALSE)
dev.off()




## Total nest survival estimates:
prod(b$estimate[1:18])  # if they start in the first week
prod(b$estimate[2:19])  # if they start in the second week
prod(b$estimate[3:20])  # if they start in the third week

prod(b$lcl[c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20)])
prod(b$ucl[c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20)])

## export results tables
write.table(CORY_AIC_table,"CORY_nest_survival_AIC_table_all_years.csv", sep=",", dec=".",row.names=F)
write.table(CORY_averaged_parameter_estimates$estimates,"CORY_nest_survival_averaged_parm_estimates_all_years.csv", sep=",",dec=".", row.names=F)



#to clean up all files created by MARK but not necessary to hold output, without asking for permission
rm(CORY.results, m2_colony, null, m_saturated, m7_colony_chamber)
cleanup(ask=FALSE)





############## simple summary ############
CORY <- read.table("Summary_NestSurvival_allyears.csv", sep=',', header=T) # this imports the data from the query to an object in R
CORY<-CORY[CORY$site!="Estrada",]
CORY$succ<-ifelse(CORY$fate==0,1,0)
CORY$egg_succ<-ifelse(CORY$lastPresent>7,1,0)
summary<-aggregate(succ~site, data=CORY, FUN='mean')
summary$breeding_success_sd<-aggregate(succ~site, data=CORY, FUN='sd')[,2]
summary$egg_success<-aggregate(egg_succ~site, data=CORY, FUN='mean')[,2]
summary$egg_success_sd<-aggregate(egg_succ~site, data=CORY, FUN='sd')[,2]
summary


