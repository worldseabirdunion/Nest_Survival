############################################################################################################
#######  NEST SURVIVAL OF MABO FROM ASCENSION AND ST HELENA  ###################################
############################################################################################################
## last updated by steffen oppel 8 July 2014

## PLEASE CITE:
##Oppel, S., Beard, A., Fox, D., Mackley, E., Leat, E., Henry, L., Clingham, E., Fowler, N., Sim, J., 
##Sommerfeld, J., Weber, N., Weber, S., Bolton, M., 2015. Foraging distribution of a tropical seabird 
##supports Ashmole's hypothesis of population regulation. Behavioral Ecology and Sociobiology 69: 915-926.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES AND CUSTOM SCRIPTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(RODBC)
library(RMark)
library(lme4)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD NEST SUCCESS DATA FROM ST HELENA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Load data from database
setwd("S:\\ConSci\\DptShare\\SteffenOppel/RSPB/UKOT/StHelena/Science/Birds/seabirds/Ringing")
setwd("C:/STEFFEN/RSPB/UKOT/StHelena/Science/Birds/seabirds/Ringing")
db <- odbcConnectAccess2007('StHelena_ringing_database_v2.accdb')
SHnests <- sqlQuery(db, "SELECT * FROM nest_success")
names(SHnests)[4:6]<-c("First_marked", "Last_active","Last_checked")
odbcClose(db)
str(SHnests)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LOAD NEST SUCCESS DATA FROM ASCENSION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Load data from database
setwd("S:\\ConSci\\DptShare\\SteffenOppel/RSPB/UKOT/Ascension/Controlled_Area/Data")
setwd("C:/STEFFEN/RSPB/UKOT/Ascension/Controlled_Area/Data")
ASInests <- read.table("MABO_nest_data.csv", header=T, sep=",")
ASInests$island<-"Ascension"
ASInests$First_marked<-as.POSIXct(ASInests$First_marked, format="%d/%m/%Y")
ASInests$Last_active<-as.POSIXct(ASInests$Last_active, format="%d/%m/%Y")
ASInests$Last_checked<-as.POSIXct(ASInests$Last_checked, format="%d/%m/%Y")
sum(ASInests$fledged)/dim(ASInests)[1]				## nesting success


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUMMARIZE NEST SUCCESS DATA FROM ST HELENA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SHnests$island<-"StHelena"
head(SHnests)
SHnests$complete<-SHnests$failed+SHnests$fledged
dim(SHnests)[1]-sum(SHnests$complete)				## number of nests excluded because not yet finished
SHnests<-SHnests[SHnests$complete==1,]
sum(SHnests$fledged)/dim(SHnests)[1]				## nesting success


### exclude nests from poor 'winter' 2013
SHnests2<-SHnests[SHnests$First_marked<as.POSIXct(as.Date("2013-09-01")),]
sum(SHnests2$fledged)/dim(SHnests2)[1]





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TEST ISLAND DIFFERENCES IN NEST SUCCESS OR NEST SURVIVAL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hist(nests$First_Marked)

###########################################################################################
### 1. MERGE DATA AND ASSIGN 'SEASONS' (= called 'Year')
###########################################################################################

head(ASInests)
nests<-data.frame(rbind(ASInests[,c(2,1,3:7)], SHnests[,c(1,3,4:6,8:9)]))
nests$island<-as.factor(nests$island)
nests$observed<-round(ifelse(as.numeric(nests$Last_checked-nests$First_marked)<164,as.numeric(nests$Last_checked-nests$First_marked),164),0)
nests$fate<-ifelse(nests$fledged>0,0,1)
nests$count<-1
nests<-nests[order(nests$First_marked),]

### EXCLUDE NESTS FROM 2011 because of low sample size ####

nests<-nests[!(nests$First_marked<as.POSIXct("2012-01-01")),]

### create different seasons, including two for 2013 (summer and winter), labelled 2014 as 2015 for convenience
nests$Year<-ifelse(nests$First_marked<as.POSIXct("2013-02-01"),2012,ifelse(nests$First_marked<as.POSIXct("2013-12-01"),2013,ifelse(nests$First_marked<as.POSIXct("2014-02-01"),2014,2015)))
nests$Year<-as.factor(nests$Year)



###########################################################################################
### 2. MAYFIELD LOGISTIC REGRESSION WITH SEASON AS RANDOM EFFECT
###########################################################################################

m1<-glmer(cbind(fate,observed)~island+(1|Year),family=binomial, data=nests, na.action=na.omit)
m0<-glmer(cbind(fate,observed)~1+(1|Year),family=binomial, data=nests, na.action=na.omit)
anova(m1,m0)

nests$predsuccess<-predict(m1, type='response')

aggregate(predsuccess~island, data=nests, FUN=mean)
aggregate(fate~island, data=nests, FUN=mean)
aggregate(observed~fate+island, data=nests, FUN=mean)


###########################################################################################
### 3. SIMILAR APPROACH BUT WITH SHAFFER LOGISTIC EXPOSURE LINK FUNCTION
###########################################################################################
## updated 3 Oct 2014
## define link function (by Ben Bolker) from: http://stackoverflow.com/questions/19012128/user-defined-link-function-for-glmer-for-known-fate-survival-modeling
## http://rpubs.com/bbolker/logregexp

## test with simulated data
#succ<-data.frame(fate=rbinom(100,1,0.65), island="ASI", Year=runif(100,2011,2014), observed=runif(100,25,164))
#succ2<-data.frame(fate=rbinom(100,1,0.35), island="STH", Year=runif(100,2011,2014), observed=runif(100,25,164))
#nests<-rbind(succ,succ2)
#nests$Year<-as.integer(nests$Year)
#nests$observed<-as.integer(nests$observed)
#head(nests)


library(MASS)
logexp <- function(exposure = 1)
{
    linkfun <- function(mu) qlogis(mu^(1/exposure))
    ## FIXME: is there some trick we can play here to allow
    ##   evaluation in the context of the 'data' argument?
    linkinv <- function(eta)  plogis(eta)^exposure
    logit_mu_eta <- function(eta) {
        ifelse(abs(eta)>30,.Machine$double.eps,
               exp(eta)/(1+exp(eta))^2)
        ## OR .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
    }
    mu.eta <- function(eta) {       
        exposure * plogis(eta)^(exposure-1) *
            logit_mu_eta(eta)
    }
    valideta <- function(eta) TRUE
    link <- paste("logexp(", deparse(substitute(exposure)), ")",
                   sep="")
    structure(list(linkfun = linkfun, linkinv = linkinv,
                   mu.eta = mu.eta, valideta = valideta, 
                   name = link),
              class = "link-glm")
}

m1<-glmer(fate~island+(1|Year),family=binomial(link=logexp(nests$observed)), data=nests, na.action=na.omit)
m0<-glmer(fate~1+(1|Year),family=binomial(link=logexp(nests$observed)), data=nests, na.action=na.omit)
anova(m1,m0)

nests$predsuccess<-predict(m1, type='response')

aggregate(predsuccess~island, data=nests, FUN=mean)
aggregate(fate~island, data=nests, FUN=mean)

###########################################################################################
### 3. SIMPLE MAYFIELD NEST SUCCESS ESTIMATE 
###########################################################################################

Mayfield_nest_success<-data.frame(site=c('Ascension','StHelena'), nest_success=0, n_nests=0, n_fledged=0)

Mayfield_nest_success$n_nests<-aggregate(count~island, data=nests, FUN=sum)[,2]
Mayfield_nest_success$n_fledged<-aggregate(fledged~island, data=nests, FUN=sum)[,2]


nests$exposure<-nests$Last_checked-nests$First_marked

for (i in c("Ascension", "StHelena")){
exposure<-as.numeric(sum(nests$exposure[nests$island==i]))
fate<-sum(nests$fate[nests$island==i])
Mayfield_dsr<-1-(fate/exposure)
Mayfield_nest_success$nest_success[Mayfield_nest_success$site==i]<-Mayfield_dsr^164    			### based on 44 days of incubation and 120 days for fledging from Nelson (1978) and Dorward (1962)
}
Mayfield_nest_success
write.table(Mayfield_nest_success, "MABO_nest_success_summary.csv",row.names=F, sep=',')





###########################################################################################
### 4. DAILY NEST SURVIVAL IN RMARK 
###########################################################################################

####  FORMATTING DATA FOR NEST SURVIVAL ANALYSIS ########

###### CONVERTING THE DATES and setting the endpoints of the season  ####################################
head(nests)
for (i in c(3:5)){
nests[,i]<-as.Date(nests[,i])
}
start<-min(nests$First_marked)-1      # the first date when a nest was initiated
nocc<-as.integer(max(nests$Last_checked)- start)    


###### CONVERT DATES TO CONTINUOUS NUMBERS (counted from day of first nest) 

nests$First_marked<-as.integer((nests$First_marked)-start)
nests$Last_active<-as.integer((nests$Last_active)-start)
nests$Last_checked<-as.integer((nests$Last_checked)-start)

#head(nests)


############## CREATE DATA FRAME IN THE ORDER MARK NEEDS IT WITH APPROPRIATE COLUMN HEADERS ####################

MABO<-data.frame(Nest_ID=nests[,1])
MABO$FirstFound<-as.integer(nests[,3])
MABO$LastPresent<-as.integer(nests[,4])
MABO$LastChecked<-as.integer(nests[,5])
MABO$Fate<-ifelse(nests$fledged>0,0,1)
MABO$ASI<-ifelse(nests$island=="Ascension",1,0)
MABO$Year<-nests$Year
MABO$EXCL<-MABO$LastChecked-MABO$LastPresent
hist(MABO$FirstFound)

MABO[MABO$EXCL>50,]


##############################################################################################
#### SET MAXIMUM EXPOSURE TIME TO 164 DAYS ########################################
##############################################################################################

## limit exposure time to 164 days for successful nests

MABO$LastChecked[MABO$Fate==0]<-ifelse(MABO$LastChecked[MABO$Fate==0]-MABO$FirstFound[MABO$Fate==0]>164,MABO$FirstFound[MABO$Fate==0]+163,MABO$LastChecked[MABO$Fate==0])
MABO$LastPresent[MABO$Fate==0]<-MABO$LastChecked[MABO$Fate==0]


## limit exposure time to 164 days for UNsuccessful nests

MABO$LastChecked[MABO$Fate==1]<-ifelse(MABO$LastChecked[MABO$Fate==1]-MABO$FirstFound[MABO$Fate==1]>164,MABO$FirstFound[MABO$Fate==1]+163,MABO$LastChecked[MABO$Fate==1])


##############################################################################################
#### ensuring data integrity for MARK analysis: LastPresent=LastChecked for nests that fledged
##############################################################################################
rm(m1,m2,m3,m4,noerror)
noerror<-"CONGRATULATIONS - no errors were found!"
#MABO$LastPresent[MABO$Fate==0] - MABO$LastChecked[MABO$Fate==0]				### THIS MUST ALL BE 0!!
error1<-(MABO$LastPresent[MABO$Fate==0] - MABO$LastChecked[MABO$Fate==0])
if (sum(error1)!=0){
check<-data.frame(Nest=MABO$Nest_ID[MABO$Fate==0], error1)
m1<-sprintf("ERROR! check the nest visits and 'hatched' column of nest %s - if the nest hatched successfully the last day it was alive must be the last day it was checked", check[error1!=0,1])
}



#MABO$LastPresent[MABO$Fate==1] < MABO$LastChecked[MABO$Fate==1]				### THIS MUST ALL BE TRUE!!
error2<-(MABO$LastPresent[MABO$Fate==1] < MABO$LastChecked[MABO$Fate==1])
if ('FALSE' %in% error2){
check<-data.frame(Nest=MABO$Nest_ID[MABO$Fate==1], error2)
m2<-sprintf("ERROR! check the nest visits and 'hatched' column of nest %s - if the nest did not hatch then the last date it was checked must be AFTER the last day it was alive", check[error2==FALSE,1])
}



#MABO$FirstFound[MABO$Fate==1] < MABO$LastChecked[MABO$Fate==1]				### THIS MUST ALL BE TRUE!!
error3<-(MABO$FirstFound[MABO$Fate==1] < MABO$LastChecked[MABO$Fate==1])
if ('FALSE' %in% error3){
check<-data.frame(Nest=MABO$Nest_ID[MABO$Fate==1], error3)
m3<-sprintf("ERROR! check the nest visits and 'hatched' column of nest %s - it is generally not possible to detect an already depredated nest", check[error3==FALSE,1])
}

#MABO$FirstFound[MABO$Fate==0] < MABO$LastChecked[MABO$Fate==0]				### THIS MUST ALL BE TRUE!!
error4<-(MABO$FirstFound[MABO$Fate==0] < MABO$LastChecked[MABO$Fate==0])
if ('FALSE' %in% error4){
check<-data.frame(Nest=MABO$Nest_ID[MABO$Fate==0], error4)
m4<-sprintf("ERROR! check the nest visits and 'hatched' column of nest %s - at least two visits are necessary to judge whether a nest was successful, so the last visit must be after the first visit", check[error4==FALSE,1])
}


### PRINT ANY EXISTING ERROR MESSAGES
all<-ls()
if('m1' %in% all){m1}else
if('m2' %in% all){m2}else
if('m3' %in% all){m3}else
if('m4' %in% all){m4}else
{noerror}





###########################################################################################
### RUNNING A SET OF CANDIDATE MODELS FOR DAILY NEST SURVIVAL ################
###########################################################################################

nocc<-max(MABO[,4])
run.MABO=function()
{

############ Environment covariate models based on research by Burns 2011 ######################
MarkPath='C:/Program Files/MARK 7.1'
# 1. Nest survival rate is constant
m1_null<-mark(MABO,nocc=nocc,model="Nest", model.parameters=list(S=list(formula=~1)))

# 2. Nest survival rate varies by island
m2_island<-mark(MABO,nocc=nocc,model="Nest", model.parameters=list(S=list(formula=~ASI)))

# 3. Nest survival rate varies by year and island
#m3_site_year<-mark(MABO,nocc=nocc,model="Nest", model.parameters=list(S=list(formula=~ASI+Year)))

# 4. Nest survival rate varies with year x island
#m4_siteXyear<-mark(MABO,nocc=nocc,model="Nest", model.parameters=list(S=list(formula=~Year)))


return(collect.models() )
}

MABO.results<-run.MABO()  
MABO_AIC_table<-print.marklist(MABO.results)

## model averaged parameter estimate for daily nest survival across all sites
## this only makes sense when the models do not assume temporal variation in nest survival
MABO_averaged_parameter_estimates<-model.average(MABO.results, vcv=T)
DAILY_SURVIVAL<-MABO_averaged_parameter_estimates$estimates[1,2:5]
DAILY_SURVIVAL^164




###########################################################################################
### SUMMARISING REAL PARAMETER ESTIMATES FROM TOP MODEL ################
###########################################################################################

## RUN THE TOP MODEL

TOP<-mark(MABO,nocc=nocc,model="Nest", model.parameters=list(S=list(formula=~PB*CAT)))


###### to estimate overall nest survival, we assume each nest has to survive 30 days ######
###### this produces nonsense with Time covariate models ######

#output<-data.frame(site=c('Man_and_Horse','ProsperousBay','Deadwood','Broad_Bottom'),  CAT=rep(c(1,0),each=4))
output<-data.frame(site=rep(c('Man_and_Horse','ProsperousBay'),2), CAT=c(1,1,0,0))

output$MH<-ifelse(output$site=='Man_and_Horse',1,0)
output$BB<-ifelse(output$site=='Broad_Bottom',1,0)
output$PB<-ifelse(output$site=='ProsperousBay',1,0)
output$DW<-ifelse(output$site=='Deadwood',1,0) 
out<-covariate.predictions(TOP,data=output,indices=1)

nest_survival<-out$estimates[,c(4,5,10:13)]
nest_survival[,3:6]<-nest_survival[,3:6]^30			## Fiona used 30 days!!
nest_survival

nest_survival$Mayfield<-Mayfield_nest_success$nest_success[Mayfield_nest_success$site %in% c('MH','PB')]
nest_survival

##### clean up all files created by MARK but not necessary to hold output, without asking for permission
rm(MABO.results, m3b_site_manage, TOP)
detach("package:verification", character.only = TRUE, unload = TRUE)
detach("package:fields", character.only = TRUE, unload = TRUE)
detach("package:spam", character.only = TRUE, unload = TRUE)
cleanup(lx=NULL,ask=FALSE, prefix="mark")


###### EXPORT RESULTS #####
#write.table(output,"MABO_daily_nest_survival_estimates.csv", sep=",", row.names=F)
write.table(nest_survival,"MABO_nest_survival_estimates.csv", sep=",", row.names=F)
write.table(MABO_AIC_table,"MABO_nest_survival_AIC_table.csv", sep=",", row.names=F)
write.table(MABO_AIC_table,"clipboard", sep="\t", row.names=F)



####################################################################################################################
######## PLOTTING THE DATA ##################################################
####################################################################################################################

## plot for paper
par(mar=c(3,4.5,0,0),oma=c(1,1,0,0))
pb<-subset(nest_survival, CAT==1)
errbar(c(1,4), pb[2:1,3], pb[2:1,5],pb[2:1,6],xlim=c(0,6), ylim=c(0,1), cex=1.5, cex.lab=1.6, pch=16, bty='n', xlab="", ylab="St Helena Plover nest survival", axes=F)
par(new=T)
mh<-subset(nest_survival, CAT==0)
errbar(c(2,5), mh[2:1,3],mh[2:1,5],mh[2:1,6], xlim=c(0,6), ylim=c(0,1), cex=1.5, cex.lab=1.6, pch=1, bty='n', xlab="", ylab="", axes=F)
axis(1, at= c(0,1.5,4.5,6), labels=c("","semi-desert","pasture",""), tck=0.015, mgp=c(3,1,0), cex.axis=1.5)
axis(2, at=seq(0,1,0.2), labels=T, tck=0.015, mgp=c(3,1,0), cex.axis=1.4, las=1)
legend("topleft",c("before cat control","after cat control"), pch=c(16,1), cex=1.6, bty='n')




### save output graphic to pdf
pdf("MABO_nest_survival_management.pdf", width=6, height=6)

par(mar=c(3,4.5,0,0),oma=c(1,1,0,0))
pb<-subset(nest_survival, CAT==1)
errbar(c(1,4), pb[2:1,3], pb[2:1,5],pb[2:1,6],xlim=c(0,6), ylim=c(0,1), cex=1.5, cex.lab=1.6, pch=16, bty='n', xlab="", ylab="St Helena Plover nest survival", axes=F)
par(new=T)
mh<-subset(nest_survival, CAT==0)
errbar(c(2,5), mh[2:1,3],mh[2:1,5],mh[2:1,6], xlim=c(0,6), ylim=c(0,1), cex=1.5, cex.lab=1.6, pch=1, bty='n', xlab="", ylab="", axes=F)
axis(1, at= c(0,1.5,4.5,6), labels=c("","semi-desert","pasture",""), tck=0.015, mgp=c(3,1,0), cex.axis=1.5)
axis(2, at=seq(0,1,0.2), labels=T, tck=0.015, mgp=c(3,1,0), cex.axis=1.4, las=1)
legend("topleft",c("before cat control","after cat control"), pch=c(16,1), cex=1.6, bty='n')

dev.off()

