rm(list=ls())

# This script contains code needed to summarize relevant results from
# GIMME-AR and GIMME-MS model fits in the lag. greater (LG) condition. It
# extracts information about estimated beta weights and fit indices, conducts
# AIC-based solution-reduction for GIMME-MS, and generates data frames
# that contain direction recall and direction precision values for 
# each simulated participant's data.


#############################################
####### 1. MAKE SUMMARY DATA FRAMES #########
#############################################
require(stringr)

#list of all similation repetitions
reps<-list.files(path = "data/")

ar.recovery.sum<-data.frame(NA,NA)
ms.recovery.sum<-data.frame(NA,NA)

colnames(ar.recovery.sum)<-c("rep","sub")
for (r in 1:length(reps)){
  tmp<-data.frame(rep(reps[r],50),seq(from=1,to=50,by=1))
  colnames(tmp)<-c("rep","sub")
  ar.recovery.sum<-rbind(ar.recovery.sum,tmp)
}
ar.recovery.sum<-ar.recovery.sum[!is.na(ar.recovery.sum$rep),]
rownames(ar.recovery.sum)<-seq(1,75000,1)


colnames(ms.recovery.sum)<-c("rep","sub")
for (r in 1:length(reps)){
  tmp<-data.frame(rep(reps[r],50),seq(from=1,to=50,by=1))
  colnames(tmp)<-c("rep","sub")
  ms.recovery.sum<-rbind(ms.recovery.sum,tmp)
}
ms.recovery.sum<-ms.recovery.sum[!is.na(ms.recovery.sum$rep),]
rownames(ms.recovery.sum)<-seq(1,75000,1)


ar.recovery.sum$tsl<-NA
ar.recovery.sum$ar<-NA

ar.recovery.sum[str_detect(ar.recovery.sum$rep,"t_50"),]$tsl<-50
ar.recovery.sum[str_detect(ar.recovery.sum$rep,"t_100"),]$tsl<-100
ar.recovery.sum[str_detect(ar.recovery.sum$rep,"t_300"),]$tsl<-300

ar.recovery.sum[str_detect(ar.recovery.sum$rep,"ar_0"),]$ar<-.0
ar.recovery.sum[str_detect(ar.recovery.sum$rep,"ar_0.1"),]$ar<-.10
ar.recovery.sum[str_detect(ar.recovery.sum$rep,"ar_0.3"),]$ar<-.30
ar.recovery.sum[str_detect(ar.recovery.sum$rep,"ar_0.5"),]$ar<-.50
ar.recovery.sum[str_detect(ar.recovery.sum$rep,"ar_0.6"),]$ar<-.60


ms.recovery.sum$tsl<-NA
ms.recovery.sum$ar<-NA

ms.recovery.sum[str_detect(ms.recovery.sum$rep,"t_50"),]$tsl<-50
ms.recovery.sum[str_detect(ms.recovery.sum$rep,"t_100"),]$tsl<-100
ms.recovery.sum[str_detect(ms.recovery.sum$rep,"t_300"),]$tsl<-300

ms.recovery.sum[str_detect(ms.recovery.sum$rep,"ar_0"),]$ar<-.0
ms.recovery.sum[str_detect(ms.recovery.sum$rep,"ar_0.1"),]$ar<-.10
ms.recovery.sum[str_detect(ms.recovery.sum$rep,"ar_0.3"),]$ar<-.30
ms.recovery.sum[str_detect(ms.recovery.sum$rep,"ar_0.5"),]$ar<-.50
ms.recovery.sum[str_detect(ms.recovery.sum$rep,"ar_0.6"),]$ar<-.60

# two different data frames for MS, one for plausible models only (PM),
# the other for all models (AM)
ms.recovery.sum.PM<-ms.recovery.sum
ms.recovery.sum.AM<-ms.recovery.sum

##################################################
######## 2. CHECK FOR IMPLAUSIBLE AR MODELS ######
##################################################

# as AR gimme does not currently record a "check psi" flag, we must go back and retrospectively
# record any models which have values >1 or <0 in the diagonal of the southeast psi matrix 

# list of psi matices
psi.con<-list()

# variable to label questionable models
ar.recovery.sum$badPsi<-NA

# pull psi values or everyone
for (r in 1:length(ar.recovery.sum$rep)){
  
  # get paths and make beta and SE lists
  psi.path<-paste("./outputAR/",
                  ar.recovery.sum$rep[r],"/individual/ind_",
                  ar.recovery.sum$sub[r],"Psi.csv",sep="")

  #in case some didn't finish
  if(file.exists(psi.path)){
  
  psi.con[[r]]<-read.csv(psi.path,row.names=1)
  
  #Is Psi bad?
  g1<-diag(as.matrix(psi.con[[r]][9:16,9:16]))>1
  l0<-diag(as.matrix(psi.con[[r]][9:16,9:16]))<0
  if (sum(c(g1,l0))>0){ar.recovery.sum$badPsi[r]<-TRUE}
  else {ar.recovery.sum$badPsi[r]<-FALSE}
  }
}

# 14 individuals (~.02%) did not converge
length(ar.recovery.sum[is.na(ar.recovery.sum$badPsi),]$badPsi)
# [1] 14
14/75000*100
#[1] 0.01866667

# 1907 individual models (2.5%) were implausible
length(ar.recovery.sum[ar.recovery.sum$badPsi==TRUE,]$badPsi)
# [1] 1907
1907/75000*100
# [1] 2.542667




###################################################################
######## 3. MULTIPLE SOLUTION REDUCTION: ALL MODELS ###############
###################################################################

# variables to label implausible or non-converged  models
ms.recovery.sum.AM$no_model<-NA
ms.recovery.sum.AM$grp_sol<-NA
ms.recovery.sum.AM$ind_sol<-NA

for (r in unique(ms.recovery.sum.AM$rep)){
  tmp.sum<-read.csv(file = paste0("./outputMS/",r,"/summaryFit.csv"))
  tmp.select<-matrix(nrow=50,ncol=length(unique(tmp.sum$grp_sol)));tmp.select<-as.data.frame(tmp.select)
  colnames(tmp.select)<-seq(from=1,to=length(unique(tmp.sum$grp_sol)),by=1)
  
  for (i in 1:50){
    for(s in 1:length(unique(tmp.sum$grp_sol))){
      if(length(tmp.sum[tmp.sum$subj==paste0("ind_",i) & tmp.sum$grp_sol==s,]$aic)>0){
        tmp.AIC<-min(tmp.sum[tmp.sum$subj==paste0("ind_",i) & tmp.sum$grp_sol==s,]$aic)
        tmp.select[i,s]<-tmp.AIC}
    }
  }  
  tmp.prop<-apply(tmp.select,2,function(x) mean(is.na(x)))
  tmp.mean<-apply(tmp.select,2,function(x) mean(x,na.rm = TRUE))
  tmp.both<-data.frame(tmp.prop,tmp.mean)
  gsol<-order(tmp.both$tmp.mean)
  gsol_out<-gsol[1];
  if(sum(tmp.prop>.25)>0){  
    tmp<-tmp.both
    tmp[tmp$tmp.prop>.25,]$tmp.mean<-Inf
    gsol_out<-order(tmp$tmp.mean); gsol_out<-gsol_out[1]} 
  if(sum(tmp.prop<=.25)==0){   
    gsol_out<-NA}
  
  ms.recovery.sum.AM[ms.recovery.sum.AM$rep==r,]$grp_sol<-gsol_out
  
  for (i in 1:50){
    if(!is.na(gsol_out)){
      isol<-tmp.sum[tmp.sum$grp_sol==gsol_out & tmp.sum$subj==paste0("ind_",i) & tmp.sum$aic==tmp.select[i,gsol_out],]$ind_sol
      ms.recovery.sum.AM[ms.recovery.sum.AM$rep==r & ms.recovery.sum.AM$sub==i,]$ind_sol<-isol[1]}
  }
  
}

# label subjects for which there is no model that converged
ms.recovery.sum.AM$no_model<-FALSE
ms.recovery.sum.AM[is.na(ms.recovery.sum.AM$ind_sol),]$no_model<-TRUE


# 13 individuals had no plausible model that converged
length(ms.recovery.sum.AM[ms.recovery.sum.AM$no_model==TRUE,]$ind_sol)
# [1] 13
13/75000*100
#[1] 0.01733333

###################################################################
######## 4. MULTIPLE SOLUTION REDUCTION: PLAUSIBLE MODELS ONLY ####
###################################################################

# varible to label questionable models
ms.recovery.sum.PM$no_model<-NA
ms.recovery.sum.PM$grp_sol<-NA
ms.recovery.sum.PM$ind_sol<-NA

for (r in unique(ms.recovery.sum.PM$rep)){
  tmp.sum<-read.csv(file = paste0("./outputMS/",r,"/summaryFit.csv"))
  tmp.select<-matrix(nrow=50,ncol=length(unique(tmp.sum$grp_sol)));tmp.select<-as.data.frame(tmp.select)
  colnames(tmp.select)<-seq(from=1,to=length(unique(tmp.sum$grp_sol)),by=1)
  
  for (i in 1:50){
    for(s in 1:length(unique(tmp.sum$grp_sol))){
      if(length(tmp.sum[tmp.sum$subj==paste0("ind_",i) & tmp.sum$grp_sol==s & tmp.sum$checkPsi==FALSE,]$aic)>0){
        tmp.AIC<-min(tmp.sum[tmp.sum$subj==paste0("ind_",i) & tmp.sum$grp_sol==s & tmp.sum$checkPsi==FALSE,]$aic)
        tmp.select[i,s]<-tmp.AIC}
    }
  }  
  tmp.prop<-apply(tmp.select,2,function(x) mean(is.na(x)))
  tmp.mean<-apply(tmp.select,2,function(x) mean(x,na.rm = TRUE))
  tmp.both<-data.frame(tmp.prop,tmp.mean)
  gsol<-order(tmp.both$tmp.mean)
  gsol_out<-gsol[1];
  if(sum(tmp.prop>.25)>0){  
    tmp<-tmp.both
    tmp[tmp$tmp.prop>.25,]$tmp.mean<-Inf
    gsol_out<-order(tmp$tmp.mean); gsol_out<-gsol_out[1]} 
  if(sum(tmp.prop<=.25)==0){   
    gsol_out<-NA}
  
  ms.recovery.sum.PM[ms.recovery.sum.PM$rep==r,]$grp_sol<-gsol_out
  
  for (i in 1:50){
    if(!is.na(gsol_out)){
      isol<-tmp.sum[tmp.sum$grp_sol==gsol_out & tmp.sum$subj==paste0("ind_",i) & tmp.sum$aic==tmp.select[i,gsol_out],]$ind_sol
      ms.recovery.sum.PM[ms.recovery.sum.PM$rep==r & ms.recovery.sum.PM$sub==i,]$ind_sol<-isol[1]}
  }
  
}

# label subjects for which there is no model that converged
ms.recovery.sum.PM$no_model<-FALSE
ms.recovery.sum.PM[is.na(ms.recovery.sum.PM$ind_sol),]$no_model<-TRUE


# 2558 individuals had no plausible model
length(ms.recovery.sum.PM[ms.recovery.sum.PM$no_model==TRUE,]$ind_sol)
# [1] 2558
2558/75000*100
#[1] 3.410667



##########################################################
######## 5. LOAD "TRUE" PARAMETERS FOR SIMULATED DATA ####
##########################################################


# pull out betas and presence/abscence for each 
# simulated data set

true.type.con<-list()
true.beta.con<-list()


for (r in 1:length(ar.recovery.sum$rep)){
  true.type.path<-paste("./levels/",
                        ar.recovery.sum$rep[r],"/ind_",
                        ar.recovery.sum$sub[r],".csv",sep="")
  true.beta.path<-paste("./true/",
                         ar.recovery.sum$rep[r],"/ind_",
                         ar.recovery.sum$sub[r],".csv",sep="")
  
  true.type.con[[r]]<-read.csv(true.type.path)
  true.beta.con[[r]]<-read.csv(true.beta.path)

}

################################################################
######## 6. LOAD BETA/SE MATRICIES FOR AR, MS.AM and MS.PM #####
################################################################

# reconstruct recovered beta matrix from indivPathEstimates.csv file 
# (to be consisitent with the MS fits)

beta.con.AR<-list()

for (r in 1:length(ar.recovery.sum$rep)){
  if(ar.recovery.sum$sub[r]==1){
    load(file = paste0("./outputAR/",ar.recovery.sum$rep[r],"/out.RData"))
    tmp<-g.out$path_est_mats
    }
  beta.con.AR[[r]]<-tmp[paste0("ind_",ar.recovery.sum$sub[r])]
}


beta.con.MS.AM<-list()

for (r in 1:length(ms.recovery.sum.AM$rep)){
  if(ms.recovery.sum.AM$sub[r]==1){
    load(file = paste0("./outputMS/",ms.recovery.sum.AM$rep[r],"/out.RData"))
    tmp<-g.out$ind_fit
    sublist<-data.frame(seq(1,50,1),rep(NA,50));colnames(sublist)<-c("a","b")
    for (e in 1:50){
      sublist[e,2]<-tmp[[1]][[e]][[1]]$subj}
  }
  if(!is.na(ms.recovery.sum.AM$ind_sol[r])){
  num<-sublist[sublist$b==paste0("ind_",ms.recovery.sum.AM$sub[r]),]$a
  beta.con.MS.AM[[r]]<-tmp[[ms.recovery.sum.AM$grp_sol[r]]][[num]][[ms.recovery.sum.AM$ind_sol[r]]]$betas
  }
  else{beta.con.MS.AM[[r]]<-NA}
}


beta.con.MS.PM<-list()

for (r in 1:length(ms.recovery.sum.PM$rep)){
  if(ms.recovery.sum.PM$sub[r]==1){
    load(file = paste0("./outputMS/",ms.recovery.sum.PM$rep[r],"/out.RData"))
    tmp<-g.out$ind_fit
    sublist<-data.frame(seq(1,50,1),rep(NA,50));colnames(sublist)<-c("a","b")
    for (e in 1:50){
      sublist[e,2]<-tmp[[1]][[e]][[1]]$subj}
  }
  if(!is.na(ms.recovery.sum.PM$ind_sol[r])){
    num<-sublist[sublist$b==paste0("ind_",ms.recovery.sum.PM$sub[r]),]$a
    beta.con.MS.PM[[r]]<-tmp[[ms.recovery.sum.PM$grp_sol[r]]][[num]][[ms.recovery.sum.PM$ind_sol[r]]]$betas
  }
  else{beta.con.MS.PM[[r]]<-NA}
}



##########################################
##### 7. COMPUTE RECOVERY STATS ##########
##########################################

# load function to calculate recovery variables for each method
source("calc_functions.R")

#run functions to calculate recovery stats
ar.recovery.sum.AM.final<-calc.rec.stats.noex(sum = ar.recovery.sum,beta = beta.con.AR, nacol = "badPsi")
ar.recovery.sum.PM.final<-calc.rec.stats(sum = ar.recovery.sum,beta = beta.con.AR, nacol = "badPsi")

ms.recovery.sum.AM.final<-calc.rec.stats(sum = ms.recovery.sum.AM,beta = beta.con.MS.AM, nacol = "no_model")
ms.recovery.sum.PM.final<-calc.rec.stats(sum = ms.recovery.sum.PM,beta = beta.con.MS.PM, nacol = "no_model")


compare.AM<-data.frame(ar.recovery.sum[,1:4])
compare.PM<-data.frame(ar.recovery.sum[,1:4])

compare.AM$dir.recall.all<-(ms.recovery.sum.AM.final$dir.recall.all-ar.recovery.sum.AM.final$dir.recall.all)
compare.AM$dir.recall.con<-(ms.recovery.sum.AM.final$dir.recall.con-ar.recovery.sum.AM.final$dir.recall.con)
compare.AM$dir.precision.all<-(ms.recovery.sum.AM.final$dir.precision.all-ar.recovery.sum.AM.final$dir.precision.all)
compare.AM$dir.precision.con<-(ms.recovery.sum.AM.final$dir.precision.con-ar.recovery.sum.AM.final$dir.precision.con)

compare.PM$dir.recall.all<-(ms.recovery.sum.PM.final$dir.recall.all-ar.recovery.sum.PM.final$dir.recall.all)
compare.PM$dir.recall.con<-(ms.recovery.sum.PM.final$dir.recall.con-ar.recovery.sum.PM.final$dir.recall.con)
compare.PM$dir.precision.all<-(ms.recovery.sum.PM.final$dir.precision.all-ar.recovery.sum.PM.final$dir.precision.all)
compare.PM$dir.precision.con<-(ms.recovery.sum.PM.final$dir.precision.con-ar.recovery.sum.PM.final$dir.precision.con)


# save it
save.image("lag_greater_recovery.RData")

