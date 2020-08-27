rm(list=ls())

###########################################
###### Values for tables 1 and 2 ##########
###########################################

# Balanced models

reps<-read.csv("reps.csv")

sum.fit.bal<-list()

for (r in 1:length(reps$rep)){
  rep<-reps$rep[r]
  sum.fit.bal[[r]]<-read.csv(file = paste0("balanced_in/outputMS/",rep,"/summaryFit.csv"))}

save(sum.fit.bal,file="ms_sum_fit_balanced.RData")


# Con greater

reps<-read.csv("reps.csv")

sum.fit.con<-list()

for (r in 1:length(reps$rep)){
  rep<-reps$rep[r]
  sum.fit.con[[r]]<-read.csv(file = paste0("con_greater_in/outputMS/",rep,"/summaryFit.csv"))}

save(sum.fit.con,file="ms_sum_fit_con.RData")


# Lag greater

reps<-read.csv("reps.csv")

sum.fit.lag<-list()

for (r in 1:length(reps$rep)){
  rep<-reps$rep[r]
  sum.fit.lag[[r]]<-read.csv(file = paste0("lag_greater_in/outputMS/",rep,"/summaryFit.csv"))}

save(sum.fit.lag,file="ms_sum_fit_lag.RData")





