rm(list=ls())

# This file contains code that generates the tables and figures for the
# manuscript.

###########################################
###### Pre-analyses #######################
###########################################

# In this section, I load RData files for the BA, LG and CG 
# conditions seperately, and then save only the relevant data frames with
# recovery stats out to be combined in a smaller "sum" RData file

#load and save list of subjects and replications for each condition
load("Balanced/balanced_recovery.RData")

subs<-data.frame(ar.recovery.sum[,1:4])
write.csv(subs,file="subs.csv",row.names = FALSE)

reps<-data.frame(ar.recovery.sum[ar.recovery.sum$sub==1,1:4])
write.csv(reps,file="reps.csv",row.names = FALSE)

rm(list=ls())

# use save_ms_summary to pull summary data
# for each group-level replication, and save in
# a single data file

# this is done in the sequence below to avoid exceeding RAM limits

# save only summary data files from each

load("Balanced/balanced_recovery.RData")

ar.recovery.sum.AM.BA<-ar.recovery.sum.AM.final
ar.recovery.sum.GM.BA<-ar.recovery.sum.GM.final
ms.recovery.sum.AM.BA<-ms.recovery.sum.AM.final
ms.recovery.sum.GM.BA<-ms.recovery.sum.GM.final

save(ar.recovery.sum.AM.BA,ar.recovery.sum.GM.BA,
     ms.recovery.sum.AM.BA,ms.recovery.sum.GM.BA,
     file="sum_recovery_all.RData");rm(list=ls())

load("Con_Greater/con_greater_recovery.RData")
load("sum_recovery_all.RData")

ar.recovery.sum.AM.CG<-ar.recovery.sum.AM.final
ar.recovery.sum.GM.CG<-ar.recovery.sum.GM.final
ms.recovery.sum.AM.CG<-ms.recovery.sum.AM.final
ms.recovery.sum.GM.CG<-ms.recovery.sum.GM.final

save(ar.recovery.sum.AM.BA,ar.recovery.sum.GM.BA,
     ms.recovery.sum.AM.BA,ms.recovery.sum.GM.BA,
     ar.recovery.sum.AM.CG,ar.recovery.sum.GM.CG,
     ms.recovery.sum.AM.CG,ms.recovery.sum.GM.CG,
     file="sum_recovery_all.RData");rm(list=ls())


load("Lag_Greater/lag_greater_recovery.RData")
load("sum_recovery_all.RData")

ar.recovery.sum.AM.LG<-ar.recovery.sum.AM.final
ar.recovery.sum.GM.LG<-ar.recovery.sum.GM.final
ms.recovery.sum.AM.LG<-ms.recovery.sum.AM.final
ms.recovery.sum.GM.LG<-ms.recovery.sum.GM.final


#make "comparison" data frames
compare.AM.BA<-data.frame(ar.recovery.sum.AM.BA[,1:4])
compare.GM.BA<-data.frame(ar.recovery.sum.AM.BA[,1:4])

compare.AM.BA$dir.recall.all<-(ms.recovery.sum.AM.BA$dir.recall.all-ar.recovery.sum.AM.BA$dir.recall.all)
compare.AM.BA$dir.recall.con<-(ms.recovery.sum.AM.BA$dir.recall.con-ar.recovery.sum.AM.BA$dir.recall.con)
compare.AM.BA$dir.precision.all<-(ms.recovery.sum.AM.BA$dir.precision.all-ar.recovery.sum.AM.BA$dir.precision.all)
compare.AM.BA$dir.precision.con<-(ms.recovery.sum.AM.BA$dir.precision.con-ar.recovery.sum.AM.BA$dir.precision.con)

compare.GM.BA$dir.recall.all<-(ms.recovery.sum.GM.BA$dir.recall.all-ar.recovery.sum.GM.BA$dir.recall.all)
compare.GM.BA$dir.recall.con<-(ms.recovery.sum.GM.BA$dir.recall.con-ar.recovery.sum.GM.BA$dir.recall.con)
compare.GM.BA$dir.precision.all<-(ms.recovery.sum.GM.BA$dir.precision.all-ar.recovery.sum.GM.BA$dir.precision.all)
compare.GM.BA$dir.precision.con<-(ms.recovery.sum.GM.BA$dir.precision.con-ar.recovery.sum.GM.BA$dir.precision.con)

compare.AM.LG<-data.frame(ar.recovery.sum.AM.BA[,1:4])
compare.GM.LG<-data.frame(ar.recovery.sum.AM.BA[,1:4])

compare.AM.LG$dir.recall.all<-(ms.recovery.sum.AM.LG$dir.recall.all-ar.recovery.sum.AM.LG$dir.recall.all)
compare.AM.LG$dir.recall.con<-(ms.recovery.sum.AM.LG$dir.recall.con-ar.recovery.sum.AM.LG$dir.recall.con)
compare.AM.LG$dir.precision.all<-(ms.recovery.sum.AM.LG$dir.precision.all-ar.recovery.sum.AM.LG$dir.precision.all)
compare.AM.LG$dir.precision.con<-(ms.recovery.sum.AM.LG$dir.precision.con-ar.recovery.sum.AM.LG$dir.precision.con)

compare.GM.LG$dir.recall.all<-(ms.recovery.sum.GM.LG$dir.recall.all-ar.recovery.sum.GM.LG$dir.recall.all)
compare.GM.LG$dir.recall.con<-(ms.recovery.sum.GM.LG$dir.recall.con-ar.recovery.sum.GM.LG$dir.recall.con)
compare.GM.LG$dir.precision.all<-(ms.recovery.sum.GM.LG$dir.precision.all-ar.recovery.sum.GM.LG$dir.precision.all)
compare.GM.LG$dir.precision.con<-(ms.recovery.sum.GM.LG$dir.precision.con-ar.recovery.sum.GM.LG$dir.precision.con)

compare.AM.CG<-data.frame(ar.recovery.sum.AM.BA[,1:4])
compare.GM.CG<-data.frame(ar.recovery.sum.AM.BA[,1:4])

compare.AM.CG$dir.recall.all<-(ms.recovery.sum.AM.CG$dir.recall.all-ar.recovery.sum.AM.CG$dir.recall.all)
compare.AM.CG$dir.recall.con<-(ms.recovery.sum.AM.CG$dir.recall.con-ar.recovery.sum.AM.CG$dir.recall.con)
compare.AM.CG$dir.precision.all<-(ms.recovery.sum.AM.CG$dir.precision.all-ar.recovery.sum.AM.CG$dir.precision.all)
compare.AM.CG$dir.precision.con<-(ms.recovery.sum.AM.CG$dir.precision.con-ar.recovery.sum.AM.CG$dir.precision.con)

compare.GM.CG$dir.recall.all<-(ms.recovery.sum.GM.CG$dir.recall.all-ar.recovery.sum.GM.CG$dir.recall.all)
compare.GM.CG$dir.recall.con<-(ms.recovery.sum.GM.CG$dir.recall.con-ar.recovery.sum.GM.CG$dir.recall.con)
compare.GM.CG$dir.precision.all<-(ms.recovery.sum.GM.CG$dir.precision.all-ar.recovery.sum.GM.CG$dir.precision.all)
compare.GM.CG$dir.precision.con<-(ms.recovery.sum.GM.CG$dir.precision.con-ar.recovery.sum.GM.CG$dir.precision.con)

# save out
save(ar.recovery.sum.AM.BA,ar.recovery.sum.GM.BA,
     ms.recovery.sum.AM.BA,ms.recovery.sum.GM.BA,
     ar.recovery.sum.AM.CG,ar.recovery.sum.GM.CG,
     ms.recovery.sum.AM.CG,ms.recovery.sum.GM.CG,
     ar.recovery.sum.AM.LG,ar.recovery.sum.GM.LG,
     ms.recovery.sum.AM.LG,ms.recovery.sum.GM.LG,
     compare.AM.BA,compare.GM.BA,
     compare.AM.LG,compare.GM.LG,
     compare.AM.CG,compare.GM.CG,
     file="sum_recovery_all.RData");rm(list=ls())

###########################################
###### Multiple Solutions #################
###########################################

# Balanced Condition
load("ms_sum_fit_balanced.RData")

sum.ms.bal<-read.csv("reps.csv")

sum.ms.bal$Gsol<-NA # number of group solutions per replication
sum.ms.bal$Isol<-NA # total number of individual-level solutions per replication
sum.ms.bal$n.subs<-50 # number of subjects per replication (always 50 in this case)

for (r in 1:length(sum.ms.bal$rep)){
  sum.ms.bal$Gsol[r]<-max(sum.fit.bal[[r]]$grp_sol)
  sum.ms.bal$Isol[r]<-length(sum.fit.bal[[r]]$grp_sol)
}

# average number of individual-level solutions per subject and group-level solution
sum.ms.bal$M.Isol.Gsol<-(sum.ms.bal$Isol/sum.ms.bal$n.subs/sum.ms.bal$Gsol)

# average number of individual-level solutions per subject across all 
# group-level solutions for that replication
sum.ms.bal$M.Isol.rep<-(sum.ms.bal$Isol/sum.ms.bal$n.subs)


# Con Greater Condition
load("ms_sum_fit_con.RData")

sum.ms.con<-read.csv("reps.csv")

sum.ms.con$Gsol<-NA # number of group solutions per replication
sum.ms.con$Isol<-NA # total number of individual-level solutions per replication
sum.ms.con$n.subs<-50 # number of subjects per replication (always 50 in this case)

for (r in 1:length(sum.ms.con$rep)){
  sum.ms.con$Gsol[r]<-max(sum.fit.con[[r]]$grp_sol)
  sum.ms.con$Isol[r]<-length(sum.fit.con[[r]]$grp_sol)
}

# average number of individual-level solutions per subject and group-level solution
sum.ms.con$M.Isol.Gsol<-(sum.ms.con$Isol/sum.ms.con$n.subs/sum.ms.con$Gsol)

# average number of individual-level solutions per subject across all 
# group-level solutions for that replication
sum.ms.con$M.Isol.rep<-(sum.ms.con$Isol/sum.ms.con$n.subs)


# Lag Greater Condition
load("ms_sum_fit_lag.RData")

sum.ms.lag<-read.csv("reps.csv")

sum.ms.lag$Gsol<-NA # number of group solutions per replication
sum.ms.lag$Isol<-NA # total number of individual-level solutions per replication
sum.ms.lag$n.subs<-50 # number of subjects per replication (always 50 in this case)

for (r in 1:length(sum.ms.lag$rep)){
  sum.ms.lag$Gsol[r]<-max(sum.fit.lag[[r]]$grp_sol)
  sum.ms.lag$Isol[r]<-length(sum.fit.lag[[r]]$grp_sol)
}

# average number of individual-level solutions per subject and group-level solution
sum.ms.lag$M.Isol.Gsol<-(sum.ms.lag$Isol/sum.ms.lag$n.subs/sum.ms.lag$Gsol)

# average number of individual-level solutions per subject across all 
# group-level solutions for that replication
sum.ms.lag$M.Isol.rep<-(sum.ms.lag$Isol/sum.ms.lag$n.subs)

#add columns for each condition and bind together
require(plyr)

sum.ms.bal$CL<-"bal"
sum.ms.con$CL<-"con"
sum.ms.lag$CL<-"lag"

sum.ms.all<-rbind(sum.ms.bal,sum.ms.con,sum.ms.lag)

# data frame of means
means.ms.all<- ddply(sum.ms.all, c("tsl", "ar", "CL"), summarise,
                     Gsol    = mean(Gsol),
                     M.Isol.Gsol    = mean(M.Isol.Gsol),
                     M.Isol.rep    = mean(M.Isol.rep)
                     )
means.ms.all$ar<-factor(means.ms.all$ar)

# plots for Figure 1 in paper
require(ggplot2)
require(grid)
require(gridExtra)

Gsol.50 <- ggplot(data = means.ms.all[means.ms.all$tsl==50,], 
                  aes(y = Gsol, x = ar, colour = CL, group=CL)) + 
                  geom_line(show.legend = FALSE,size=1.5) + ylim(0,12.5) + theme_bw() +
                                  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  labs(title="50 Time Points",x="autoregressive strength",y="Gsol/rep")+theme(plot.title = element_text(hjust = 0.5))


Gsol.100 <- ggplot(data = means.ms.all[means.ms.all$tsl==100,], 
                  aes(y = Gsol, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,12.5) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  labs(title="100 Time Points",x="autoregressive strength",y="")+theme(plot.title = element_text(hjust = 0.5))


Gsol.300 <- ggplot(data = means.ms.all[means.ms.all$tsl==300,], 
                   aes(y = Gsol, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,12.5) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  labs(title="300 Time Points",x="autoregressive strength",y="")+theme(plot.title = element_text(hjust = 0.5))


Isol.Gsol.50 <- ggplot(data = means.ms.all[means.ms.all$tsl==50,], 
                  aes(y = M.Isol.Gsol, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,6) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  labs(title="",x="autoregressive strength",y=expression(paste(mu,"Isol/Gsol")))+theme(plot.title = element_text(hjust = 0.5))


Isol.Gsol.100 <- ggplot(data = means.ms.all[means.ms.all$tsl==100,], 
                   aes(y = M.Isol.Gsol, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,6) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  labs(title="",x="autoregressive strength",y="")+theme(plot.title = element_text(hjust = 0.5))


Isol.Gsol.300 <- ggplot(data = means.ms.all[means.ms.all$tsl==300,], 
                   aes(y = M.Isol.Gsol, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,6) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  labs(title="",x="autoregressive strength",y="")+theme(plot.title = element_text(hjust = 0.5))


Isol.rep.50 <- ggplot(data = means.ms.all[means.ms.all$tsl==50,], 
                       aes(y = M.Isol.rep, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,20) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  labs(title="",x="autoregressive strength",y=expression(paste(mu,"Isol/rep")))+theme(plot.title = element_text(hjust = 0.5))


Isol.rep.100 <- ggplot(data = means.ms.all[means.ms.all$tsl==100,], 
                        aes(y = M.Isol.rep, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,20) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  labs(title="",x="autoregressive strength",y="")+theme(plot.title = element_text(hjust = 0.5))


Isol.rep.300 <- ggplot(data = means.ms.all[means.ms.all$tsl==300,], 
                        aes(y = M.Isol.rep, x = ar, colour = CL, group=CL)) + 
  geom_line(show.legend = FALSE,size=1.5) + ylim(0,20) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  labs(title="",x="autoregressive strength",y="")+theme(plot.title = element_text(hjust = 0.5))


grid.arrange(Gsol.50,Gsol.100,Gsol.300,
             Isol.Gsol.50,Isol.Gsol.100,Isol.Gsol.300,
             Isol.rep.50,Isol.rep.100,Isol.rep.300,
             ncol=3)
g<-arrangeGrob(Gsol.50,Gsol.100,Gsol.300,
               Isol.Gsol.50,Isol.Gsol.100,Isol.Gsol.300,
               Isol.rep.50,Isol.rep.100,Isol.rep.300,
               ncol=3)
ggsave(filename = "Figure_1.tif",g,units = "in",dpi=300,
       width = 10,height = 10,device="tiff")

##############################################
########  Direction Recovery Tables ##########
##############################################

#load data frames with recovery summary information
load("sum_recovery_all.RData")

#make overall summary data frame
require(plyr)

# summary data frames 
con.ar.sums.AM.BA<- ddply(ar.recovery.sum.AM.BA, c("tsl", "ar"), summarise,
                     hits.con    = sum(dir.hits.con,na.rm = TRUE),
                     misses.con    = sum(dir.misses.con,na.rm = TRUE),
                     cr.con    = sum(dir.cr.con,na.rm = TRUE),
                     fp.con    = sum(dir.fp.con,na.rm = TRUE),
                     dir.rec = mean(dir.recall.con,na.rm=TRUE),
                     dir.prec = mean(dir.precision.con,na.rm=TRUE)
)
con.ar.sums.AM.BA$recall<-con.ar.sums.AM.BA$hits.con/(con.ar.sums.AM.BA$hits.con+con.ar.sums.AM.BA$misses.con)
con.ar.sums.AM.BA$precision<-con.ar.sums.AM.BA$hits.con/(con.ar.sums.AM.BA$hits.con+con.ar.sums.AM.BA$fp.con)

con.ar.sums.AM.LG<- ddply(ar.recovery.sum.AM.LG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ar.sums.AM.LG$recall<-con.ar.sums.AM.LG$hits.con/(con.ar.sums.AM.LG$hits.con+con.ar.sums.AM.LG$misses.con)
con.ar.sums.AM.LG$precision<-con.ar.sums.AM.LG$hits.con/(con.ar.sums.AM.LG$hits.con+con.ar.sums.AM.LG$fp.con)

con.ar.sums.AM.CG<- ddply(ar.recovery.sum.AM.CG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ar.sums.AM.CG$recall<-con.ar.sums.AM.CG$hits.con/(con.ar.sums.AM.CG$hits.con+con.ar.sums.AM.CG$misses.con)
con.ar.sums.AM.CG$precision<-con.ar.sums.AM.CG$hits.con/(con.ar.sums.AM.CG$hits.con+con.ar.sums.AM.CG$fp.con)


con.ms.sums.AM.BA<- ddply(ms.recovery.sum.AM.BA, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE),
                          dir.rec = mean(dir.recall.con,na.rm=TRUE),
                          dir.prec = mean(dir.precision.con,na.rm=TRUE)
)
con.ms.sums.AM.BA$recall<-con.ms.sums.AM.BA$hits.con/(con.ms.sums.AM.BA$hits.con+con.ms.sums.AM.BA$misses.con)
con.ms.sums.AM.BA$precision<-con.ms.sums.AM.BA$hits.con/(con.ms.sums.AM.BA$hits.con+con.ms.sums.AM.BA$fp.con)

con.ms.sums.AM.LG<- ddply(ms.recovery.sum.AM.LG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ms.sums.AM.LG$recall<-con.ms.sums.AM.LG$hits.con/(con.ms.sums.AM.LG$hits.con+con.ms.sums.AM.LG$misses.con)
con.ms.sums.AM.LG$precision<-con.ms.sums.AM.LG$hits.con/(con.ms.sums.AM.LG$hits.con+con.ms.sums.AM.LG$fp.con)

con.ms.sums.AM.CG<- ddply(ms.recovery.sum.AM.CG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ms.sums.AM.CG$recall<-con.ms.sums.AM.CG$hits.con/(con.ms.sums.AM.CG$hits.con+con.ms.sums.AM.CG$misses.con)
con.ms.sums.AM.CG$precision<-con.ms.sums.AM.CG$hits.con/(con.ms.sums.AM.CG$hits.con+con.ms.sums.AM.CG$fp.con)

# output for Table 1:

table1<-data.frame(c(con.ar.sums.AM.BA$ar,con.ar.sums.AM.BA$ar),
                   c(con.ar.sums.AM.BA$recall,con.ar.sums.AM.BA$precision),
                   c(con.ar.sums.AM.LG$recall,con.ar.sums.AM.LG$precision),
                   c(con.ar.sums.AM.CG$recall,con.ar.sums.AM.CG$precision),
                   c(con.ms.sums.AM.BA$recall,con.ms.sums.AM.BA$precision),
                   c(con.ms.sums.AM.LG$recall,con.ms.sums.AM.LG$precision),
                   c(con.ms.sums.AM.CG$recall,con.ms.sums.AM.CG$precision)
                   ) 

colnames(table1)<-c("auto.con","ar.BA","ar.LG","ar.CG","ms.BA","ms.LG","ms.CG")

table1$diff.BA<-table1$ms.BA-table1$ar.BA
table1$diff.LG<-table1$ms.LG-table1$ar.LG
table1$diff.CG<-table1$ms.CG-table1$ar.CG

# write out to paste into table:
write.csv(table1,row.names = FALSE,file="table_1_output.csv")


# Same for lagged relations

# summary data frames 
lag.ar.sums.AM.BA<- ddply(ar.recovery.sum.AM.BA, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ar.sums.AM.BA$recall<-lag.ar.sums.AM.BA$hits.lag/(lag.ar.sums.AM.BA$hits.lag+lag.ar.sums.AM.BA$misses.lag)
lag.ar.sums.AM.BA$precision<-lag.ar.sums.AM.BA$hits.lag/(lag.ar.sums.AM.BA$hits.lag+lag.ar.sums.AM.BA$fp.lag)

lag.ar.sums.AM.LG<- ddply(ar.recovery.sum.AM.LG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ar.sums.AM.LG$recall<-lag.ar.sums.AM.LG$hits.lag/(lag.ar.sums.AM.LG$hits.lag+lag.ar.sums.AM.LG$misses.lag)
lag.ar.sums.AM.LG$precision<-lag.ar.sums.AM.LG$hits.lag/(lag.ar.sums.AM.LG$hits.lag+lag.ar.sums.AM.LG$fp.lag)

lag.ar.sums.AM.CG<- ddply(ar.recovery.sum.AM.CG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ar.sums.AM.CG$recall<-lag.ar.sums.AM.CG$hits.lag/(lag.ar.sums.AM.CG$hits.lag+lag.ar.sums.AM.CG$misses.lag)
lag.ar.sums.AM.CG$precision<-lag.ar.sums.AM.CG$hits.lag/(lag.ar.sums.AM.CG$hits.lag+lag.ar.sums.AM.CG$fp.lag)


lag.ms.sums.AM.BA<- ddply(ms.recovery.sum.AM.BA, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ms.sums.AM.BA$recall<-lag.ms.sums.AM.BA$hits.lag/(lag.ms.sums.AM.BA$hits.lag+lag.ms.sums.AM.BA$misses.lag)
lag.ms.sums.AM.BA$precision<-lag.ms.sums.AM.BA$hits.lag/(lag.ms.sums.AM.BA$hits.lag+lag.ms.sums.AM.BA$fp.lag)

lag.ms.sums.AM.LG<- ddply(ms.recovery.sum.AM.LG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ms.sums.AM.LG$recall<-lag.ms.sums.AM.LG$hits.lag/(lag.ms.sums.AM.LG$hits.lag+lag.ms.sums.AM.LG$misses.lag)
lag.ms.sums.AM.LG$precision<-lag.ms.sums.AM.LG$hits.lag/(lag.ms.sums.AM.LG$hits.lag+lag.ms.sums.AM.LG$fp.lag)

lag.ms.sums.AM.CG<- ddply(ms.recovery.sum.AM.CG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ms.sums.AM.CG$recall<-lag.ms.sums.AM.CG$hits.lag/(lag.ms.sums.AM.CG$hits.lag+lag.ms.sums.AM.CG$misses.lag)
lag.ms.sums.AM.CG$precision<-lag.ms.sums.AM.CG$hits.lag/(lag.ms.sums.AM.CG$hits.lag+lag.ms.sums.AM.CG$fp.lag)

# output for Supplemental Table of Lagged relations:

table.lagged<-data.frame(c(lag.ar.sums.AM.BA$ar,lag.ar.sums.AM.BA$ar),
                   c(lag.ar.sums.AM.BA$recall,lag.ar.sums.AM.BA$precision),
                   c(lag.ar.sums.AM.LG$recall,lag.ar.sums.AM.LG$precision),
                   c(lag.ar.sums.AM.CG$recall,lag.ar.sums.AM.CG$precision),
                   c(lag.ms.sums.AM.BA$recall,lag.ms.sums.AM.BA$precision),
                   c(lag.ms.sums.AM.LG$recall,lag.ms.sums.AM.LG$precision),
                   c(lag.ms.sums.AM.CG$recall,lag.ms.sums.AM.CG$precision)
) 

colnames(table.lagged)<-c("auto.lag","ar.BA","ar.LG","ar.CG","ms.BA","ms.LG","ms.CG")

table.lagged$diff.BA<-table.lagged$ms.BA-table.lagged$ar.BA
table.lagged$diff.LG<-table.lagged$ms.LG-table.lagged$ar.LG
table.lagged$diff.CG<-table.lagged$ms.CG-table.lagged$ar.CG

# write out to paste into table:
write.csv(table.lagged,row.names = FALSE,file="Supp_lagged_table_output.csv")



# Table 2 to report effects of relative con/lag strength and autoregressive strength
table2<-data.frame(con.ar.sums.AM.BA$tsl,con.ar.sums.AM.BA$ar)
colnames(table2)<-c("tsl","ar")
table2$AR.recall<-(con.ar.sums.AM.LG$recall-con.ar.sums.AM.CG$recall)
table2$MS.recall<-(con.ms.sums.AM.LG$recall-con.ms.sums.AM.CG$recall)
table2$AR.prec<-(con.ar.sums.AM.LG$precision-con.ar.sums.AM.CG$precision)
table2$MS.prec<-(con.ms.sums.AM.LG$precision-con.ms.sums.AM.CG$precision)
write.csv(table2,row.names = FALSE,file="Table2_output.csv")

########################################################################
######## Recovery Tables: Plausible Models Only (Supplemental) #########
########################################################################


# summary data frames 
con.ar.sums.GM.BA<- ddply(ar.recovery.sum.GM.BA, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ar.sums.GM.BA$recall<-con.ar.sums.GM.BA$hits.con/(con.ar.sums.GM.BA$hits.con+con.ar.sums.GM.BA$misses.con)
con.ar.sums.GM.BA$precision<-con.ar.sums.GM.BA$hits.con/(con.ar.sums.GM.BA$hits.con+con.ar.sums.GM.BA$fp.con)

con.ar.sums.GM.LG<- ddply(ar.recovery.sum.GM.LG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ar.sums.GM.LG$recall<-con.ar.sums.GM.LG$hits.con/(con.ar.sums.GM.LG$hits.con+con.ar.sums.GM.LG$misses.con)
con.ar.sums.GM.LG$precision<-con.ar.sums.GM.LG$hits.con/(con.ar.sums.GM.LG$hits.con+con.ar.sums.GM.LG$fp.con)

con.ar.sums.GM.CG<- ddply(ar.recovery.sum.GM.CG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ar.sums.GM.CG$recall<-con.ar.sums.GM.CG$hits.con/(con.ar.sums.GM.CG$hits.con+con.ar.sums.GM.CG$misses.con)
con.ar.sums.GM.CG$precision<-con.ar.sums.GM.CG$hits.con/(con.ar.sums.GM.CG$hits.con+con.ar.sums.GM.CG$fp.con)


con.ms.sums.GM.BA<- ddply(ms.recovery.sum.GM.BA, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ms.sums.GM.BA$recall<-con.ms.sums.GM.BA$hits.con/(con.ms.sums.GM.BA$hits.con+con.ms.sums.GM.BA$misses.con)
con.ms.sums.GM.BA$precision<-con.ms.sums.GM.BA$hits.con/(con.ms.sums.GM.BA$hits.con+con.ms.sums.GM.BA$fp.con)

con.ms.sums.GM.LG<- ddply(ms.recovery.sum.GM.LG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ms.sums.GM.LG$recall<-con.ms.sums.GM.LG$hits.con/(con.ms.sums.GM.LG$hits.con+con.ms.sums.GM.LG$misses.con)
con.ms.sums.GM.LG$precision<-con.ms.sums.GM.LG$hits.con/(con.ms.sums.GM.LG$hits.con+con.ms.sums.GM.LG$fp.con)

con.ms.sums.GM.CG<- ddply(ms.recovery.sum.GM.CG, c("tsl", "ar"), summarise,
                          hits.con    = sum(dir.hits.con,na.rm = TRUE),
                          misses.con    = sum(dir.misses.con,na.rm = TRUE),
                          cr.con    = sum(dir.cr.con,na.rm = TRUE),
                          fp.con    = sum(dir.fp.con,na.rm = TRUE)
)
con.ms.sums.GM.CG$recall<-con.ms.sums.GM.CG$hits.con/(con.ms.sums.GM.CG$hits.con+con.ms.sums.GM.CG$misses.con)
con.ms.sums.GM.CG$precision<-con.ms.sums.GM.CG$hits.con/(con.ms.sums.GM.CG$hits.con+con.ms.sums.GM.CG$fp.con)

# output for contemp. relations considering only plausible models:

table.con.plausible<-data.frame(c(con.ar.sums.GM.BA$ar,con.ar.sums.GM.BA$ar),
                   c(con.ar.sums.GM.BA$recall,con.ar.sums.GM.BA$precision),
                   c(con.ar.sums.GM.LG$recall,con.ar.sums.GM.LG$precision),
                   c(con.ar.sums.GM.CG$recall,con.ar.sums.GM.CG$precision),
                   c(con.ms.sums.GM.BA$recall,con.ms.sums.GM.BA$precision),
                   c(con.ms.sums.GM.LG$recall,con.ms.sums.GM.LG$precision),
                   c(con.ms.sums.GM.CG$recall,con.ms.sums.GM.CG$precision)
) 

colnames(table.con.plausible)<-c("auto.con","ar.BA","ar.LG","ar.CG","ms.BA","ms.LG","ms.CG")

table.con.plausible$diff.BA<-table.con.plausible$ms.BA-table.con.plausible$ar.BA
table.con.plausible$diff.LG<-table.con.plausible$ms.LG-table.con.plausible$ar.LG
table.con.plausible$diff.CG<-table.con.plausible$ms.CG-table.con.plausible$ar.CG

# write out to paste into table:
write.csv(table.con.plausible,row.names = FALSE,file="table_con_plausible_output.csv")



# summary data frames 
lag.ar.sums.GM.BA<- ddply(ar.recovery.sum.GM.BA, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ar.sums.GM.BA$recall<-lag.ar.sums.GM.BA$hits.lag/(lag.ar.sums.GM.BA$hits.lag+lag.ar.sums.GM.BA$misses.lag)
lag.ar.sums.GM.BA$precision<-lag.ar.sums.GM.BA$hits.lag/(lag.ar.sums.GM.BA$hits.lag+lag.ar.sums.GM.BA$fp.lag)

lag.ar.sums.GM.LG<- ddply(ar.recovery.sum.GM.LG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ar.sums.GM.LG$recall<-lag.ar.sums.GM.LG$hits.lag/(lag.ar.sums.GM.LG$hits.lag+lag.ar.sums.GM.LG$misses.lag)
lag.ar.sums.GM.LG$precision<-lag.ar.sums.GM.LG$hits.lag/(lag.ar.sums.GM.LG$hits.lag+lag.ar.sums.GM.LG$fp.lag)

lag.ar.sums.GM.CG<- ddply(ar.recovery.sum.GM.CG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ar.sums.GM.CG$recall<-lag.ar.sums.GM.CG$hits.lag/(lag.ar.sums.GM.CG$hits.lag+lag.ar.sums.GM.CG$misses.lag)
lag.ar.sums.GM.CG$precision<-lag.ar.sums.GM.CG$hits.lag/(lag.ar.sums.GM.CG$hits.lag+lag.ar.sums.GM.CG$fp.lag)


lag.ms.sums.GM.BA<- ddply(ms.recovery.sum.GM.BA, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ms.sums.GM.BA$recall<-lag.ms.sums.GM.BA$hits.lag/(lag.ms.sums.GM.BA$hits.lag+lag.ms.sums.GM.BA$misses.lag)
lag.ms.sums.GM.BA$precision<-lag.ms.sums.GM.BA$hits.lag/(lag.ms.sums.GM.BA$hits.lag+lag.ms.sums.GM.BA$fp.lag)

lag.ms.sums.GM.LG<- ddply(ms.recovery.sum.GM.LG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ms.sums.GM.LG$recall<-lag.ms.sums.GM.LG$hits.lag/(lag.ms.sums.GM.LG$hits.lag+lag.ms.sums.GM.LG$misses.lag)
lag.ms.sums.GM.LG$precision<-lag.ms.sums.GM.LG$hits.lag/(lag.ms.sums.GM.LG$hits.lag+lag.ms.sums.GM.LG$fp.lag)

lag.ms.sums.GM.CG<- ddply(ms.recovery.sum.GM.CG, c("tsl", "ar"), summarise,
                          hits.lag    = sum(dir.hits.lag,na.rm = TRUE),
                          misses.lag    = sum(dir.misses.lag,na.rm = TRUE),
                          cr.lag    = sum(dir.cr.lag,na.rm = TRUE),
                          fp.lag    = sum(dir.fp.lag,na.rm = TRUE)
)
lag.ms.sums.GM.CG$recall<-lag.ms.sums.GM.CG$hits.lag/(lag.ms.sums.GM.CG$hits.lag+lag.ms.sums.GM.CG$misses.lag)
lag.ms.sums.GM.CG$precision<-lag.ms.sums.GM.CG$hits.lag/(lag.ms.sums.GM.CG$hits.lag+lag.ms.sums.GM.CG$fp.lag)

# output for Supplemental Table of Lagged relations:

table.lagged.plausible<-data.frame(c(lag.ar.sums.GM.BA$ar,lag.ar.sums.GM.BA$ar),
                         c(lag.ar.sums.GM.BA$recall,lag.ar.sums.GM.BA$precision),
                         c(lag.ar.sums.GM.LG$recall,lag.ar.sums.GM.LG$precision),
                         c(lag.ar.sums.GM.CG$recall,lag.ar.sums.GM.CG$precision),
                         c(lag.ms.sums.GM.BA$recall,lag.ms.sums.GM.BA$precision),
                         c(lag.ms.sums.GM.LG$recall,lag.ms.sums.GM.LG$precision),
                         c(lag.ms.sums.GM.CG$recall,lag.ms.sums.GM.CG$precision)
) 

colnames(table.lagged.plausible)<-c("auto.lag","ar.BA","ar.LG","ar.CG","ms.BA","ms.LG","ms.CG")

table.lagged.plausible$diff.BA<-table.lagged.plausible$ms.BA-table.lagged.plausible$ar.BA
table.lagged.plausible$diff.LG<-table.lagged.plausible$ms.LG-table.lagged.plausible$ar.LG
table.lagged.plausible$diff.CG<-table.lagged.plausible$ms.CG-table.lagged.plausible$ar.CG

# write out to paste into table:
write.csv(table.lagged.plausible,row.names = FALSE,file="table_lagged_plausible_output.csv")

# Same as CG/LG comparison table as above but with plausible models only (for Supplemental)
table2.plausible<-data.frame(con.ar.sums.GM.BA$tsl,con.ar.sums.GM.BA$ar)
colnames(table2.plausible)<-c("tsl","ar")
table2.plausible$AR.recall<-(con.ar.sums.GM.LG$recall-con.ar.sums.GM.CG$recall)
table2.plausible$MS.recall<-(con.ms.sums.GM.LG$recall-con.ms.sums.GM.CG$recall)
table2.plausible$AR.prec<-(con.ar.sums.GM.LG$precision-con.ar.sums.GM.CG$precision)
table2.plausible$MS.prec<-(con.ms.sums.GM.LG$precision-con.ms.sums.GM.CG$precision)
write.csv(table2.plausible,row.names = FALSE,file="Table2_plausible_output.csv")


##############################################################
########  Direction Recovery Plots: Main Manuscript ##########
##############################################################

require(ggplot2)
require(grid)
require(gridExtra)

# Contemp. recall in the balanced condition (Figure 2)
tmp.d<-ar.recovery.sum.AM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of directions recalled in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.AM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of directions recalled in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Figure_2.tif",g,units = "in",dpi=300,
       width = 10,height = 8,device="tiff")



# Contemp. precision in the balanced condition (Figure 4)
tmp.d<-ar.recovery.sum.AM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of true directions in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.AM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of true directions in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Figure_4.tif",g,units = "in",dpi=300,
       width = 10,height = 8,device="tiff")



# Compare MS and AR contemp recall in all conditions (Figure 3)

tmp.d<-compare.AM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",y="MS-AR proportion recalled")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5))) 


tmp.d<-compare.AM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion recalled")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 

tmp.d<-compare.AM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot7<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion recalled")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot8<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot9<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 


g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               plot7,plot8,plot9,
               left=textGrob("Con. Greater                            Balanced                            Lag. Greater",
                             rot=90,gp = gpar(fontsize=22), vjust=.3), ncol=3)
ggsave(filename = "Figure_3.tif",g,units = "in",dpi=300,
       width = 11,height = 12,device="tiff")


# Compare MS and AR contemp. precision in all conditions (Figure 5)
tmp.d<-compare.AM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",y="MS-AR proportion true directions")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5))) 


tmp.d<-compare.AM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion true directions")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 

tmp.d<-compare.AM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot7<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion true directions")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot8<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot9<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 


g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               plot7,plot8,plot9,
               left=textGrob("Con. Greater                            Balanced                            Lag. Greater",
                             rot=90,gp = gpar(fontsize=22), vjust=.3), ncol=3)
ggsave(filename = "Figure_5.tif",g,units = "in",dpi=300,
       width = 11,height = 12,device="tiff")






#######################################################################
########  Direction Recovery Plots: LG and CG (Supplemental) ##########
#######################################################################

# Contemp. recall in the LG condition 
tmp.d<-ar.recovery.sum.AM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of directions recalled in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.AM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of directions recalled in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Recall_Con_LG.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")



# Contemp. precision in the LG condition 
tmp.d<-ar.recovery.sum.AM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of true directions in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.AM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of true directions in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Precision_Con_LG.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")



# Contemp. recall in the CG condition 
tmp.d<-ar.recovery.sum.AM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of directions recalled in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.AM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of directions recalled in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Recall_Con_CG.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")



# Contemp. precision in the CG condition 
tmp.d<-ar.recovery.sum.AM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of true directions in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.AM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of true directions in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Precision_Con_CG.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")


############################################################################
########  Direction Recovery Plots: Plausible Models Only (Supplemental) ###
############################################################################


# Contemp. recall in the balanced condition, plausible models only
tmp.d<-ar.recovery.sum.GM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of directions recalled in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.GM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of directions recalled in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Recall_Con_Balanced_plausible.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")



# Contemp. precision in the balanced condition, plausible models only
tmp.d<-ar.recovery.sum.GM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of true directions in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.GM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of true directions in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Precision_Con_Balanced_plausible.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")





# Contemp. recall in the LG condition, plausible models only
tmp.d<-ar.recovery.sum.GM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of directions recalled in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.GM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of directions recalled in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Recall_Con_LG_plausible.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")



# Contemp. precision in the LG condition, plausible models only
tmp.d<-ar.recovery.sum.GM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of true directions in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.GM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of true directions in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Precision_Con_LG_plausible.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")



# Contemp. recall in the CG condition, plausible models only
tmp.d<-ar.recovery.sum.GM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of directions recalled in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.GM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of directions recalled in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Recall_Con_CG_plausible.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")



# Contemp. precision in the CG condition, plausible models only
tmp.d<-ar.recovery.sum.GM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",
                                                 y="Proportion of true directions in AR model")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Greens")+theme_bw()+ylim(0,1)+
          theme(plot.title = element_text(hjust = 0.5))) 
tmp.d<-ms.recovery.sum.GM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",
                                                 y="Proportion of true directions in MS model")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Blues")+theme_bw()+ylim(0,1)) 

grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               ncol=3)
ggsave(filename = "Precision_Con_CG_plausible.jpg",g,units = "in",dpi=300,
       width = 10,height = 8,device="jpeg")







# Compare MS and AR contemp recall in all conditions, plausible models only

tmp.d<-compare.GM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",y="MS-AR proportion recalled")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5))) 


tmp.d<-compare.GM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion recalled")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 

tmp.d<-compare.GM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot7<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion recalled")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot8<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot9<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.recall.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 

# grid.arrange(plot1,plot2,plot3,
#              plot4,plot5,plot6,
#              plot7,plot8,plot9,
#              ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               plot7,plot8,plot9,
               left=textGrob("Con. Greater                            Balanced                            Lag. Greater",
                             rot=90,gp = gpar(fontsize=22), vjust=.3), ncol=3)
ggsave(filename = "Recall_con_differences_plausible.jpg",g,units = "in",dpi=300,
       width = 11,height = 12,device="jpeg")


# Compare MS and AR contemp precision in all conditions, plausible models only

tmp.d<-compare.GM.LG
tmp.d$ar<-as.factor(tmp.d$ar)
plot1<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="50 Time Points",x="autoregressive strength",y="MS-AR proportion true directions")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))
plot2<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="100 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5)))    
plot3<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="300 Time Points",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)+
          theme(plot.title = element_text(hjust = 0.5))) 


tmp.d<-compare.GM.BA
tmp.d$ar<-as.factor(tmp.d$ar)
plot4<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion true directions")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot5<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot6<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 

tmp.d<-compare.GM.CG
tmp.d$ar<-as.factor(tmp.d$ar)
plot7<-(ggplot(data=tmp.d[tmp.d$tsl==50,],aes(group=ar,x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="MS-AR proportion true directions")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))
plot8<-(ggplot(data=tmp.d[tmp.d$tsl==100,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8))    
plot9<-(ggplot(data=tmp.d[tmp.d$tsl==300,],aes(x=ar,y=dir.precision.con,fill=ar))+
          geom_boxplot(show.legend = FALSE)+labs(title="",x="autoregressive strength",y="")+
          scale_fill_brewer(palette="Reds")+theme_bw()+ylim(-.8,.8)) 

# grid.arrange(plot1,plot2,plot3,
#              plot4,plot5,plot6,
#              plot7,plot8,plot9,
#              ncol=3)
g<-arrangeGrob(plot1,plot2,plot3,
               plot4,plot5,plot6,
               plot7,plot8,plot9,
               left=textGrob("Con. Greater                            Balanced                            Lag. Greater",
                             rot=90,gp = gpar(fontsize=22), vjust=.3), ncol=3)
ggsave(filename = "Precision_con_differences_plausible.jpg",g,units = "in",dpi=300,
       width = 11,height = 12,device="jpeg")





