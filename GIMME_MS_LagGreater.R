rm(list=ls())

#load relevant packages
library('gimme')
library('doMC')
library('foreach')

#list directories
reps<-list.files(path = "/Sim_Study/lag_greater_in/data/")

#Set up number of cores for parallel processing 

registerDoMC(cores = 22)

foreach (r=reps,.errorhandling="remove") %dopar%
{

# Run 
g.out<-gimmeSEM(data = paste("/Sim_Study/lag_greater_in/data/",r,sep=""),
         out = paste("/Sim_Study/lag_greater_in/outputMS/",r,sep=""),
         sep = ",", #comma delimited
         header = TRUE,
         ar = FALSE,
         ms_allow = TRUE,
         plot = FALSE,
         subgroup = FALSE,
         paths = NULL,
         sub_feature = "lag & contemp",
         confirm_subgroup = FALSE,
         conv_vars = NULL,
         conv_length = 16,
         conv_interval = 2,
         mult_vars = NULL,
         mean_center_mult = FALSE,
         standardize = TRUE,
         groupcutoff = .75,
         subcutoff = .5,
         diagnos = FALSE)

save(g.out,file=paste("/Sim_Study/lag_greater_in/outputMS/",r,"/out.RData",sep=""))

}

 

