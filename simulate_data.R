rm(list=ls())

# This script, large parts of which were adapted from scripts
# used in the previous Lane et al. (2019)  simulation study published
# in Psychological Methods (preprint openly available at: https://psyarxiv.com/d3nk9/), 
# generates simulated uSEM model data for all conditions in the current study

#######################################################
#### 1. load functions needed to simulate the data ####
#######################################################

# This function generates group-level relations

# input parameters:
# p.con = proportion of paths that are contemporaneous (vs. lagged)
# nvar = number of variables in the time series
# AR = average strength of autoregressive realtions (SD=.10)
# dens = density of network (not including autoregressive relations)
# p.group = proportion of relations at the group (vs. individual) level
# con.b = average strength of contemporaneous relations (SD=.10)
# lag.b = average strength of cross-lagged relations (SD=.10)


mat.generate.asw <- function(p.con, nvar,AR,dens,p.group,con.b,lag.b){
  repeat{
  A   <- matrix(0, ncol = nvar, nrow = nvar) # contemporaneous
  Phi <- matrix(0, ncol = nvar, nrow = nvar) # lagged (AR in diag)
  pos.all <- 2*(nrow(A)*ncol(A)) - 2*nrow(A) # all available (non-ar) spots
  cnt.all <- dens*p.group*pos.all # number of group paths given density and proportion of group paths
  indices <- which(Phi == 0, arr.ind = TRUE)
  indices <- indices[which(indices[,1] != indices[,2]), ]# kick out diag
  row.col <- sample(1:nrow(indices), cnt.all, replace = F)# sample n=cnt.all paths 
  # group paths (in case there are an odd number, randomly chose whether con or lag gets one extra)
  n.p.1<-row.col[1:round(p.con*length(row.col))]
  n.p.2<-row.col[(length(n.p.1)+1):length(row.col)]
  n.p<-list(n.p.1,n.p.2)
  rand<-sample(1:2,size=2,replace=FALSE)
  grp.con <- n.p[[ rand[1] ]]
  grp.lag <- n.p[[ rand[2] ]]

  Phi[indices[grp.lag,]] <- lag.b
  A[indices[grp.con,]]   <- con.b
  diag(Phi)              <- AR # Insert AR terms here!
  
  # disallow bidirectional contemporaneous paths and A matrices with max eigenvalues >1
  A.test<-1*(A!=0)
  A.test<-(A.test+t(A.test)) 
  if ( (max(A.test)!=2) & (max(abs(eigen(A, only.values = FALSE)$values))<1) ) break
  
  }
  
  all <- cbind(Phi, A)
  ind.pres <- which(all != 0, arr.ind = T)
  ind.pres <- ind.pres[which(ind.pres[,1] != ind.pres[,2]), ]
  
  level <- "grp"
  
  all.lvl           <- matrix(NA, ncol = ncol(all), nrow = nrow(all))
  all.lvl[ind.pres] <- level
  diag(all.lvl)     <- "grp"
  
  all_sub1 <- all
  all_lvl1 <- all.lvl
  
  res <- list(sub1 = all_sub1, 
              lvl1 = all_lvl1)
  return(res)
}



# This function generates group-level relations adds individual-level
# relations to the matrix, adds noise, and simulates the time series

# input parameters:
# mat = matrix of group-level relations
# lvl = matric of input relation levels (e.g., "group")
# p.con = proportion of paths that are contemporaneous (vs. lagged)
# p.group = proportion of relations at the group (vs. individual) level
# con.b = average strength of contemporaneous relations (SD=.10)
# lag.b = average strength of cross-lagged relations (SD=.10)

ts.generate.asw <- function (mat, lvl, t,dens,p.group,con.b,lag.b,p.con) {
  repeat {
    
    repeat{
    v <- ncol(mat)/2 #calc nvars from matrix
    Phi <- mat[, 1:v] # pull out Phi
    A   <- mat[, (v+1):(v*2)] # pull out A
    A_ind        <- matrix(0, ncol=v, nrow=v)
    ## finds indices of zero elements in matrix
    indices.A      <- which(A==0,arr.ind=T)
    indices.Phi      <- which(Phi==0,arr.ind=T)
    ## removes indices of diagonal elements from consideration
    indices.A      <- indices.A[which(indices.A[,1]!=indices.A[,2]),]
    indices.Phi      <- indices.Phi[which(indices.Phi [,1]!=indices.Phi[,2]),]
    
    ## determines the number of individual paths to add
    # (in case there are an odd number, randomly chose whether con or lag gets one extra)
    pos.all<-(v*v-v)*2
    cnt.all <- dens*p.group*pos.all  
    rand<-sample(c(round(cnt.all/2),(round(cnt.all)-round(cnt.all/2))),size=2,replace=FALSE)
    ## randomly selects row/col combinations for individual-level paths
    row.col.A      <- sample(1:nrow(indices.A), rand[1], replace = F) 
    row.col.Phi      <- sample(1:nrow(indices.Phi), rand[2], replace = F) 
    
    # Betas for lagged and contemporaneous paths
    Phi[indices.Phi[row.col.Phi,]] <- lag.b
    A[indices.A[row.col.A,]]     <- con.b
    
    # add noise to A betas, SD =.1
    noise.inds      <- which(A != 0, arr.ind = TRUE) 
    A[noise.inds]   <- A[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .1)
    
    # add noise to Phi betas, excluding ar terms, SD = .1
    noise.inds      <- which(Phi != 0, arr.ind = TRUE) 
    noise.inds      <- noise.inds[which(noise.inds[,1] != noise.inds[,2]), ]# kick out diag
    Phi[noise.inds] <- Phi[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .1)
    
    # add noise to AR terms, SD = .1 (if you want)
    noise.inds      <- which(Phi !=Inf, arr.ind = TRUE)    
    noise.inds      <- noise.inds[which(noise.inds[,1] == noise.inds[,2]), ]# include diag
    Phi[noise.inds] <- Phi[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .1)
    
    #disallow bidirectional contemporaneous paths and and A matrices with max eigenvalues >1
    A.test<-1*(A!=0)
    A.test<-(A.test+t(A.test)) 
    if( (max(A.test)!=2) & (max(abs(eigen(A, only.values = FALSE)$values))<1) ) break
    }
    
          
    st <- (t+50) #Alex added, now robust to any t 
    noise <- matrix(rnorm(v*st,0,1),v) #
    I     <- diag(v) # identity matrix
    time  <- matrix(0,nrow=v, ncol=(st+1))
    time1 <- matrix(0,nrow=v, ncol=st)
    
    # simulate data points for each time step
    for (i in 1:st){
      time1[,i]  <- solve(I-A)%*%(Phi%*%time[,i] + noise[,i])
      time[,i+1] <- time1[,i]
    }               
    time1  <- time1[,(51:(50+t))] # Fixed this, was initially 50:, causing ts to be one too long
    series <- t(time1)
    paths  <- cbind(Phi, A)
  if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 
      & abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list   <- list("series"  = series,
                 "paths"   = paths,
                 "levels"  = lvl)
  return(list)
}



#######################################################
#### 2. Simulate Data: Balanced Condition #############
#######################################################

# enter simulation parameters
v             <- c(8) # Number of variables
n             <- c(50) # number of individuals
t             <- c(50, 100, 300) # Number of time points
rep           <- seq(1:100) # replications per condition 
ar            <-c(.0,.1,.3,.5,.6) # ar paths to try
conditions    <- expand.grid(t, n, v, ar,rep)
all           <- rbind(conditions)
colnames(all) <- c("t", "n", "v", "ar", "rep")

# This loop generates a folder name for each condition and replication
for (i in 1:nrow(all)){
  all$folder[i] <- paste("t",all$t[i],
                         "ar",all$ar[i],"rep",all$rep[i],sep="_")
}

rownames(all) <- NULL
colnames(all) <- c("t","n","v","ar","rep","folder")
all$t         <- as.numeric(as.character(all$t))
all$n         <- as.numeric(as.character(all$n))
all$ar         <- as.numeric(as.character(all$ar))


# name directories to place simulated data in (Change)
dir.create('/Sim_Study/balanced_in')
data.path <- '/Sim_Study/balanced_in/data'
true.path <- '/Sim_Study/balanced_in/true'
level.path <- '/Sim_Study/balanced_in/levels'

dir.create(data.path)
dir.create(true.path)
dir.create(level.path)

# Creates actual folders for each iteration
folders <- all$folder
for (i in 1:nrow(all)){
  data <- file.path(data.path, folders[i])
  dir.create(data)
  true <- file.path(true.path, folders[i])
  dir.create(true)
  level <- file.path(level.path, folders[i])
  dir.create(level)
}

# Does the simulations
for (i in 1:nrow(all)){
  if (!length(list.files(file.path(data.path,all$folder[i]))) %in% c(50)) {
    
    # generate group matrix for each simulated data set
    res <- mat.generate.asw(p.con = .50, 
                            nvar = all$v[i], AR=all$ar[i],
                            p.group = .50,dens = .20,
                            con.b = .3,lag.b = -.3)
    
    #for each individual,generate matrix and ts
    for (a in 1:all$n[i]){
      
      out <- ts.generate.asw(mat = res$sub1,
                             lvl = res$lvl1,
                             t   = all$t[i],
                             p.group = .50,dens = .20,
                             con.b = .3,lag.b = -.3,
                             p.con = .50)
      
      out$series <- round(out$series,digits=5)
      

      write.csv(out$series,
                file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)

      write.csv(out$paths,
                file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)

      write.csv(out$levels,
                file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
    }
    
  }
}


# make files for fit outputs
dir.create('/Sim_Study/balanced_in/outputAR')
dir.create('/Sim_Study/balanced_in/outputMS')

ms.path <- '/Sim_Study/balanced_in/outputAR'
ar.path <- '/Sim_Study/balanced_in/outputMS'

folders <- all$folder
for (i in 1:nrow(all)){
  ms <- file.path(ms.path, folders[i])
  dir.create(ms)
  ar <- file.path(ar.path, folders[i])
  dir.create(ar)
}



#######################################################
#### 3. Simulate Data: Con. Greater (CG) Condition ####
#######################################################

# enter simulation parameters
v             <- c(8) # Number of variables
n             <- c(50) # number of individuals
t             <- c(50, 100, 300) # Number of time points
rep           <- seq(1:100) # replications per condition 
ar            <-c(.0,.1,.3,.5,.6) # ar paths to try
conditions    <- expand.grid(t, n, v, ar,rep)
all           <- rbind(conditions)
colnames(all) <- c("t", "n", "v", "ar", "rep")

# This loop generates a folder name for each condition and replication
for (i in 1:nrow(all)){
  all$folder[i] <- paste("t",all$t[i],
                         "ar",all$ar[i],"rep",all$rep[i],sep="_")
}

rownames(all) <- NULL
colnames(all) <- c("t","n","v","ar","rep","folder")
all$t         <- as.numeric(as.character(all$t))
all$n         <- as.numeric(as.character(all$n))
all$ar         <- as.numeric(as.character(all$ar))


# name directories to place simulated data in
dir.create('/Sim_Study/con_greater_in')
data.path <- '/Sim_Study/con_greater_in/data'
true.path <- '/Sim_Study/con_greater_in/true'
level.path <- '/Sim_Study/con_greater_in/levels'

dir.create(data.path)
dir.create(true.path)
dir.create(level.path)

# Creates actual folders for each iteration
folders <- all$folder
for (i in 1:nrow(all)){
  data <- file.path(data.path, folders[i])
  dir.create(data)
  true <- file.path(true.path, folders[i])
  dir.create(true)
  level <- file.path(level.path, folders[i])
  dir.create(level)
}


# Does the simulations
for (i in 1:nrow(all)){
  if (!length(list.files(file.path(data.path,all$folder[i]))) %in% c(50)) {
    
    # generate group matrix for each simulated data set
    res <- mat.generate.asw(p.con = .50, 
                            nvar = all$v[i], AR=all$ar[i],
                            p.group = .50,dens = .20,
                            con.b = .6,lag.b = -.3)
    
    #for each individual,generate matrix and ts
    for (a in 1:all$n[i]){
      
      out <- ts.generate.asw(mat = res$sub1,
                             lvl = res$lvl1,
                             t   = all$t[i],
                             p.group = .50,dens = .20,
                             con.b = .6,lag.b = -.3,
                             p.con = .50)
      
      out$series <- round(out$series,digits=5)
      
      
      write.csv(out$series,
                file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$paths,
                file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$levels,
                file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
    }
    
  }
}

# make files for fit outputs
dir.create('/Sim_Study/con_greater_in/outputAR')
dir.create('/Sim_Study/con_greater_in/outputMS')

ms.path <- '/Sim_Study/con_greater_in/outputAR'
ar.path <- '/Sim_Study/con_greater_in/outputMS'

folders <- all$folder
for (i in 1:nrow(all)){
  ms <- file.path(ms.path, folders[i])
  dir.create(ms)
  ar <- file.path(ar.path, folders[i])
  dir.create(ar)
}


#######################################################
#### 4. Simulate Data: Lag. Greater (LG) Condition ####
#######################################################

# enter simulation parameters
v             <- c(8) # Number of variables
n             <- c(50) # number of individuals
t             <- c(50, 100, 300) # Number of time points
rep           <- seq(1:100) # replications per condition 
ar            <-c(.0,.1,.3,.5,.6) # ar paths to try
conditions    <- expand.grid(t, n, v, ar,rep)
all           <- rbind(conditions)
colnames(all) <- c("t", "n", "v", "ar", "rep")

# This loop generates a folder name for each condition and replication
for (i in 1:nrow(all)){
  all$folder[i] <- paste("t",all$t[i],
                         "ar",all$ar[i],"rep",all$rep[i],sep="_")
}

rownames(all) <- NULL
colnames(all) <- c("t","n","v","ar","rep","folder")
all$t         <- as.numeric(as.character(all$t))
all$n         <- as.numeric(as.character(all$n))
all$ar         <- as.numeric(as.character(all$ar))


# name directories to place simulated data in
dir.create('/Sim_Study/lag_greater_in')
data.path <- '/Sim_Study/lag_greater_in/data'
true.path <- '/Sim_Study/lag_greater_in/true'
level.path <- '/Sim_Study/lag_greater_in/levels'

dir.create(data.path)
dir.create(true.path)
dir.create(level.path)

# Creates actual folders for each iteration
folders <- all$folder
for (i in 1:nrow(all)){
  data <- file.path(data.path, folders[i])
  dir.create(data)
  true <- file.path(true.path, folders[i])
  dir.create(true)
  level <- file.path(level.path, folders[i])
  dir.create(level)
}

# Does the simulations
for (i in 1:nrow(all)){
  if (!length(list.files(file.path(data.path,all$folder[i]))) %in% c(50)) {
    
    # generate group matrix for each simulated data set
    res <- mat.generate.asw(p.con = .50, 
                            nvar = all$v[i], AR=all$ar[i],
                            p.group = .50,dens = .20,
                            con.b = .3,lag.b = -.6)
    
    #for each individual,generate matrix and ts
    for (a in 1:all$n[i]){
      
      out <- ts.generate.asw(mat = res$sub1,
                             lvl = res$lvl1,
                             t   = all$t[i],
                             p.group = .50,dens = .20,
                             con.b = .3,lag.b = -.6,
                             p.con = .50)
      
      out$series <- round(out$series,digits=5)
      
      
      write.csv(out$series,
                file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$paths,
                file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$levels,
                file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
    }
    
  }
}


# make files for fit outputs
dir.create('/Sim_Study/lag_greater_in/outputAR')
dir.create('/Sim_Study/lag_greater_in/outputMS')

ms.path <- '/Sim_Study/lag_greater_in/outputAR'
ar.path <- '/Sim_Study/lag_greater_in/outputMS'

folders <- all$folder
for (i in 1:nrow(all)){
  ms <- file.path(ms.path, folders[i])
  dir.create(ms)
  ar <- file.path(ar.path, folders[i])
  dir.create(ar)
}


