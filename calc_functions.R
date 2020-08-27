# These functions allow you to quickly generate summary 
# direction recovery statistics. They are used in the following
# scripts: 

# analyze_recovery_balanced.R
# analyze_recovery_con_greater.R
# analyze_recovery_lag_greater.R

# The functions below produce identical output, except the
# ".noex" function does not exclude cases with impluasible 
# parameter estimates

calc.rec.stats<-function(sum,beta,nacol){
  
  #name recall and precision variables for each
  #simulated participant
  
  sum$dir.hits.all<-NA
  sum$dir.misses.all<-NA
  sum$dir.cr.all<-NA
  sum$dir.fp.all<-NA
  
  sum$dir.hits.con<-NA
  sum$dir.misses.con<-NA
  sum$dir.cr.con<-NA
  sum$dir.fp.con<-NA
  
  sum$dir.hits.lag<-NA
  sum$dir.misses.lag<-NA
  sum$dir.cr.lag<-NA
  sum$dir.fp.lag<-NA
  
  sum$dir.hits.grp<-NA
  sum$dir.misses.grp<-NA
  
  sum$dir.hits.ind<-NA
  sum$dir.misses.ind<-NA
  
  # now, generate the relevant stats
  for (r in 1:length(sum$rep)){
    
    if((!is.na(sum[r,nacol]) && (sum[r,nacol]==FALSE))){
      

      #take simulated beta matrix, make it binary
      tmp.true.all<-as.matrix(true.type.con[[r]])
      tmp.true.all[is.na(tmp.true.all)]<-0
      tmp.true.all[tmp.true.all!=0]<-1
      class(tmp.true.all)<-"numeric"
      # make all AR paths 1 because we are not measuring their recovery
      diag(tmp.true.all)<-1
      
      #binary simulated matrix of individual paths only 
      tmp.true.ind<-as.matrix(true.type.con[[r]])
      tmp.true.ind[is.na(tmp.true.ind) | tmp.true.ind=="grp"]<-0
      tmp.true.ind[tmp.true.ind=="ind"]<-1
      class(tmp.true.ind)<-"numeric"
      
      #binary simulated matrix of group paths only 
      tmp.true.grp<-as.matrix(true.type.con[[r]])
      tmp.true.grp[is.na(tmp.true.grp) | tmp.true.grp=="ind"]<-0
      tmp.true.grp[tmp.true.grp=="grp"]<-1
      class(tmp.true.grp)<-"numeric"
      # make all AR paths 1 because we are not measuring their recovery
      diag(tmp.true.grp)<-1
      
      #binary matrix of recovered paths
      tmp.rec<-as.matrix(as.data.frame(beta[[r]]))
      tmp.rec[tmp.rec!=0]<-1
      # fill in all of the AR paths even if they weren't recovered, 
      # as we aren't measuring their recovery
      diag(tmp.rec[,])<-1
      class(tmp.rec)<-"numeric"
      
      
      #overall direction recall/precision
      diff.mat<-tmp.true.all-(tmp.rec*10)
      sum$dir.hits.all[r]=(sum(diff.mat==-9)-8)
      sum$dir.misses.all[r]=(sum(diff.mat==1))
      sum$dir.cr.all[r]=(sum(diff.mat==0)-8)
      sum$dir.fp.all[r]=(sum(diff.mat==-10))
      
      #direction recall for group paths
      diff.mat<-tmp.true.grp-(tmp.rec*10)
      sum$dir.hits.grp[r]=(sum(diff.mat==-9)-8)
      sum$dir.misses.grp[r]=(sum(diff.mat==1))
      
      #direction recall for individual paths 
      diff.mat<-tmp.true.ind-(tmp.rec*10)
      sum$dir.hits.ind[r]=(sum(diff.mat==-9))
      sum$dir.misses.ind[r]=(sum(diff.mat==1))
      
      #direction recall/precision for contemporaneous paths
      diff.mat<-tmp.true.all-(tmp.rec*10)
      diff.mat<-diff.mat[,9:16]
      sum$dir.hits.con[r]<-(sum(diff.mat==-9))
      sum$dir.misses.con[r]<-(sum(diff.mat==1))
      sum$dir.cr.con[r]<-(sum(diff.mat==0)-8)
      sum$dir.fp.con[r]<-(sum(diff.mat==-10))
      
      #direction recall/precision for lagged paths
      diff.mat<-tmp.true.all-(tmp.rec*10)
      diff.mat<-diff.mat[,1:8]
      sum$dir.hits.lag[r]<-(sum(diff.mat==-9)-8)
      sum$dir.misses.lag[r]<-(sum(diff.mat==1))
      sum$dir.cr.lag[r]<-(sum(diff.mat==0))
      sum$dir.fp.lag[r]<-(sum(diff.mat==-10))
      
    }
  }
  
  
  # recovery indices for each person
  sum$dir.recall.all<-sum$dir.hits.all/(sum$dir.hits.all+sum$dir.misses.all)
  sum$dir.precision.all<-sum$dir.hits.all/(sum$dir.hits.all+sum$dir.fp.all)
  
  sum$dir.recall.grp<-sum$dir.hits.grp/(sum$dir.hits.grp+sum$dir.misses.grp)
  sum$dir.recall.ind<-sum$dir.hits.ind/(sum$dir.hits.ind+sum$dir.misses.ind)
  
  sum$dir.recall.con<-sum$dir.hits.con/(sum$dir.hits.con+sum$dir.misses.con)
  sum$dir.precision.con<-sum$dir.hits.con/(sum$dir.hits.con+sum$dir.fp.con)
  
  sum$dir.recall.lag<-sum$dir.hits.lag/(sum$dir.hits.lag+sum$dir.misses.lag)
  sum$dir.precision.lag<-sum$dir.hits.lag/(sum$dir.hits.lag+sum$dir.fp.lag)
  
  out<-sum
  
  out
  
}

calc.rec.stats.noex<-function(sum,beta,nacol){
  
  
  #name recall and precision variables for each
  #simulated participant
  
  sum$dir.hits.all<-NA
  sum$dir.misses.all<-NA
  sum$dir.cr.all<-NA
  sum$dir.fp.all<-NA
  
  sum$dir.hits.con<-NA
  sum$dir.misses.con<-NA
  sum$dir.cr.con<-NA
  sum$dir.fp.con<-NA
  
  sum$dir.hits.lag<-NA
  sum$dir.misses.lag<-NA
  sum$dir.cr.lag<-NA
  sum$dir.fp.lag<-NA
  
  sum$dir.hits.grp<-NA
  sum$dir.misses.grp<-NA
  
  sum$dir.hits.ind<-NA
  sum$dir.misses.ind<-NA
  
  # now, generate the relevant stats 
  for (r in 1:length(sum$rep)){
    
    if((!is.na(sum[r,nacol]))){
      
      
      #take simulated beta matrix, make it binary
      tmp.true.all<-as.matrix(true.type.con[[r]])
      tmp.true.all[is.na(tmp.true.all)]<-0
      tmp.true.all[tmp.true.all!=0]<-1
      class(tmp.true.all)<-"numeric"
      # make all AR paths 1 because we are not measuring their recovery
      diag(tmp.true.all)<-1
      
      #binary simulated matrix of individual paths only 
      tmp.true.ind<-as.matrix(true.type.con[[r]])
      tmp.true.ind[is.na(tmp.true.ind) | tmp.true.ind=="grp"]<-0
      tmp.true.ind[tmp.true.ind=="ind"]<-1
      class(tmp.true.ind)<-"numeric"
      
      #binary simulated matrix of group paths only 
      tmp.true.grp<-as.matrix(true.type.con[[r]])
      tmp.true.grp[is.na(tmp.true.grp) | tmp.true.grp=="ind"]<-0
      tmp.true.grp[tmp.true.grp=="grp"]<-1
      class(tmp.true.grp)<-"numeric"
      # make all AR paths 1 because we are not measuring their recovery
      diag(tmp.true.grp)<-1
      
      #binary matrix of recovered paths
      tmp.rec<-as.matrix(as.data.frame(beta[[r]]))
      tmp.rec[tmp.rec!=0]<-1
      # fill in all of the AR paths even if they weren't recovered, 
      # as we aren't measuring their recovery
      diag(tmp.rec[,])<-1
      class(tmp.rec)<-"numeric"
      
      
      #overall direction recall/precision
      diff.mat<-tmp.true.all-(tmp.rec*10)
      sum$dir.hits.all[r]=(sum(diff.mat==-9)-8)
      sum$dir.misses.all[r]=(sum(diff.mat==1))
      sum$dir.cr.all[r]=(sum(diff.mat==0)-8)
      sum$dir.fp.all[r]=(sum(diff.mat==-10))
      
      #direction recall for group paths
      diff.mat<-tmp.true.grp-(tmp.rec*10)
      sum$dir.hits.grp[r]=(sum(diff.mat==-9)-8)
      sum$dir.misses.grp[r]=(sum(diff.mat==1))
      
      #direction recall for individual paths 
      diff.mat<-tmp.true.ind-(tmp.rec*10)
      sum$dir.hits.ind[r]=(sum(diff.mat==-9))
      sum$dir.misses.ind[r]=(sum(diff.mat==1))
      
      #direction recall/precision for contemporaneous paths
      diff.mat<-tmp.true.all-(tmp.rec*10)
      diff.mat<-diff.mat[,9:16]
      sum$dir.hits.con[r]<-(sum(diff.mat==-9))
      sum$dir.misses.con[r]<-(sum(diff.mat==1))
      sum$dir.cr.con[r]<-(sum(diff.mat==0)-8)
      sum$dir.fp.con[r]<-(sum(diff.mat==-10))
      
      #direction recall/precision for lagged paths
      diff.mat<-tmp.true.all-(tmp.rec*10)
      diff.mat<-diff.mat[,1:8]
      sum$dir.hits.lag[r]<-(sum(diff.mat==-9)-8)
      sum$dir.misses.lag[r]<-(sum(diff.mat==1))
      sum$dir.cr.lag[r]<-(sum(diff.mat==0))
      sum$dir.fp.lag[r]<-(sum(diff.mat==-10))
      
    }
  }
  
  
  # recovery indices for each person
  sum$dir.recall.all<-sum$dir.hits.all/(sum$dir.hits.all+sum$dir.misses.all)
  sum$dir.precision.all<-sum$dir.hits.all/(sum$dir.hits.all+sum$dir.fp.all)
  
  sum$dir.recall.grp<-sum$dir.hits.grp/(sum$dir.hits.grp+sum$dir.misses.grp)
  sum$dir.recall.ind<-sum$dir.hits.ind/(sum$dir.hits.ind+sum$dir.misses.ind)
  
  sum$dir.recall.con<-sum$dir.hits.con/(sum$dir.hits.con+sum$dir.misses.con)
  sum$dir.precision.con<-sum$dir.hits.con/(sum$dir.hits.con+sum$dir.fp.con)
  
  sum$dir.recall.lag<-sum$dir.hits.lag/(sum$dir.hits.lag+sum$dir.misses.lag)
  sum$dir.precision.lag<-sum$dir.hits.lag/(sum$dir.hits.lag+sum$dir.fp.lag)
  
  out<-sum
  
  out
  
}