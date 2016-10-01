migrate_harvester<-function(wd,n=3,models,multilocus=F){ #wd<-"/Users/eric/Datasets/simulations/ibdsim/migrate_pareto2_mtdna_results"; n=3; models<-c("panmixia","5island","5stepping.stone","5stepping.stone1","10island","10stepping.stone","10stepping.stone1")
  
  setwd(wd)
  likelists<-list() #initialize an empty list
  
  for(r in list.files(wd)){
    wd1<-file.path(wd,r)
    setwd(wd1)
    marglike<-data.frame(model=character(0),thermodynamic=numeric(0),bezier.corrected=numeric(0),harmonic.mean=numeric(0),stringsAsFactors=F) #initialize a data frame to take the values
    l=1 #initialize l
    for(m in models){ #i<-"stepping.stone"
      wd2<-file.path(wd1,m)
      print(wd2)
      if(!file.exists(wd2)){next}
      setwd(wd2)
      outfile<-scan(file="outfile",what="character",sep="\n") #scan in the outfile, separating at each newline
      
      if(multilocus=F){
      #get the result from thermodynamic integration
      thermoline<-grep("(1a)",outfile,value=T) #find the line with the thermodynamic likelihood on it
      if(length(thermoline)==0){next}
      thermoline<-strsplit(thermoline,split="=",fixed=T) #split it up
      thermo<-as.numeric(substr(thermoline[[1]][2],start=1,stop=12)) #grab the thermodynamic likelihood
      bezistart<-grep("\\(",strsplit(thermoline[[1]][2],split="")[[1]])+1
      bezier<-as.numeric(substr(thermoline[[1]][2],start=bezistart,stop=bezistart+11)) #and the bezier-corrected value
      #get the harmonic mean
      harmoline<-grep("\\(2\\) H",outfile,value=T) #find the line with harmonic mean likelihood on it
      harmoline<-strsplit(harmoline,split="=",fixed=T) #grab the harmonic mean
      harmo<-as.numeric(substr(harmoline[[1]][2],start=1,stop=12))
      marglike[l,]<-c(m,thermo,bezier,harmo) #add this as a row to the data frame
      l=l+1
      }
      
      if(multilocus=T){
        all.locus.line<-outfile[grep("\\[Scaling",outfile,value=F)-1] #find the line with all three values on it (it comes right before the line that has the scaling factor on it)
        all.locus.line<-strsplit(all.locus.line," +")
        thermo<-as.numeric(all.locus.line[[1]][3])
        bezier<-as.numeric(all.locus.line[[1]][4])
        harmo<-as.numeric(all.locus.line[[1]][5])
        marglike[l,]<-c(m,thermo,bezier,harmo) #add this as a row to the data frame
        l=l+1
      }
    }
    likelists[[r]]<-marglike
  }
  setwd(wd)
  return(likelists)
}


bfcalcs<-function(df,ml="bezier.corrected"){
  df$thermodynamic<-as.numeric(df$thermodynamic)
  df$bezier.corrected<-as.numeric(df$bezier.corrected)
  df$harmonic<-as.numeric(df$harmonic)
  mlcol<-df[,ml] 
  bmvalue<-mlcol[which.max(mlcol)]
  lbf<-2*(mlcol-bmvalue)
  choice<-rank(-mlcol)
  modelprob<-exp(lbf/2)/sum(exp(lbf/2))
  dfall<-cbind(df,lbf,choice,modelprob)
  return(dfall)
}


getmodels<-function(dfr){
  model1<-dfr$model[which(dfr$choice==1)]
  model2<-dfr$model[which(dfr$choice==2)]
  modelprob1<-dfr$modelprob[which(dfr$choice==1)]
  modelprob2<-dfr$modelprob[which(dfr$choice==2)]
  
  print(model1)
  print(model2)
  print(modelprob1)
  c(model1,modelprob1,model2,modelprob2)
}

#modeltables must be a list of bfcalcs output
#ml_type = "bezier.corrected","thermodynamic", or "harmonic.mean"
bfcalcs_reps<-function(modeltables,models,ml_type="bezier.corrected"){
  
  require(perm)
  output<-list()
  
  likes<-rbind(cbind(modeltables[[1]],rep=1), cbind(modeltables[[2]],rep=2), cbind(modeltables[[3]], rep=3))
  likes$model<-factor(likes$model, models) #factor to turn missing models into NA
  likes<-likes[!(is.na(likes$model)),] #remove NAs 
  likes$model<-factor(likes$model, models) #and refactor
  mean.ml<-as.vector(by(likes[[ml_type]],likes$model,mean)) #take the means for each model
  sd.ml<-as.vector(by(likes[[ml_type]],likes$model,sd))
    
  bmvalue<-mean.ml[which.max(mean.ml)] #get the best mean marg like
  lbf<-2*(mean.ml-bmvalue) #get the log bayes factors
  choice<-rank(-mean.ml) # make a choice column
  modelprob<-exp(lbf/2)/sum(exp(lbf/2)) # relative model probability
 
  dfall<-data.frame(models,mean.ml,sd.ml,lbf,choice,modelprob) #put it all together
  
  #do permutation tests
  kmeans.p<-permKS(x = likes[[ml_type]], g = likes$model, method="exact.mc")$p.value
  ttest.p<-permTS(x=likes[[ml_type]][which(likes$model==dfall$models[which(dfall$choice==1)])], 
                y=likes[[ml_type]][which(likes$model==dfall$models[which(dfall$choice==2)])], 
                alternative="greater", method="exact.mc")$p.value #non-parametric t-test
  
  output[["table"]]<-dfall
  output["perm_anova.p"]<-kmeans.p
  output["perm_t.test.p"]<-ttest.p
  return(output)
  
}

