migrate_harvester<-function(wd,n=3,models){ #wd<-"/Users/eric/Datasets/simulations/ibdsim/migrate_pareto2_mtdna_results"; n=3; models<-c("panmixia","5island","5stepping.stone","5stepping.stone1","10island","10stepping.stone","10stepping.stone1")
  
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
