#setwd('/pub1/data/mg_projects/projects/web_script/R/')
library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/8c1de46ae43736a08386ee1df0e7f42d/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/dae42f7681eca28355c6ba9105a3b5ff',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()
#stx=rbind()
tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "GEO Data press"))
  #logs=c(logs,paste0('geting data:',paste0(paste0(names(opt),'=',opt),collapse = ',')))
  logs=c(logs,paste0('run lasso_cox.R-',basename(opt$outfile)))
  #library("rjson")
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 100)
  exp_path=unlist(data$exp_path)
  samples=unlist(data$samples)
  times=unlist(data$times)
  events=unlist(data$events)
  lambda.min=unlist(data$lambda)
  nfolds=unlist(data$nfolds)
  oldFolder=unlist(data$oldFolder)
  
  set.seed(123456789)
  library("glmnet") 
  library('survival')
  
  if(oldFolder!=''){
    file.copy(opt$infile,paste0(oldFolder,'/clini.json'),overwrite = T)
    load(file =  paste0(oldFolder,'/fit1_cv.RData'))
    logs=c(logs,'load fit1_cv.RData')
    load(file =  paste0(oldFolder,'/fit.RData'))
    logs=c(logs,'load fit.RData')
    load(file =  paste0(oldFolder,'/mgdata.RData'))
    logs=c(logs,'load mgdata.RData')
    logs=c(logs,paste0('reset lambda=',lambda.min))
  }else{
    oldFolder=opt$outfile
    logs=c(logs,paste0('run init',oldFolder))
    if(!file.exists(paste0(oldFolder,'/clini.json'))){
      write.table('',file = paste0(oldFolder,'/clini.json'),quote = F,row.names = F,sep = '\t')
    }
    logs=c(logs,paste0('run readdata'))
  dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F
                        ,na.strings="NA",data.table = F)
  rNames=unique(dat[,1])
  rNames=rNames[which(rNames!='')]
  dat=dat[match(rNames,dat[,1]),]
  row.names(dat)=dat[,1]
  dat=dat[,-1]
  logs=c(logs,paste0('read data,row=',nrow(dat),',col=',ncol(dat)))
  uSamples=unique(samples)
  t_inds=match(uSamples,samples)
  
  dat=dat[,match(uSamples,colnames(dat))]
  times=times[t_inds]
  events=events[t_inds]
  logs=c(logs,paste0('clean data,row=',nrow(dat),',col=',ncol(dat)))
  #sum(apply(dat1, 2,sd)==0)
  dat1=t(dat)
  colnames(dat1)=paste0('C',1:ncol(dat1))
  time=as.numeric(times)
  event=as.numeric(events)
  #t.ind=which(time>0)
  y=Surv(time,event)
  #nfolds=1
  #nfolds=3
  #dim(dat1)
  #sum(time>0)
  logs=c(logs,'starting fit cv')
  fit1_cv = cv.glmnet(as.matrix(dat1), y, family = "cox", nfolds=nfolds)
  logs=c(logs,'end fit cv,start fit')
  
  fit<-glmnet(dat1, y, family = "cox")
  logs=c(logs,'end fit')
  
  logs=c(logs,'output PLD data')
  tgc=data.frame(lambda=fit1_cv$lambda,cvm=fit1_cv$cvm,cvup=fit1_cv$cvup,cvlo=fit1_cv$cvlo,cvsd=fit1_cv$cvsd)
  
  write.table(tgc,file = paste0(opt$outfile,'/batchPLD.txt')
              ,row.names = F,col.names = T,quote = F,sep = '\t')  
  
  #cv_fit$lambda.min
  #cv_fit$lambda.1se
  
  #ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B')
  logs=c(logs,'output Coef data')
  fit.coef=fit$beta[(apply(fit$beta,1,function(x){
    return(sum(x!=0))
  })>0),]
  fit.coef=as.matrix(fit.coef)
  row.names(fit.coef)=row.names(dat)[match(row.names(fit.coef),colnames(dat1))]
  mtx1=rbind(fit$lambda,fit.coef)
  write.table(cbind(Tag=row.names(mtx1),mtx1),file = paste0(opt$outfile,'/batchCoef.txt')
              ,row.names = F,col.names = F,quote = F,sep = '\t')  
  logs=c(logs,'save RData')
  #head(mtx1)
  #dim(mtx1)
  save(fit1_cv,file =  paste0(opt$outfile,'/fit1_cv.RData'))
  save(fit,file =  paste0(opt$outfile,'/fit.RData'))
  mgdata<-cbind(times,events,t(dat))
  save(mgdata,file =  paste0(opt$outfile,'/mgdata.RData'))
  }
  #fit.coef[1:10,1:10]
  logs=c(logs,paste0('lambda.min=',fit1_cv$lambda.min))
  logs=c(logs,paste0('lambda.1se=',fit1_cv$lambda.1se))
  
  logs=c(logs,'set lambda.min')
  if(lambda.min==1){
    lambda=fit1_cv$lambda.min
  }else{
    lambda=lambda.min
  }
  coefficients<-coef(fit,s=lambda)
  Active.Index<-which(coefficients[,1]!=0)
  if(length(Active.Index)==0){
    logs=c(logs,paste0('not found active genes'))
    #write.table(cbind('',''),file = paste0(oldFolder,'/coef.txt')
    #            ,row.names = F,col.names = T,quote = F,sep = '\t')  
  }else{
    genes=row.names(coefficients)[Active.Index]
    Active.coefficients<-coefficients[Active.Index]  
    lst.genes.exp=mgdata[,as.numeric(gsub('C','',genes))+2]
    logs=c(logs,paste0('found genes:',paste0(colnames(lst.genes.exp),collapse = '#,#')))
    logs=c(logs,paste0('gene coefficients:',paste0(Active.coefficients,collapse = '#,#')))
    
    coef_dat=cbind(Tag=colnames(lst.genes.exp),Coef=Active.coefficients)
    write.table(coef_dat,file = paste0(oldFolder,'/coef.txt')
                ,row.names = F,col.names = T,quote = F,sep = '\t')  
    #dim(lst.genes.exp)
    #head(lst.genes.exp)
    riskscore=lst.genes.exp%*%Active.coefficients
    risk_OS=cbind(Time=mgdata[,1],Status=mgdata[,2],RiskScore=riskscore[,1])
    logs=c(logs,paste0('output riskcore'))
    write.table(cbind(Tag=row.names(risk_OS),risk_OS),file = paste0(oldFolder,'/riskscore.txt')
                ,row.names = F,col.names = T,quote = F,sep = '\t')  
    logs=c(logs,paste0('output filter gene exp'))
    write.table(cbind(Tag=row.names(lst.genes.exp),lst.genes.exp),file = paste0(oldFolder,'/filterGeneExp.txt')
                ,row.names = F,col.names = T,quote = F,sep = '\t')
    logs=c(logs,paste0('all runed'))
    #rel=cbind(Time=mgdata[,1],Status=mgdata[,2],RiskScore=riskscore[,1],lst.genes.exp)
    #rel=rbind(c(NA,NA,NA,Active.coefficients),rel)
    #head(rel)
    #row.names(dat)[match(genes,colnames(dat1))]
  }
},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  if(oldFolder!=opt$outfile) write.table(logs,file = paste0(oldFolder,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})


