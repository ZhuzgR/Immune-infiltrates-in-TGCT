#setwd('/pub1/data/mg_projects/projects/web_script/R/')
library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/36b7244f705a53bef20bca4092c674c8/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/dae42f7681eca28355c6ba9105a3b5ff',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()
#stx=rbind()
tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "GEO Data press"))
  #logs=c(logs,paste0('geting data:',paste0(paste0(names(opt),'=',opt),collapse = ',')))
  logs=c(logs,paste0('run batch_cox.R-',basename(opt$outfile)))
  #library("rjson")
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 100)
  exp_path=unlist(data$exp_path)
  samples=unlist(data$samples)
  times=unlist(data$times)
  events=unlist(data$events)
  
  dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F
                        ,na.strings="NA",data.table = F)
  row.names(dat)=dat[,1]
  dat=dat[,-1]
  logs=c(logs,paste0('read data,row=',nrow(dat),',col=',ncol(dat)))
  uSamples=unique(samples)
  t_inds=match(uSamples,samples)
  
  dat=dat[,match(uSamples,colnames(dat))]
  times=times[t_inds]
  events=events[t_inds]
  logs=c(logs,paste0('clean data,row=',nrow(dat),',col=',ncol(dat)))
  library(survival)
  t(apply(dat, 1, function(x){
    gp=as.numeric(x)
    dt=c(NA,NA,NA,NA,NA,NA)
    tryCatch({
      nDat1=data.frame(time=times,status=events,group=gp)
      cx=coxph(Surv(time, status) ~ group, data=nDat1)
      cxs=summary(cx)
      hr=cxs$conf.int[1,c(1,3,4)]
      dt=c(hr,cxs$logtest[3],cxs$sctest[3],cxs$waldtest[3])
    },error = function(e) {
    }, finally = {
      return(dt)
    })
  }))->all_cx
  logs=c(logs,paste0('runed batch cox',nrow(all_cx),',col=',ncol(all_cx)))
  colnames(all_cx)=c('HR','Lower','Upper','Likelihood','logrank','Wald')
  all_cx=all_cx[order(as.numeric(all_cx[,5])),]
  write.table(cbind(Tag=row.names(all_cx),all_cx),file = paste0(opt$outfile,'/batchCox.txt')
              ,row.names = F,col.names = T,quote = F,sep = '\t')  
  logs=c(logs,paste0('outputed batch cox'))
  
  dat1=cbind(Time=times,Status=events,t(dat[match(row.names(all_cx),row.names(dat)),]))
  row.names(dat1)=colnames(dat)
  logs=c(logs,paste0('merged exp'))
  write.table(cbind(Tag=colnames(dat1),t(dat1)),file = paste0(opt$outfile,'/batchCoxExp.txt')
              ,row.names = F,col.names = T,quote = F,sep = '\t')  
  logs=c(logs,paste0('outputed exp'))
},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})


