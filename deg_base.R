#setwd('/pub1/data/mg_projects/projects/web_script/R/')
library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_file_catche/runing/tool_runing/c055430ed9b4f109b0ae12591bd48044/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/dae42f7681eca28355c6ba9105a3b5ff',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()
tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "GEO Data press"))
  logs=c(logs,paste0('run deg_base.R-',basename(opt$outfile)))
  #library("rjson")
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 100)
  
  if(unlist(data$method)=='filter'){
    outpath=unlist(data$outpath)
    pvalue=as.numeric(as.character(unlist(data$pvalue)))
    lfc=as.numeric(as.character(unlist(data$lfc)))
    isFDR=as.numeric(as.character(unlist(data$isFDR)))
    
    
    data1<-jsonlite::stream_in(file(paste0(outpath,'/input.json')),pagesize = 100)
    
    compare=unlist(data1$compare)
    groups=unlist(data1$groups)
    types=unlist(data1$types)
    method=data1$method
    cmp=unique(compare)
    file_dat=rbind()
    for(u in cmp){
      #u=cmp[1]
      #paste0(opt$outfile,'/',u,'_exp.mtx')
      deg_path=paste0(outpath,'/',u,'_deg.mtx')
      inds=which(compare==u)
      group=groups[inds]
      tp=types[inds]
      ulab=(group[which(tp==1)][1])
      dlab=(group[which(tp==0)][1])
      
      dat=data.table::fread(deg_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F,na.strings="NA",data.table = F)
      lfcs=as.numeric(as.character(dat[,2]))
      if(isFDR==1){
        if(method=='limma'|method=='limmaVoom'){
          ps=as.numeric(as.character(dat[,6]))  
        }else if(method=='rankTest'|method=='TTest'|method=='deseq2'|method=='edgR'){
          ps=as.numeric(as.character(dat[,4]))  
        }
      }else{
        if(method=='limma'|method=='limmaVoom'){
          ps=as.numeric(as.character(dat[,5]))  
        }else if(method=='rankTest'|method=='TTest'|method=='deseq2'|method=='edgR'){
          ps=as.numeric(as.character(dat[,3]))  
        }
      }
      t.inds=which(ps<pvalue&abs(lfcs)>lfc)
      uCnt=length(which(ps<pvalue&lfcs>lfc))
      dCnt=length(which(ps<pvalue&lfcs<lfc*-1))
      write.table(dat[t.inds,],file = paste0(opt$outfile,'/',u,'_filter.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
      file_dat=rbind(file_dat,c(paste0(ulab,'-vs-',dlab)
                                ,paste0(opt$outfile,'/',u,'_filter.mtx'),uCnt,dCnt))
    }
    write.table(file_dat
                ,file = paste0(opt$outfile,'/stat.mtx'),row.names = F,col.names = F,quote = F,sep = '\t')
  }else{
    
  exp_path=unlist(data$exp_path)
  samples=unlist(data$samples)
  groups=unlist(data$groups)
  types=unlist(data$types)
  compare=unlist(data$compare)
  
  #compare=rep(0,length(types))
  dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F,na.strings="NA",data.table = F)
  #head(dat)
  logs=c(logs,paste0('read data nrow=',nrow(dat),',ncol=',ncol(dat)))
  unm=unique(dat[,1])
  dat=dat[match(unm,dat[,1]),]
  row.names(dat)=dat[,1]
  
  dat=dat[,-1]
  
  limmaDEG=function(exp,group,ulab,dlab){
    library(limma)
    ind1=which(group==ulab)
    ind2=which(group==dlab)
    sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2))) 
    eset=exp[,c(ind1,ind2)]
    fl <- as.factor(sml)
    design <- model.matrix(~fl+0)
    colnames(design) <- levels(fl)
    cont.matrix<-makeContrasts(contrasts='G1-G0',levels=design)
    fit<-lmFit (eset,design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
    eset=eset[match(row.names(tT),row.names(eset)),]
    regulated=ifelse(tT$logFC>0,'Up','Down')
    lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
    all.deg.cnt=cbind()
    for(lfc in lfcs){
      deg1=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.05)]
      deg2=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.01)]
      deg3=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.05)]
      deg4=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.01)]
      all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                      ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                      ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                      ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
    }
    row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
    colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
    return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=all.deg.cnt))
  }
  
  testDEG=function(exp,group,ulab,dlab,rank=1){
    ind1=which(group==ulab)
    ind2=which(group==dlab)
    
    gmn=min(exp,na.rm = T)
    if(gmn<0){
      cv_exp=apply(exp, 2, function(x){
        #print((x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))
        return ((x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))
      })
      fc=apply(cv_exp, 1, function(x){
        return(mean(x[ind1],na.rm = T)/mean(x[ind2],na.rm = T))
      })
    }else{
      fc=apply(exp, 1, function(x){
        return(mean(x[ind1],na.rm = T)/mean(x[ind2],na.rm = T))
      })
    }
    pv=apply(exp, 1, function(x){
      #x=as.numeric(exp[1,])
      x[ind1]->x1
      x[ind2]->x2
      x1=x1[!is.na(x1)]
      x2=x2[!is.na(x2)]
      p=1
      if(length(x1)>2&length(x2)>2){
      if(rank==1){
        wilcox.test(x1,x2)->tt
        p=tt$p.value
      }else{
        t.test(x1,x2)->tt
        p=tt$p.value
      }
      }
      return(p)
    })
    fdr=p.adjust(pv)
    fc=log2(fc)
    eset=exp[,c(ind1,ind2)]
    #sum(fc>0.1)
    regulated=ifelse(fc>0,'Up','Down')
    lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
    all.deg.cnt=cbind()
    for(lfc in lfcs){
      deg1=regulated[which(abs(fc)>lfc&pv<0.05)]
      deg2=regulated[which(abs(fc)>lfc&pv<0.01)]
      deg3=regulated[which(abs(fc)>lfc&fdr<0.05)]
      deg4=regulated[which(abs(fc)>lfc&fdr<0.01)]
      all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                      ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                      ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                      ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
    }
    row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
    colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
    tT=cbind(lfc=fc,pvalue=pv,FDR=fdr)
    row.names(tT)=row.names(eset)
    return(list(Exp=eset[order(pv),],Group=group[c(ind1,ind2)],DEG=tT[order(pv),],Summary=all.deg.cnt))
  }
  
  countDEG=function(lfcValues,pvalue,fdr){
    lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
    regulated=ifelse(lfcValues>0,'Up','Down')
    all.deg.cnt=cbind()
    for(lfc in lfcs){
      deg1=regulated[which(abs(lfcValues)>lfc&pvalue<0.05)]
      deg2=regulated[which(abs(lfcValues)>lfc&pvalue<0.01)]
      deg3=regulated[which(abs(lfcValues)>lfc&fdr<0.05)]
      deg4=regulated[which(abs(lfcValues)>lfc&fdr<0.01)]
      all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                      ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                      ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                      ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
    }
    row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
    colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
    return(all.deg.cnt)
  }
  
  edgRDEG=function(exp,group,ulab,dlab){
    library(edgeR)
    ind1=which(group==ulab)
    ind2=which(group==dlab)
    eset=exp[,c(ind1,ind2)]
    group_list <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))
    y <- DGEList(counts=eset,group=group_list) 
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y) 
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y) 
    tT=data.frame(logFC=et$table$logFC,pvalue=et$table$PValue,FDR=p.adjust(et$table$PValue),logCPM=et$table$logCPM)
    row.names(tT)=row.names(et$table)
    tT=tT[order(tT$FDR),]
    eset <- cpm(y, log=F, prior.count=2)
    eset=eset[match(row.names(tT),row.names(eset)),]
    return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=countDEG(tT$logFC,tT$pvalue,tT$FDR)))
  }
  
  DESeq2DEG=function(exp,group,ulab,dlab){
    library(DESeq2)
    ind1=which(group==ulab)
    ind2=which(group==dlab)
    eset=exp[,c(ind1,ind2)]
    group_list <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))
    
    dds <- DESeqDataSetFromMatrix(countData = eset,
                                  colData = data.frame(row.names=factor(colnames(eset)), group_list=group_list),
                                  design = ~ group_list)
    dds2 <- DESeq(dds)
    resultsNames(dds2)
    res <-  results(dds2, contrast=c("group_list","G1","G0"))
    resOrdered <- res[order(res$padj),]
    tT=data.frame(logFC=res$log2FoldChange,pvalue=res$pvalue,FDR=res$padj)
    row.names(tT)=row.names(res)
    tT=tT[order(tT$FDR),]
    eset=counts(dds2, normalized=TRUE)
    eset=eset[match(row.names(tT),row.names(eset)),]
    return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=countDEG(tT$logFC,tT$pvalue,tT$FDR)))
  }
  
  limmaVoom=function(exp,group,ulab,dlab){
    library(limma)
    ind1=which(group==ulab)
    ind2=which(group==dlab)
    eset=exp[,c(ind1,ind2)]
    sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2))) 
    fl <- as.factor(sml)
    design <- model.matrix(~fl+0)
    colnames(design) <- levels(fl)
    v <- voom(eset, design, plot=F, normalize="quantile")
    fit <- lmFit(v, design)
    cont.matrix<-makeContrasts(contrasts='G1-G0',levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
    tT=tT[order(tT$adj.P.Val),]
    eset=v$E
    eset=eset[match(row.names(tT),row.names(eset)),]
    return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=countDEG(tT$logFC,tT$P.Value,tT$adj.P.Val)))
  }
  logs=c(logs,paste0('run ',data$method,' nrow=',nrow(dat),',ncol=',ncol(dat)))
  if(unlist(data$method)=='limma'){
    if(data$log==1){
      dat=log2(dat)
    }
    cmp=unique(compare)
    all_stat=rbind()
    for(u in cmp){
      #u=cmp[1]
      inds=which(compare==u)
      group=groups[inds]
      exp=dat[,match(samples[inds],colnames(dat))]
      tp=types[inds]
      ulab=(group[which(tp==1)][1])
      dlab=(group[which(tp==0)][1])
      deg=limmaDEG(exp,group,ulab,dlab)
      all_stat=rbind(all_stat,
      cbind(Group=paste0(ulab,'-vs-',dlab),deg$Summary))
      write.table(cbind(Tag=row.names(deg$Exp),deg$Exp)
                  ,file = paste0(opt$outfile,'/',u,'_exp.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
      write.table(cbind(Tag=row.names(deg$DEG),deg$DEG)
                  ,file = paste0(opt$outfile,'/',u,'_deg.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
      logs=c(logs,paste0('output deg:',ulab,'-vs-',dlab))
    }
    write.table(cbind(Tag=row.names(all_stat),all_stat)
                ,file = paste0(opt$outfile,'/stat.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
    logs=c(logs,paste0('succ'))
    
  }else if(unlist(data$method)=='rankTest'|unlist(data$method)=='TTest'){
    
    if(unlist(data$method)=='rankTest'){
      rank=1
    }else{
      rank=0
    }
    if(data$log==1){
      dat=log2(dat)
    }
    cmp=unique(compare)
    all_stat=rbind()
    for(u in cmp){
      #u=cmp[1]
      inds=which(compare==u)
      group=groups[inds]
      exp=dat[,match(samples[inds],colnames(dat))]
      tp=types[inds]
      ulab=(group[which(tp==1)][1])
      dlab=(group[which(tp==0)][1])
      deg=testDEG(exp,group,ulab,dlab,rank)
      all_stat=rbind(all_stat,cbind(Group=paste0(ulab,'-vs-',dlab),deg$Summary))
      write.table(cbind(Tag=row.names(deg$Exp),deg$Exp)
                  ,file = paste0(opt$outfile,'/',u,'_exp.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
      write.table(cbind(Tag=row.names(deg$DEG),deg$DEG)
                  ,file = paste0(opt$outfile,'/',u,'_deg.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
      logs=c(logs,paste0('output deg:',ulab,'-vs-',dlab))
    }
    write.table(cbind(Tag=row.names(all_stat),all_stat)
                ,file = paste0(opt$outfile,'/stat.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
    logs=c(logs,paste0('succ'))
  }else if(unlist(data$method)%in%c('limmaVoom','deseq2','edgR')){
      ctr=as.numeric(unlist(data$cutRowN))/100
      #head(dat)
      dat=floor(dat)
      ctr.ps=apply(dat, 1, function(x){
        return(sum(x>0,na.rm = T))
      })/ncol(dat)
      t.inds=which(ctr.ps>ctr)
      if(length(t.inds)>2){
        dat=dat[t.inds,]
        cmp=unique(compare)
        all_stat=rbind()
        for(u in cmp){
          inds=which(compare==u)
          group=groups[inds]
          exp=dat[,match(samples[inds],colnames(dat))]
          tp=types[inds]
          ulab=(group[which(tp==1)][1])
          dlab=(group[which(tp==0)][1])
          deg=NULL
          if(unlist(data$method)=='limmaVoom'){
            deg=limmaVoom(exp,group,ulab,dlab)
          }else if(unlist(data$method)=='deseq2'){
            deg=DESeq2DEG(exp,group,ulab,dlab)
          }else if(unlist(data$method)=='edgR'){
            deg=edgRDEG(exp,group,ulab,dlab)
          }
          all_stat=rbind(all_stat,cbind(Group=paste0(ulab,'-vs-',dlab),deg$Summary))
          write.table(cbind(Tag=row.names(deg$Exp),deg$Exp)
                      ,file = paste0(opt$outfile,'/',u,'_exp.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
          write.table(cbind(Tag=row.names(deg$DEG),deg$DEG)
                      ,file = paste0(opt$outfile,'/',u,'_deg.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
          logs=c(logs,paste0('output deg:',ulab,'-vs-',dlab))
        }
        write.table(cbind(Tag=row.names(all_stat),all_stat)
                    ,file = paste0(opt$outfile,'/stat.mtx'),row.names = F,col.names = T,quote = F,sep = '\t')
        logs=c(logs,paste0('succ'))
      }else{
        logs=c(logs,paste0('error:','filter count last<3'))
      }
  }
  }
},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})


