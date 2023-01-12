library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/689a26b32a93e8cc5345b11f59bfdef9/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/b5d509f2ddaf74dd2ca07303a86e46a3',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()
#stx=rbind()
tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "Data press"))
  #logs=c(logs,paste0('geting data:',paste0(paste0(names(opt),'=',opt),collapse = ',')))
  logs=c(logs,paste0('run gsea.R-',basename(opt$outfile)))
  #library("rjson")
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 1000)
  
  exp_path=unlist(data$exp_path)
  
  #exp_path='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/GSE73452.txt'
  #dbPath='/pub1/data/mg_projects/projects/web_script/source'
  #dbName='c2.cp.kegg.v7.4.symbols.gmt'
  #####Example 1######
  #exp_path='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/GSE73452.txt'
  #dbPath='/pub1/data/mg_projects/projects/web_script/source'
  #dbName='c2.cp.kegg.v7.4.symbols.gmt'
  #outFolder='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/GSEA'
  ###########
  samples=unlist(data$samples)
  groups=unlist(data$groups)
  dbName=unlist(data$dbMode)
  dbPath=unlist(data$dbPath)
  gmtNames=unlist(data$gmtNames)
  gmtGenes=unlist(data$gmtGenes)
  method=unlist(data$method)#
  gene=unlist(data$gene)
  geneSplit=unlist(data$geneSplit)
  outFolder=opt$outfile
  command=NULL
  
  #paste0(round(groups,2),collapse = "','")
  #paste0(samples,collapse = "','")
  #samples=row.names(dat)
  #groups=apply(dat, 1,function(x){
  #  return (median(x[which(egfr.exp>4.970198)])/median(x[which(egfr.exp<=4.970198)]))
  #})
  #method='rank'
  #samples=samples[1:1000]
  #groups=groups[1:1000]
  #samples1=samples[order(groups)]
  #gmtGenes=c(samples1[sample(1:100,30)],
  #           samples1[sample(1:200,40)],
  #           samples1[sample(900:1000,50)])
  #gmtNames=c(rep('Signature 1',30),rep('Signature 2',40),rep('Signature 3',50))
  #dbName='Custom'
  #paste0(gmtNames,collapse = "','")
  
  nxt=TRUE
  if(method!='rank'){
    dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F
                          ,na.strings="NA",data.table = F,skip = 0,fill=T)
    #head(dat)
    #dim(dat)
    uName=unique(dat[,1])
    dat=dat[match(uName,dat[,1]),]
    row.names(dat)=dat[,1]
    dat=dat[,-1]  
    all.genes=row.names(dat)
    if(method=='gene'){
      geneSplit=ceiling(as.numeric(geneSplit))
      g_ind=which(row.names(dat)==gene)
      if(length(g_ind)==1){
        g_exp=as.numeric(dat[g_ind,])
        g_q=quantile(g_exp,seq(0,1,0.01))
        groups=ifelse(g_exp>g_q[geneSplit+1],'H','L')
        samples=colnames(dat)
      }else{
        nxt=FALSE
        logs=c(logs,paste0('#Stop:not found gene ',gene,'!'))
      }
    }
    
    if(nxt){
      smp.cm=intersect(samples,colnames(dat))
      if(length(smp.cm)<3){
        nxt=FALSE
        logs=c(logs,'#Stop:intersect Sample <3!')
      }else{
        t_ind=which(samples%in%smp.cm)
        samples=samples[t_ind]
        groups=groups[t_ind]
        if(length(table(groups))!=2){
          nxt=FALSE
          logs=c(logs,'#Stop:intersect Sample Group !=2 !')
        }else{
          if(sum(table(groups)>2)!=2){
            nxt=FALSE
            logs=c(logs,'#Stop:intersect Sample Group <2 !')
          }
        }
      }
    }
    
  }else{
    logs=c(logs,'output rank file!')
    dat=cbind(Gene=samples,Value=groups)
    write.table(dat,file = paste0(outFolder,'/exp_data.txt')
                ,quote = F,row.names = F,sep = '\t')
    logs=c(logs,'outputed rank file!')
    all.genes=unique(samples)
  }
  
  #egfr.exp=as.numeric(dat[which(row.names(dat)=='EGFR') ,])
  #groups=ifelse(egfr.exp>median(egfr.exp),'H','L')
  #samples=colnames(dat)
  
  jar_path=paste0(dbPath,'/MG_GSEA.jar')
  gmt_path=paste0(dbPath,'/',dbName)
  
  #outFolder='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/GSEA'
  
  if(nxt){
    if(dbName=='Custom'){
      logs=c(logs,'outputing custom gmt')
      all.list=list()
      gmt_path=paste0(outFolder,'/',dbName,'.gmt')
      for(gn in unique(gmtNames)){
        t_inds=which(gmtNames==gn)
        gens=intersect(unique(gmtGenes[t_inds]),all.genes)
        if(length(gens)>5){
          gs=GSEABase::GeneSet(setName=gn, setIdentifier=paste0("101")
                               ,geneIds=gens
                               ,GSEABase::SymbolIdentifier()) 
          all.list=c(all.list,list(gs))
        }else{
          logs=c(logs,paste0('remove ',gn,' in gmt,gene<5'))
        }
      }
      if(length(all.list)>0){
        gsc <- GSEABase::GeneSetCollection(all.list)
        GSEABase::toGmt(gsc, gmt_path)
        logs=c(logs,'outputed Custom.gmt')
      }else{
        nxt=FALSE
        logs=c(logs,'#Gene not found!')
        logs=c(logs,'#Stop:remove all Custom GMT!')
      }
    }else{
      gmt=clusterProfiler::read.gmt(gmt_path)
      gmt=gmt[gmt[,2]%in%all.genes,]
      if(sum(table(gmt[,1])>5)==0){
        nxt=FALSE
        logs=c(logs,'#Gene not found!')
        logs=c(logs,'#Stop:remove all GMT!')
      }
    }
  }
  
  if(nxt){
  if(method!='rank'){
    logs=c(logs,'maping group!')
    exp_dat=dat[,match(samples,colnames(dat))]
    sample_dat=cbind(Sample=samples,Group=groups)
    logs=c(logs,paste0('run compare[',unique(groups)[1],'#-vs-#',unique(groups)[2],']'));
    write.table(sample_dat,file = paste0(outFolder,'/exp_data_sample.txt'),quote = F,row.names = F,sep = '\t')
    write.table(cbind(Tag=row.names(exp_dat),exp_dat),file = paste0(outFolder,'/exp_data.txt')
                ,quote = F,row.names = F,sep = '\t')
    logs=c(logs,'outputed group!')
    
    command=paste0('java -jar ',jar_path,' exp_group '
                   ,paste0(outFolder,'/exp_data.txt'),' ',paste0(outFolder,'/exp_data_sample.txt'),' '
                   ,outFolder,' ',gmt_path,' ','false',' ',3,' ',5,' ',5000) 
  }else{
    command=paste0('java -jar ',jar_path,' rank '
                   ,paste0(outFolder,'/exp_data.txt'),' ',1,' ',outFolder,' ',gmt_path,' ','false',' ','3',' ',5,' ',5000) 
  }
  logs=c(logs,'runing GSEA!')
  glogs=system(command, intern = T, 
              ignore.stdout = FALSE, ignore.stderr = FALSE, 
              wait = TRUE, input = NULL, show.output.on.console = TRUE, 
              minimized = FALSE, invisible = TRUE)
  logs=c(logs,'runed GSEA!')
  logs=c(logs,glogs)
  getFile=function(fod){
    fl=c()
    for(d in dir(fod)){
      if(fs::is_dir(paste0(fod,'/',d))){
        fl1=getFile(paste0(fod,'/',d))
        fl=c(fl,paste0(d,'/',fl1))
      }else{
        fl=c(fl,d)
      }
    }
    return(fl)
  }
  logs=c(logs,'select file!')
  for(d in dir(outFolder)){
    #print(d)
    #d='7364ec348ae3428c9d9058df1a224ba9.cls'
    #logs=c(logs,d)
    if(fs::is_dir(paste0(outFolder,'/',d))){
      logs=c(logs,paste0('#GSEAed# ',d))
      logs=c(logs,paste0('GF:',getFile(paste0(outFolder,'/',d))))
    }
  }
  logs=c(logs,'selected file!')
  
}
},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})

#library(ggstatsplot)
