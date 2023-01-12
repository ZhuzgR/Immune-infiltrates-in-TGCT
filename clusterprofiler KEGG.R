library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/48d62680c879ba86522cfdfa2c641491/input.json',
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
  logs=c(logs,paste0('run clusterprofiler.R-',basename(opt$outfile)))
  #library("rjson")
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 1000)
  
  dbName=unlist(data$dbMode)
  genes=unlist(data$genes)
  genefc=unlist(data$genefc)
  gmtNames=unlist(data$gmtNames)
  gmtGenes=unlist(data$gmtGenes)
  
  pAdjustMethod=unlist(data$pAdjustMethod)#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
  outFolder=opt$outfile
  
  dbPath='/pub1/data/mg_projects/projects/web_script/source'
  logs=c(logs,paste0('input gene length=',length(genes)))

  library(clusterProfiler)
  gmt2tab<-function(dbName){
    gmt_path=paste0(dbPath,'/',dbName)
    prex='';
    if(dbName=='c2.cp.biocarta.v7.4.symbols.gmt') prex='BIOCARTA_';
    if(dbName=='c2.cp.kegg.v7.4.symbols.gmt') prex='KEGG_';
    if(dbName=='c2.cp.pid.v7.4.symbols.gmt') prex='PID_';
    if(dbName=='c2.cp.reactome.v7.4.symbols.gmt') prex='REACTOME_';
    if(dbName=='c2.cp.wikipathways.v7.4.symbols.gmt') prex='WP_';
    if(dbName=='c4.cm.v7.4.symbols.gmt') prex='MODULE_';
    if(dbName=='c5.go.bp.v7.4.symbols.gmt') prex='GOBP_';
    if(dbName=='c5.go.cc.v7.4.symbols.gmt') prex='GOCC_';
    if(dbName=='c5.go.mf.v7.4.symbols.gmt') prex='GOMF_';
    if(dbName=='c5.hpo.v7.4.symbols.gmt') prex='HP_';
    if(dbName=='h.all.v7.4.symbols.gmt') prex='HALLMARK_';
    if(prex!=''){
      cgeneset=GSEABase::getGmt(gmt_path)
      all.tb=rbind()
      for(i in 1:length(cgeneset)){
        all.tb=rbind(all.tb,cbind(cgeneset[[i]]@setName,cgeneset[[i]]@geneIds))
      }
      all.tb[,1]=gsub(paste0('^',prex),'',all.tb[,1])
      TERM2NAME=data.frame(term =paste0('ID',1:length(unique(all.tb[,1]))),name=unique(all.tb[,1]))
      TERM2GENE=data.frame(term =as.character(TERM2NAME[match(all.tb[,1],as.character(TERM2NAME[,2])),1]),gene=all.tb[,2])
      return(list(TERM2NAME=TERM2NAME,TERM2GENE=TERM2GENE,subMap=NULL))
    }else{
      if(dbName=='KEGG'){
        ft = tidyfst::import_fst(paste0(dbPath,'/KEGG_map.fst'),as.data.table = F)
        TERM2NAME=unique(data.frame(term =ft[,2],name=ft[,3]))
        TERM2GENE=data.frame(term =ft[,2],gene=ft[,1])
        subMap=unique(data.frame(term =ft[,2],name=ft[,3],ft[,4:5]))
        return(list(TERM2NAME=TERM2NAME,TERM2GENE=TERM2GENE,subMap=subMap))
      }else if(dbName=='GO'){
        ft = tidyfst::import_fst(paste0(dbPath,'/GO_map.fst'),as.data.table = F)
        TERM2NAME=unique(data.frame(term =ft[,2],name=ft[,3]))
        TERM2GENE=data.frame(term =ft[,2],gene=ft[,1])
        subMap=unique(data.frame(term =ft[,2],name=ft[,3],ft[,4]))
        return(list(TERM2NAME=TERM2NAME,TERM2GENE=TERM2GENE,subMap=subMap))
      }else if(dbName=='Custom'){
        TERM2NAME=data.frame(term =paste0('ID',1:length(unique(gmtNames))),name=unique(gmtNames))
        TERM2GENE=data.frame(term =as.character(TERM2NAME[match(gmtNames,as.character(TERM2NAME[,2])),1]),gene=gmtGenes)
        return(list(TERM2NAME=TERM2NAME,TERM2GENE=TERM2GENE,subMap=NULL))
      }
    }
    return(NULL)
  }
  logs=c(logs,paste0('geting ',dbName,' DB map'))
  
  #dbName='c2.cp.kegg.v7.4.symbols.gmt'
  #  dbTab=gmt2tab('KEGG')
  
  #genes=as.character(gene2kegg[,1][1:100])
  #pAdjustMethod=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")[1]
  
  logs=c(logs,paste0('geting ',dbName,' DB map'))
  dbTab=gmt2tab(dbName)
  cmp=intersect(dbTab$TERM2GENE[,2],genes)
  logs=c(logs,paste0('intersect DB Gene length=',length(cmp)))
  if(length(cmp)>3){
    genes=unique(genes)
    enrich.tab=rbind()
    if(dbName=='GO'){
      #dbTab=gmt2tab("GO")
      enrich1=clusterProfiler::enricher(genes, pvalueCutoff = 1
                                       , pAdjustMethod = pAdjustMethod
                                       ,minGSSize = 3, maxGSSize = 5000
                                       , qvalueCutoff = 1
                                       ,TERM2GENE=dbTab$TERM2GENE[dbTab$TERM2GENE[,1]%in%dbTab$subMap[which(dbTab$subMap[,3]=='BP'),1],],
                                       TERM2NAME = dbTab$TERM2NAME)
      logs=c(logs,paste0('succ enriched GO_BP length=',nrow(enrich1)))
      if(!is.null(enrich1)){
        enrich.tab=rbind(enrich.tab,cbind(enrich1@result,ONT='BP'))
      }
      enrich2=clusterProfiler::enricher(genes, pvalueCutoff = 1
                                        , pAdjustMethod = pAdjustMethod
                                        ,minGSSize = 3, maxGSSize = 5000
                                        , qvalueCutoff = 1
                                        ,TERM2GENE=dbTab$TERM2GENE[dbTab$TERM2GENE[,1]%in%dbTab$subMap[which(dbTab$subMap[,3]=='CC'),1],],
                                        TERM2NAME = dbTab$TERM2NAME)
      logs=c(logs,paste0('succ enriched GO_CC length=',nrow(enrich2)))
      if(!is.null(enrich2)){
        enrich.tab=rbind(enrich.tab,cbind(enrich2@result,ONT='CC'))
      }
      enrich3=clusterProfiler::enricher(genes, pvalueCutoff = 1
                                        , pAdjustMethod = pAdjustMethod
                                        ,minGSSize = 3, maxGSSize = 5000
                                        , qvalueCutoff = 1
                                        ,TERM2GENE=dbTab$TERM2GENE[dbTab$TERM2GENE[,1]%in%dbTab$subMap[which(dbTab$subMap[,3]=='MF'),1],],
                                        TERM2NAME = dbTab$TERM2NAME)
      logs=c(logs,paste0('succ enriched GO_MF length=',nrow(enrich3)))
      if(!is.null(enrich3)){
        enrich.tab=rbind(enrich.tab,cbind(enrich3@result,ONT='MF'))
      }
    }else{
      enrich=clusterProfiler::enricher(genes, pvalueCutoff = 1
                                            , pAdjustMethod = pAdjustMethod
                                            ,minGSSize = 3, maxGSSize = 5000
                                            , qvalueCutoff = 1
                                            ,TERM2GENE=dbTab$TERM2GENE,
                                            TERM2NAME = dbTab$TERM2NAME)
      logs=c(logs,paste0('succ enriched length=',nrow(enrich)))
      if(!is.null(enrich)){
        enrich.tab=enrich@result
        if(dbName!='KEGG'){
          enrich.tab=enrich.tab[,-1]
        }else{
          enrich.tab=cbind(enrich.tab,dbTab$subMap[match(enrich.tab[,1],dbTab$subMap[,1]),3:4])
        }
      }
    }
    
    if(!is.null(enrich.tab)){
      logs=c(logs,'outputing enrich!')
      enrich.tab=enrich.tab[order(enrich.tab$p.adjust),]
      write.table(enrich.tab,file = paste0(opt$outfile,'/enrichResult.txt'),quote = F,row.names = F,col.names = T,sep = '\t')
      logs=c(logs,'output succed enrich!')
      #if(dbName=='GO'){
      #  logs=c(logs,'GOSemSim starting!') 
      #}
    }else{
      logs=c(logs,paste0('intersect DB Gene length<3'))
    }
    #colnames(enrich.tab)
  }else{
    #go1 <- c("GO:0004022", "GO:0004024", "GO:0004023","GO:0009055", "GO:0020037")
    #d <- new("GOSemSimDATA",ont = 'CC', 
    #         metadata = metadata(org.Hs.eg.db))
    #bps=c('GO:0001775','GO:0002252')
    #orangeCorr=GOSemSim::termSim(bps, bps, d, method = c("Wang", "Resnik", "Rel", "Jiang", "Lin")[1])
    #orangeClust <- hclust(dist(orangeCorr, method="euclidean"), method="complete")
    #dynamicCut <- dynamicTreeCut::cutreeDynamic(orangeClust, minClusterSize=1, method="hybrid"
    #                            , distM=as.matrix(dist(orangeCorr, method="euclidean"))
    #                            , deepSplit=4, maxCoreScatter=NULL, minGap=NULL
    #                            , maxAbsCoreScatter=NULL, minAbsGap=NULL)
    logs=c(logs,paste0('intersect DB Gene length<3'))
}

},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})

#library(ggstatsplot)
