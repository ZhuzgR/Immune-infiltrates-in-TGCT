library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/732700f6bf741bdd6ab5a9a518e5d002/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/dae42f7681eca28355c6ba9105a3b5ff',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()

tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "Data press"))
  #logs=c(logs,paste0('geting data:',paste0(paste0(names(opt),'=',opt),collapse = ',')))
  logs=c(logs,paste0('run ConsensusClusterPlus.R-',basename(opt$outfile)))
  
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 1000)
  outFolder=opt$outfile
  exp_path=data$exp_path
  maxK=data$maxK
  reps=data$reps#10
  pItem=data$pItem#0.8
  clusterAlg=data$clusterAlg#'hc'#pam,km,kmdist
  innerLinkage='average'
  finalLinkage='average'
  distance=data$distance#'pearson'# 'spearman','euclidean','binary','maximum','canberra','minkowski"
  
  logs=c(logs,'开始读取表达谱')
  dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F
                        ,na.strings="NA",data.table = F,fill = T)  
  dat=dat[match(unique(dat[,1]),dat[,1]),]
  row.names(dat)=dat[,1]
  exp=dat[,-1]
  logs=c(logs,paste0('表达谱行=',nrow(exp),',列=',ncol(exp)))
  
  triangle=function (m, mode = 1) 
  {
    n = dim(m)[1]
    nm = matrix(0, ncol = n, nrow = n)
    fm = m
    nm[upper.tri(nm)] = m[upper.tri(m)]
    fm = t(nm) + nm
    diag(fm) = diag(m)
    nm = fm
    nm[upper.tri(nm)] = NA
    diag(nm) = NA
    vm = m[lower.tri(nm)]
    if (mode == 1) {
      return(vm)
    }
    else if (mode == 3) {
      return(fm)
    }
    else if (mode == 2) {
      return(nm)
    }
  }
  run_cdf=function (ml, breaks = 100) 
  {
    k = length(ml)
    this_colors = rainbow(k - 1)
    areaK = c()
    all_line=rbind()
    line_gap=c()
    for (i in 1:length(ml)) {
      v = triangle(ml[[i]], mode = 1)
      h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1/breaks))
      h$counts = cumsum(h$counts)/sum(h$counts)
      thisArea = 0
      for (bi in 1:(length(h$breaks) - 1)) {
        thisArea = thisArea + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi])
        bi = bi + 1
      }
      areaK = c(areaK, thisArea)
      h_lower=quantile(h$counts,seq(0,1,0.01))['5%']
      h_upper=quantile(h$counts,seq(0,1,0.01))['95%']
      line_gap=c(line_gap,h_upper-h_lower)
      all_line=rbind(all_line,cbind(K=i+1,X=h$mids,Y=h$counts))
    }
    deltaK = areaK[1]
    for (i in 2:(length(areaK))) {
      deltaK = c(deltaK, (areaK[i] - areaK[i - 1])/areaK[i - 1])
    }
    return(list(CDF=all_line,Delta=cbind(K=1 + (1:length(deltaK)),Y=deltaK,lineGap=line_gap)))
  }
  getTreeNode=function(hc_obj){
    dend=as.dendrogram(hc_obj)
    capture.output(str(dend))->stx2
    gsub('^ \\s+','',stx2)->stx2
    nodeMap=rbind()
    leafs=rbind()
    all_node_obj=rbind()
    for(i in length(stx2):1){
      if(length(grep('--leaf',stx2[i]))==0){
        h=unlist(stringr::str_split(stx2[i],'members at h = '))[2]
        h=gsub(']','',h)
        lefNode=unlist(stringr::str_split(stx2[i],' members at h = '))[1]
        lefNode=unlist(stringr::str_split(lefNode,'branches and '))[2]
        sub1=paste0('C',i+1)
        all_node_obj[which(all_node_obj[,1]==paste0('C',i+1)),4]=paste0('C',i)
        for(j in 1:nrow(all_node_obj)){
          sub_mx=all_node_obj[which(all_node_obj[,1]==paste0('C',i+j+1)),]
          if(sub_mx[4]==''){
            sub2=sub_mx[1]
            break()
          }
        }
        all_node_obj[which(all_node_obj[,1]==sub2),4]=paste0('C',i)
        all_node_obj=rbind(all_node_obj,c(paste0('C',i),sub1,sub2,'',lefNode))
        nodeMap=rbind(nodeMap,c('',paste0('C',i),paste0('C',i+1),sub2,h))
      }else{
        node=gsub(' ','',unlist(stringr::str_split(stx2[i],'--leaf'))[2])
        node=gsub('^ \\s+','',node)
        node=gsub('\\s+ $','',node)
        node=gsub('^"','',node)
        node=gsub('"$','',node)
        leafs=rbind(leafs,c(node,paste0('C',i)))
        all_node_obj=rbind(all_node_obj,c(paste0('C',i),'','','',1))
        nodeMap=rbind(nodeMap,c(node,paste0('C',i),'','',0))
      }
    }
    od=rep(nrow(nodeMap),nrow(nodeMap))
    od[match(leafs[,1],nodeMap[,1])]=nrow(leafs):1
    nodeMap=nodeMap[order(od),]
    return(nodeMap)
  }
  
  library(ConsensusClusterPlus)
  Kvec = 2:maxK
  logs=c(logs,'开始对行进行中心化处理')
  #print(head(exp))
  exp=as.matrix(exp)
  d = sweep(exp,1, apply(exp,1,median,na.rm=T))
  colnames(d)=paste0('C',1:ncol(d))
  
  logs=c(logs,'开始一致性聚类')
  rcc = ConsensusClusterPlus(d,maxK=maxK,reps=reps,pItem=pItem,pFeature=1,title="t",distance=distance,clusterAlg=clusterAlg
                             ,innerLinkage=innerLinkage,seed = 123456)
  logs=c(logs,'聚类完成，开始计算样本聚类一致性')
  resICL = calcICL(rcc,title="t",plot = NULL)
  logs=c(logs,'保存一致性聚类结果')
  save(exp,d,rcc,resICL,file = paste0(outFolder,'/AllCluster.RData'))
  
  #cmt=rcc[[2]]$consensusMatrix
  zip_consensusMatrix=function(cmt){
    m_cmt=data.table::melt(cmt)
    cmt_dt=c()
    for(i in 1:(ncol(cmt)-1)){
      t_m_cmt=m_cmt[m_cmt[,2]==i,]
      cmt_dt=c(cmt_dt,t_m_cmt[match((i+1):ncol(cmt),t_m_cmt[,1]),3])
    }
    #which(cmt_dt-cmt_dt2!=0)
    str_mt=rep('',length(cmt_dt))
    str_mt[which(cmt_dt==0)]='-'
    str_mt[which(cmt_dt==1)]='+'
    str_mt[which(cmt_dt!=0&cmt_dt!=1)]=paste0('.',floor(cmt_dt[which(cmt_dt!=0&cmt_dt!=1)]*1000))
    return(paste0(str_mt,collapse = ''))
  }
  #length(cmt_dt)
  #dim(rcc[[i]]$consensusMatrix)
  all_ml=list()
  clusterConsensus=c()
  all_cluster_tree=rbind()
  all_consensusMatrix=rbind()
  all_consensusMatrix_zip=rbind()
  for(i in Kvec){
    #all_class=rbind(all_class,rcc[[i]]$consensusClass)
    all_ml=c(all_ml,list(rcc[[i]]$ml))
    
    all_consensusMatrix=rbind(all_consensusMatrix,cbind(K=i,rcc[[i]]$consensusMatrix))
    #i=2
    all_consensusMatrix_zip=rbind(all_consensusMatrix_zip,c(i,zip_consensusMatrix(rcc[[i]]$consensusMatrix)))
    
    hc=rcc[[i]]$consensusTree
    treeNode=getTreeNode(hc)
    all_cluster_tree=rbind(all_cluster_tree,cbind(K=i,treeNode))
    clusterConsensus=c(clusterConsensus,mean(resICL$clusterConsensus[resICL$clusterConsensus[,1]==i,3]))
  }
  #dim(rcc[[i]]$consensusMatrix)
  #length(all_ml)
  cdf_value=run_cdf(all_ml)
  delta_vl0=cbind(1,cdf_value$CDF)
  delta_vl=cbind(cdf_value$Delta,clusterConsensus)
  colnames(delta_vl0)=paste0('C',1:4)
  colnames(delta_vl)=paste0('C',1:4)
  delta_v=rbind(delta_vl0,delta_vl)
  clusterConsensus=cbind(0,resICL$clusterConsensus)
  colnames(clusterConsensus)=paste0('C',1:4)
  delta_v=rbind(delta_v,clusterConsensus)
  
  #head(delta_v)
  itemConsensus=resICL$itemConsensus[which(resICL$itemConsensus[,4]!=0),]
  itemConsensus[,3]=match(itemConsensus[,3],colnames(d))
  #exp=d
  itemConsensus=cbind(0,'',itemConsensus)
  colnames(itemConsensus)=paste0('C',1:6)
  colnames(all_cluster_tree)=paste0('C',1:6)
  all_cluster_tree_out=rbind(all_cluster_tree,as.matrix(itemConsensus))
  #head(all_cluster_tree)
  #all_cluster_tree[nrow(all_cluster_tree):(nrow(all_cluster_tree)-10),]
  #length(rcc)
  #rcc[[1]][,1:10]
  cor_clust=unique(as.character(t(rcc[[1]])))
  #setdiff(cor_clust,rcc[[1]][6,])
  consensusClusterAll=apply(rcc[[1]], 2,function(x){
    return(match(x,cor_clust))
  })
  tmp1=rbind(colnames(exp),as.matrix(consensusClusterAll))
  #dim(tmp1)
  #length(cor_clust)
  #unique(rcc[[1]][9,])
  #head(tmp1[,1:10])
  consensusClusterAll_out=cbind(cor_clust,tmp1)
  #consensusClusterAll_out[,1:10]
  #head(all_cluster_tree_out)
  #head(delta_v)
  #head(all_consensusMatrix)
  #head(consensusClusterAll_out)
  logs=c(logs,'输出聚类结果')
  write.table(all_cluster_tree_out,file = paste0(outFolder,'/cluster_tree.txt'),quote = F,row.names = F,col.names = F,sep = '\t')
  write.table(delta_v,file = paste0(outFolder,'/cdf_delta.txt'),quote = F,row.names = F,col.names = F,sep = '\t')
  #dim(all_consensusMatrix)
  #511*9
  #all_consensusMatrix[1,]
  logs=c(logs,'zipConsensusMatrix')
  #all_consensusMatrix_zip[1,]
  #dim(all_consensusMatrix_zip)
  #apply(array, margin, ...)
  write.table(all_consensusMatrix,file = paste0(outFolder,'/consensusMatrix.txt'),quote = F,row.names = F,col.names = F,sep = '\t')
  write.table(all_consensusMatrix_zip,file = paste0(outFolder,'/consensusMatrixZip.txt'),quote = F,row.names = F,col.names = F,sep = '\t')
  
  write.table(consensusClusterAll_out,file = paste0(outFolder,'/cluster_all.txt'),quote = F,row.names = F,col.names = F,sep = '\t')
  logs=c(logs,'run succ')
},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})



