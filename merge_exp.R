COLORS = c("red","green3","blue","cyan","magenta","yellow","gray","black",
           "orange","darkred","green","darkblue","darkcyan","darkmagenta",
           "darkorchid1","darkgoldenrod3","aquamarine","antiquewhite",
           "darkolivegreen3");

# create color map
#-------------------------------------------------------------------------------
makeColorMap = function(eset, label)
{
  colMap = list();
  vec = unique(as.vector(pData(eset)[,label]));
  for(i in 1:length(vec)) { colMap[[ vec[i] ]] = COLORS[i]; }
  return(colMap);
}  

# create color vector for plots
#-------------------------------------------------------------------------------
makeColorVec = function(eset, label, colMap)
{
  labels = as.vector(pData(eset)[,label]);
  return(as.vector(unlist(sapply(labels, function(x) { colMap[x]; }))));
}

# apply geneBMCESet on list of eSets
#-------------------------------------------------------------------------------
mergeBMC = function(esets)
{
  esets = lapply(esets, geneBMCESet);
  
  eset = mergeNONE(esets);
  return(eset);
}

#-------------------------------------------------------------------------------
geneBMCESet = function(eset)
{
  # SM: Changed / in - since exprs already on log scale
  # This corresponds to shifting mean to 0
  exprs(eset) = exprs(eset) - rowMeans(exprs(eset));
  return(eset);
}
# Apply Combat on list of eSets
#
# Code taken from: http://jlab.byu.edu/ComBat/Download_files/ComBat.R
# (slightly adapted to work with ExpressionSets)
#-------------------------------------------------------------------------------
mergeCOMBAT = function(esets)
{
  raw_merged = mergeNONE(esets)
  batchInfo = NULL;
  for(i in 1:length(esets))
  {
    batchInfo = c(batchInfo, rep(i,ncol(esets[[i]])));
  }
  
  saminfo = cbind(rownames(pData(raw_merged)),
                  rownames(pData(raw_merged)),
                  batchInfo)
  colnames(saminfo) = c("Array name", "Sample name", "Batch")
  
  dat = exprs(raw_merged)
  
  design <- design.mat(saminfo)
  
  batches <- list.batch(saminfo)
  n.batch <- length(batches)
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  
  ##Standardize Data across genes
  B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
  
  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
  if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
  
  ##Get regression batch effect parameters
  batch.design <- design[,1:n.batch]
  gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
  delta.hat <- NULL
  for (i in batches)
  {
    delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
  }
  
  ##Find Priors
  gamma.bar <- apply(gamma.hat, 1, function(x){return(mean(x,na.rm=T))})
  t2 <- apply(gamma.hat, 1,  function(x){return(var(x,na.rm=T))})
  a.prior <- apply(delta.hat, 1, function(x)
  {
    m=mean(x,na.rm=T); s2=var(x,na.rm=T); return((2*s2+m^2)/s2)
  })
  b.prior <- apply(delta.hat, 1, function(x){
    m=mean(x,na.rm=T); s2=var(x,na.rm=T); return((m*s2+m^3)/s2)
  })
  
  ##Find EB batch adjustments
  gamma.star <- delta.star <- NULL
  for (i in 1:n.batch)
  {
    temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
    gamma.star <- rbind(gamma.star,temp[1,])
    delta.star <- rbind(delta.star,temp[2,])
  }
  
  ### Normalize the Data ###
  bayesdata <- s.data
  j <- 1
  for (i in batches)
  {
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
  
  eset=raw_merged
  exprs(eset)=bayesdata
  return(eset)	
}

#Helper functions
#-------------------------------------------------------------------------------
build.design <- function(vec, des=NULL, start=2)
{
  tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
  for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
  cbind(des,tmp)
}

#-------------------------------------------------------------------------------
design.mat <- function(saminfo)
{
  tmp <- which(colnames(saminfo) == 'Batch')
  tmp1 <- as.factor(saminfo[,tmp])
  msg("  => Found",nlevels(tmp1),'batches');
  design <- build.design(tmp1,start=1)
  ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
  msg("  => Found",ncov,'covariate(s)');
  if(ncov>0){
    for (j in 1:ncov){
      tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
      design <- build.design(tmp1,des=design)
    }
  }
  design
}

#-------------------------------------------------------------------------------
list.batch <- function(saminfo)
{
  tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
  batches <- NULL
  for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
  batches
}

#-------------------------------------------------------------------------------
aprior <- function(gamma.hat)
{
  m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2
}

#-------------------------------------------------------------------------------
it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001)
{
  n <- apply(!is.na(sdat),1,sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
    d.new <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old,na.rm = T)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

#-------------------------------------------------------------------------------
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}
# apply geneNormESet on list of eSets
#-------------------------------------------------------------------------------
mergeGENENORM = function(esets)
{
  esets = lapply(esets, geneNormESet);
  eset = mergeNONE(esets);  
}

#-------------------------------------------------------------------------------
geneNormESet = function(eset)
{
  #eset=eset1
  exp = apply(exprs(eset), 1, function(x){x-mean(x, na.rm=TRUE)});
  exp = apply(exp, 1, function(x){x/sd(x, na.rm=TRUE)});
  exprs(eset)=exp
  #eset1=ExpressionSet(assayData =as.matrix(cbind(d1, d2)))
  return(eset);
}
# Concatenate different esets. No data transformation is applied in this
# function.
#-------------------------------------------------------------------------------
mergeNONE = function(esets)
{
  #esets=set_all
  eset1 = esets[[1]];
  annot1 = annotation(eset1)
  
  for(i in 2:length(esets))
  {
    eset2 = esets[[i]];
    d1 = exprs(eset1);
    d2 = exprs(eset2);
    
    #-----------------------------------------------------
    # Rebuild fData
    #-----------------------------------------------------
    cg = sort(intersect(rownames(d1), rownames(d2))); 
    # If too few overlapping genes
    if(length(cg) < (min(dim(d1)[1],dim(d2)[1])/100))
    {
      msg(" ! WARNING ! Number of common genes < 1%")
    }
    fData = fData(eset1)[cg,]
    
    #-----------------------------------------------------
    # Rebuild pData
    #-----------------------------------------------------
    
    # Create new pData: first common annotations, then unique ones
    p1 = pData(eset1)
    p2 = pData(eset2)
    cp = sort(intersect(colnames(p1),colnames(p2)));
    tp = sort(unique(union(colnames(p1),colnames(p2))))
    
    sp1 = setdiff(colnames(p1),cp)
    sp2 = setdiff(colnames(p2),cp)
    
    pheno = matrix(NA,ncol=length(tp),nrow=nrow(p1)+nrow(p2))
    rownames(pheno) = c(rownames(p1),rownames(p2))
    colnames(pheno) = tp;
    
    if(length(cp)!=0)
    {
      pheno[1:nrow(p1),cp] = as.matrix(p1[,cp])
      pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)),cp] = as.matrix(p2[,cp])
    }
    
    if(length(sp1)!=0)
    {
      pheno[1:nrow(p1),sp1] = as.matrix(p1[,sp1])
    }
    if(length(sp2)!=0)
    {
      pheno[(nrow(p1)+1):(nrow(p1)+nrow(p2)), sp2] = as.matrix(p2[,sp2])
    }
    
    pData = as.data.frame(pheno)
    
    #-----------------------------------------------------
    # Rebuild rest eset
    #-----------------------------------------------------
    
    # drop = FALSE prevents to loose the column name in case there is only 1 sample in one exprs
    d1 = d1[cg, , drop = FALSE];
    d2 = d2[cg, , drop = FALSE];
    
    
    #eset1 = new("ExpressionSet");    
    #exprs(eset1) = cbind(d1, d2);
    eset1=ExpressionSet(assayData =as.matrix(cbind(d1, d2)))
    pData(eset1) = pData;
    fData(eset1) = fData;
    
    annot1 = c(annot1, annotation(eset2));
    
  }
  
  annotation(eset1) = unique(annot1);
  return(eset1);
}
#-------------------------------------------------------------------------------
merge = function(esets, method= "NONE")
{
  
  if(!is.element(method, MERGING_METHODS))
  {
    err("Merging Method: ",method," is not supported.");
  }
  
  
  if(method == "NONE") {
    msg("Run with no additional merging technique...");
    mergedEset = mergeNONE(esets);
  } 	else if(method=="BMC") 
  { 
    msg("Run BMC...");
    mergedEset = mergeBMC(esets);      
  } else if(method=="COMBAT")   
  { 
    msg("Run COMBAT...");
    mergedEset = mergeCOMBAT(esets);
  } else if(method=="GENENORM")
  { 	
    msg("Run GENENORM..."); 
    mergedEset = mergeGENENORM(esets); 
  } else if(method=="XPN")
  { 
    msg("Run XPN...");
    mergedEset = mergeXPN(esets);
  }
  
  return(matchExprsPheno(mergedEset));
  
}
# apply xpn on list of eSets
#-------------------------------------------------------------------------------
mergeXPN = function(esets)
{
  aux = function(x,y)
  {
    res = xpn(x,y);
    return(mergeNONE(list(res[[1]], res[[2]])));
  }
  
  return(combineByTwo(esets,aux,rand=TRUE));
}
# plot gene-wise box
#-------------------------------------------------------------------------------
plotGeneWiseBoxPlot = function(eset, colLabel, batchLabel, gene=NULL, 
                               legend=TRUE, file=NULL, ...)
{
  if(! is.null(file)) { pdf(file, width=12, height=7); }
  
  if(is.null(gene)) 
  { 
    s = sample(1:nrow(exprs(eset)), 1); 
    gene = rownames(exprs(eset))[s];
  }
  
  colMap = makeColorMap(eset, colLabel); 
  
  #--Add margin to the right for the legend 
  tmp = par()$mar;
  if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }
  
  #--Get values per class and batch for boxplot
  values = list();
  names = NULL;
  colVec = NULL;
  
  batches = sort(unique(pData(eset)[,batchLabel]));
  for(batch in batches)
  {  
    idx1 = pData(eset)[,batchLabel]==batch;
    classes = sort(unique(pData(eset)[idx1,colLabel]));
    for(class in classes)
    {
      idx2 = pData(eset)[,colLabel]==class;
      idx = idx1 & idx2;
      values[[length(values)+1]] = as.vector(exprs(eset)[gene,idx]);
      names = c(names,batch);
      colVec = c(colVec,colMap[[class]]);
    }
    # Mimick empty bar between batches
    values[[length(values)+1]] = NA
    names = c(names,"");
    colVec = c(colVec,"white");
  } 
  
  #--Calculate plot boundaries 
  min_x = 0;
  max_x = length(values);
  min_y = min(unlist(values), na.rm=TRUE);
  max_y = max(unlist(values), na.rm=TRUE);
  
  #--Create dummy empty plot and set axis
  plot(1,
       xlab="",
       ylab="",
       xaxt="n",
       yaxt="n",
       xlim=c(min_x,max_x),
       ylim=c(min_y,max_y));
  
  #--Add boxplots
  boxplot(values, 
          col=colVec,
          outline=FALSE, 
          add = TRUE,
          range=0,
          xaxt="n",
          panel.first={U = par("usr"); rect(U[1],U[3],U[2],U[4],
                                            col="azure2",
                                            border="black",
                                            lwd=3);},
          ...,
          main=gene);
  
  #--Reset X-axis
  axis(1,at=min_x+1:max_x,labels=names, tck=0, las=2); 
  
  #--Add legend
  if(legend)
  {
    x = max_x + (max_x - min_x) * 0.1;
    y = max_y - (max_y - min_y) * 0.1; 
    
    legend(x,y,
           legend = names(colMap),
           pt.lwd=2,
           pch = 19,
           col = unlist(colMap),
           box.lwd=3,
           bg="azure2");
    
    #--Reset margin
    par(xpd=F, mar=tmp);
  }
  
  if(! is.null(file)) { dev.off(); }
}
# plot colored MDS plot
#-------------------------------------------------------------------------------
plotMDS = function(eset, colLabel, symLabel, legend=TRUE, file=NULL, ...)
{
  if(! is.null(file)) { pdf(file, width=12, height=7); }
  
  mds = cmdscale(dist(t(exprs(eset))), eig=TRUE);
  
  colMap = makeColorMap(eset, colLabel); 
  colVec = makeColorVec(eset, colLabel, colMap);
  
  #--Add margin to the right for the legend 
  tmp = par()$mar;
  if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }
  
  range_x = range(mds$points[,1]);
  range_y = range(mds$points[,2]);
  
  plot(mds$points,
       col=colVec,
       pch=as.numeric(pData(eset)[,symLabel]),
       panel.first={ U = par("usr");
       rect(U[1],U[3],U[2],U[4],
            col="azure2",
            border="black",
            lwd=3)},
       lwd=2,
       xlab="",
       ylab="",
       xlim=range_x,
       ylim=range_y,
       ...);
  
  if(legend)
  {
    x = range_x[2] + (range_x[2]-range_x[1])*0.1;
    y = range_y[2] - (range_y[2]-range_y[1])*0.1;
    
    syms = unique(pData(eset)[,symLabel])
    legend(x,y,
           legend = syms,
           pt.lwd=2,
           pch = as.numeric(syms),
           box.lwd=3,
           bg="azure2");
    
    legend(x,y-(length(syms)*(range_y[2]-range_y[1])*0.1),
           legend = names(colMap),
           pt.lwd=2,
           pch=19,
           col = unlist(colMap),
           box.lwd=3,
           bg="azure2");
    
    #--Reset margin
    par(xpd=F, mar=tmp)
  }
  
  if(! is.null(file)) { dev.off(); }
}
# plot colored RLE plot
#-------------------------------------------------------------------------------
plotRLE = function(eset, colLabel, legend=TRUE, file=NULL, ...)
{
  if(! is.null(file)) { pdf(file, width=12, height=7); }
  
  colMap = makeColorMap(eset, colLabel); 
  colVec = makeColorVec(eset, colLabel, colMap);
  
  deviations = exprs(eset) - rowMedians(exprs(eset));
  
  #--Add margin to the right for the legend 
  tmp = par()$mar;
  if(legend) { par(xpd=T, mar=par()$mar+c(0,0,0,8)); }
  
  min_x = 1;
  max_x = ncol(eset);
  min_y = min(deviations);
  max_y = max(deviations);
  #min_y = -2;
  #max_y = 2;
  
  #--Create dummy empty plot and set axis
  plot(1,
       xlab="",
       ylab="",
       xaxt="n",
       xlim=c(1,ncol(eset)),
       ylim=c(min_y,max_y));
  
  #--Add layout
  U = par("usr");
  lines(c(U[1],U[2]),c(0,0),
        lwd=1,
        panel.first={rect(U[1],U[3],U[2],U[4],
                          col="azure2",
                          border="black",
                          lwd=3);});
  
  #--Add boxplots
  boxplot.matrix(deviations, 
                 col=colVec,
                 outline=FALSE, 
                 add = TRUE,
                 #names=rep("",ncol(eset)),
                 #range=1.5,
                 range=0,
                 xaxt="n",
                 ...);
  
  #--Add legend
  if(legend)
  {
    x = max_x + (max_x - min_x) * 0.1;
    y = max_y - (max_y - min_y) * 0.1; 
    
    legend(x,y,
           legend = names(colMap),
           pt.lwd=2,
           pch = 19,
           col = unlist(colMap),
           box.lwd=3,
           bg="azure2");
    
    #--Reset margin
    par(xpd=F, mar=tmp);
  }
  
  if(! is.null(file)) { dev.off(); }
}
.test <- function() BiocGenerics:::testPackage("inSilicoMerging")
library(inSilicoDb);

#-------------------------------------------------------------------------------
pp  = function(...) { paste(...,sep=""); }
err = function(...) { stop(...,call.=FALSE); }
msg = function(...) { message(paste("  INSILICOMERGING:",...)); } 
war = function(...) { msg(" ! WARNING ! ",...); }

#-------------------------------------------------------------------------------
MERGING_METHODS = c("BMC", "COMBAT", "NONE", "GENENORM", "XPN");

#-------------------------------------------------------------------------------
isOdd = function(x) { as.logical(x%%2) };

# Combines recursively all elements of a list pairwise in the following manner:
# [ A ; B ; C ; D ; E ]
#   => [ E ; AB ; CD ]
#   => [ CD ; EAB ]
#   => [ CDEAB ]
#-------------------------------------------------------------------------------
combineByTwo = function(x, fun, rand=FALSE, ...)
{
  if(rand) { x = sample(x); }
  
  n = length(x);
  if(n==1) { return(x[[1]]); }
  
  res = NULL;
  if(isOdd(n)) { res = c(res,x[[n]]); }
  for(i in seq(1,n-1,by=2))
  {
    res = c(res,fun(x[[i]],x[[i+1]],...));
  }
  combineByTwo(res,fun);
}

#-------------------------------------------------------------------------------
commonGenes = function(lst)
{
  common_genes = identify_common_genes(lst);
  lst = lapply(lst,function(x){return(x[rownames(exprs(x))%in%common_genes,])});
  return(lst);
}

#-------------------------------------------------------------------------------
identify_common_genes = function(lst)
{
  temp = rownames(exprs(lst[[1]]))
  for(eset in lst)
  {
    temp = intersect(rownames(exprs(eset)),temp);
  }
  return(temp);
}

#-------------------------------------------------------------------------------
matchExprsPheno = function(eset) {
  #make sure the samples are ordered in the same way
  sampleNames = colnames(exprs(eset));  
  pData(eset) = pData(eset)[sampleNames, , drop = FALSE];
  return(eset);
}

#-------------------------------------------------------------------------------

##--------------------------------------------------------------------
# xpn(eset1,eset2)
# accum_array
##--------------------------------------------------------------------
# Implemenation is based on XPN implementation in Matlab accompanying 
# the paper: A.A. Shabalin et al., Merging two gene-expression studies
# via cross-platform normalization.Bioinformatics, Vol. 24, no. 9, 
# 2008, pages 1154-1160
#
# Code based on matlab code taken from: https://genome.unc.edu/xpn/
##--------------------------------------------------------------------

repmat = function(a,n,m) {kronecker(matrix(1,n,m),a)}

#---------------------------------------------------------------------
xpn=function(eset1,eset2,geneClusters=25,sampleClusters=5,nRep=32)
{
  #Find common genes
  res = commonGenes(list(eset1,eset2));
  eset1 = res[[1]];
  eset2 = res[[2]];
  
  x1 = exprs(eset1);
  x2 = exprs(eset2);
  
  #Size check
  p=nrow(x1)
  n1=ncol(x1)
  n2=ncol(x2)
  if(nrow(x1)!=nrow(x2))
  {
    err("Different number of genes!")
  }
  
  #Check for NaNs
  if(any(is.nan(c(x1))) || any(is.nan(c(x2))))
  {
    err("Missing values in input data!")
  }
  
  #Check for geneClusters argument
  if(nargs()<3)
  {
    geneClusters=min(25,max(floor(p/50),1))
  }
  
  #Check for sampleClusters argument
  if(nargs()<4)
  {
    sampleClusters=min(5,floor(min(n1,n2)/4))
  }
  
  if(geneClusters <= 0)
  {
    err("geneClusters is too small!")
  }
  
  if(sampleClusters <= 0)
  {
    err("sampleClusters is too small!")
  }
  
  if(nRep <= 0)
  {
    err("nRep is too small!")
  }
  
  if(nRep < 15)
  {
    war("nRep is small, should be at least 15!")
  }
  
  #Check for constant rows
  constrows<-sd(t(cbind(x1,x2)))==0
  if(any(constrows))
  {
    if(all(constrows))
    {
      war("Constant input data! Will be removed.")
    }
    
    x1fix=x1[as.numeric(!constrows)*1:n1,]
    x2fix=x2[as.numeric(!constrows)*1:n2,]
    
    exprs(eset1) = x1fix;
    exprs(eset2) = x2fix;
    
    XPN_result=xpn(eset1,
                   eset2,
                   geneClusters,
                   sampleClusters,
                   nRep)
    y1=XPN_result[[1]]
    y2=XPN_result[[2]]
    
    z1=matrix(0,p,n1)
    z1[as.numeric(!constrows)*1:n1,]=y1
    z1[as.numeric(constrows),]=x1[as.numeric(constrows),]
    
    z2=matrix(0,p,n2)
    z2[as.numeric(!constrows)*1:n2,]=y2
    z2[as.numeric(constrows),]=x2[as.numeric(constrows),]	
    return(list(z1,z2))
  }
  
  #Check skewness
  z=c(x1)
  m=mean(z)
  s=sd(z)
  skewness1=mean((z-m)^3)/s^3
  
  z=c(x2)
  m=mean(z)
  s=sd(z)
  skewness2=mean((z-m)^3)/s^3
  
  if(abs(skewness1)>8 || abs(skewness2)>8)
  {
    war("Input data is excessively skewed")
    war("The data may require log-transformation")
  }
  
  #Transform platforms to have 0 mean
  x1med=matrix(0,p,1)
  x2med=matrix(0,p,1)
  for(i in 1:p)
  {
    x1med[i]=median(x1[i,])
    x2med[i]=median(x2[i,])
  }
  
  x1primeMed=x1-repmat(x1med,1,n1)
  x2primeMed=x2-repmat(x2med,1,n2)
  xxprimeMed=cbind(x1primeMed, x2primeMed) 
  
  TheSum1=matrix(0,nrow(x1),ncol(x1))
  TheSum2=matrix(0,nrow(x2),ncol(x2))
  
  #Begin Cycle
  for(cycle in 1:nRep)
  {
    IDR=kmeans(xxprimeMed,geneClusters,1000,1)[[1]]
    repeatje=TRUE
    counter = 0;  # Count the number of troublesome clusters
    maxCounter = 200;
    while(repeatje)
    {
      if(counter > maxCounter) 
      { err("Unable to find clusters. Number of samples too small..."); }
      IDC=kmeans(t(xxprimeMed),sampleClusters,1000,1)[[1]]
      IDC1=IDC[1:n1];
      IDC2=IDC[(n1+1):(n1+n2)];
      countL1 = matrix(table(IDC1),1); 
      countL2 = matrix(table(IDC2),1);
      if (length(countL1) < sampleClusters || length(countL2) < sampleClusters)
      { 
        #msg("Troublesome cluster:",counter,"/",maxCounter); 
        counter = counter + 1; 
        repeatje=TRUE; 
      }
      else
      { 
        if((cycle && 8)==0) { msg("Iteration:",cycle,"/",nRep); } 
        repeatje=FALSE; 
      }
    }
    
    #GLP means and GP variance G -gene, L -colCluster, P- study
    sumGL1=accum_array(x1,IDC1,sampleClusters)
    sumGL2=accum_array(x2,IDC2,sampleClusters)
    
    sumGL1sq=accum_array(x1^2,IDC1,sampleClusters)
    sumGL2sq=accum_array(x2^2,IDC2,sampleClusters)
    
    xbar1=sumGL1 / repmat(countL1,p,1)
    xbar2=sumGL2 / repmat(countL2,p,1)
    
    sigmaGL1sq=sumGL1sq / repmat(countL1,p,1) - xbar1^2
    sigmaGL2sq=sumGL2sq / repmat(countL2,p,1) - xbar2^2
    
    #Apply MLE to each gene cluster
    A1=matrix(0,geneClusters,sampleClusters)
    A2=matrix(0,geneClusters,sampleClusters)
    b1=matrix(0,p,1)
    b2=matrix(0,p,1)
    c1=matrix(0,p,1)
    c2=matrix(0,p,1)
    s1=matrix(0,p,1)
    s2=matrix(0,p,1)
    
    #print("Start core function XPN_MLE")
    for(i in 1:geneClusters)
    {
      xpn_result=XPN_MLE(xbar1[IDR==i,],sigmaGL1sq[IDR==i,],countL1)
      A1[i,]=xpn_result[[1]]
      b1[IDR==i] = xpn_result[[2]]
      c1[IDR==i] = xpn_result[[3]]
      s1[IDR==i] = xpn_result[[4]]
      
      xpn_result=XPN_MLE(xbar2[IDR==i,],sigmaGL2sq[IDR==i,],countL2)
      A2[i,]=xpn_result[[1]]; 
      b2[IDR==i] = xpn_result[[2]]; 
      c2[IDR==i] = xpn_result[[3]];
      s2[IDR==i] = xpn_result[[4]]; 
      
    }
    
    #Average the estimates accross platforms
    A = (A1 * repmat(countL1,geneClusters,1) + A2 * repmat(countL2,geneClusters,1))/repmat(countL1+countL2,geneClusters,1)
    b = (b1 %*% n1 + b2 %*% n2)/(n1+n2)
    c = (c1 %*% n1 + c2 %*% n2)/(n1+n2)
    s = (s1 %*% n1 + s2 %*% n2)/(n1+n2)
    
    #KsiP residual of the model - P - study
    Ksi1 = (x1 - A1[IDR,IDC1] * repmat(b1,1,n1) - repmat(c1,1,n1))/repmat(sqrt(s1),1,n1)
    Ksi2 = (x2 - A2[IDR,IDC2] * repmat(b2,1,n2) - repmat(c2,1,n2))/repmat(sqrt(s2),1,n2)
    
    #xPstar transformed data - P - study
    x1star = A[IDR,IDC1] * repmat(b,1,n1) + repmat(c,1,n1) + Ksi1 * repmat(sqrt(s),1,n1)
    x2star = A[IDR,IDC2] * repmat(b,1,n2) + repmat(c,1,n2) + Ksi2 * repmat(sqrt(s),1,n2)
    
    #Check whether NaNs are created by XPN code
    if(any(any(is.nan(x1star))) || any(any(is.nan(x2star))))
    {
      err("NaN values are produced by the XPN code!")
    }
    
    #Aggregate results of several tries in TheSumP
    TheSum1=TheSum1+x1star
    TheSum2=TheSum2+x2star
    
  }#end cycle 1:nRep
  
  #Take the average over several tries
  z1=TheSum1/nRep
  z2=TheSum2/nRep
  
  exprs(eset1) = z1;
  exprs(eset2) = z2;
  return(list(eset1,eset2))
}

#Internal MLE procedure for XPN
XPN_MLE=function(xbar, sigma2,nj)
{
  repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}
  #xbar matrix of averages of X within column clusters (given gene, platform)
  #sigma2 - variance of the elements within column clusters (given gene,platform)
  #nj - number of elements in the column clusters
  
  #Remark: the platform/gene cluster is fixed
  
  #n = number of samples in the platform
  n=sum(nj)
  
  #I - number of genes in the gene cluster
  #J - number of sample clusters
  
  ##--> UPDATE 30Nov2009 JONATAN
  
  if(!is.matrix(xbar))
  {
    I = 1   
    J = length(xbar)
  }
  else
  {
    I=nrow(xbar)
    J=ncol(xbar)
  }
  
  #I=nrow(xbar)
  #J=ncol(xbar)
  
  ##--> UPDATE 30Nov2009 JONATAN
  
  #A,b,c,s2-parameters of the model
  A=matrix(0,1,J)
  b=matrix(1,I,1)
  c=matrix(0,I,1)
  s2=matrix(1,I,1)
  
  #previous values
  old=cbind(A,t(b),t(c),t(s2))*0
  current=cbind(A,t(b),t(c),t(s2))
  
  while(sum((current-old)^2)>1e-16 * max(old))
  {
    old=current
    
    #iteratively update the values
    c = rowSums((xbar - repmat(b,1,J)*repmat(A,I,1)) * repmat(nj,I,1))/n
    
    #fix sign of b
    if(sum(b)<0) { b = -b; }
    
    A = colSums(repmat(b,1,J) * (xbar - repmat(c,1,J)) / repmat(s2,1,J),1)/sum(b^2/s2)
    
    #Enforce constraints on A
    A=A-mean(A)
    A = A * sqrt(J/sum(A^2))
    
    b = rowSums(repmat(t(A),I,1) * (xbar - repmat(c,1,J)) * repmat(nj,I,1),2)/sum(A^2*nj)
    s2=rowSums(((xbar - repmat(c,1,J) - repmat(t(A),I,1) * repmat(b,1,J))^2 + sigma2)*repmat(nj,I,1))/n
    
    ##--> UPDATE 30Nov2009 JONATAN
    
    #s2(s2==0) = realmin('double')
    #s2[s2==0]=1e-16
    s2[s2<1e-16] = 1e-16
    
    ##--> UPDATE 30Nov2009 JONATAN
    
    A=t(A)
    current=cbind(A,t(b),t(c),t(s2))
  }
  
  #Return something
  return(list(A,b,c,s2))
}

#-------------------------------------------------------------------------------

accum_array = function(my_data, clustering, nr_clusters)
{
  accumulated=matrix(0,nrow(my_data),nr_clusters)
  for(j in 1:nr_clusters)
  {
    accumulated[,j]=rowSums(my_data[,clustering==j,drop=F])
  }
  return(accumulated)
}



#setwd('/pub1/data/mg_projects/projects/web_script/R/')
library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/b827ca09ff34750b619c027b77f444d4/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/c21c7bf994ef50564605c70a902c0fa',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()
mgBox=cbind()
mgDen=rbind()
umap.dt=cbind()
tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "GEO Data press"))
  logs=c(logs,paste0('run merge_exp.R-',basename(opt$outfile)))
  #logs=c(logs,paste0('geting data:',paste0(paste0(names(opt),'=',opt),collapse = ',')))
  #library("rjson")
  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 100)
  md=c('XPN','COMBAT','BMC','GENENORM','limma')
  paths=unlist(data$paths)
  method=unlist(data$method)#method='limma'
  logs=c(logs,'数据矩阵：',paths,'方法：',method)
  #paths=c(path1,path2,path3)
  #path1='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/merged_input.txt'
  #path2='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/GSE63127.txt'
  #path3='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/limma_input.txt'
  #basename(path1) 

  
  set_all=list()
  set_name=c()
  set_sName=c()
  all_genes=c()
  for(pth in paths){
    fn=sub('\\..*$', '', gsub('^\\.','',basename(pth)))
    dat=data.table::fread(pth, sep = "\t",header = T,stringsAsFactors = F,check.names = F,na.strings="NA",data.table = F)
    unms=unique(dat[,1])
    unms=unms[which(unms!='')]
    dat=dat[match(unms,dat[,1]),]
    row.names(dat)=dat[,1]
    dat=dat[,-1]
    logs=c(logs,paste0('读取数据',pth,',行：',nrow(dat),'列：',ncol(dat)))
    eset1=ExpressionSet(assayData = as.matrix(dat))
    set_all=c(set_all,list(eset1))
    set_name=c(set_name,rep(fn,ncol(dat)))
    set_sName=c(set_sName,colnames(dat))
    all_genes=rbind(all_genes,cbind(row.names(dat),fn))
  }
  #head(all_genes)
  colnames(all_genes)=c('V1','V2')
  all_genes=as.data.frame(all_genes)
  #head(all_genes[,2])
  all_gene_tab=data.table::dcast(all_genes,V1~V2)
  #head(all_gene_tab)
  row.names(all_gene_tab)=all_gene_tab[,1]
  all_gene_tab=all_gene_tab[,-1]
  #head(all_gene_tab)
  all_gene_tab=apply(all_gene_tab,2,function(x){
    return(ifelse(is.na(x),0,1))
  })
  all_gene_tab=all_gene_tab[,match(unique(set_name),colnames(all_gene_tab))]
  all_gene_tab.ven=apply(all_gene_tab, 1, function(x){return(paste0(x,collapse = '/'))})
  all_gene_tab.ven=cbind(names(all_gene_tab.ven),all_gene_tab.ven)
  colnames(all_gene_tab.ven)=c('Tag',paste0(colnames(all_gene_tab),collapse = '/'))
  logs=c(logs,paste0('统计数据集交集并输出韦恩图'))
  write.table(all_gene_tab.ven
              ,file = paste0(opt$outfile,'/data_venn.txt'),row.names = F,col.names = T,quote = F,sep = '\t')
  mgBox=cbind(set_sName,match(set_name,unique(set_name)))
  umap.dt=cbind(set_sName,match(set_name,unique(set_name)))
  
  logs=c(logs,paste0('合并数据'))
  mgt=merge(set_all,method="NONE")
  #str(set_all[[1]])
  
  #dim(set_all[[1]]@assayData$exprs)
  #dim(set_all[[2]]@assayData$exprs)
  #set_all[[1]]@assayData$exprs[c(30,90),c(90)]
  
  mgt.xp=exprs(mgt)
  logs=c(logs,paste0('输出合并数据'))
  write.table(cbind(Tag=row.names(mgt.xp),mgt.xp)
              ,file = paste0(opt$outfile,'/before_merge.txt'),row.names = F,col.names = T,quote = F,sep = '\t')
  
  library(umap)
  nn=ceiling(ncol(mgt.xp)/3)
  if(nn>15) nn=15
  if(nn<2) nn=2
  logs=c(logs,paste0('设置UMAP最近邻个数为',nn))
  getUmap=function(exp){
    exp <- exp[!duplicated(exp), ]
    iris.umap = umap(t(exp), n_neighbors = nn, random_state = 123)
    umaply=iris.umap$layout
    umaply=umaply[match(colnames(exp),row.names(umaply)),]
    #row.names(umaply)=set_name[match(row.names(umaply),colnames(mtr))]
    #write.table(cbind(row.names(umaply),umaply),file = paste0(opt$outfile,'/umap.tab')
    #            ,quote = F,row.names = F,col.names = F,sep = '\t')
    return(umaply)
  }
  getDensty=function(exp,group){
    all_den=rbind()
    for(s in unique(group)){
      dt=density(as.numeric(exp[,which(group==s)]))
      all_den=rbind(all_den,cbind(s,dt$x,dt$y))
    }
    all_den[,1]=match(all_den[,1],unique(group))
    return(all_den)
  }
  
  before.box=t(apply(mgt.xp, 2, function(x){
    quantile(x,seq(0,1,0.25), na.rm=TRUE)
  }))
  mgBox=cbind(mgBox,before.box)
  before.den=getDensty(mgt.xp,set_name)
  mgDen=rbind(mgDen,cbind(0,before.den))
  umap.dt=cbind(umap.dt,getUmap(mgt.xp))

  mtr=NULL
  logs=c(logs,paste0('去除批次效应'))
  if(stringr::str_to_lower(method)=='limma'){
    mtr=limma::removeBatchEffect(exprs(mgt),batch = set_name)
  }else if(stringr::str_to_upper(method)%in%md){
    mtr=merge(set_all,method=stringr::str_to_upper(method))
    mtr=exprs(mtr)
  }
  if(!is.null(mtr)){
    logs=c(logs,paste0('输出去除批次效应的数据集'))
    write.table(cbind(Tag=row.names(mtr),mtr)
                ,file = paste0(opt$outfile,'/after_merge.txt'),row.names = F,col.names = T,quote = F,sep = '\t')
    myt.box=t(apply(mtr, 2, function(x){
      quantile(x,seq(0,1,0.25))
    }))
    mgBox=cbind(mgBox,myt.box)
    myt.den=getDensty(mtr,set_name)
    mgDen=rbind(mgDen,cbind(1,myt.den))
    umap.dt=cbind(umap.dt,getUmap(mtr))
    logs=c(logs,'succ all')
  }

},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(mgBox
              ,file = paste0(opt$outfile,'/boxplot.mtx'),row.names = F,col.names = F,quote = F,sep = '\t')
  write.table(mgDen
              ,file = paste0(opt$outfile,'/densty.mtx'),row.names = F,col.names = F,quote = F,sep = '\t')
  write.table(umap.dt
              ,file = paste0(opt$outfile,'/umap.mtx'),row.names = F,col.names = F,quote = F,sep = '\t')
  
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})

