#' Title
#'
#' @param fsdata
#' @param r
#'
#' @return
#' @export
#'
#' @examples
featureSelectionPlot = function(fsdata,r=0.6){

  counts = getVarsFromFS(fsdata,r,T)
  names(counts) = prettyPredictorName(names(counts))
  fsin = getVarsMetaDataFromFS(fsdata,r,F)$meaneffects
  fsin = as.data.frame(fsin,stringsAsFactors=F)
  fsin$meaneffect = -1*as.numeric(fsin$meaneffect)
  fsin$meaneffect = 100*fsin$meaneffect/max(fsin$meaneffect)
  fsin = fsin[order(fsin$meaneffect,decreasing=T),]
  enrichedfeatures = fsin$feature[1:10]
  enrichedfeffect = fsin$meaneffect[1:10]
  last = nrow(fsin)
  depletedfeatures = fsin$feature[last:last-5]
  depletedfeffect = -1*fsin$meaneffect[last:last-5]

  colors=rep("green",length(fsin$meaneffect))
  xmin = min(fsin$meaneffect) + min(fsin$meaneffect)*0.5
  xmax = 2*max(fsin$meaneffect)
  colors[fsin$meaneffect < 0] = "red"

  #Now we have to order based on categories
  fsin$feature = prettyPredictorName(fsin$feature)
  lineheight = NULL
  #First exact
  appearance = NULL
  alloccs = NULL
  labels = c("Genetic constraints", "Co-exp. MM", "Co-exp. Adjacency",
             "Expression","Gene complexity")
  usedLabels = NULL
  i = 1
  for(pattern in c("gnomad|RVIS|LoFTool","^MM","^Adjacency","^Expression")){
    occurrences = grep(pattern,fsin$feature)
    if(length(occurrences) > 0){
      #Then complexity
      appearance = c(appearance,occurrences)
      usedLabels = c(usedLabels,labels[i])
      if(length(lineheight))
        lineheight = c(lineheight,tail(lineheight,n=1) + length(occurrences))
      else
        lineheight = length(occurrences)
    }
    i = i +1
  }
  if(length(which(!(1:length(fsin$feature) %in% appearance)))){
    appearance = c(appearance,which(!(1:length(fsin$feature) %in% appearance)))
    #lineheight = c(lineheight,tail(lineheight,n=1) + length(occurrences))
    usedLabels = c(usedLabels,labels[i])
  }

  fsin$meaneffect = fsin$meaneffect[appearance]
  fsin$feature = fsin$feature[appearance]
  colors = colors[appearance]

  pointsizes = counts[match(fsin$feature,names(counts))]
  pointsizes = 2*pointsizes/max(pointsizes)
  plot(x=fsin$meaneffect,y=1:length(fsin$meaneffect),
       main="Features selected for ML",
       xlim=c(xmin,xmax),
       type="p",col=colors,pch=19,cex=pointsizes,
       xlab=paste0("regression coefficient"),
       ylab="Predictor categories")

  abline(v=0,lty="dotted",col="blue")
  for(lh in lineheight){
    abline(h=lh+0.5,lty="dotted",col="blue")
  }

  text(x=rep(xmin+xmin*0.1,5),y=c(0,lineheight[1:(length(lineheight))])+1,pos=4,
       labels=usedLabels,cex=0.5,col="blue")


  text(x=fsin$meaneffect,y=1:length(fsin$meaneffect),
       cex=0.5,
       labels=prettyPredictorName(fsin$feature))

  legend("topright",fill=c("red","green"),
         title="Meaning of direction",
         col=c("red","green"),
         legend=c("depleted for disease","enriched for disease"),cex=0.6)

}


umapPlot = function(fsdata,r=0.6,ensemble){

  vars = getVarsFromFS(fsdata,r,F)
  ctrltype = ensemble[[1]]$controls
  panel = ensemble[[1]]$panel
  controlgenes = NULL

  if(ctrltype == "allghosh" || ctrltype == "allgenome"){
    controlgenes = ensemble[[1]]$controlgenes
  }
  stopifnot(!is.null(controlgenes))

  ncontrols = length(controlgenes)
  diseasegenes = unique(unlist(lapply(ensemble,function(x){return(x$genes[!(x$genes %in% x$controlgenes)])})))
  alldata = fromGenes2MLData(genes=diseasegenes,
                             which.controls="allgenome")

  dataforumap = alldata[,vars]

  umapresult = umap(scale(dataforumap))
  naxes = 2
  allpcadata = as.data.frame(cbind(umapresult$layout,
                                   alldata$gene,
                                   alldata$condition),
                             stringsAsFactors=F)
  colnames(allpcadata) = c(paste0("Axis",1:2),"gene","condition")
  for(i in 1:naxes)
    allpcadata[,i] = as.numeric(allpcadata[,i])
  alldata = allpcadata
  cols = rep("lightgrey",nrow(alldata))
  cols[alldata$gene %in% diseasegenes] = "orange"
  cols[alldata$gene %in% ensemble$eval$finalpreds] = "red"
  theorder = order(cols)
  allpcadata = allpcadata[theorder,]

  axis1=1
  axis2=2
  xlab="1st UMAP Axis"
  ylab="2nd UMAP Axis"

  title = panel
  plot(allpcadata[,axis1],allpcadata[,axis2],col=cols[theorder],
       cex=0.5,type="p",
       xlab=xlab,
       ylab=ylab,
       main=paste0(title," : ",sum(cols=="red")," ",sum(cols=="orange")))
  legend("bottomright",fill=c("grey","orange","red"),
         title="Gene types",
         col=c("grey","orange","red"),
         legend=c("genome","disease","predictions"),cex=0.6)
}

pcaPlot = function(fsdata,r=0.6,ensemble,bestPCAs=F){

  vars = getVarsFromFS(fsdata,r,F)
  ctrltype = ensemble[[1]]$controls
  panel = ensemble[[1]]$panel
  controlgenes = NULL

  if(ctrltype == "allghosh" || ctrltype == "allgenome"){
    controlgenes = ensemble[[1]]$controlgenes
  }
  stopifnot(!is.null(controlgenes))

  ncontrols = length(controlgenes)
  diseasegenes = unique(unlist(lapply(ensemble,function(x){return(x$genes[!(x$genes %in% x$controlgenes)])})))
  alldata = fromGenes2MLData(genes=diseasegenes,
                             which.controls="allgenome")

  dataforpca = alldata[,vars]
  dataforpca = prcomp(dataforpca,scale=T,retx=T)
  naxes = which(cumsum(dataforpca$sdev^2/sum(dataforpca$sdev^2)) >= 0.8)[1]
  #And now the plot
  allpcadata = as.data.frame(cbind(dataforpca$x[,1:naxes],
                                   alldata$gene,alldata$condition),
                             stringsAsFactors=F)
  colnames(allpcadata) = c(paste0("PCA",1:naxes),"gene","condition")
  for(i in 1:naxes)
    allpcadata[,i] = as.numeric(allpcadata[,i])
  alldata = allpcadata
  cols = rep("lightgrey",nrow(alldata))
  cols[alldata$gene %in% diseasegenes] = "orange"
  cols[alldata$gene %in% ensemble$eval$finalpreds] = "red"
  theorder = order(cols)
  allpcadata = allpcadata[theorder,]

  #Which axes to plot?
  if(bestPCAs){
    condition = vector(mode="numeric",length=nrow(alldata))
    condition[] = 0
    condition[alldata$condition == "Disease"] = 1
    pvals = NULL
    cors = NULL
    for(i in 1:naxes){
      singletest = cor.test(condition,allpcadata[,i])
      pvals = c(pvals,singletest$p.value)
      cors = c(cors,singletest$estimate)
    }
    cors = abs(cors)
    axis1 = order(cors,decreasing=T)[1]
    axis2 = order(cors,decreasing=T)[2]
    xlab=paste0("PCA",axis1," P < ",signif(pvals[axis1],2))
    ylab=paste0("PCA ",axis2," P < ",signif(pvals[axis2],2))
  }else{
    axis1=1
    axis2=2
    xlab="1st PCA"
    ylab="2nd PCA"
  }

  title = panel
  oldpar = par()
  plot(allpcadata[,axis1],allpcadata[,axis2],col=cols[theorder],cex=0.5,type="p",
       xlab=xlab,
       ylab=ylab,
       main=paste0(title," : ",sum(cols=="red")," ",sum(cols=="orange")))
  legend("bottomright",fill=c("grey","orange","red"),
         title="Gene types",
         col=c("grey","orange","red"),
         legend=c("genome","disease","predictions"),cex=0.6)
}


annotateWithAmelie = function(ensemble,
                              phenotype="HP:0001300",
                              getNullDistribution=F,
                              nNull=100,
                              only.journals=T){
  panel = ensemble[[1]]$panel
  genes = ensemble$eval$finalpreds
  allresults = NULL
  cat("Working with",panel,"\n")
  cat("Working with",phenotype,"\n")

  if(getNullDistribution){
    genes = getCodingGenome()
    gsize = length(ensemble$eval$finalpreds)
    for(i in 1:nNull){
      cat("Random query:",i,"\n")
      rgenes = genes[sample(1:length(genes),gsize)]
      allresults[[i]] = fromJSON(postForm('https://amelie.stanford.edu/api/',verify=F,
                                          genes=paste0(rgenes,collapse=","),
                                          phenotypes=phenotype))
    }
    return(allresults)
  }

  cat("Calling Amelie now\n")
  results = fromJSON(postForm('https://amelie.stanford.edu/api/',verify=F,
                              genes=paste0(genes,collapse=","),
                              phenotypes=phenotype))
  #results = readRDS("~/tmp/results.rds")
  cat("Done!\n")

  if(length(results)){
    for(i in 1:length(results)){
      gene = unlist(results[[i]][[1]])
      partial = results[[i]][[2]]
      if(length(partial) >0){
        partial = t(sapply(partial,function(x){ return(c(x[[1]],x[[2]]))}))[,c(1,2),drop=FALSE]
        #print(partial)
        allresults = rbind(allresults,cbind(rep(panel,nrow(partial)),
                                            rep(gene,nrow(partial)),
                                            rep(phenotype,nrow(partial)),
                                            partial))
      }

      else{
        #cat("Nothing for",gene,"\n")
        allresults = rbind(allresults,c(panel,gene,phenotype,NA,NA))
      }
    }
  }


  #Now we get the title and journal
  by=200
  n = length(unique(na.omit(allresults[,5])))
  indexes = seq(by,n,by)
  if(!(n %in% indexes))
    indexes = c(indexes,n)
  lastindex = 1
  titles = NULL
  jnames = NULL
  allids = unique(na.omit(allresults[,5]))
  journals = NULL
  result = NULL
  for(index in indexes){
    #cat("From",index,"\n")
    ids = paste0(allids[lastindex:index],collapse=",")
    result = c(result,fromJSON(postForm("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",db="pubmed",id=ids,retmode="json"))$result)
    setupids = allids[lastindex:index]
    lastindex = index
  }
  for(i in 1:nrow(allresults)){
    if(!is.na(allresults[i,5])){
      titles = c(titles,result[[as.character(allresults[i,5])]]$title)
      jnames = c(jnames,result[[as.character(allresults[i,5])]]$fulljournalname)
    }else{
      titles = c(titles,"not applicable")
      jnames = c(jnames,"not applicable")

    }
  }

  allresults = cbind(allresults,titles,jnames)
  colnames(allresults) = c("panel","gene","phenotype","confidence","pubmedid","title","journal")
  return(as.data.frame(allresults,stringsAsFactors=F))
}

amelieStudy = function(rndfile = "~/Dropbox/KCL/talks/nih2019/pdAmelieRandom.rds",
                       ameliefile="~/tmp/results.rds",
                       doplot=F,panel="PD",
                       ensemble=NULL,
                       phenotype=NULL,
                       nNull=100){
  genecount = NULL
  brutecount = NULL

  stopifnot(!is.null(rndfile) | !is.null(ensemble))

  if(is.null(rndfile))
    result = annotateWithAmelie(ensemble=ensemble,
                                phenotype=phenotype,
                                getNullDistribution=T,
                                nNull=nNull)
  else
    result = readRDS(rndfile)

  if(is.null(ameliefile))
    result[[length(result) + 1]] = annotateWithAmelie(ensemble=ensemble,
                                                      phenotype=phenotype,
                                                      getNullDistribution=F)
  else
    result[[length(result) + 1]] = readRDS(ameliefile)

  for(j in 1:length(result)){
    results = result[[j]]
    gcount = 0
    bcount = 0
    for(i in 1:length(results)){
      gene = unlist(results[[i]][[1]])
      partial = results[[i]][[2]]
      if(length(partial) > 0){
        partial = t(sapply(partial,function(x){ return(c(x[[1]],x[[2]]))}))[,c(1,2),drop=FALSE]
        #print(partial)
        rndresults = rbind(rndresults,cbind(rep(j,nrow(partial)),partial))
        gcount = gcount + nrow(partial)
        bcount = bcount + sum(as.numeric(partial[,1]))
      }
    }
    genecount = c(genecount,gcount)
    brutecount = c(brutecount,bcount)
  }

  if(doplot){
    nrandom = length(genecount) - 1
    par(mfrow=c(1,2))
    xlim=c(min(genecount),max(genecount))
    plot(density(genecount[1:nrandom]),xlab="Hits (gene, phenotype) found",
         main=paste0("Do we get more ",panel," hits?"),
         cex=0.8,xlim=xlim)
    abline(v=genecount[nrandom + 1],col="red")

    plot(density(brutecount[1:nrandom]/genecount[1:nrandom]),
         xlab="Mean score per hit",
         main=paste0("Are ",panel," hits of better quality?"),
         cex=0.8)
    abline(v=brutecount[nrandom + 1]/genecount[nrandom + 1],col="red")

  }
  par(mfrow=c(1,1))
  return(list(gcount=genecount,scorecount=brutecount,result=result))
}

seedGenesVsPredictions = function(model,
                                  panel,
                                  whereToCompare="rea"){


  preds = model$eval$finalpreds
  genes = unique(unlist(lapply(model,function(x){ return(x$genes[x$condition == "Disease"])})))
  preds = bitr(preds,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  genes = bitr(genes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  mycluster = NULL
  mycluster$Panel = genes$ENTREZID
  mycluster$Predictions = preds$ENTREZID
  #for(domain in c("BP","CC","MF")){
  cat("Going for panel ",panel,"\n")

  if(whereToCompare == "rea"){
    compFunc = "enrichPathway"
    library(ReactomePA)
    enr = compareCluster(geneCluster = mycluster,
                         fun = compFunc,
                         pAdjustMethod ="BH",
                         pvalueCutoff=0.05,
                         readable=T)
  }else if(whereToCompare == "keg"){
    compFunc = "enrichKEGG"
    enr = compareCluster(geneCluster = mycluster,
                         fun = compFunc,
                         pAdjustMethod ="BH",
                         pvalueCutoff=0.05)
  }else if(whereToCompare %in% c("BP","MF","CC")){
    compFunc = "enrichGO"
    universe = getCodingGenome()
    universe = bitr(universe,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    enr = compareCluster(geneCluster = mycluster,
                         fun = compFunc,
                         OrgDb = org.Hs.eg.db,
                         ont=whereToCompare,
                         universe=universe,
                         pAdjustMethod ="BH",
                         pvalueCutoff=0.05,
                         readable=T)

  }else stop(paste0(whereToCompare," is not something we can use\n"))



  npreds = sum(enr@compareClusterResult$Cluster == "Predictions")
  npanel = sum(enr@compareClusterResult$Cluster == "Panel")
  common = table(enr@compareClusterResult$ID)
  common = names(common)[common == 2]
  npredsr = length(common)/npreds
  npanelr = length(common)/npanel

  if(length(common) > 5){
    enr@compareClusterResult = enr@compareClusterResult[enr@compareClusterResult$ID %in% common,]
  }

  ncats = 20
  if(nrow(enr@compareClusterResult) < ncats)
    ncats = nrow(enr@compareClusterResult)
  mask = order(enr@compareClusterResult$p.adjust,decreasing=F)
  enr@compareClusterResult = enr@compareClusterResult[mask,]


  dotplot(enr,showCategory=ncats,title=paste0(panel,", ",
                                              length(common)," matchs,",
                                              "(",100*signif(npanelr,2),"% of ",npanel,"), (",
                                              100*signif(npredsr,2),"% of ",npreds,")"),
          font.size=6)
}

tsneStudy = function(fsdata,r=0.6,ensemble,outpath="~/Downloads//",altPanel=NULL){

  alldata = fromGenes2MLData(genes=NULL,which.controls="allgenome",
                                     filter=c("DPI","DSI","ESTcount","constitutiveexons"))
  allgenes = alldata$gene
  panel = ensemble[[1]]$panel
  vars = getVarsFromFS(fsdata,r,F)
  localalldata = scale(alldata[,vars])
  perplexity=50
  if(file.exists("~/tmp/result.rds"))
    result = readRDS("~/tmp/result.rds")
  else{
    result = Rtsne(localalldata, dims = 2,
                   perplexity=perplexity,
                   verbose=TRUE,
                   max_iter=300,
                   check_duplicates=F)
    saveRDS(result,"~/tmp/result.rds")
  }

  if(!is.null(altPanel))
    panel = altPanel

  colors = rep("grey",nrow(result$Y))
  genes = unique(unlist(lapply(ensemble,function(x){return(x$genes[!(x$genes %in% x$controlgenes)])})))
  colors[which(allgenes %in% genes)] = "red"
  mask = match(genes,allgenes)

  #Distance
  buddies = NULL
  for(gene in mask){
    g2all = apply(result$Y,1,function(x,y){
      sqrt(sum((x - y) ^ 2))
    },y=result$Y[gene,])
    buddies[[allgenes[gene]]] = allgenes[order(g2all,decreasing=F)][2:20]
  }

  querybuddies = NULL
  for(gene in (names(buddies))){
    querybuddies[[gene]] = c(buddies[[gene]],gene)
  }
  if(!file.exists(paste0(outpath,"/",panel,"tsnGO.csv"))){
    go = gprofiler(query=querybuddies,
                   src_filter=c("KEGG","REAC","HP","OMIM"))
    if(length(go) > 0)
      write.csv(go,paste0(outpath,"/",panel,"tsnGO.csv"))
  }

  mask2 = allgenes %in% unlist(buddies)
  predictions = ensemble$eval$finalpreds
  mask3 = allgenes %in% predictions


  #Test
  maskdisease = allgenes %in% genes
  accdist = 0
  for(gene in which(mask3 == T)){
    g2all = apply(result$Y[maskdisease,],1,function(x,y){
      sqrt(sum((x - y) ^ 2))
    },y=result$Y[gene,])
    accdist = accdist + min(g2all)
  }
  accdist = accdist/length(mask3)
  for(i in 1:100){
    rnddists = NULL
    rndmask = sample(1:length(allgenes),sum(mask3))
    singlernddist = 0
    for(gene in rndmask){
      g2all = apply(result$Y[maskdisease,],1,function(x,y){
        sqrt(sum((x - y) ^ 2))
      },y=result$Y[gene,])
      singlernddist = singlernddist + min(g2all)
    }
    singlernddist = singlernddist/sum(mask3)
    rnddists = c(rnddists,singlernddist)
  }
  foldchange = signif(mean(rnddists)/accdist,3)


  pdf(paste0(outpath,"/",panel,"tsneA.pdf"))
  plot(result$Y[colors == "grey",],t='p',xlab="1st axis",ylab="2nd axis",
       main=paste0("tsne of ",panel),
       col="grey",pch=19,cex=0.5)
  dev.off()

  pdf(paste0(outpath,"/",panel,"tsneB.pdf"))
  plot(result$Y[colors == "grey",],t='p',xlab="1st axis",ylab="2nd axis",
       main=paste0("tsne of ",panel),
       col="grey",pch=19,cex=0.5)

  lines(result$Y[colors == "red",],t='p',main="tsne",col="red",pch=19,cex=0.2)
  text(x=result$Y[mask,1],y=result$Y[mask,2],labels=genes,pos=2,cex=0.7)
  dev.off()

  pdf(paste0(outpath,"/",panel,"tsneC.pdf"))
  plot(result$Y[colors == "grey",],t='p',xlab="1st axis",ylab="2nd axis",
       main=paste0("tsne of ",panel),
       col="grey",pch=19,cex=0.5)
  lines(result$Y[mask2,],t='p',col="orange",pch=19,cex=0.2)

  lines(result$Y[colors == "red",],t='p',main="tsne",col="red",pch=19,cex=0.2)
  text(x=result$Y[mask,1],y=result$Y[mask,2],labels=genes,pos=2,cex=0.7)
  dev.off()

  pdf(paste0(outpath,"/",panel,"tsneD.pdf"))
  plot(result$Y[colors == "grey",],t='p',xlab="1st axis",ylab="2nd axis",
       main=paste0("tsne of ",panel,", rnd distance fold ",foldchange),
       col="grey",pch=19,cex=0.5)
  #lines(result$Y[mask2,],t='p',col="orange",pch=19,cex=0.2)
  lines(result$Y[mask3,],t='p',col="pink",pch=19,cex=0.2)

  lines(result$Y[colors == "red",],t='p',main="tsne",col="red",pch=19,cex=0.2)
  text(x=result$Y[mask,1],y=result$Y[mask,2],labels=genes,pos=2,cex=0.7)
  dev.off()


  print(accdist)
  print(rnddists)
}
