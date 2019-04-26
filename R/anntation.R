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
  legend("bottomright",fill=c("blue","orange","red"),
         title="Gene types",
         col=c("blue","orange","red"),
         legend=c("genome","disease","predictions"),cex=0.6)
}
