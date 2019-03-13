

#' Title
#'
#' @param genes
#' @param seed
#' @param useSMOTE
#' @param sizes
#' @param k
#' @param trnProp
#' @param repeats
#' @param gacontrols
#'
#' @return
#' @export
#'
#' @examples
featureSelection = function(genes=NULL,
                            seed=12345,
                            sizes = c(5,10,20),
                            k=5,
                            trnProp=0.9,
                            repeats=10,
                            gacontrols=-1){

  library(DMwR)
  allgenes = genes
  cat("We'll work with",length(allgenes),"disease genes\n")
  mldata = fromGenes2MLData(genes=allgenes,
                            which.controls="allgenome")
  mask = !(colnames(mldata) %in% c("gene","condition"))

  mldata[,mask] = scale(mldata[,mask])
  set.seed(seed)
  fsruns = NULL
  gcondition = mldata$condition

  for(run in 1:repeats){

      the_n = floor(length(genes)*trnProp)
      lmldata = mldata[c(sample(which(mldata$condition == "Disease"),the_n),
                         sample(which(mldata$condition == "Nondisease"),the_n)),]
      condition = lmldata$condition
      lmldata$gene = NULL
      lmldata$condition = NULL



    fcondition = as.factor(condition)
    ctrl = caret::rfeControl(functions=caret::rfFuncs,
                             method="cv",
                             number=k,
                             verbose=T,
                             saveDetails=T)

    lmProfile = caret::rfe(x=lmldata,
                           y = fcondition,
                           sizes = sizes,
                           rfeControl=ctrl)

    mask = !unlist(apply(lmldata,2,function(x){ length(unique(x)) == TRUE}))
    if(sum(!mask) > 0){
      cat("Removing",paste0(colnames(lmldata)[!mask],collapse=", "),
          "from data before going for glm\n")
      lmldata = lmldata[,mask]
    }
    condition[condition == "Disease"] = 0
    condition[condition == "Nondisease"] = 1
    condition = as.numeric(condition)

    lmProfile$optVariables = lmProfile$optVariables[lmProfile$optVariables %in% colnames(lmldata)]
    fit = glm(formula=condition ~ . ,data=lmldata[,lmProfile$optVariables])
    lmProfile$myFit = fit
    fsruns[[run]] = lmProfile
  }
  return(fsruns)

}

#' Title
#'
#' @param fsdata
#' @param r
#'
#' @return
#' @export
#'
#' @examples
getVarsFromFS = function(fsdata,r=0.6,counts=F){
  availfiles = 0
  vars = NULL
  k = length(fsdata)
  for(i in 1:k){
    fsdatal = fsdata[[i]]
    if(is.null(fsdatal$optVariables)){
      vars = c(vars,names(sort(table(fsdatal$variables$var[fsdatal$variables$Variables == fsdatal$bestSubset]),
                               decreasing=T)))
    }else
      vars = c(vars,fsdatal$optVariables[fsdatal$optVariables %in% names(na.omit(fsdata[[i]]$myFit$coefficients))])
  }
  minAppearance = floor(r*k)

  fcounts = table(vars)
  fcounts = fcounts[fcounts >= minAppearance]
  fcounts = sort(fcounts,decreasing=T)
  if(counts)
    return(fcounts)
  return(unique(names(fcounts)))
}



getVarsMetaDataFromFS = function(allfsdata,r=r,panel="Unknown"){

  allfeatures = NULL
  allsigns = NULL
  alleffects = NULL
  allmeaneffect = NULL
  alleffectstab = NULL
  allmeaneffectstab = NULL

  availfiles = 0
  cat("Working now with",panel,"\n")
  features = NULL
  meaneffect = NULL
  folds = length(allfsdata)

  selected = getVarsFromFS(allfsdata,r=r)
  fcounts = getVarsFromFS(allfsdata,r=r,T)
  allfeatures = fcounts
  signs = NULL
  effects = NULL
  models = NULL
  folds = length(allfsdata)
  for(k in 1:folds){
    fsdata = allfsdata[[k]]
    if(!is.null(fsdata$ga)){
      models[[as.character(k)]] = fsdata$fit
      fit = fsdata$fit$finalModel
    }else{
      models[[as.character(k)]] = fsdata$myFit
      fit = fsdata$myFit
    }

    for(sel in selected){
      if(sel %in% names(fit$coefficients)){
        mask = names(fit$coefficients) == sel
        if(is.null(signs[[sel]])){
          signs[[sel]] = NULL
          if(!is.na(fit$coefficients[mask]))
            signs[[sel]] = as.character(sign(fit$coefficients[mask]))
          else
            signs[[sel]] = "0"
          signs = as.list(signs)
        }else{
          older = signs[[sel]]
          if(!is.na(fit$coefficients[mask]))
            signs[[sel]] = paste0(older,":",as.character(sign(fit$coefficients[mask])))
          else
            signs[[sel]] = paste0(older,":0")
        }

        mask = names(fit$effects) == sel
        if(sum(mask)){
          if(is.null(effects[[sel]])){
            effects[[sel]] = NULL
            if(!is.na(fit$effects[mask])){
              effects[[sel]] = as.character(signif(fit$effects[mask],3))
              meaneffect[[sel]] = c(meaneffect[[sel]],fit$effects[mask])
            }

            else
              effects[[sel]] = "0"
            effects = as.list(effects)
            meaneffect = as.list(meaneffect)
          }else{
            older = effects[[sel]]
            if(!is.na(fit$effects[mask])){
              effects[[sel]] = paste0(older,":",as.character(signif(fit$effects[mask],3)))
              meaneffect[[sel]] = c(meaneffect[[sel]],fit$effects[mask])
            }

            else
              effects[[sel]] = paste0(older,":0")

          }
        }


      }
    }
  }

  allsigns = unlist(signs)
  alleffects = unlist(effects)
  mask = order(abs(unlist(lapply(meaneffect,mean))),decreasing=T)
  allmeaneffect = unlist(lapply(meaneffect,mean))[mask]

  localdata = as.matrix(unlist(allmeaneffect))
  localdata = cbind(rownames(localdata),localdata)
  allmeaneffectstab = rbind(allmeaneffectstab,
                            cbind(rep(panel,
                                      nrow(localdata)),
                                  localdata,
                                  rep(folds,nrow(localdata))))

  #print(allmeaneffectstab)
  localdata = as.matrix(unlist(alleffects))
  localdata = cbind(rownames(localdata),localdata)
  localdata = cbind(rep(panel,nrow(localdata)),
                    localdata)

  alleffectstab = rbind(alleffectstab,localdata)

  #print(alleffectstab)

  colnames(alleffectstab) = c("panel","predictor","foldeffects")
  colnames(allmeaneffectstab) = c("panel","feature","meaneffect","coverage")

  return(list(features=allfeatures,signs=allsigns,effects=alleffects,meaneffects=allmeaneffectstab))
}

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

