
#for i in `cat panelsApr2018`;
#do echo "cd /SAN/neuroscience/WT_BRAINEAC/ml/tournamentnewDB/; Rscript -e \"library(caret);library(G2PML);panel=\\\"$i\\\";saveRDS(featureSelection(genes=getGelGenes(panel=panel),k=5,repeats=40),\\\"$i.fs.rds\\\")\"" |  qsub  -S /bin/bash -N $i -l h_rt=48:0:0 -l tmem=8G,h_vmem=8G -o
#/SAN/neuroscience/WT_BRAINEAC/ml/tournamentnewDB/$i.fs.o -e /SAN/neuroscience/WT_BRAINEAC/ml/tournamentnewDB/$i.fs.e ; done

#for i in `cat panelsApr2018`; do echo
#"cd /SAN/neuroscience/WT_BRAINEAC/ml/tournamentnewDB/; Rscript -e \"library(caret);library(G2PML);panel=\\\"$i\\\";saveRDS(ensembleLearnKFold(panel=\\\"$i\\\",fs=featureSelection(genes=getGelGenes(panel=panel),k=3,repeats=20),nboot=20,auto=T,k=5,maxTrials=5),\\\"$i.Apr2018.rds\\\")\"" |  qsub -S /bin/bash -N $i -l h_rt=72:0:0 -l tmem=8G,h_vmem=8G
#-o /SAN/neuroscience/WT_BRAINEAC/ml/tournamentnewDB/$i.apr.o -e /SAN/neuroscience/WT_BRAINEAC/ml/tournamentnewDB/$i.apr.e ; done

#' Which features are most important for your gene list?
#'
#' \code{featureSelection} computes the best features that discriminate between
#' your list of disease genes and control genes. Uses bootstrapping to form
#' balanced sets of disease and non-disease genes then selects the best
#' features based on a random forest algorithm implemented through the
#' \code{\link{caret::rfe}} function.
#'
#' @param genes chr vector. Gene symbols - can be returned from
#'   \code{\link{getGenesFromPanelApp}}.
#' @param seed num scalar. Random seed for reproducibility.
#' @param useSMOTE
#' @param sizes int vector. Sizes to be used in the recursive feature
#'   elimination \code{\link{caret::rfe}}.
#' @param k int scalar. Factor by which to split training set for k-fold cross
#'   validation.
#' @param trnProp num scalar. Between 0-1 - proportion of disease genes to keep
#'   when bootstrapping.
#' @param repeats int scalar. Number of times you want to bootstrap/iterate. For
#'   each iteration, \code{featureSelection} will compute rfe on an random
#'   proportion (trnProp) of disease genes and a random size-matched set of
#'   controls.
#' @param gacontrols
#'
#' @return list of length repeats. Each element contains an rfe class object
#'   fitted for a set of randomly sampled disease and control genes.
#' @export
#'
#' @examples
#' genes <- getGenesFromPanelApp(disorder="Neurology and neurodevelopmental disorders",
#'   panel="Parkinson Disease and Complex Parkinsonism", color = "green")
#' featureSelection(genes, controls = "allgenome")
#'
#' @seealso For more details on rfe:
#'   \url{http://topepo.github.io/caret/recursive-feature-elimination.html}
featureSelection = function(genes=NULL,
                            seed=12345,
                            sizes = c(5,10,20),
                            k=5,
                            controls="allghosh",
                            trnProp=0.9,
                            repeats=10,
                            gacontrols=-1){

  library(DMwR)
  allgenes = genes
  cat("We'll work with",length(allgenes),"disease genes\n")
  mldata = fromGenes2MLData(genes=allgenes,
                            which.controls=controls)
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

#' Which features are most important for your gene list?
#'
#' \code{getVarsFromFS} converts the output of \code{\link{featureSelection}}
#' function to a list of features that were most important across bootstrapping
#'
#' @param fsdata rfe object. Output of \code{\link{featureSelection}}.
#' @param r num scalar. Between 0-1. Minimum proportion of bootstraps that must
#'   have selected that feature as important.
#' @param counts lgl scalar. If TRUE, will return counts each "important"
#'   feature appeared amongst bootstrap iterations. If FALSE, the names of the
#'   features in the same order.
#'
#' @return Either counts or names of the most important features determined by
#'   rfe.
#' @export
#'
#' @example
#' getVarsFromFS(allfsdata, r = 0.6, counts = F)
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

#'Which direction and to what extent do your important features affect the
#'gene's likelihood of being disease causing?
#'
#'\code{getVarsFromFS} extracts from the output of
#'\code{\link{featureSelection}} the most important features based a voting
#'system across bootstraps. Then for the most-voted features obtain the sign and
#'effect sizes of each feature.
#'
#'@inheritParams getVarsFromFS
#'
#'@param allfsdata list of rfe objects. Output of \code{\link{featureSelection}}.
#'@param panel chr scalar. The GE panel name from which input genes were
#'  extracted.
#'
#'@return list with all counts, signs, effect sizes and mean effect sizes of all
#'  important features across all bootstraps
#'
#'@export
#'
#'@examples
#'getVarsMetaDataFromFS(allfsdata, r = 0.6, panel = "Unknown")
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


