

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
                            useSMOTE=F,
                            sizes = c(5,10,20),
                            k=5,
                            trnProp=0.9,
                            repeats=10,
                            gacontrols=-1){

  library(DMwR)
  allgenes = genes
  cat("We'll work with",length(allgenes),"disease genes\n")
  mldata = fromGenes2MLData(genes=allgenes,
                            which.controls="allgenome",
                            filter=c("DPI","DSI","ESTcount","constitutiveexons"))
  set.seed(seed)
  fsruns = NULL

  for(run in 1:repeats){
    if(!useSMOTE){
      the_n = floor(length(genes)*trnProp)
      lmldata = mldata[c(sample(which(mldata$condition == "Disease"),the_n),
                         sample(which(mldata$condition == "Nondisease"),the_n)),]
      condition = lmldata$condition
      lmldata$gene = NULL
      lmldata$condition = NULL

    }else{
      if(gacontrols >  1){
        lmldata = mldata[c(which(mldata$condition == "Disease"),
                           sample(which(mldata$condition == "Nondisease"),gacontrols)),]
        cat("Some of the controls are",
            paste0(lmldata$gene[lmldata$condition == "Nondisease"][1:10],collapse=","),"\n")
      }
      lmldata$gene = NULL
      cat("Before applying SMOTE to conventional feature selection\n")
      print(table(lmldata$condition))
      lmldata$condition = as.factor(lmldata$condition)
      lmldata = SMOTE(condition ~.,data=lmldata,
                      perc.over=gacontrols,
                      perc.under=100)
      condition = as.character(lmldata$condition)
      lmldata$condition = NULL
      cat("After applying SMOTE\n")
      print(table(condition))

    }

    condition = as.factor(condition)
    ctrl = caret::rfeControl(functions=caret::rfFuncs,
                      method="cv",
                      number=k,
                      verbose=T,
                      saveDetails=T)

    lmProfile = caret::rfe(x=lmldata,
                    y = condition,
                    sizes = sizes,
                    rfeControl=ctrl)

    mldatasc = scale(lmldata[,lmProfile$optVariables])
    mask = !unlist(apply(mldatasc,2,function(x){ length(unique(x)) == TRUE}))
    if(sum(!mask) > 0){
      cat("Removing",paste0(colnames(mldatasc)[!mask],collapse=", "),
          "from data before going for glm\n")
      mldatasc = mldatasc[,mask]
    }

    fit = glm(condition ~.,data=data.frame(mldatasc),
              family=binomial(link = "logit"))
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
getVarsFromFS = function(fsdata,r=0.6){
  availfiles = 0
  vars = NULL
  k = length(fsdata)
  for(i in 1:k){
      fsdatal = fsdata[[i]]
      if(is.null(fsdatal$optVariables)){
        vars = c(vars,names(sort(table(fsdatal$variables$var[fsdatal$variables$Variables == fsdatal$bestSubset]),
                                 decreasing=T)))
      }else
        vars = c(vars,fsdatal$optVariables)
  }
  minAppearance = floor(r*k)
  #if(k <= 3)
  #  minAppearance = k
  #else
  #  minAppearance = k - 1

  fcounts = table(vars)
  fcounts = fcounts[fcounts >= minAppearance]
  fcounts = sort(fcounts,decreasing=T)
  unique(names(fcounts))
}


fsGetVars = function(file,k=10,exp="fsrfs"){
  vars = NULL

  if(file.exists(file)){
    fsdata = readRDS(file)
    cat("We will use non GA feature selection on vars\n")
    if(!is.null(fsdata$variables)){
      vars = names(sort(table(fsdata$variables$var[fsdata$variables$Variables == fsdata$bestSubset]),
                        decreasing=T))
    }
  }else{
    availfiles = 0
    vars = NULL
    for(i in 1:k){
      localfile = paste0(file,"_",exp,"_",i,".rds")
      if(file.exists(localfile)){
        cat("We will use GA feature selection on vars, reading fold\n")
        availfiles = availfiles + 1
        fsdata = readRDS(localfile)
        if(is.null(fsdata$optVariables)){
          vars = c(vars,names(sort(table(fsdata$variables$var[fsdata$variables$Variables == fsdata$bestSubset]),
                                   decreasing=T)))
        }else
          vars = c(vars,fsdata$optVariables)
      }else
        cat("File",localfile,"is not there\n")
    }
    if(availfiles <= 3)
      minAppearance = availfiles
    else
      minAppearance = availfiles - 1

    fcounts = table(vars)
    fcounts = fcounts[fcounts >= minAppearance]
    fcounts = sort(fcounts,decreasing=T)
    vars = names(fcounts)

  }
  return(vars)
}
