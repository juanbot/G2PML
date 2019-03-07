
#' Title
#'
#' @param methods
#' @param k
#' @param seed
#' @param nboot
#' @param tuneLength
#' @param panel
#' @param auto
#' @param maxTrials
#' @param controls
#' @param fs
#' @param fsThreshold
#' @param qmeasure
#' @param nsamps
#'
#' @return
#' @export
#'
#' @examples
ensembleLearnKFold = function(panel="Congenital_myopathy",
                              methods=c("rpart",
                                        "C5.0Tree",
                                        "svmRadial",
                                        "rf",
                                        "sparseLDA",
                                        "kknn","naive_bayes"),
                              k=5,
                              #The control genes to use
                              seed=12345,
                              #Number of simple models to create within the ensemble
                              nboot=200,
                              #This is useful to generate a null distribution of models to
                              #Number of different values to test for each ML algorithm
                              #hyperparameter
                              tuneLength=5,
                              #Number of repeats to be done by Caret
                              nsamps=5,
                              auto=F,
                              maxTrials=5,
                              controls = "allghosh",
                              fs="own",
                              fsThreshold=0.6,
                              qmeasure="kappa"){

  panelf = paste0(system.file("g2pml/gel/geMarch2017/", "", package = "G2PML"),
                  panel,".csv")

  stopifnot(file.exists(panelf))

  methods = unlist(lapply(methods,
                          function(x){ tryCatch({ checkInstall(getModelInfo(x)$library); return(x)},
                                                error = function(e){ print(e); NULL },0)}))

  cat(paste0("Available methods to use are ",paste0(methods,collapse=", ")),"\n")

  #Check which models are installed to try them


  cat("Working now with genes from",panelf,"\n")
  finalEnsembles = NULL
  #for(method in methods){
  method = methods[1]
  exprStr = ""
  expid = paste0(c(panel,method),collapse="_")
  cat("Working with",method,"\n")
  cat("Reading disease genes from",panelf,"\n")
  genes = read.csv(panelf,stringsAsFactors=F)
  genes = genes$GeneSymbol[genes$LevelOfConfidence == "HighEvidence"]
  genes = na.omit(G2PML::fromSymbol2Hugo(genes))

  if(controls == "allghosh")
    ctrlpath=system.file("g2pml/controlgenes/allghosh/", "", package = "G2PML")
  else
    stop(paste0("Control type not found:",controls))

  #Getting control genes
  ctrlfile = paste0(ctrlpath,"/",panel,"_controls.tsv")
  cat("Reading controls from",ctrlfile,"\n")
  ctrlgenes = read.delim(ctrlfile,header=F,stringsAsFactors=F)$V1
  ctrlgenes = na.omit(G2PML::fromSymbol2Hugo(ctrlgenes))
  cat("We get",length(ctrlgenes),"control genes\n")


  #Configuration parameters
  if(length(genes) < 50){
    leaveoutpos = as.integer(length(genes)*0.3)
  }else{
    leaveoutpos = as.integer(length(genes)/k)
  }
  indexespos = matrix(nrow=k,ncol=length(genes)-leaveoutpos)
  evalindexespos = matrix(nrow=k,ncol=leaveoutpos)
  leaveoutneg = as.integer(length(ctrlgenes)/k)
  indexesneg = matrix(nrow=k,ncol=length(ctrlgenes)-leaveoutneg)
  evalindexesneg = matrix(nrow=k,ncol=leaveoutneg)

  if(typeof(fs) == "character"){
    if(fs == "own")
      pathfs = system.file("g2pml/fs/", "", package = "G2PML")
    else
      stop(paste0("Feature selection strategy not found:",fs))
    vars=NULL
    fsfile = paste0(pathfs,"/featureSelection",panel)
    if(!file.exists(paste0(fsfile,"_fsrfs_1.rds")) &
       !file.exists(paste0(fsfile,"_fsrfs_2.rds")) &
       !file.exists(paste0(fsfile,"_fsrfs_3.rds")) &
       !file.exists(paste0(fsfile,"_fsrfs_4.rds")) &
       !file.exists(paste0(fsfile,"_fsrfs_5.rds")))
      fsfile=NULL
  }else
    vars = getVarsFromFS(fs,r=fsThreshold)



  for(fold in 1:k){
    sequence = 1:length(genes)
    indexespos[fold,] = sample(x=sequence,length(genes) -  leaveoutpos)
    evalindexespos[fold,] = sequence[!(sequence %in% indexespos[fold,])]
    lexpid = paste0(expid,"_kfold","_",fold)

    sequence = 1:length(ctrlgenes)
    indexesneg[fold,] = sample(x=sequence,length(ctrlgenes) -  leaveoutneg)
    evalindexesneg[fold,] = sequence[!(sequence %in% indexesneg[fold,])]

    genespos = genes[indexespos[fold,]]
    genespos = G2PML::fromSymbol2Hugo(genespos)
    genesneg = ctrlgenes[indexesneg[fold,]]
    genesneg = G2PML::fromSymbol2Hugo(genesneg)
    genestogo = c(genespos,genesneg)
    condition = c(rep("Disease",length(genespos)),
                  rep("Nondisease",length(genesneg)))

    if(auto)
      finalEnsembles[[fold]] = ensembleLearnAutonomous(panel=panel,
                                                       methods=methods,
                                                       condition=condition,
                                                       genes=genestogo,
                                                       nboot=nboot,
                                                       nsamps=nsamps,
                                                       tuneLength=tuneLength,
                                                       fsfile=fsfile,
                                                       maxTrials=maxTrials,
                                                       qmeasure=qmeasure,
                                                       vars=vars,
                                                       evalctrl=ctrlgenes[evalindexesneg[fold,]],
                                                       evaldisease=genes[evalindexespos[fold,]])
    else
      finalEnsembles[[fold]] = ensembleLearn(expID=lexpid,
                                             panel=panel,
                                             condition=condition,
                                             useSeedGene=useSeedGene,
                                             useSmote=useSmote,
                                             method=method,
                                             genes=genestogo,
                                             nboot=nboot,
                                             nsamps=nsamps,
                                             tuneLength=tuneLength,
                                             fsfile=fsfile,
                                             vars=vars,
                                             evalctrl=ctrlgenes[evalindexesneg[fold,]],
                                             evaldisease=genes[evalindexespos[fold,]])
  }

  finalEnsembles$eval = evalEnsemblesOneShot(finalEnsembles)
  return(finalEnsembles)
}

ensembleLearnAutonomous  = function(genes,
                                    panel=NULL,
                                    #If condition is NULL then all genes are disease
                                    condition=NULL,
                                    #The ID for referring to the results
                                    expID = "BootstrapEnsemble",
                                    #The caret Method to use
                                    methods=c("rpart",
                                              "C5.0Tree",
                                              "svmRadial",
                                              "rf",
                                              "sparseLDA",
                                              "kknn"),
                                    #If we do feature selection, these are the vars we use
                                    fsfile=NULL,
                                    vars=NULL,
                                    #The control genes to use
                                    controls="ghosh",
                                    #We don't want these predictors to be used for learning
                                    filter=c("DSI","DPI","ESTcount","constitutiveexons"),
                                    seed=12345,
                                    #Number of simple models to create within the ensemble
                                    nboot=200,
                                    #Number of different values to test for each ML algorithm
                                    #hyperparameter
                                    tuneLength=5,
                                    #Number of repeats to be done by Caret
                                    nsamps=5,
                                    evalStep = 1,
                                    maxTrials=10,
                                    minKappa = 0.1,
                                    minPPV = 0.6,
                                    cutoff=0.9,
                                    alpha=0.7,
                                    #Where to save stuff
                                    out.folder="tmp/",
                                    debug=F,
                                    usereplacement=F,
                                    qmeasure="sens",
                                    useppv=F,
                                    evaldisease=NULL,
                                    evalctrl=NULL,
                                    ...)
{
  cat("Calling ensembleLearn with\n")
  genes = na.omit(fromSymbol2Hugo(genes))
  cat("Genes(",length(genes),"):",paste0(genes[1:3],collapse=","),"...\n")
  cat("We'll use",nboot,"sub-models\n")
  cat("Proportion between disease and non disease\n")
  print(table(condition))
  cat("Method list:",paste0(methods,collapse=","),"\n")
  cat("Controls:",controls,"\n")
  cat("expID:",expID,"\n")

  methods = unlist(lapply(methods,
                          function(x){ tryCatch({ checkInstall(getModelInfo(x)$library); return(x)},
                                                error = function(e){ print(e); NULL },0)}))

  cat(paste0("Available methods to use are ",paste0(methods,collapse=", ")),"\n")

  if(is.null(vars)){
    if(!is.null(fsfile)){
      vars = fsGetVars(file=fsfile)
      cat("We will use feature selection on vars\n")
      print(vars)
    }
  }

  modelEvals = NULL
  modelHits = NULL

  set.seed(seed)
  evaluation = NULL
  prefix = NULL
  ensemble = NULL

  alldata.in = fromGenes2MLData(genes=genes,
                                addcontrols=F,
                                vars=vars,
                                condition=condition,
                                filter=filter)

  if(debug){
    alldata.in = alldata.in[,c(1:10,ncol(alldata.in))]
  }

  set.seed(seed)
  evaluation = NULL
  prefix = NULL
  ensemble = NULL

  #Separate controls from disease
  controlsd = alldata.in[alldata.in$condition == "Nondisease",]
  data.in = alldata.in[alldata.in$condition == "Disease",]
  ncontrols = nrow(data.in)
  allresults = NULL
  fits = NULL

  isMinimumQuality = function(qmeasure,value){
    if(qmeasure == "costerror"){
      return(value < 0.9)
    }
    if(qmeasure == "kappa" | qmeasure == "mccc")
      return(value > 0.1)
    if(qmeasure == "bacc" | qmeasure == "sens")
      return(value > 0.3)
  }

  isOptimalValue = function(qmeasure,value){
    if(qmeasure == "costerror"){
      return(value == 0)
    }
    return(value == 1)
  }

  isProgressMade = function(qmeasure,value,lastValue){
    if(qmeasure == "costerror"){
      return(value <= lastValue)
    }
    return(value >= lastValue)
  }

  i = 1
  if(qmeasure == "costerror"){
    lastKappa = 1
    minKappa=0.9
  }
  else
    lastKappa = -1
  nulltrials = 0
  maxTrials = maxTrials*length(methods)
  methodIndex = 1
  optInfo = list()
  nonOptimal=T
  while(i <= nboot & nulltrials < maxTrials & nonOptimal){
    cat("Starting with bootstrap",i,"\n")
    cat("We will select controls sampling from the whole set\n")
    ctrlmask = sample(1:nrow(controlsd),nrow(data.in))
    if(usereplacement)
      localdata.in = rbind(data.in[sample(1:nrow(data.in),
                                          nrow(data.in),
                                          replace=T),],
                           controlsd[ctrlmask,])
    else{
      then = floor(nrow(data.in)*0.9)
      ctrlmask = sample(1:nrow(controlsd),then)
      localdata.in = rbind(data.in[sample(1:nrow(data.in),then),],
                           controlsd[ctrlmask,])

    }


    method = methods[methodIndex]
    cat("Learning with method",method,"\n")
    result = caretLearn(tuneLength=tuneLength,
                        nsamps=nsamps,
                        in.file=localdata.in,
                        model.with.all = T,
                        method=method)


    ctrlgenes = localdata.in$gene[localdata.in$condition == "Nondisease"]

    ensemble[[i]] = list(model=result$model$finalModel,
                         attsused=result$attsused,
                         method=method,
                         controls=ctrlgenes)


    toreturn = result$results[,c("ROC","Sens","Spec","ROCSD","SensSD","SpecSD")]
    toreturn = as.vector(toreturn)
    cat("Done with bootstrap",i,"\n")

    allresults = rbind(allresults,cbind(rep(i,nrow(toreturn)),toreturn))
    cat("Now let's evaluate\n")

    #Do we get an improvement???
    model.indexes = as.numeric(names(sort(tapply(allresults$ROC,allresults[,1],
                                                 function(x){ max(x)}),
                                          decreasing=T)))
    testModel = list(genes=genes,
                     method=method,
                     panel=panel,
                     controls=controls,
                     controlgenes=controlsd$gene,
                     vars=vars,
                     condition=condition,
                     nboot=length(ensemble),
                     model=ensemble,
                     eval=allresults,
                     modelindexes=model.indexes,
                     modelprefix=prefix)

    preds = ensemblePredictAllGenomev2(ensemble=testModel,
                                       n=testModel$nboot,
                                       cutoff=cutoff,
                                       vars=vars)

    testModel$preds = preds
    testModel$evaldisease = evaldisease
    testModel$evalctrl = evalctrl

    testModelEvaluation = evalEnsembleOneShot(testModel,cutoff=cutoff,alpha=alpha)

    #Should we keep evalating or not??
    #We access the kappa value
    if(useppv){
      kappa = mean(testModelEvaluation$ppv)
      cat("Using PPV\n")
      minThreshold = minPPV
    }else{
      cat("Using as kappa,",qmeasure,"\n")
      #kappa = mean(testModelEvaluation$kappa)
      kappa = mean(unlist(testModelEvaluation[qmeasure]))
      minThreshold = minKappa
    }

    print(kappa)
    if(!is.nan(kappa)){

      if(isOptimalValue(qmeasure,kappa))
        nonOptimal=F
      else{
        if(isMinimumQuality(qmeasure,kappa)){
          cat("The ensemble is good enough with the new model, should we keep it?\n")

          if(i >= evalStep){
            cat("We now test whether the model should we kept it?\n")

            if(isProgressMade(qmeasure,kappa,lastKappa)){
              cat("We improve kappa from",lastKappa,"to",kappa,"at iteration",i,"using method",methods[methodIndex],"\n")
              lastKappa = kappa
              optInfo[[length(optInfo) + 1]] = list(kappa=kappa,method=methods[methodIndex],success=T)
              nulltrials = 0
              i = i + 1
              methodIndex = (methodIndex %% length(methods)) + 1

              singleeval = evalEnsembleOneShot(ensemble=testModel,cutoff=0.5)
              modelEvals = rbind(modelEvals,c(panel,i,singleeval["auc"],singleeval["ppv"],
                                              singleeval["kappa"],singleeval["mcc"],
                                              singleeval["bacc"],singleeval["sens"],
                                              singleeval["spec"],singleeval["costerror"]))
              modelHits = rbind(modelHits,singleeval$hits)
              #print(modelEvals)
              #print(modelHits)

            }else{
              cat("We can't improve kappa from",lastKappa,"to",kappa,"at iteration",i,", doing backtrack\n")
              ensemble =  ensemble[1:(i-1)]
              optInfo[[length(optInfo) + 1]] = list(kappa=kappa,method=methods[methodIndex],success=F)
              nulltrials = nulltrials + 1
              allresults = allresults[allresults[,1] != i,]
              cat("Still",maxTrials - nulltrials,"left\n")
              methodIndex = (methodIndex %% length(methods)) + 1

              cat("We'll try now with method",methods[methodIndex],"\n")
            }
            if(useppv)
              cat("PPVs so far",paste0(unlist(lapply(optInfo,function(x){return(x["kappa"])})),collapse=","),"\n")
            else

              cat("Kappas so far",paste0(unlist(lapply(optInfo,function(x){return(x["kappa"])})),collapse=","),"\n")

          }else{ #We simply add the model and wait until we have enough of them
            cat("Yes, we keep it. Not enough models in the ensemble yet\n")
            optInfo[[length(optInfo) + 1]] = list(kappa=kappa,method=methods[methodIndex],success=T)
            i = i + 1
            methodIndex = (methodIndex %% length(methods)) + 1
            lastKappa = kappa
          }
        }else{ #At this point, the whole ensemble is bad, we drop the last model
          cat("We have to drop the last model, no minimum quality achieved\n")
          ensemble[[i]] =  NULL
          allresults = allresults[allresults[,1] != i,]
          methodIndex = (methodIndex %% length(methods)) + 1

          cat("We'll try now with method",methods[methodIndex],"\n")
        }
      }



    }else{
      cat("This iteration is NULL, kappa is not valid\n")
      nulltrials = nulltrials + 1
      allresults = allresults[allresults[,1] != i,]
      cat("Still",maxTrials - nulltrials,"left\n")
      methodIndex = (methodIndex %% length(methods)) + 1
      cat("We'll try now with method",methods[methodIndex],"\n")
    }

  }
  cat("Leaving the loop and finishing up\n")
  genes = c(data.in$gene,controlsd$gene)

  finalModel = list(genes=genes,
                    panel=panel,
                    controls=controls,
                    controlgenes=controlsd$gene,
                    vars=vars,
                    condition=condition,
                    nboot=length(ensemble),
                    model=ensemble,
                    eval=allresults,
                    modelindexes=1:length(ensemble),
                    modelprefix=prefix,
                    optinfo=optInfo)

  preds = ensemblePredictAllGenomev2(ensemble=finalModel,
                                     n=finalModel$nboot,
                                     cutoff=cutoff,
                                     vars=vars)
  finalModel$preds = preds
  finalModel$cutoff = cutoff

  #finalModel$metadata = mllearn.genEnsembleMetaData(finalModel,panel)

  finalModel$method = "multimethod"
  finalModel$evaldisease = evaldisease
  finalModel$evalctrl = evalctrl
  finalModel$modelEvals = modelEvals
  finalModel$modelHits = modelHits
  return(finalModel)
}


#' Title
#'
#' @param genes
#' @param panel
#' @param condition
#' @param expID
#' @param method
#' @param fsfile
#' @param controls
#' @param filter
#' @param seed
#' @param nboot
#' @param tuneLength
#' @param nsamps
#' @param out.folder
#' @param debug
#' @param usereplacement
#' @param useSeedGene
#' @param useSmote
#' @param evalEach
#' @param experiment
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ensembleLearn  = function(genes,
                          panel=NULL,
                          #If condition is NULL then all genes are disease
                          condition=NULL,
                          #The ID for referring to the results
                          expID = "BootstrapEnsemble",
                          #The caret Method to use
                          method="J48",
                          #If we do feature selection, these are the vars we use
                          fsfile=NULL,
                          vars=NULL,
                          #The control genes to use
                          controls="ghosh",
                          #We don't want these predictors to be used for learning
                          filter=c("DSI","DPI","ESTcount","constitutiveexons"),
                          seed=12345,
                          #Number of simple models to create within the ensemble
                          nboot=200,
                          #Number of different values to test for each ML algorithm
                          #hyperparameter
                          tuneLength=5,
                          #Number of repeats to be done by Caret
                          nsamps=5,
                          evalctrl,
                          evaldisease,
                          debug=F,
                          usereplacement=T,
                          evalEach=seq(from=1,to=nboot,by=3),
                          experiment="DRN",
                          ...)
{
  cat("Calling ensembleLearn with\n")
  cat("Genes(",length(genes),"):",paste0(genes[1:3],collapse=","),"...\n")
  cat("We'll use",nboot,"sub-models\n")
  cat("Proportion between disease and non disease\n")
  print(table(condition))
  cat("Method:",method,"\n")
  cat("Controls:",controls,"\n")
  cat("expID:",expID,"\n")

  if(is.null(vars)){
    if(!is.null(fsfile)){
      vars = fsGetVars(file=fsfile)
      cat("We will use feature selection on vars\n")
      print(vars)
    }
  }


  set.seed(seed)
  evaluation = NULL
  prefix = NULL
  ensemble = NULL

  alldata.in = fromGenes2MLData(genes=genes,
                                addcontrols=F,
                                vars=vars,
                                condition=condition,
                                filter=filter,...)

  if(debug){
    alldata.in = alldata.in[,c(1:10,ncol(alldata.in))]
  }

  #Separate controls from disease
  controlsd = alldata.in[alldata.in$condition == "Nondisease",]
  data.in = alldata.in[alldata.in$condition == "Disease",]
  ncontrols = nrow(data.in)
  allresults = NULL
  fits = NULL

  modelEvals = NULL
  modelHits = NULL
  for(i in 1:nboot){
    cat("Starting with bootstrap",i,"\n")

    cat("We will select controls sampling from the whole set\n")
    ctrlmask = sample(1:nrow(controlsd),nrow(data.in))
    if(usereplacement)
      localdata.in = rbind(data.in[sample(1:nrow(data.in),nrow(data.in),replace=T),],
                           controlsd[ctrlmask,])
    else
      localdata.in = rbind(data.in,controlsd[ctrlmask,])
    result = caretLearn(tuneLength=tuneLength,
                        nsamps=nsamps,
                        in.file=localdata.in,
                        model.with.all = T,
                        method=method)

    ctrlgenes = localdata.in$gene[localdata.in$condition == "Nondisease"]
    ensemble[[i]] = list(model=result$model$finalModel,
                         attsused=result$attsused,
                         method=method,
                         controls=ctrlgenes)


    toreturn = as.vector(result$results)
    cat("Done with bootstrap",i,"\n")
    allresults = rbind(allresults,cbind(rep(i,nrow(toreturn)),toreturn))

    if(i %in% evalEach){

      model.indexes = as.numeric(names(sort(tapply(allresults$ROC,allresults[,1],
                                                   function(x){ max(x)}),
                                            decreasing=T)))
      #print("Temporal model indexes order")
      #print(model.indexes)
      #if(controls == "clustering"){
      genes = c(data.in$gene,controlsd$gene)

      tmpModel = list(genes=genes,
                      panel=panel,
                      method=method,
                      controls=controls,
                      controlgenes=controlsd$gene,
                      vars=vars,
                      condition=condition,
                      nboot=i,
                      model=ensemble,
                      eval=allresults,
                      modelindexes=model.indexes,
                      modelprefix=prefix)
      #if(genPredictions){
      preds = ensemblePredictAllGenomev2(ensemble=tmpModel,
                                         n=i,
                                         cutoff=0.5,
                                         vars=vars,
                                         method=method)

      tmpModel$preds = preds
      #We'll get now the genes leaved out
      tmpModel$evaldisease = evaldisease
      tmpModel$evalctrl = evalctrl

      singleeval = evalEnsembleOneShot(ensemble=tmpModel,cutoff=0.5)
      modelEvals = rbind(modelEvals,c(panel,i,singleeval["auc"],singleeval["ppv"],
                                      singleeval["kappa"],singleeval["mcc"],
                                      singleeval["bacc"],singleeval["sens"],
                                      singleeval["spec"],singleeval["costerror"]))
      modelHits = rbind(modelHits,singleeval$hits)
      #print(modelEvals)
      #print(modelHits)
    }

  }
  colnames(modelEvals) = c("panel","boot","auc","ppv","kappa","mcc",
                           "bacc","sens","spec","costerror")

  #write.table(modelEvals,paste0(out.model,"_evals.tsv"),col.names=T,row.names=F,sep="\t",quote=F)
  colnames(modelHits) = c("panel","method","q","models","predictedgenes","evalgenes","hits","enrichment")
  #write.table(modelHits,paste0(out.model,"_hits.tsv"),col.names=T,row.names=F,sep="\t",quote=F)
  colnames(allresults)[1] = "bootstrap"
  model.indexes = as.numeric(names(sort(tapply(allresults$ROC,allresults$bootstrap,
                                               function(x){ max(x)}),
                                        decreasing=T)))
  #print("Model indexes order")
  #print(model.indexes)
  #if(controls == "clustering"){
  genes = c(data.in$gene,controlsd$gene)
  finalModel = list(genes=genes,
                    panel=panel,
                    method=method,
                    controls=controls,
                    controlgenes=controlsd$gene,
                    vars=vars,
                    condition=condition,
                    nboot=nboot,
                    model=ensemble,
                    eval=allresults,
                    modelindexes=model.indexes,
                    modelprefix=prefix)
  #if(genPredictions){
  preds = ensemblePredictAllGenomev2(ensemble=finalModel,
                                     n=nboot,
                                     cutoff=0.5,
                                     vars=vars,
                                     method=method)

  finalModel$preds = preds
  finalModel$modelHits = modelHits
  finalModel$allresults = allresults
  finalModel$evalctrl = evalctrl
  finalModel$evaldisease = evaldisease



  return(finalModel)
}

caretLearn = function(in.file,
                      genes=NULL,
                      method="Rborist",
                      out.file=NULL,
                      tuneLength=10,
                      nsamps=10,
                      trainingProportion=0.8,
                      tuneGrid=NULL,
                      model.with.all=F)
{
  library(caret)

  if(is.null(genes)){
    if(typeof(in.file) == "character"){
      if(length(grep(".rds$",in.file)) > 0)
        data.in = readRDS(in.file)
      else
        data.in = read.delim(in.file,sep=",")
    }else
      data.in = in.file
  }else{
    data.in = fromGenes2MLData(genes)
  }


  data.in = data.in[,colnames(data.in) != "gene"]
  numcolumns = which(colnames(data.in) != "condition")
  mask = which(startsWith(colnames(data.in),"RankedMM"))
  mask = c(mask,which(startsWith(colnames(data.in),"ExprSpecific")))
  mask = c(mask,which(startsWith(colnames(data.in),"AdjSpecific")))


  for(i in mask){
    data.in[,i] = as.factor(data.in[,i])
  }

  #Removing factors with less than 2 levels
  mask = NULL
  for(i in 1:ncol(data.in)){
    if(is.factor(data.in[,i])){
      if(length(levels(data.in[,i])) >= 2){
        mask = c(mask,TRUE)
      }else
        mask = c(mask,FALSE)
    }else if(typeof(data.in[,i]) == "character" & i != which(colnames(data.in) == "condition")){
      if(length(unique(data.in[,i])) <= 1 | sum(data.in[,i] == "1") < 5 | sum(data.in[,i] == "0") < 5)
        mask = c(mask,FALSE)
      else mask = c(mask,TRUE)
    } else
      mask = c(mask,TRUE)
  }

  if(sum(!mask) > 0)
    cat("After preparing the data for learning, we quick out",sum(!mask)," atributes\n")

  #Keep only those which are informative
  data.in = data.in[,mask]
  cat("After preparing the data for learning, we still have",ncol(data.in)," atributes\n")


  #Do a repeated cross-validation evaluation with class probabilities and two class summaries
  ctrl = trainControl(method="repeatedcv",
                      repeats=nsamps,
                      classProbs=T,
                      summaryFunction=twoClassSummary)

  #model.with.all is going to be true if we call this function from a bootstrapping model
  #learning function
  if(model.with.all){


    cat("***********Calling Caret train method with",nrow(data.in),
        "samples and",ncol(data.in),"attributes (all data)\n")

    fit = train(condition ~.,
                data = data.in,
                method=method,
                trControl=ctrl,
                metric="ROC",
                tuneLength=tuneLength)


    plsClasses = predict(fit, newdata = data.in)
    cfm = confusionMatrix(data = plsClasses,
                          reference=as.factor(data.in$condition),
                          positive="Disease")

    data.in$condition = NULL
    cat("Finished, returning the model\n")
    return(list(model=fit,results=fit$results,cfm=cfm,attsused=colnames(data.in)))
  }


  cat("**************************************Calling Caret train method with",nrow(data.in),
      "samples and",ncol(data.in),"attributes\n")
  #We should be here when we call this function in the conventional way
  inTrain = createDataPartition(y=data.in$condition,p=trainingProportion,list=F)
  training = data.in[inTrain,]
  testing = data.in[-inTrain,]

  if(!is.null(tuneGrid))
    fit = train(condition ~ .,data=training,
                method=method,
                trControl=ctrl,
                tuneGrid=tuneGrid,
                metric="ROC")
  else
    fit = train(condition ~ .,data=training,
                method=method,
                trControl=ctrl,
                tuneLength=tuneLength,
                metric="ROC")

  #Evaluate with testing data
  plsClasses = predict(fit, newdata = testing)
  cfm = confusionMatrix(data = plsClasses,
                        reference=testing$condition,
                        positive="Disease")
  return(list(model=fit,results=fit$results,cfm=cfm))
}

ensemblePredictAllGenomev2 = function(ensemble,
                                      n = 200,
                                      cutoff=0.9,
                                      vars=NULL,
                                      method=NULL){
  cat("Predicting All Genome v2\n")
  allgenes = getCodingGenome()
  return(ensemblePredict(genes=allgenes,
                         ensemble=ensemble$model,
                         n=n,
                         vars=vars,
                         cutoff=cutoff,
                         model.indexes=ensemble$modelindexes[1:n],
                         method=method))
}



ensemblePredict = function(genes,
                           ensemble,
                           n = 200,
                           cutoff=0.9,
                           silent=T,
                           vars=NULL,
                           model.indexes=NULL,
                           method=NULL){

  alldata.in = fromGenes2MLData(genes=genes,addcontrols=F,vars=vars)
  genes = alldata.in$gene
  alldata.in$gene = NULL
  alldata.in$condition = NULL

  preds = NULL
  if(is.null(model.indexes))
    model.indexes = 1:n
  for(i in model.indexes){
    data.in = alldata.in
    model = ensemble[[i]]$model
    method = ensemble[[i]]$method

    attsused = ensemble[[i]]$attsused
    data.in = data.in[,attsused]
    cat("Predicting with simple model",i,"and method",method,"\n")

    mask = c(grep("ExprSpecific",colnames(data.in)),
             grep("AdjSpecific",colnames(data.in)),
             grep("RankedMMSpecific",colnames(data.in)))
    colnames(data.in)[mask] = paste0(colnames(data.in[mask]),"1")

    if(nrow(data.in) == 0)
      return(NULL)

    #if(!silent)
    #print(model$finalModel)
    if(method %in% c("xgbTree","xgbLinear","xgbDART")){
      #print(colnames(as.matrix(data.in)))
      #print(model$xNames)
      numpred = predict(model,newdata=as.matrix(data.in),type="response")
      localpred = vector(mode="character",length=length(numpred))
      localpred[numpred > 0.5] = "Disease"
      localpred[numpred <= 0.5] = "Nondisease"
    }else if(method %in% c("C5.0Tree","rpart","C5.0","naive_bayes","rf","J48")){
      localpred = as.character(predict(model,newdata=data.in,type="class"))
    }else if (method %in% c("kknn","svmRadial")){
      library(kernlab)
      localpred = as.character(predict(model,newdata=data.in))
    }else if (method %in% c("lda","sparseLDA")){
      localpred = as.character(predict(model,newdata=data.in)$class)
    }else{
      numpred = predict(model,newdata=data.in,type="response")
      localpred = vector(mode="character",length=length(numpred))
      localpred[numpred <= 0.5] = "Disease"
      localpred[numpred > 0.5] = "Nondisease"
    }
    #print(table(localpred))

    preds$preds[[i]] = localpred
  }

  preds$genes = genes
  preds$predictions = g2pmlpredict(preds,ratio=cutoff)
  #print(table(preds$predictions$prediction))

  return(list(genes=preds$genes,
              unknowngenes=genes[!(genes %in% preds$genes)],
              prediction=preds$predictions$prediction,
              quality=preds$predictions$quality,
              rawpredictions=preds$preds))
}

g2pmlpredict = function(preds,ratio=0.5,useDontKnow=F){

  if(is.null(preds))
    return(preds)

  genes = preds$genes
  ngenes = length(preds$genes)
  disease = vector(mode="numeric",length=ngenes)
  disease[] = 0
  nondisease = disease

  preds = preds$preds
  npreds = length(preds)
  for(i in 1:npreds){
    mask = preds[[i]] == "Disease"
    disease[mask] = disease[mask] + 1
    nondisease[!mask] = nondisease[!mask] + 1
  }
  disease = disease/npreds
  nondisease = nondisease/npreds

  if(useDontKnow){
    outlabel = rep("DontKnow",ngenes)
    outlabel[disease >= ratio] = "Disease"
    outlabel[nondisease >= ratio] = "Nondisease"
  }else{
    outlabel = rep("Nondisease",ngenes)
    outlabel[disease >= ratio] = "Disease"
  }


  confidence = disease
  #confidence[nondisease >= ratio] = nondisease[nondisease >= ratio]
  names(confidence) = outlabel

  return(list(prediction=outlabel,quality=confidence,genes=genes))
}

evalEnsemblesOneShot = function(ensembles,
                                useppv=T,
                                cutoff=0.7,
                                qmeasure="mcc",
                                doGO=F){

  k = length(ensembles)
  globalallhits = NULL
  globalhitsperfold = NULL

  allpredictions = NULL
  alldiseasegenes = NULL
  allhits = NULL
  hitsperfold = NULL

  results = NULL
  kappa = 0
  modelcount = 0
  used = 0
  evalgenes = NULL
  for(i in 1:k){
    untilTree = ensembles[[i]]$nboot
    #for(untilTree in 1:ensembles[[1]]$nboot){
    print(untilTree)

    cat("First we evaluate without predictions\n")
    evaluation = evalEnsembleOneShot(ensemble = ensembles[[i]],
                                     cutoff=cutoff,
                                     trees=untilTree)
    ensembles[[i]]$evaluation = evaluation
  }

  for(q in c(0.5,0.8,0.9,0.99)){
    diseasegenes = NULL
    for(i in 1:k){
      cat("Now we generate predictions\n")
      if(useppv){
        kappa = kappa + ensembles[[i]]$evaluation$ppv
      }else{
        kappa = kappa + unlist(ensembles[[i]]$evaluation[qmeasure])
      }
      milked = milkEnsemble(ensemble = ensembles[[i]],
                            cutoff=q,
                            trees=untilTree,
                            remove=c(ensembles[[i]]$genes[ensembles[[i]]$condition == "Disease"])) #,ensembles[[i]]$evaldisease))
      hits = getHits(panel=ensembles[[i]]$panel,
                     genes=milked,
                     brandnew=ensembles[[i]]$evaldisease)

      hitsperfold = rbind(hitsperfold,c(ensembles[[i]]$panel,
                                        ensembles[[i]]$method,
                                        i,q,untilTree,
                                        length(milked),
                                        hits$newgenes,
                                        hits$hits,
                                        hits$fold))
      diseasegenes = c(diseasegenes, milked)

    }

    if(length(diseasegenes) > 0){
      diseasegenes = table(diseasegenes)
      #diseasegenes = as.numeric(diseasegenes/used)
      diseasegenes = sort(diseasegenes,decreasing=T)
      genes = names(diseasegenes)
      diseasegenes = as.vector(diseasegenes)
      diseasegenes = diseasegenes/k

      finaltable = as.data.frame(cbind(genes=genes,quality=diseasegenes),stringsAsFactors=F)
      colnames(finaltable) = c("gene","quality")
      allpredictions = rbind(allpredictions,cbind(rep(ensembles[[i]]$panel,nrow(finaltable)),finaltable))
      for(qkfold in c(0,0.3,0.5)){
        localfinaltable = finaltable
        localfinaltable = localfinaltable[as.numeric(localfinaltable[,2]) > qkfold,]

        #Now the hits
        hits = getHits(panel=ensembles[[i]]$panel,genes=localfinaltable[,1])
        if(!is.null(hits)){
          allhits = rbind(allhits,c(ensembles[[i]]$panel,
                                    ensembles[[i]]$method,
                                    q,
                                    qkfold,
                                    nrow(localfinaltable),
                                    hits$newgenes,
                                    hits$hits,
                                    hits$fold))
          #print(allhits)
        }

        cat("Trying gProfilerR query for gene quality",q,"and fold quality",qkfold,
            "and",nrow(localfinaltable),"genes\n")
        if(nrow(localfinaltable) > 0 & nrow(localfinaltable) < 4000 & doGO){
          cat("Doing gProfilerR query for gene quality",q,"and fold quality",qkfold,
              ",yields",nrow(localfinaltable),"genes\n")
          res = gProfileR::gprofiler(query=localfinaltable[,1],
                                     src_filter=c("GO","KEGG","REAC","HP","OMIM"))
          #res=matrix(nrow=0,ncol=0)
          cat("We got",nrow(res),"enrichment terms\n")
          if(nrow(res) > 0){
            res$intersection = NULL
            write.table(res,paste0(out.path,"/",localsaveprefix,"_predsGO_",q,"_",qkfold,".tsv"),
                        row.names=F,quote=F,sep="\t")

          }

        }
      }
    }else{
      cat("No genes were predicted at any cutoff\n")
    }

    if(!is.null(allhits)){
      colnames(allhits) = c("panel","method","ensembleq","kfoldq","predictions","newgenes","hits","enrichment")
      globalallhits = rbind(globalallhits,allhits)
    }
    if(!is.null(hitsperfold)){
      c(i,q,untilTree,length(milked),hits$newgenes,hits$hits,hits$fold)
      colnames(hitsperfold) = c("panel","method","fold","q","trees","predictions","evalgenes","hits","enrichment")
      globalhitsperfold = rbind(globalhitsperfold,hitsperfold)
    }

  }
  return(list(hits=globalallhits,hitsperfold=globalhitsperfold))
}

evalEnsembleOneShot = function(ensemble,
                               cutoff=0.9,
                               balance=50,
                               alpha=0.9,
                               trees=-1){

  if(length(trees) == 1){
    if(trees < 0 | trees > length(ensemble$model))
      trees = length(ensemble$model)
  }

  if(trees < length(ensemble$model)){
    allgenes = getCodingGenome()
    ensemble$preds = ensemblePredict(genes=allgenes,
                                     ensemble=ensemble$model,
                                     n=n,
                                     vars=ensemble$vars,
                                     cutoff=cutoff,
                                     model.indexes=ensemble$modelindexes[1:trees],
                                     method=ensemble$method)
  }

  #Let us eval the disease genes
  gindexes = which(ensemble$preds$genes %in% ensemble$evaldisease)
  allcgindexes = which(ensemble$preds$genes %in% ensemble$evalctrl)
  preds = ensemble$preds$rawpredictions

  alloutpreds = NULL
  allstats = NULL
  cat("Trying with",trees,"trees\n")
  outpreds = NULL
  localstats = NULL
  done = 0
  for(b in 1:balance){
    cgindexes = sample(allcgindexes,length(gindexes))
    gtruth = c(rep("Disease",
                   length(gindexes)),
               rep("Nondisease",
                   length(cgindexes)))

    allindexestoeval = c(gindexes,cgindexes)
    numpreds = rep("Disease",length(allindexestoeval))
    numpreds[ensemble$preds$quality[allindexestoeval] < cutoff] = "Nondisease"
    itstats = stats(numpreds,gtruth,alpha)
    #if(sum(is.nan(unlist(itstats))) == 0){
    localstats = rbind(localstats,c(itstats$ppv,itstats$npv,
                                    itstats$kappa,itstats$mcc,
                                    itstats$bacc,itstats$f1,itstats$auc,
                                    itstats$sens,
                                    itstats$spec,itstats$costerror))
    # done = done + 1
    #}
  }
  hitsperfold=NULL
  for(q in c(0.5,0.8,0.9,0.99)){
    diseasegenes = NULL
    milked = milkEnsemble(ensemble = ensemble,
                          cutoff=q,
                          trees=trees,
                          remove=c(ensemble$genes[ensemble$condition == "Disease"])) #,ensemble$evaldisease))

    hits = getHits(panel=ensemble$panel,
                   genes=milked,
                   brandnew=ensemble$evaldisease)

    hitsperfold = rbind(hitsperfold,c(ensemble$panel,
                                      ensemble$method,
                                      q,
                                      trees,
                                      length(milked),
                                      hits$newgenes,
                                      hits$hits,
                                      hits$fold))


  }


  #print(localstats)
  localstats = na.omit(localstats)
  return(list(auc=mean(localstats[,7]),
              ppv=mean(localstats[,1]),
              kappa=mean(localstats[,3]),
              mcc=mean(localstats[,4]),
              bacc=mean(localstats[,5]),
              sens=mean(localstats[,8]),
              spec=mean(localstats[,9]),
              costerror=mean(localstats[,10]),hits=hitsperfold))
}

milkEnsemble = function(ensemble,
                        cutoff=0.9,
                        remove=NULL,
                        trees=-1){
  if(!is.null(ensemble$model)){
    if(trees < 0 | trees > length(ensemble$model))
      trees = length(ensemble$model)

    if(trees < length(ensemble$model)){
      allgenes = getCodingGenome()
      ensemble$preds = ensemblePredict(genes=allgenes,
                                       ensemble=ensemble$model,
                                       n=n,
                                       vars=ensemble$vars,
                                       cutoff=cutoff,
                                       model.indexes=ensemble$modelindexes[1:trees],
                                       method=ensemble$method)
    }

  }
  #Let us eval the disease genes
  if(!is.null(remove)){
    genemask = !(ensemble$preds$genes %in% remove)
    ensemble$preds$genes = ensemble$preds$genes[genemask]
    ensemble$preds$quality = ensemble$preds$quality[genemask]
  }
  return(ensemble$preds$genes[ensemble$preds$quality >= cutoff])
}

getHits = function(panel,genes,
                   newgenespath=paste0(system.file("g2pml/gel/", "", package = "G2PML"),
                                       "geOctober2018/"),
                   oldgenespath=paste0(system.file("g2pml/gel/", "", package = "G2PML"),
                                       "geApril2018/"),
                   evidence = c("HighEvidence","ModerateEvidence"),
                   brandnew=NULL,
                   backgroundgenes=18000)
{
  cat("Getting hits for panel",panel,"\n")
  if(is.null(brandnew)){
    newgenes = read.csv(paste0(newgenespath,"/",panel,".csv"),stringsAsFactors=F)
    newgenes = newgenes$GeneSymbol[newgenes$LevelOfConfidence %in% evidence]
    oldgenes = read.csv(paste0(oldgenespath,"/",panel,".csv"),stringsAsFactors=F)
    oldgenes = oldgenes$GeneSymbol[oldgenes$LevelOfConfidence %in% "HighEvidence"]
    brandnew = setdiff(newgenes,oldgenes)
    brandnew = na.omit(fromSymbol2Hugo(brandnew))
    print("new genes from newer panel are")
    print(brandnew)
  }


  if(length(brandnew) > 0){
    bgrndprob = length(genes)/backgroundgenes
    localhits = sum(brandnew %in% genes)
    cat("We got",localhits,"hits\n")
    fold = (localhits/length(brandnew))/bgrndprob
    return(list(hits=localhits,fold=fold,newgenes=length(brandnew),predictions=length(genes)))
  }else
    cat("No new genes to get hits from\n")
  return(NULL)
}


stats = function(preds,gtruth,alpha=0.5){
  library(pROC)
  tp = sum(preds == "Disease" & gtruth == "Disease")
  fp = sum(preds == "Disease" & gtruth == "Nondisease")
  tn = sum(preds == "Nondisease" & gtruth == "Nondisease")
  fn = sum(preds == "Nondisease" & gtruth == "Disease")
  pos = sum(gtruth == "Disease")
  posuseful = sum(gtruth == "Disease" & preds != "DontKnow")
  neg = sum(gtruth == "Nondisease")
  neguseful = sum(gtruth == "Nondisease" & preds != "DontKnow")

  dnratio = sum(preds == "DontKnow")/length(preds)
  diseaseratio = sum(preds == "Disease")/length(preds)
  ndiseaseratio = sum(preds == "Nondisease")/length(preds)
  derror = fp/(tp + fp)

  fprate = fp/(tp+fp)
  fnrate = fn/(tn+fn)
  costerror = (alpha*fprate + (1-alpha)*fnrate)/(neg + pos)
  if(is.nan(costerror))
    costerror = 0

  sensitivity = tp/(tp + fn)
  specificity = tn/(fp + tn)
  ppv = tp/(tp + fp)
  npv = tn/(tn + fn)
  if(is.nan(npv))
    npv = 0
  acc = (tp + tn)/sum(preds == "Disease" | preds == "Nondisease")
  bacc = (tp/pos + tn/neg)/2
  f1 = 2*tp/(2*tp + fp + fn)
  mcc = (tp*tn - fp*fn)/sqrt(as.numeric(tp+fp)*as.numeric(tp+fn)*as.numeric(tn+fp)*as.numeric(tn+fn))
  if(is.nan(mcc))
    mcc = -1
  #cat("mcc:",mcc,"\n")
  randomacc = ((tn+fp)*(tn+fn) + (fn+tp)*(fp+tp))/(tp+tn+fp+fn)^2
  kappa = (((tp+tn)/(tp+tn+fp+fn)) - randomacc) / (1 - randomacc)
  #if(length(unique(preds)) == 2)
  #  auc = auc(preds,gtruth,plotit=F)
  #else
  auc = 0.4
  return(list(auc=auc,sens=sensitivity,spec=specificity,
              ppv=ppv,npv=npv,acc=acc,bacc=bacc,f1=f1,
              mcc=mcc,kappa=kappa,
              diseaseratio=diseaseratio,
              ndiseaseratio=ndiseaseratio,
              dnratio=dnratio,
              costerror=costerror,
              derror=derror))
}

