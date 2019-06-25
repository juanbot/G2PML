

#' Title
#'
#' @param atts
#'
#' @return
#' @export
#'
#' @examples
prettyPredictorName = function(atts){
  atts = gsub("ExprSpecific","Expression",atts)
  atts = gsub("RankedMMSpecificRank","",atts)
  atts = gsub("AdjSpecificAdj","Adjacency",atts)
  atts = gsub("1$","",atts)
  return(atts)

}


#' Title
#'
#' @param out.path
#' @param outfile
#'
#' @return
#' @export
#'
#' @examples
generateMLDB = function(handlers=dataHandlerMethods()){

  allgenes = read.table(paste0(system.file("g2pml/", "", package = "G2PML"),
                               "workinggenes.txt"),
                        stringsAsFactors=F,header=F)

  allgenes = fromSymbol2Hugo(allgenes$V1)
  allgenes = na.omit(allgenes)
  dataall = generateMLDBByGene(casecontrolset=getSingleConditionSet(genes=allgenes),
                               handlers = handlers)
  dataft = dataall[,!(colnames(dataall) %in%
                        c("alt3.5EST","ExACpLI","ExACpRec","ExACpNull","ExACpMiss","DPI",
                          "DSI","ESTcount","constitutiveexons"))]
  return(dataft)
}

generateMLDBByGene = function(casecontrolset=getCaseControlSet(which.ones="ge_neurogenes"),
                              handlers=dataHandlerMethods(),
                              filter=NULL){


  #We keep only those that have an ensemble ID
  hugo.genes = fromSymbol2Hugo(casecontrolset$genes)
  mask = !is.na(hugo.genes)
  casecontrolset = casecontrolset[mask,]
  cat("When generating Specificity Data Set keeping",length(mask),"genes out of",length(hugo.genes),
      "because of no translation to HGCN\n")

  hugo.genes = hugo.genes[mask]

  valgenes = casecontrolset$genes
  #myens <<- fromGeneName2EnsemblBM(valgenes)
  names(valgenes) = myens
  #names(valgenes) = fromGeneName2EnsemblBM(valgenes)

  mldata = getAllGeneValues(genes=valgenes,
                            queries = handlers)

  mldata = as.data.frame(cbind(mldata,casecontrolset$condition))
  colnames(mldata)[ncol(mldata)] = "condition"


  #Deal with NAs how?
  #String
  mldata = dealWithNA(mldata,which(colnames(mldata) == "String"),"mean")
  #ExACpLI
  mldata = dealWithNA(mldata,which(colnames(mldata) == "ExACpLI"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "ExACpRec"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "ExACpNull"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "ExACpMiss"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "gnomadpLI"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "gnomadpRec"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "gnomadpNull"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "gnomadpMiss"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "gnomadOELoF"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "gnomadOEMiss"),"mean")

  mldata = dealWithNA(mldata,which(colnames(mldata) == "LoFTool"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "RVIS"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "EvoTol"),"mean")

  mldata = dealWithNA(mldata,which(colnames(mldata) == "pAD"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "pAR"),"mean")


  #Biomart on gene complexity
  mldata = dealWithNA(mldata,which(colnames(mldata) == "GeneLength"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "TranscriptCount"),"mean")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "GCcontent"),"mean")

  mldata = dealWithNA(mldata,which(colnames(mldata) == "CountsOverlap"),"meani")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "CountsProtCodOverlap"),"meani")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "NumJunctions"),"meani")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "IntronicLength"),"meani")

  #Disgenet
  mldata = dealWithNA(mldata,which(colnames(mldata) == "DSI"),"stzero")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "DPI"),"stzero")

  #HEXEvent attributes
  mldata = dealWithNA(mldata,which(colnames(mldata) == "constitutiveexons"),"meani")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "alt3EST"),"meani")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "alt5EST"),"meani")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "alt3.5EST"),"meani")
  mldata = dealWithNA(mldata,which(colnames(mldata) == "ESTcount"),"meani")

  if(!is.null(filter)){
    cat("Removing columns for attributes",paste0(filter,collapse=" and "),"\n")
    mldata = mldata[,!(colnames(mldata) %in% filter)]
  }

  return(mldata)

}

dataHandlerMethods = function(){
  return(c(getGeneValuesGeneLength,
           getGeneValuesTranscriptCount,
           getGeneValuesCountsOverlap,
           getGeneValuesNumJunctions,
           getGeneValuesIntronicLength,
           getGeneValuesGCcontent,
           getGeneValuesString,
           getGeneValuesLoFTool,
           getGeneValuesEvoTol,
           getGeneValuesRVIS,
           getGeneValuespAD,
           getGeneValuespAR,
           getGeneValuesgnomadpLI,
           getGeneValuesgnomadpRec,
           getGeneValuesgnomadpNull,
           getGeneValuesgnomadpMiss,
           getGeneValuesgnomadoeLoF,
           getGeneValuesgnomadoeMiss,
           getGeneValuesRankedMMSpecificity,
           getGeneValuesExpressionSpecificity,
           getGeneValuesAdjacencySpecificity))
}

#' Title
#'
#' @param genes
#' @param queries
#'
#' @return
#' @export
#'
#' @examples
getAllGeneValues = function(genes,
                            queries = dataHandlerMethods()){

  toreturn = data.frame(list(gene=genes),stringsAsFactors=F)

  for(query in queries){
    partialquery = getGeneValues(genes=genes,query=query)
    if(sum(genes == partialquery$gene) == length(genes)){
      toreturn = cbind(toreturn,partialquery[,2:ncol(partialquery),drop=FALSE])

    }else{
      fquery = NULL
      mask = match(genes,partialquery$gene)
      dummyna = rep(NA,ncol(partialquery) - 1)
      for(i in 1:length(genes)){
        if(is.na(mask[i]))
          fquery = rbind(fquery,dummyna)
        else
          fquery = rbind(fquery,partialquery[mask[i],2:ncol(partialquery)])
      }
      colnames(fquery) = colnames(partialquery)[2:ncol(partialquery)]
      toreturn = cbind(toreturn,fquery)
    }

  }
  return(toreturn)
}

#' Title
#'
#' @param data.in
#' @param colindex
#' @param how
#'
#' @return
#' @export
#'
#' @examples
dealWithNA = function(data.in,colindex,how=c("stzero","mean","min")){
  if(how == "stzero"){
    for(x in colindex){
      cat("Dealing with NA in",colnames(data.in)[x],"with stzero\n")
      cat("% of NA is",100*sum(is.na(data.in[,x]))/length(data.in[,x]),"\n")
      data.in[is.na(data.in[,x]),x] = 0
    }
  }else if(how == "mean"){
    for(x in colindex){
      cat("Dealing with NA in",colnames(data.in)[x],"with mean\n")
      cat("% of NA is",100*sum(is.na(data.in[,x]))/length(data.in[,x]),"\n")
      the.mean = mean(na.omit(data.in[,x]))
      data.in[is.na(data.in[,x]),x] = the.mean
    }
  }else if(how == "meani"){
    for(x in colindex){
      cat("Dealing with NA in",colnames(data.in)[x],"with mean\n")
      cat("% of NA is",100*sum(is.na(data.in[,x]))/length(data.in[,x]),"\n")
      the.mean = as.integer(mean(na.omit(data.in[,x])))
      data.in[is.na(data.in[,x]),x] = the.mean
    }
  }else if(how == "min"){
    for(x in colindex){
      cat("Dealing with NA in",colnames(data.in)[x],"with min\n")
      cat("% of NA is",100*sum(is.na(data.in[,x]))/length(data.in[,x]),"\n")
      the.min = min(na.omit(data.in[,x]))
      data.in[is.na(data.in[,x]),x] = the.min
    }
  }else if(how == "mode"){
    for(x in colindex){
      cat("Dealing with NA in",colnames(data.in)[x],"with mode\n")
      cat("% of NA is",100*sum(is.na(data.in[,x]))/length(data.in[,x]),"\n")
      vals = unique(na.omit(data.in[,x]))
      tabs = table(na.omit(data.in[,x]))
      the.mode = as.integer(names(tabs)[which.max(tabs)])
      data.in[is.na(data.in[,x]),x] = the.mode
      #print(data.in[,2])
    }
  }
  return(data.in)
}

getGeneValues = function(genes,query,...){
  return(query(genes,...))
}

getGeneValuesGeneLength = function(genes,drop.na=F){
  cat("Generating GeneLength data for",length(genes),"genes\n")
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)
  fileout = paste0(system.file("g2pml/", "", package = "G2PML"),"/genelength.txt.zip")
  datain = read.table(unz(fileout,"genelength.txt"),stringsAsFactors=F,header=T)
  mask = match(ens.genes,datain$ensembl_gene_id)
  dataout = NULL
  j = 1
  for(i in mask){
    if(is.na(i)){
      dataout = rbind(dataout,c(genes[j],NA))
    }else
      dataout = rbind(dataout,c(genes[j],abs(datain$end_position[i]-datain$start_position[i])))
    j = j + 1
  }
  colnames(dataout) = c("gene","GeneLength")
  dataout = as.data.frame(dataout,stringsAsFactors=F)
  dataout$GeneLength = as.numeric(dataout$GeneLength)
  dataout$gene = genes
  return(dataout)

}

getGeneValuesTranscriptCount = function(genes,drop.na=F){
  cat("Generating TranscriptCount data for",length(genes),"genes\n")
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)
  fileout = paste0(system.file("g2pml/", "", package = "G2PML"),
                   "transcriptcount.txt.zip")
  datain = read.table(unz(fileout,"transcriptcount.txt"),stringsAsFactors=F,header=T)
  mask = match(ens.genes,datain$ensembl_gene_id)
  dataout = NULL

  j = 1
  for(i in mask){
    if(is.na(i)){
      dataout = rbind(dataout,c(genes[j],NA))
    }else
      dataout = rbind(dataout,c(genes[j],datain$transcript_count[i]))
    j = j + 1
  }
  colnames(dataout) = c("gene","TranscriptCount")
  dataout = as.data.frame(dataout,stringsAsFactors=F)
  dataout$TranscriptCount = as.numeric(dataout$TranscriptCount)
  dataout$gene = genes
  return(dataout)

}

getGeneValuesGCcontent = function(genes,drop.na=F){
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)
  fileout = paste0(system.file("g2pml/", "", package = "G2PML"),
                   "gccontent.txt.zip")
  datain = read.table(unz(fileout,"gccontent.txt"),stringsAsFactors=F,header=T)
  mask = match(ens.genes,datain$ensembl_gene_id)
  dataout = NULL

  j = 1
  for(i in mask){
    if(is.na(i)){
      dataout = rbind(dataout,c(genes[j],NA))
    }else
      dataout = rbind(dataout,c(genes[j],datain$percentage_gene_gc_content[i]))
    j = j + 1
  }
  colnames(dataout) = c("gene","GCcontent")
  dataout = as.data.frame(dataout,stringsAsFactors=F)
  dataout$GCcontent = as.numeric(dataout$GCcontent)
  dataout$gene = genes
  return(dataout)
}

getGeneValuesCountsOverlap = function(genes,drop.na=F){
  cat("Generating CountsOverlap data for",length(genes),"genes\n")
  filein = paste0(system.file("g2pml/", "", package = "G2PML"),"annotation.zip")
  datain = read.table(unz(filein,"annotation.txt"),stringsAsFactors=F,sep="\t",header=T)
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)
  mask = match(ens.genes,datain$gname)
  dataout = NULL

  j = 1
  for(i in mask){
    if(is.na(i)){
      dataout = rbind(dataout,c(genes[j],NA))
    }else
      dataout = rbind(dataout,c(genes[j],datain$countsOverlap[i]))
    j = j + 1
  }
  if(drop.na)
    dataout = dataout[!is.na(dataout[,2]),,drop=F]
  colnames(dataout) = c("gene","CountsOverlap")
  dataout = as.data.frame(dataout,stringsAsFactors=F)
  dataout$CountsOverlap = as.numeric(dataout$CountsOverlap)
  return(dataout)
}

getGeneValuesCountsProtCodOverlap = function(genes,drop.na=F){
  cat("Generating CountsProtCodOverlap data for",length(genes),"genes\n")
  filein = paste0(system.file("g2pml/", "", package = "G2PML"),"annotation.zip")
  datain = read.table(unz(filein,"annotation.txt"),stringsAsFactors=F,sep="\t",header=T)
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)
  mask = match(ens.genes,datain$gname)
  dataout = NULL

  j = 1
  for(i in mask){
    if(is.na(i)){
      dataout = rbind(dataout,c(genes[j],NA))
    }else
      dataout = rbind(dataout,c(genes[j],datain$countsProCod[i]))
    j = j + 1
  }
  if(drop.na)
    dataout = dataout[!is.na(dataout[,2]),,drop=F]
  colnames(dataout) = c("gene","CountsProtCodOverlap")
  dataout = as.data.frame(dataout,stringsAsFactors=F)
  dataout$CountsProtCodOverlap = as.numeric(dataout$CountsProtCodOverlap)
  return(dataout)
}

#This returns the number of unique junctions across the gene
#In principle, it is expected this is correlated with number
#of exons, not necessarily with number of transcripts
getGeneValuesNumJunctions = function(genes,drop.na=F){
  cat("Generating NumJunctions data for",length(genes),"genes\n")
  filein = paste0(system.file("g2pml/", "", package = "G2PML"),"annotation.zip")
  datain = read.table(unz(filein,"annotation.txt"),stringsAsFactors=F,sep="\t",header=T)
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)
  mask = match(ens.genes,datain$gname)
  dataout = NULL

  j = 1
  for(i in mask){
    if(is.na(i)){
      dataout = rbind(dataout,c(genes[j],NA))
    }else
      dataout = rbind(dataout,c(genes[j],datain$numJunctions[i]))
    j = j + 1
  }
  if(drop.na)
    dataout = dataout[!is.na(dataout[,2]),,drop=F]
  colnames(dataout) = c("gene","NumJunctions")
  dataout = as.data.frame(dataout,stringsAsFactors=F)
  dataout$NumJunctions = as.numeric(dataout$NumJunctions)
  return(dataout)
}

#This returns the sum of lengths of each intron within the gene
#If NA, the gene only has an exon, thus NA is equal to 0
getGeneValuesIntronicLength = function(genes,drop.na=F){
  cat("Generating IntronicLength data for",length(genes),"genes\n")
  filein = paste0(system.file("g2pml/", "", package = "G2PML"),"annotation.zip")
  datain = read.table(unz(filein,"annotation.txt"),stringsAsFactors=F,sep="\t",header=T)
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)
  mask = match(ens.genes,datain$gname)
  dataout = NULL

  j = 1
  for(i in mask){
    if(is.na(i)){
      dataout = rbind(dataout,c(genes[j],NA))
    }else
      dataout = rbind(dataout,c(genes[j],datain$intronicLength[i]))
    j = j + 1
  }
  if(drop.na)
    dataout = dataout[!is.na(dataout[,2]),,drop=F]
  colnames(dataout) = c("gene","IntronicLength")
  dataout = as.data.frame(dataout,stringsAsFactors=F)
  dataout$IntronicLength = as.numeric(dataout$IntronicLength)
  return(dataout)
}

getGeneValuesDSI = function(genes,drop.na=F){
  cat("Generating Disease specificity index data for",length(genes),"genes\n")
  disgenet = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),"disgenet_allgeneleveldata.tsv"),
                        stringsAsFactors=F,sep="\t")
  if(drop.na)
    genes = genes[!is.na(match(genes,disgenet$Symbol))]
  toreturn = cbind(genes,disgenet$DSI[match(genes,disgenet$Symbol)])
  colnames(toreturn) = c("gene","DSI")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$DSI = as.numeric(toreturn$DSI)
  return(toreturn)
}

getGeneValuesDPI = function(genes,drop.na=F){
  cat("Generating Disease pleiotropy index data for",length(genes),"genes\n")
  disgenet = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),"disgenet_allgeneleveldata.tsv"),
                        stringsAsFactors=F,sep="\t")
  if(drop.na)
    genes = genes[!is.na(match(genes,disgenet$Symbol))]
  toreturn = cbind(genes,disgenet$DPI[match(genes,disgenet$Symbol)])
  colnames(toreturn) = c("gene","DPI")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$DPI = as.numeric(toreturn$DPI)
  return(toreturn)
}

getGeneValuesLoFTool = function(genes,drop.na=F){
  cat("Generating LoFTool scores data for",length(genes),"genes\n")
  loftool = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                             "genesScoresNBarahona.txt.zip"),
                           "genesScoresNBarahona.txt"),
                      stringsAsFactors=F,sep=" ")
  pli = tapply(loftool$LoFtool,loftool$Gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","LoFTool")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$LoFTool = as.numeric(toreturn$LoFTool)
  return(toreturn)
}

getGeneValuesEvoTol = function(genes,drop.na=F){
  cat("Generating EvoTol scores data for",length(genes),"genes\n")
  evotol = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),"EvoTolWebQuery13March2019.txt.zip"),
                          "EvoTolWebQuery13March2019.txt"),
                       stringsAsFactors=F,sep="\t")
  pli = tapply(evotol$evotol,evotol$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","EvoTol")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$EvoTol = as.numeric(toreturn$EvoTol)
  return(toreturn)
}

#https://academic.oup.com/bioinformatics/article/33/4/471/2525582
getGeneValuesRVIS = function(genes,drop.na=F){
  cat("Generating RVIS scores data for",length(genes),"genes\n")
  rvis = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                              "RVIS_Unpublished_ExACv2_March2017.zip"),"RVIS_Unpublished_ExACv2_March2017.txt"),
                       stringsAsFactors=F,sep="\t")
  pli = tapply(rvis$OEratio,rvis$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","RVIS")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$RVIS = as.numeric(toreturn$RVIS)
  return(toreturn)
}

#http://genetics.bwh.harvard.edu/genescores/selection.html
getGeneValuespAD = function(genes,drop.na=F){
  cat("Generating probability of AD Model of inheritance scores data for",length(genes),"genes\n")
  pad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                           "moi.zip"),"moi.txt"),
                    stringsAsFactors=F,sep="\t")
  pli = tapply(pad$p_AD,pad$gene_symbol,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","pAD")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$pAD = as.numeric(toreturn$pAD)
  return(toreturn)
}

#http://genetics.bwh.harvard.edu/genescores/selection.html
getGeneValuespAR = function(genes,drop.na=F){
  cat("Generating probability of AR Model of inheritance scores data for",length(genes),"genes\n")
  pad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                          "moi.zip"),"moi.txt"),
                   stringsAsFactors=F,sep="\t")
  pli = tapply(pad$p_AR,pad$gene_symbol,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","pAR")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$pAR = as.numeric(toreturn$pAR)
  return(toreturn)
}


getGeneValuesgnomadpLI = function(genes,drop.na=F){
  cat("Generating gnomADpLI data for",length(genes),"genes\n")
  gnomad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                             "gnomADrelease_2.1_ht_constraint_constraint.zip"),
                          "gnomADrelease_2.1_ht_constraint_constraint.txt"),
                      stringsAsFactors=F,sep="\t")
  pli = tapply(gnomad$pLI,gnomad$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","gnomadpLI")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$gnomadpLI = as.numeric(toreturn$gnomadpLI)
  return(toreturn)
}

getGeneValuesgnomadpRec = function(genes,drop.na=F){
  cat("Generating gnomADpRec data for",length(genes),"genes\n")
  gnomad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                                 "gnomADrelease_2.1_ht_constraint_constraint.zip"),
                          "gnomADrelease_2.1_ht_constraint_constraint.txt"),
                      stringsAsFactors=F,sep="\t")
  pRec = tapply(gnomad$pRec,gnomad$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pRec)))]
  toreturn = cbind(genes,pRec[match(genes,names(pRec))])
  colnames(toreturn) = c("gene","gnomadpRec")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$gnomadpRec = as.numeric(toreturn$gnomadpRec)
  return(toreturn)
}

getGeneValuesgnomadpMiss = function(genes,drop.na=F){
  cat("Generating gnomADpMiss data for",length(genes),"genes\n")
  gnomad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                                 "gnomADrelease_2.1_ht_constraint_constraint.zip"),
                          "gnomADrelease_2.1_ht_constraint_constraint.txt"),
                      stringsAsFactors=F,sep="\t")
  pMiss = tapply(gnomad$mis_z,gnomad$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pMiss)))]
  toreturn = cbind(genes,pMiss[match(genes,names(pMiss))])
  colnames(toreturn) = c("gene","gnomadpMiss")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$gnomadpMiss = as.numeric(toreturn$gnomadpMiss)
  return(toreturn)
}

getGeneValuesgnomadoeLoF = function(genes,drop.na=F){
  cat("Generating gnomADoeLoF data for",length(genes),"genes\n")
  gnomad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                                 "gnomADrelease_2.1_ht_constraint_constraint.zip"),
                          "gnomADrelease_2.1_ht_constraint_constraint.txt"),
                      stringsAsFactors=F,sep="\t")
  pMiss = tapply(gnomad$oe_lof,gnomad$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pMiss)))]
  toreturn = cbind(genes,pMiss[match(genes,names(pMiss))])
  colnames(toreturn) = c("gene","gnomadOELoF")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$gnomadOELoF = as.numeric(toreturn$gnomadOELoF)
  return(toreturn)
}

getGeneValuesgnomadoeMiss = function(genes,drop.na=F){
  cat("Generating gnomADoeMiss data for",length(genes),"genes\n")
  gnomad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                                 "gnomADrelease_2.1_ht_constraint_constraint.zip"),
                          "gnomADrelease_2.1_ht_constraint_constraint.txt"),
                      stringsAsFactors=F,sep="\t")
  pMiss = tapply(gnomad$oe_mis,gnomad$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pMiss)))]
  toreturn = cbind(genes,pMiss[match(genes,names(pMiss))])
  colnames(toreturn) = c("gene","gnomadOEMiss")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$gnomadOEMiss = as.numeric(toreturn$gnomadOEMiss)
  return(toreturn)
}
getGeneValuesgnomadpNull = function(genes,drop.na=F){
  cat("Generating gnomADpNull data for",length(genes),"genes\n")
  gnomad = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                                 "gnomADrelease_2.1_ht_constraint_constraint.zip"),
                          "gnomADrelease_2.1_ht_constraint_constraint.txt"),
                      stringsAsFactors=F,sep="\t")
  pNull = tapply(gnomad$pNull,gnomad$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pNull)))]
  toreturn = cbind(genes,pNull[match(genes,names(pNull))])
  colnames(toreturn) = c("gene","gnomadpNull")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$gnomadpNull = as.numeric(toreturn$gnomadpNull)
  return(toreturn)
}

getGeneValuesExACpLI = function(genes,drop.na=F){
  cat("Generating ExACpLI data for",length(genes),"genes\n")
  exac = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),"fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"),
                    stringsAsFactors=F,sep="\t")
  pli = tapply(exac$pLI,exac$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","ExACpLI")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$ExACpLI = as.numeric(toreturn$ExACpLI)
  return(toreturn)
}

getGeneValuesExACpRec = function(genes,drop.na=F){
  cat("Generating ExACpRec data for",length(genes),"genes\n")
  exac = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),"fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"),
                    stringsAsFactors=F,sep="\t")
  pli = tapply(exac$pRec,exac$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","ExACpRec")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$ExACpRec = as.numeric(toreturn$ExACpRec)
  return(toreturn)
}

getGeneValuesExACpNull = function(genes,drop.na=F){
  cat("Generating ExACpNull data for",length(genes),"genes\n")
  exac = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),"fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"),
                    stringsAsFactors=F,sep="\t")
  pli = tapply(exac$pNull,exac$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","ExACpNull")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$ExACpNull = as.numeric(toreturn$ExACpNull)
  return(toreturn)
}

getGeneValuesExACpMiss = function(genes,drop.na=F){
  cat("Generating ExACpMiss data for",length(genes),"genes\n")
  exac = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),"fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"),
                    stringsAsFactors=F,sep="\t")
  #pli = tapply(exac$mu_mis,exac$gene,max)
  pli = tapply(exac$mis_z,exac$gene,max)
  if(drop.na)
    genes = genes[!is.na(match(genes,names(pli)))]
  toreturn = cbind(genes,pli[match(genes,names(pli))])
  colnames(toreturn) = c("gene","ExACpMiss")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$ExACpMiss = as.numeric(toreturn$ExACpMiss)
  return(toreturn)
}



getGeneValuesHexConstExons = function(genes,drop.na=F){
  hexevent = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),
                               "hexeventdbFull.txt"),
                        stringsAsFactors=F,sep="\t")
  exons = table(hexevent$genename)
  exons = exons[names(exons) != "onlyEST"]

  if(drop.na)
    genes = genes[!is.na(match(genes,names(exons)))]
  toreturn = cbind(genes,exons[match(genes,names(exons))])
  colnames(toreturn) = c("gene","constitutiveexons")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$constitutiveexons = as.numeric(toreturn$constitutiveexons)
  return(toreturn)
}

getGeneValuesHexESTCount = function(genes,drop.na=F){
  hexevent = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),
                               "hexeventdbFull.txt"),
                        stringsAsFactors=F,sep="\t")
  exons = tapply(hexevent$count,as.factor(hexevent$genename),sum)
  exons = exons[names(exons) != "onlyEST"]

  if(drop.na)
    genes = genes[!is.na(match(genes,names(exons)))]
  toreturn = cbind(genes,exons[match(genes,names(exons))])
  colnames(toreturn) = c("gene","ESTcount")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$ESTcount = as.numeric(toreturn$ESTcount)
  return(toreturn)
}

getGeneValuesHexAlt3EST = function(genes,drop.na=F){
  cat("Generating Alt 3' EST data for",length(genes),"genes\n")
  hexevent = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),
                               "hexeventdbFull.txt"),
                        stringsAsFactors=F,sep="\t")
  exons = tapply(hexevent$alt3,as.factor(hexevent$genename),sum)
  exons = exons[names(exons) != "onlyEST"]

  if(drop.na)
    genes = genes[!is.na(match(genes,names(exons)))]
  toreturn = cbind(genes,exons[match(genes,names(exons))])
  colnames(toreturn) = c("gene","alt3EST")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$alt3EST = as.numeric(toreturn$alt3EST)
  return(toreturn)

}

getGeneValuesHexAlt5EST = function(genes,drop.na=F){
  cat("Generating Alt 5' EST data for",length(genes),"genes\n")
  hexevent = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),
                               "hexeventdbFull.txt"),
                        stringsAsFactors=F,sep="\t")
  exons = tapply(hexevent$alt5,as.factor(hexevent$genename),sum)
  exons = exons[names(exons) != "onlyEST"]

  if(drop.na)
    genes = genes[!is.na(match(genes,names(exons)))]
  toreturn = cbind(genes,exons[match(genes,names(exons))])
  colnames(toreturn) = c("gene","alt5EST")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$alt5EST = as.numeric(toreturn$alt5EST)
  return(toreturn)
}

getGeneValuesString = function(genes,drop.na=F){
  cat("Generating String experimental score data for",length(genes),"genes\n")
  domino = read.delim(unz(paste0(system.file("g2pml/", "", package = "G2PML"),
                             "domino_score_all_final_03.04.17.txt.zip"),"domino_score_all_final_03.04.17.txt"),
                      stringsAsFactors=F,sep="\t")
  domino$Gene = fromSymbol2Hugo(domino$Gene)
  #pli = tapply(exac$mu_mis,exac$gene,max)
  stringsc = domino$STRING.combinedscore[match(genes,domino$Gene)]
  if(drop.na){
    genes = genes[is.na(stringsc)]
    stringsc = na.omit(stringsc)
  }
  toreturn = cbind(genes,stringsc)
  colnames(toreturn) = c("gene","String")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$String = as.numeric(toreturn$String)
  return(toreturn)
}

getGeneValuesHexAlt3.5EST = function(genes,drop.na=F){
  cat("Generating Alt 3' & 5' EST data for",length(genes),"genes\n")
  hexevent = read.delim(paste0(system.file("g2pml/", "", package = "G2PML"),
                               "hexeventdbFull.txt"),
                        stringsAsFactors=F,sep="\t")
  exons = tapply(hexevent$alt3.5,as.factor(hexevent$genename),sum)
  exons = exons[names(exons) != "onlyEST"]

  if(drop.na)
    genes = genes[!is.na(match(genes,names(exons)))]
  toreturn = cbind(genes,exons[match(genes,names(exons))])
  colnames(toreturn) = c("gene","alt3.5EST")
  toreturn = as.data.frame(toreturn,stringsAsFactors=F)
  toreturn$alt3.5EST = as.numeric(toreturn$alt3.5EST)
  return(toreturn)
}

availableSpecificityData = function(ens.genes){
  expr = readRDS(paste0(system.file("g2pml/", "", package = "G2PML"),
                        "adjmatrix.rds"))
  return(ens.genes %in% colnames(expr))
}

getSpecificityMatrix = function(which.one,fold=3,visibility=5){
  if(which.one == "Expression")
    return(readRDS(paste0(system.file("g2pml/", "", package = "G2PML"),
                          "exprmatrix.",fold,".",visibility,".rds")))
  if(which.one == "Adjacency")
    return(readRDS(paste0(system.file("g2pml/", "", package = "G2PML"),
                          "adjmatrix.",fold,".",visibility,".rds")))
  if(which.one == "RankedMM")
    return(readRDS(paste0(system.file("g2pml/", "", package = "G2PML"),
                          "rankedmmmatrix.",fold,".",visibility,".rds")))
  return(NULL)
}
getGeneValuesExpressionSpecificity = function(genes,use.cluster=F,ens.correction=NULL,
                                              fold=3,visibility=5){
  cat("Generating expression data specificity for",length(genes),"genes\n")
  expr = readRDS(paste0(system.file("g2pml/", "", package = "G2PML"),
                        "exprmatrix.",fold,".",visibility,".rds"))
  rownames(expr) = unlist(lapply(rownames(expr),fromLong2ShortTissue))

  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)

  if(!is.null(ens.correction)){
    for(i in seq(from=1,to=length(ens.correction),by=2)){
      mask = which(ens.genes == ens.correction[i])
      if(sum(mask) > 0){
        ens.genes[mask] = ens.correction[i+1]
        cat("Substituting", ens.correction[i],"by",ens.correction[i+1],"\n")
      }
    }
  }

  datamatrix = matrix(ncol=nrow(expr),nrow=length(genes))
  datamatrix[,] = 0
  colnames(datamatrix) = rownames(expr)
  rownames(datamatrix) = genes
  maskexpressedgenes = ens.genes %in% colnames(expr)
  #else
  #	ens.genes = genes
  datamatrix[maskexpressedgenes,] = t(expr[,match(ens.genes[maskexpressedgenes],colnames(expr))])
  expr = as.data.frame(cbind(genes,datamatrix),stringsAsFactors=F)
  colnames(expr)[1] = "gene"
  colnames(expr)[2:ncol(expr)] = paste0("ExprSpecific",colnames(expr)[2:ncol(expr)])
  return(expr)
}

getGeneValuesAdjacencySpecificity = function(genes,use.cluster=F,ens.correction=NULL,
                                             fold=3,visibility=5){
  cat("Generating adjacency data specificity for",length(genes),"genes\n")
  expr = readRDS(paste0(system.file("g2pml/", "", package = "G2PML"),"adjmatrix.",fold,".",visibility,".rds"))
  rownames(expr) = unlist(lapply(rownames(expr),fromLong2ShortTissue))
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)

  if(!is.null(ens.correction)){
    for(i in seq(from=1,to=length(ens.correction),by=2)){
      mask = which(ens.genes == ens.correction[i])
      if(sum(mask) > 0){
        ens.genes[mask] = ens.correction[i+1]
        cat("Substituting", ens.correction[i],"by",ens.correction[i+1],"\n")
      }
    }
  }
  datamatrix = matrix(ncol=nrow(expr),nrow=length(genes))
  datamatrix[,] = 0
  colnames(datamatrix) = rownames(expr)
  rownames(datamatrix) = genes

  maskexpressedgenes = ens.genes %in% colnames(expr)
  #else
  #	ens.genes = genes
  datamatrix[maskexpressedgenes,] = t(expr[,match(ens.genes[maskexpressedgenes],colnames(expr))])
  expr = as.data.frame(cbind(genes,datamatrix),stringsAsFactors=F)
  colnames(expr)[1] = "gene"
  colnames(expr)[2:ncol(expr)] = paste0("AdjSpecific",colnames(expr)[2:ncol(expr)])
  return(expr)

}

getGeneValuesRankedMMSpecificity = function(genes,use.cluster=F,
                                            ens.correction=NULL,
                                            fold=3,visibility=5){

  cat("Generating Ranked MM data specificity for",length(genes),"genes\n")
  expr = readRDS(paste0(system.file("g2pml/", "", package = "G2PML"),
                        "rankedmmmatrix.",fold,".",visibility,".rds"))
  rownames(expr) = unlist(lapply(rownames(expr),fromLong2ShortTissue))
  if(is.null(names(genes)))
    ens.genes = fromGeneName2EnsemblBM(genes)
  else
    ens.genes = names(genes)

  if(!is.null(ens.correction)){
    for(i in seq(from=1,to=length(ens.correction),by=2)){
      mask = which(ens.genes == ens.correction[i])
      if(sum(mask) > 0){
        ens.genes[mask] = ens.correction[i+1]
        cat("Substituting", ens.correction[i],"by",ens.correction[i+1],"\n")
      }
    }
  }
  datamatrix = matrix(ncol=nrow(expr),nrow=length(genes))
  datamatrix[,] = 0
  colnames(datamatrix) = rownames(expr)
  rownames(datamatrix) = genes

  maskexpressedgenes = ens.genes %in% colnames(expr)
  #else
  #	ens.genes = genes
  datamatrix[maskexpressedgenes,] = t(expr[,match(ens.genes[maskexpressedgenes],colnames(expr))])
  expr = as.data.frame(cbind(genes,datamatrix),stringsAsFactors=F)
  colnames(expr)[1] = "gene"
  colnames(expr)[2:ncol(expr)] = paste0("RankedMMSpecific",colnames(expr)[2:ncol(expr)])
  return(expr)
}

#' Title
#'
#' @param panel
#' @param when
#' @param evidence
#'
#' @return
#' @export
#'
#' @examples
getGelGenes = function(panel,when="geOctober2018",evidence="HighEvidence"){
  gpath = paste0(system.file("g2pml/gel/", "", package = "G2PML"),when,"/")
  g = read.csv(paste0(gpath,"/",panel,".csv"),stringsAsFactors=F)
  return(g$GeneSymbol[g$LevelOfConfidence %in% evidence])
}

#' Obtaining ML features for your genes of interest
#'
#' \code{fromGenes2MLData} obtains the genetic properties (transcriptomic,
#' coexpression, genetic constraint.. etc) for a given set of gene symbols
#'
#' @param genes chr vector. Gene symbols of your disease genes - can be returned
#'   from \code{\link{getGenesFromPanelApp}}.
#' @param addcontrols lgl scalar. Do you want to add a set of control genes?
#' @param which.controls chr scalar. One of "allghosh", "allgenome",
#'   "clustering", "gauss", "gausskfold" specifying the set of control genes you
#'   would like to use.
#' @param condition chr vector. Vector of length genes describing which are
#'   "Disease" and "Nondisease".
#' @param vars chr vector. Names of features you would like to include.
#' @param filter chr vector. Names of features you would like to exclude.
#' @param ... additional arguments for clustering, only used when which.controls
#'   is "clustering".
#'
#' @return df with features of input genes formatted for ML.
#' @export
#'
#' @examples
#' genes = getGenesFromPanelApp(disorder="Neurology and neurodevelopmental disorders",
#'   panel="Parkinson Disease and Complex Parkinsonism", color = "green")
#' genedata = fromGenes2MLData(genes=genes, which.controls="allgenome")
fromGenes2MLData = function(genes,
                            addcontrols=T,
                            which.controls="allghosh",
                            condition=NULL,
                            vars=NULL,
                            filter=c("DPI","DSI","ESTcount","constitutiveexons"),
                            ...){
  if(addcontrols){
    if(which.controls == "allghosh"){
      ndgenes = read.csv(paste0(system.file("g2pml/", "", package = "G2PML"),
                                "Ghoshpaperdiseasenondisease.csv"),stringsAsFactors=F)
      #Keep only those that are non-disease ones, and convert from Ensembl to gene ID
      ndgenes = fromEnsembl2GeneNameBM(ndgenes$gene[ndgenes$type == "non-disease"])
      #Detect whether some of the control genes are disease too and remove
      mask = ndgenes %in% genes
      ndgenes = ndgenes[!mask]
      if(sum(mask) > 0)
        cat("We have to remove",sum(mask),"disease genes from the controls\n")
      #Now we create the ML datafile with control and disease genes
      #We don't want DPI, DSI, ESTCount and constitutive exons so we filter them
      #And write everything into tmp/myfile.txt
      disease = getSingleConditionSet(genes=genes)
      nondisease = getSingleConditionSet(genes=ndgenes,condition="Nondisease")
      data = genLearningDataSet(casecontrolset=rbind(disease,nondisease),
                                filter=filter)
      if(!is.null(vars))
        data = data[,c("gene","condition",vars)]
      return(data)
    }else if(which.controls == "allgenome"){
      allgenes = getCodingGenome()
      mask = allgenes %in% genes
      allgenes = allgenes[!mask]
      if(sum(mask) > 0)
        cat("ALLGENOME mode: We have to remove",sum(mask),"disease genes from the all genome controls\n")
      #Now we create the ML datafile with control and disease genes
      #We don't want DPI, DSI, ESTCount and constitutive exons so we filter them
      #And write everything into tmp/myfile.txt
      disease = getSingleConditionSet(genes=genes)
      nondisease = getSingleConditionSet(genes=allgenes,condition="Nondisease")
      data = genLearningDataSet(casecontrolset=rbind(disease,nondisease),
                                filter=filter)
      if(!is.null(vars))
        data = data[,c("gene","condition",vars)]
      return(data)
    }else if(which.controls == "clustering"){
      allgenes = doLearningByClusteringv2(genes=genes,vars=vars,plots=F,...)$controls
      mask = allgenes %in% genes
      allgenes = allgenes[!mask]
      if(sum(mask) > 0)
        cat("We have to remove",sum(mask),"disease genes from clustering controls\n")
      disease = getSingleConditionSet(genes=genes)
      nondisease = getSingleConditionSet(genes=allgenes,condition="Nondisease")
      data = genLearningDataSet(casecontrolset=rbind(disease,nondisease),
                                filter=filter)
      if(!is.null(vars))
        data = data[,c("gene","condition",vars)]
      return(data)
    }else if(which.controls == "gauss"){
      allgenes = doLearningByClusteringv3(genes=genes,
                                          panel="Unknown",
                                          plots=F)$controls
      mask = allgenes %in% genes
      allgenes = allgenes[!mask]
      if(sum(mask) > 0)
        cat("We have to remove",sum(mask),"disease genes from clustering controls\n")
      disease = getSingleConditionSet(genes=genes)
      nondisease = getSingleConditionSet(genes=allgenes,condition="Nondisease")
      data = genLearningDataSet(casecontrolset=rbind(disease,nondisease),
                                filter=filter)
      if(!is.null(vars))
        data = data[,c("gene","condition",vars)]
      return(data)
    }else if(which.controls == "gausskfold"){
      allgenes = doLearningByclusteringv3KFold(genes=genes,
                                               panel="Unknown",
                                               vars=vars,
                                               folds=10)
      mask = allgenes %in% genes
      allgenes = allgenes[!mask]
      if(sum(mask) > 0)
        cat("We have to remove",sum(mask),"disease genes from clustering controls\n")
      disease = getSingleConditionSet(genes=genes)
      nondisease = getSingleConditionSet(genes=allgenes,condition="Nondisease")
      data = genLearningDataSet(casecontrolset=rbind(disease,nondisease),
                                filter=filter)
      if(!is.null(vars))
        data = data[,c("gene","condition",vars)]
      return(data)
    }else {
      cat("No controls ",which.controls,"recognized in fromGenes2MLData\n")
      return(NULL)
    }

  }
  if(!is.null(condition)){
    dgenes = genes[condition == "Disease"]
    ndgenes = genes[condition == "Nondisease"]
    disease = getSingleConditionSet(genes=dgenes)
    nondisease = getSingleConditionSet(genes=ndgenes,condition="Nondisease")
    data = genLearningDataSet(casecontrolset=rbind(disease,nondisease),
                              filter=filter)
    if(!is.null(vars))
      data = data[,c("gene","condition",vars)]
    return(data)
  }
  disease = getSingleConditionSet(genes=genes)
  data = genLearningDataSet(casecontrolset=disease,
                            filter=filter)
  if(!is.null(vars))
    data = data[,c("gene","condition",vars)]
  return(data)

}


genLearningDataSet = function(casecontrolset=getCaseControlSet(which.ones="ge_neurogenes"),
                              sep="\t",
                              filter=NULL,
                              newdata=T){

  f.in = paste0(system.file("g2pml/", "", package = "G2PML"),
                "mlDBMarch2019.txt.zip")

  genes = casecontrolset$gene
  mldata = read.delim(unz(f.in,"mlDBMarch2019.txt"),stringsAsFactors=F,sep="\t")
  mask = match(genes,mldata$gene)
  nulls = sum(is.na(mask))
  if(nulls > 0){
    toprint = genes[is.na(mask)]
    if(length(toprint) > 5)
      toprint = toprint[1:5]
    cat("We have to remove",nulls,"genes",toprint,
        " because we don't know them\n")
  }
  mldata = mldata[na.omit(mask),]
  mldata$condition = casecontrolset$condition[match(mldata$gene,casecontrolset$gene)]



  if(!is.null(filter)){
    cat("Removing columns for attributes",paste0(filter,collapse=" and "),"\n")
    mldata = mldata[,!(colnames(mldata) %in% filter)]
  }
  return(mldata)
}

getSingleConditionSet = function(genes=NULL,condition="Disease"){
  out.data = cbind(genes,rep(condition,length(genes)))
  colnames(out.data) = c("genes","condition")
  return(as.data.frame(out.data,stringsAsFactors=F))
}

fromSymbol2Hugo = function(genes,debug=F){

  if(debug){
    if(exists("hugodb")) rm(hugodb,envir=globalenv())
    if(exists("hugoalias")) rm(hugoalias,envir=globalenv())
    if(exists("prevsymbol")) rm(prevsymbol,envir=globalenv())
  }

  if(!exists("hugodb")){
    hugodb <<- read.delim(paste0(system.file("g2pml", "", package = "G2PML"),
                                 "/gene_with_protein_product.txt"),stringsAsFactors=F)
    hugoalias <<- NULL
    prevsymbol <<- NULL
    for(i in 1:nrow(hugodb)){
      aliases = hugodb[i,"alias_symbol"]
      if(aliases != ""){
        symbol = hugodb[i,"symbol"]
        aliases = strsplit(aliases,"\\|")
        lapply(aliases[[1]],function(x){ hugoalias[[x]] <<- symbol })
      }

      previous = hugodb[i,"prev_symbol"]
      if(previous != ""){
        symbol = hugodb[i,"symbol"]
        previous = strsplit(previous,"\\|")
        lapply(previous[[1]],function(x){ prevsymbol[[x]] <<- symbol })
        #prevsymbol[[previous]] <<- symbol
      }
    }
  }

  #Where are genes?
  symbolgenes = genes[genes %in% hugodb$symbol]
  notfound = genes[!(genes %in% symbolgenes)]
  foundasalias = notfound[notfound %in% names(hugoalias)]
  stillnotfound = notfound[!(notfound %in% foundasalias)]
  foundasprev = stillnotfound[stillnotfound %in% names(prevsymbol)]
  stillnotfound = stillnotfound[!(stillnotfound %in% foundasprev)]

  if(length(foundasalias) > 0){
    #cat(foundasalias,"#")
    notfoundsymbol = unlist(lapply(foundasalias,function(x){ return(hugoalias[[x]])}))
    mask = match(symbolgenes,genes)
    genes[mask] = symbolgenes
    mask = match(foundasalias,genes)
    genes[mask] = notfoundsymbol
  }

  if(length(foundasprev) > 0){
    symbol = unlist(lapply(foundasprev,function(x){ return(prevsymbol[[x]]) }))
    mask = match(foundasprev,genes)
    genes[mask] = symbol
  }

  mask = match(stillnotfound,genes)
  genes[mask] = NA

  return(genes)
}

getCodingGenome = function(){
  genes = unique(read.table(paste0(system.file("g2pml/", "", package = "G2PML"),
                                   "codinggenome.txt"),
                            header=F,stringsAsFactors=F)$V1)

  genes = fromSymbol2Hugo(genes)
  return(na.omit(genes))
}

fromEnsembl2GeneNameBM = function(genes,use38=T){
  if(use38){
    ensembl <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_name"
  }else{
    ensembl <- biomaRt::useMart(host="jun2013.archive.ensembl.org",
                       biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_id"
  }

  attributes <- c("ensembl_gene_id",external.gene.att)
  genes.with.name = biomaRt::getBM(attributes=attributes, filters="ensembl_gene_id", values=genes,mart=ensembl)
  cat("From",length(genes),"Ensembl IDs we got",nrow(genes.with.name),"genes with external gene name\n")
  thematch = match(genes,genes.with.name$ensembl_gene_id)
  outgenes = genes.with.name[,2][thematch]
  if(sum(is.na(thematch)) == 0)
    return(outgenes)
  cat("Couldn't conver",sum(is.na(thematch)),"genes\n")
  outgenes[is.na(thematch)] = genes[is.na(thematch)]
  return(outgenes)


  thematch = match(genes,genes.with.name$external_gene_name)
  outgenes = genes.with.name[,2][thematch]
  outgenes[is.na(thematch)] = genes[is.na(thematch)]
  cat("fromGeneName2EnsemblBM, couldn't convert",sum(is.na(thematch)),"genes\n")
  return(outgenes)
}

fromGeneName2EnsemblBM = function(genes,use38=T){

  if(use38){
    ensembl <- biomaRt::useMart(host="www.ensembl.org",
                       biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_name"
  }else{
    ensembl <- biomaRt::useMart(host="jun2013.archive.ensembl.org",
                       biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
    external.gene.att = "external_gene_id"
  }

  attributes <- c(external.gene.att,"ensembl_gene_id")
  genes.with.name = biomaRt::getBM(attributes=attributes, filters="hgnc_symbol", values=genes,mart=ensembl)
  cat("From",length(genes),"gene IDs we got",nrow(genes.with.name),"genes with Ensemble name\n")
  #if(nrow(genes.with.name) >= length(genes))
  thematch = match(genes,genes.with.name$external_gene_name)
  outgenes = genes.with.name[,2][thematch]
  outgenes[is.na(thematch)] = genes[is.na(thematch)]
  cat("fromGeneName2EnsemblBM, couldn't convert",sum(is.na(thematch)),"genes\n")
  return(outgenes)
}

fromLong2ShortTissue = function(tissue){

  if(tissue == "Adipose - Subcutaneous")
    return("AdiposeSub")
  if(tissue == "Adipose - Visceral (Omentum)")
    return("AdiposeVisceral")
  if(tissue == "Adrenal Gland")
    return("AdrenalGland")
  if(tissue == "Artery - Aorta")
    return("ArteryAorta")
  if(tissue == "Artery - Coronary")
    return("ArteryCoronary")
  if(tissue == "Artery - Tibial")
    return("ArteryTibial")
  if(tissue == "Brain - Amygdala")
    return("Amygdala")
  if(tissue == "Brain - Substantia nigra")
    return("Substantianigra")
  if(tissue == "Brain - Spinal cord (cervical c-1)")
    return("Spinalcord")
  if(tissue == "Brain - Anterior cingulate cortex (BA24)")
    return("AntCingCortex")
  if(tissue == "Brain - Caudate (basal ganglia)")
    return("Caudate")
  if(tissue == "Brain - Cerebellar Hemisphere")
    return("CerebHemisphere")
  if(tissue == "Brain - Cerebellum")
    return("Cerebellum")
  if(tissue == "Brain - Cortex")
    return("Cortex")
  if(tissue == "Brain - Frontal Cortex (BA9)")
    return("FCortex")
  if(tissue == "Brain - Hippocampus")
    return("Hippocampus")
  if(tissue == "Brain - Hypothalamus")
    return("Hypothalamus")
  if(tissue == "Brain - Nucleus accumbens (basal ganglia)")
    return("NucAccumbens")
  if(tissue == "Brain - Putamen (basal ganglia)")
    return("Putamen")

  if(tissue == "Breast - Mammary Tissue")
    return("Breast")

  if(tissue == "Cells - EBV-transformed lymphocytes")
    return("CellsLymphocytes")
  if(tissue == "Cells - Leukemia cell line (CML)")
    return("CellsLeukemiaCL")

  if(tissue == "Cells - Transformed fibroblasts")
    return("CellsFirbroblasts")

  if(tissue == "Cells - Transformed firbroblasts")
    return("CellsFirbroblasts")

  if(tissue == "Colon - Sigmoid")
    return("ColonSigmoid")
  if(tissue == "Colon - Transverse")
    return("ColonTransverse")
  if(tissue == "Esophagus - Gastroesophageal Junction")
    return("EsophGastJunction")
  if(tissue == "Esophagus - Mucosa")
    return("EsophMucosa")
  if(tissue == "Esophagus - Muscularis")
    return("EsophMuscularis")



  if(tissue == "Heart - Atrial Appendage")
    return("HeartAtrialApp")
  if(tissue == "Heart - Left Ventricle")
    return("HeartLeftVent")
  if(tissue == "Muscle - Skeletal")
    return("MuscleSkeletal")
  if(tissue == "Nerve - Tibial")
    return("NerveTibial")
  if(tissue == "Small Intestine - Terminal Ileum")
    return("SmallIntestine")

  if(tissue == "Skin - Not Sun Exposed (Suprapubic)")
    return("SkinSuprapubic")
  if(tissue == "Skin - Sun Exposed (Lower leg)")
    return("SkinLowerLeg")
  if(tissue == "Whole Blood")
    return("WholeBlood")
  return(tissue)

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
getPanelsFromPanelApp = function(){
  req <- curl::curl_fetch_memory("panelapp.genomicsengland.co.uk/WebServices/list_panels/")
  myjson = jsonlite::fromJSON(rawToChar(req$content))
  myjson$results$Relevant_disorders = NULL
  return(myjson$result)
}

#' Title
#'
#' @param disorder
#' @param panel
#'
#' @return
#' @export
#'
#' @examples
getPanelFromPanelApp = function(disorder="Neurology and neurodevelopmental disorders",
                                panel="Parkinson Disease and Complex Parkinsonism")
{
  panels = getPanelsFromPanelApp()
  if(!(panel %in% panels$Name)){
    cat("Panel",panel,"not found at the Panel App, available panels are\n")
    print(panels$Name)
  }
  stopifnot(panel %in% panels$Name)
  panelweb = gsub(" ","%20",panel)
  req <- curl::curl_fetch_memory(paste0("panelapp.genomicsengland.co.uk/WebServices/get_panel/",
                                    panelweb,"/"))
    myjson = jsonlite::fromJSON(rawToChar(req$content))
    out = cbind(myjson$result$Genes$GeneSymbol,myjson$result$Genes$LevelOfConfidence,
                myjson$result$Genes$ModeOfInheritance,myjson$result$Genes$ModeOfPathogenicity)
    colnames(out) = c("GeneSymbol","LevelOfConfidence","ModeOfInheritance","ModeOfPathogenicity")
  return(as.data.frame(out,stringsAsFactors=F))
}

#' Title
#'
#' @param disorder
#' @param panel
#' @param color
#'
#' @return
#' @export
#'
#' @examples
getGenesFromPanelApp = function(disorder="Neurology and neurodevelopmental disorders",
                     panel="Parkinson Disease and Complex Parkinsonism",
                     color){

  evidence = NULL
  if("red" %in% color)
    evidence = c(evidence,"LowEvidence")
  if(color == "amber")
    evidence = c(evidence,"ModerateEvidence")
  if("green" %in% color)
    evidence = c(evidence,"HighEvidence")
  genes = getPanelFromPanelApp(disorder,panel)
  return(genes$GeneSymbol[genes$LevelOfConfidence %in% evidence])
}


