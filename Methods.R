#######################################################################################
#                                                                                     #
# List of the Differentially Expression (DE) methods used for scRNA-seq data analysis #
#                                                                                     #
#######################################################################################
setwd(....) ########set working directory
countData <- read.table(file = "Islam.txt", row.names = 1, header = T, sep = "\t") ###Example data

######## Filtering of data to remove lowly expressed genes
#' @param countData UMI data matrix with genes as rows and columns as cells
#' @param threshold minimum of number of non-zero counts for cells to be retained
#' @param Cellsize Minimum library sizes of the cells to be retained in the data
FilterData <- function (countData, threshold, Cellsize){
  countData <- countData[rowSums(countData)!=0, ]
  f0 <- function(y) length(y) - sum(y==0)
  non.0 <- apply(countData, 1, f0)
  id <- which( non.0 > threshold)
  if(is.null(Cellsize))
    countData <- countData[id, ]
  if(!is.null(Cellsize))
    countData <- countData[id, colSums(countData) > Cellsize]
  return(countData)
}

##########processed data
countData <- FilterData (countData = countData, threshold, Cellsize) 
###enter the values for threshold and Cellsize 

###########Program for methods
#' @param CountData countData UMI data matrix with genes as rows and columns as cells
#' @param group vector of 1 and -1, represents the group memberships of the cells
#' 
##### 1 DESeqNB method (ver 1.39.0)
RunDESeqNB = function(CountData , group){
  require("DESeq")
  gr1 <- length(which(group==1))
  gr2 <- length(which(group==2))
  #dataPack = data.frame(row.names= colnames(FilteredData), condition = c(rep("GR1",gr1) ,rep("GR2",gr2)))
  conds <- factor (c(rep("GR1",gr1) ,rep("GR2",gr2)))
  # Initializing a CountDataSet (DESeq data structure for count matrix):
  cds <- newCountDataSet(CountData, conds)
  # DESeq normaliztion:
  cds <- estimateSizeFactors (cds)
  #print (sizeFactors(cds))
  #print (head(counts(cds, normalized = TRUE)))
  # Variance estimation:
  cds <- estimateDispersions (cds)
  #print(str(fitInfo(cds)))
  # Negative binomial test:
  DESeq_results = nbinomTest(cds, "GR1", "GR2")
  return(DESeq_results)
}

##### 2 DESeqLRT (version 1.28.1)
runDEseqLRT <- function(CountData, group){
  group <- as.factor(group)
  CountData <- as.matrix(CountData)
  require(DESeq2)
  dds <- DESeqDataSetFromMatrix(CountData, DataFrame(group), ~ group)
  ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
  resLRT <- results(ddsLRT)
  return(resLRT)
}

##### 3 DEGSeq (ver 1.42.0)
runDEGSeq <- function(CountData, group){
  require(DEGseq)
  outputDir <- file.path(getwd(), "GSE75790")
  CountGroup1 <- cbind(rownames(CountData), CountData[, which(group==1)])
  CountGroup2 <- cbind(rownames(CountData), CountData[, which(group==2)])
  res <- DEGexp(geneExpMatrix1=CountGroup1, geneCol1=1, expCol1=c(2:length(which(group==1))+1),
                geneExpMatrix2=CountGroup2, geneCol2=1, expCol2=c(2:length(which(group==2))+1),
                method="LRT", outputDir = outputDir)
  result<-read.table(file=file.path(outputDir,"output_score.txt"),header = T,sep = "\t")
  colnames(result)<-c("GeneNames", "value1", "value2", "log2_Fold_change", "log2_Fold_change_normalized","pval","qval1","qval2","Signature")
  return(result)
}

##### 4 LIMMA (version 3.44.1)
runLIMMA <- function(CountData, group){
  library(edgeR)
  require(limma)
  group <- as.factor(group)
  d0 <- DGEList(CountData)
  mm <- model.matrix(~0 + group)
  y <- voom(d0, mm, plot = F)
  fit <- lmFit(y, mm)
  contr <- makeContrasts(group1 - group2, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  pval <- tmp$p.value
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  return(top.table)
}

##### 5 ROTS (version 1.16.0)

TMMnormalization <- function(countTable){
  ## TMM normalization based on edgeR package:
  require("edgeR")
  nf=calcNormFactors(countTable ,method= "TMM")
  nf= colSums(countTable)*nf
  scalingFactors = nf/mean(nf)
  countTableTMM <- t(t(countTable)/scalingFactors)
  return(countTableTMM)
}

#####Running ROTS:
runROTS = function(CountData, group){
  require(ROTS)
  # countData is the input raw count table
  # B is the number of bootstraps
  # k is the number of top ranked genes
  # FDR is the fdr threshold for the final detections
  # First the genes with 0 value for all the cells were filtered out:
  # TMM normalization over the raw count matrix:
  TMM_Filtered_Data = TMMnormalization(CountData)
  # Running ROTS over the filtered and normalized data:
  ROTS_filtered_TMM = ROTS(data = TMM_Filtered_Data, groups = group , B = 1000, K = 6000 )
  reslt <- cbind(ROTS_filtered_TMM$logfc, ROTS_filtered_TMM$pvalue, ROTS_filtered_TMM$FDR)
  colnames(reslt) <- c("LogFC", "pvalue", "FDR")
  # Rows/genes with FDR smaller than the desired threshold:
  return(reslt)
}

##### 6 edgeR (ver 3.30.3)
runedgeRLRT <- function (CountData, group) {
  require(edgeR)
  group <- as.factor(group)
  y <- DGEList(counts = CountData, group = group)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, model.matrix(~group))
  fit <- glmFit(y,design = model.matrix(~group))
  lrt <- glmLRT(fit, coef=2:length(levels(group)))
  return(lrt$table)
}

##### 7 EdgeR with quasilikelihood function
runedgeRQLF <- function (CountData, group) {
  require(edgeR)
  group <- as.factor(group)
  y <- DGEList(counts = CountData, group = group)
  y <- calcNormFactors(y)
  # Fit the NB GLMs with QL methods
  y <- estimateDisp(y, model.matrix(~group))
  fit <- glmQLFit(y, design = model.matrix(~group), robust=TRUE)
  results <- glmQLFTest(fit)
}


##### 8 Monocle (version 2.16.0)
runMonocle <- function (CountData, group, UMI = T, n.cores = 2) {
  require(monocle)
  CountData <- round(as.matrix(CountData))
  storage.mode(CountData) <- "integer"
  pd <- new("AnnotatedDataFrame",
            data = data.frame(CellType = as.factor(group), row.names = colnames(CountData)))
  fd <- new("AnnotatedDataFrame",
            data = data.frame(gene_short_name = rownames(CountData), row.names = rownames(CountData)))
  # IF UMI data
  if (UMI){
    MNC <- newCellDataSet(CountData,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily=negbinomial.size())
  } else {
    # IF non-UMI data, first convert to TPM outside this function
    MNC <- newCellDataSet(CountData,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit=0.1,
                          expressionFamily=tobit(Lower=0.1))
    rpc_matrix <- relative2abs(MNC) # Census
    MNC <- newCellDataSet(rpc_matrix,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit=0.5,
                          expressionFamily=negbinomial.size())
  }
  MNC <- estimateSizeFactors(MNC)
  mnc.de <- differentialGeneTest(MNC, fullModelFormulaStr="~CellType", cores = n.cores, verbose = T)
  return(mnc.de)
}

##### 9 MAST (ver 1.14.0)

runMAST <- function (CountData, group, UMI = T) {
  require(MAST)
  group <- as.factor(group)
  if(UMI){
    expr <- log2(edgeR::cpm(countData)+1)
  }
  else { # else put tpm in
    expr <- log(CountData+1)
  }
  sca <- FromMatrix(exprsArray = expr)
  cdr2 <-colSums(assay(sca)>0)
  colData(sca)$cngeneson <- scale(cdr2)
  colData(sca)$condition <- group
  zlmCond <- zlm(~condition + cngeneson, sca)
  summaryCond <- summary(zlmCond, doLRT='condition2')
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='condition2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='condition2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  return(fcHurdle)
}

##### 10 scDD (ver 1.12.0)
runscDD <- function(CountData, group, UMI = T){
  require(scDD)
  group <- as.factor(group)
  if(UMI){
    require(edgeR)
    expr <- log2(edgeR::cpm(CountData)+1)
  }
  else { # else put tpm in
    expr <- log(CountData+1)
  }
  sca=SingleCellExperiment(assays = list(normcounts = as.matrix(expr)))
  #sca <- FromMatrix(exprsArray = expr)
  colData(sca)$condition <- group
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  scDatExSim <- scDD(sca, prior_param=prior_param, testZeroes=FALSE)
  RES <- results(scDatExSim)
  return(RES)
}

#####  11 DEsingle (ver 1.8.2)

runDEsingle <- function(CountData, group){
  group <- as.factor(group)
  CountData <- as.matrix(CountData)
  require(DEsingle)
  require(stats)
  res <- DEsingle(counts = CountData, group = group, parallel = FALSE)
  #out <- res; row.names(out) <- row.names(res)
  #remove(res)
  return(res)
}

##### 12 DECENT (ver 1.10)
runDECENT <- function(CountData, group){
  require(DECENT)
  group <- as.factor(group)
  #X <- model.matrix(~group)
  CountData <- as.matrix(CountData)
  out <- decent(data.obs = CountData, X=~group, use.spikes = F, CE.range = c(0.02, 0.1))
  return(out)
}

##### 13 BPSC (ver 0.99.2)

runBPSC <- function(CountData, group){
  require(BPSC)
  CountData <- as.matrix(CountData)
  controlIds <- which(group==1)
  group <- as.factor(group)
  design <- model.matrix(~group)
  res <- BPglm(data=CountData, controlIds=controlIds, design=design, coef=2, estIntPar=FALSE)
  res.pval <- res$PVAL; names(res.pval) <- row.names(CountData)
  return(res.pval)
}

##### 14 NODES (ver 0.0.0.9010)

runNODES <- function(CountData, group){
  require(preprocessCore)
  require(NODES)
  require(MetaDE)
  cell <- colnames(CountData)
  CountData <- as.matrix(CountData)
  CountData <- pQ(CountData)
  idd <- match(colnames(CountData), cell)
  res <- NODES(data = CountData, group[idd], r = 20, smooth_points = 10000, zper = 0.5)
  return(res)
}

##### 15 EMDomics (ver 2.18.0)
runEMDomics <- function(CountData, group, n){
  require(EMDomics)
  group <- c(rep("1", length(which(group==1))), rep("2", length(which(group==2))), rep("3", n))
  CountData <- as.matrix(CountData)
  group <- as.vector(group)
  names(group) <- colnames(CountData)
  results <- calculate_cvm(data = CountData, outcomes = group, nperm = 100, pairwise.p = TRUE, parallel=FALSE)
  
  res1 <- cbind(results$pairwise.cvm.table, results$pairwise.q.table)
  return(res1)
}

###### 16 T-test (stats ver 4.0.2)
runTtest <- function(CountData, group){
  CountData <- TMMnormalization(CountData)
  id1 <- which(group==1)
  id2 <- which(group==2)
  result <- matrix(0, nrow(CountData), 2)
  for (i in 1:nrow(CountData)) {
    res <- t.test(as.vector(CountData[i, id1]), as.vector(CountData[i, id2]),  alternative = "two.sided",
                  paired = FALSE, var.equal = FALSE)
    result[i, ] <- c(res$p.value, res$statistic)
  }
  rownames(result) <- rownames(countData)
  colnames(result) <- c("Pvalue", "t-statistic")
  return(result)
}

###### 17 Wilcoxon test (stats 4.0.2)
runWilcox <- function (CountData, group){
  CountData <- TMMnormalization(CountData)
  id1 <- which(group==1)
  id2 <- which(group==2)
  result <- matrix(0, nrow(CountData), 2)
  for (i in 1:nrow(CountData)) {
    res <- wilcox.test(as.vector(CountData[i, id1]), as.vector(CountData[i, id2]),  alternative = "two.sided",
                       paired = FALSE, var.equal = FALSE)
    result[i, ] <- c(res$p.value, res$statistic)
  }
  rownames(result) <- rownames(countData)
  colnames(result) <- c("Pvalue", "Statistic")
  return(result)
}

###### 18 NBPSeq (ver 0.3.0)

runNBSeq <- function (CountData, group){
  require(NBPSeq)
  CountData <- as.matrix(CountData)
  reslt <- nbp.test(counts = CountData, grp.ids = group, grp1 = 1, grp2 = 2, norm.factors = rep(1, dim(CountData)[2]),
                    model.disp = "NBQ", lib.sizes = colSums(CountData))
  reslt1 <- cbind(reslt$log.fc, reslt$p.values, reslt$q.values)
  rownames(reslt1) <- rownames(CountData)
  colnames(reslt1) = c("logFc", "p-value", "q-value")
  return(reslt1)
}

##### 19 EBSeq (ver 1.28.0)

runEBSeq <- function(CountData, group){
  CountData <- as.matrix(CountData)
  group <- as.factor(group)
  require(EBSeq)
  Sizes = MedianNorm(CountData)
  EBOut = EBTest(Data = CountData,
                 Conditions = group,
                 sizeFactors = Sizes, maxround = 50)
  PP = GetPPMat(EBOut)
  return(PP)
}

###############Running the codes
setwd(............) #####set working directory for outputs
resDESeqNB <- RunDESeqNB(CountData = countData, group = group)
resDEseqLRT <- runDEseqLRT(CountData = countData, group = group)
resDEGSeq <- runDEGSeq(CountData = countData, group = group)
resLIMMA <- runLIMMA(CountData = countData, group = group)
resROTS <- RunROTS (CountData = countData, group = group)
resEdgeRLRT <- runedgeRLRT (CountData = countData, group = group)
resEdgRQLF <- runedgeRQLF (CountData = countData, group = group)
resMonocle <- runMonocle(CountData = countData, group = group, UMI = T, n.cores = 2)
resMAST <- runMAST (CountData = countData, group = group, UMI = T)
resSCD <- runscDD (CountData = countData, group = group, UMI = T)
resDEsingle <- runDEsingle (CountData = countData, group = group)
resDECENT <- runDECENT (CountData = countData, group = group)
resBPSC <- runBPSC (CountData = countData, group = group)
resNOD <- runNODES (CountData = countData, group = group)
resEMD <- runEMDomics (CountData = countData, group = group, n=4)
resT <- runTtest (CountData = countData, group = group)
resWil <- runWilcox (CountData = countData, group = group)
resNBseq <- runNBSeq (CountData = countData, group = group)
resEBSeq <- runEBSeq (CountData = countData, group = group)
