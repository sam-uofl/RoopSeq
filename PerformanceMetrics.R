##############################################################################
#                                                                            #
#  Computation of performance metrics                                        #
#                                                                            #
##############################################################################
#' @param genes Detected genes by the DE method
#' @param ref List of reference genes
#' @space Whole list of genes in the data
#' 
calcMeasures <- function(genes, ref, space) {
  TP <- length(intersect(ref, genes))
  #if(intersect(ref, genes) == NULL) TP = 0
  FP <- length(setdiff(genes, ref))
  negref <- setdiff(space, ref)
  neggen <- setdiff(space, genes)
  TN <- length(intersect(neggen, negref))
  FN <- length(setdiff( ref, genes))
  
  TPR <- TP/(TP + FN)
  FDR <- FP/(TP+FP)
  FPR <- FP /(FP + TN)
  PPR <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  ACC <- (TP + TN) / (TP + TN + FP + FN) 
  #MCC <- ((TP * TN) - (FP * FN)) / ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ^ 0.5
  F1 <- 2* TP / (2 * TP + FP + FN)
  #OUT <- list(TP = TP, FP = FP, TN = TN, FN = FN, TPR = TPR, FPR = FPR, FDR = FDR, PPR = PPR,
  #           NPV = NPV, ACC = ACC, MCC = MCC, F1 = F1)
  OUT <- c(TP, FP, TN, FN, TPR, FPR, FDR, PPR,  NPV,  ACC, F1)
  names(OUT) <- c("TP", "FP", "TN", "FN", "TPR", "FPR", "FDR", "PPR", "NPV", "ACC", "F1")
  
  return(OUT)
}

######## Computation of FPR (= 1-Specificity)
FPR <- function(genes, ref) {
  neg <- setdiff(genes, ref); return(sapply(1:length(genes), function(x) mean(neg %in% genes[1:x])))
}

#########Computation of TPR (=Sensitivity)
TPR <- function(genes, ref) {
  return(sapply(1:length(genes), function(x) mean(ref %in% genes[1:x])))
}
###############Computation of FDRA
FDR <- function(genes, ref) {
  TP <- length(intersect(ref, genes))
  #if(intersect(ref, genes) == NULL) TP = 0
  FP <- length(setdiff(genes, ref))
  FDR <- FP/(TP+FP)
  return(FDR)
}


############### MCDM - TOPSIS Analysis
#' @param decision matrix with rows as methods and columns as decision criteria
#' @param weights weights for each criteria
#' @impact vector of "+" and "-" shows the impact of the criterion on the performance of the methods

TOPSIS <- function (decision, weights, impacts) {
  
  weights <- weights/sum(weights)
  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for (i in 1:nrow(decision)) {
    for (j in 1:ncol(decision)) {
      N[i, j] <- decision[i, j]/sqrt(sum(decision[, j]^2))
    }
  }
  W = diag(weights)
  V = N %*% W
  u <- as.integer(impacts == "+") * apply(V, 2, max) + 
    as.integer(impacts == "-") * apply(V, 2, min)
  l <- as.integer(impacts == "-") * apply(V, 2, max) + 
    as.integer(impacts == "+") * apply(V, 2, min)
  distance_u = function(x) {
    sqrt(sum((x - u)^2))
  }
  distance_l = function(x) {
    sqrt(sum((x - l)^2))
  }
  du <- apply(V, 1, distance_u)
  dl <- apply(V, 1, distance_l)
  score <- dl/(dl + du)
  return(data.frame(alt.row = 1:nrow(decision), score = score, 
                    rank = rank(-score)))
}