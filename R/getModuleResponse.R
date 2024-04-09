#' @title Calculate the drug response levels of the modules
#' @description This function predicts the clinical drug response of patients with a specific subtype
#'     based on the expression profile of the module, and calculates the number of drugs whose
#'     sensitivity score is affected by the expression of the module.
#'
#' @param moduleDF A data frame storing expression values of the module, with rows representing proteins and columns representing samples.
#'
#' @return The number of significantly affected and unaffected drugs.
#' @import utils oncoPredict ConsensusClusterPlus
#' @importFrom stats phyper wilcox.test
#' @export
#'
#' @examples
#' \dontrun{
#'   sig_count <- getModuleResponse(moduleDF)
#' }
#'
getModuleResponse <- function(moduleDF){

  data(GDSC2_Expr)
  data(GDSC2_Res)

  # Predict patients' clinical responses based on module expressions.
  testExpr<- as.matrix(moduleDF)
  calcPhenotype(trainingExprData = GDSC2_Expr,
                trainingPtype = GDSC2_Res,
                testExprData = testExpr,
                batchCorrect = 'eb',
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 10,
                printOutput = TRUE,
                removeLowVaringGenesFrom = 'rawData',
                cc = TRUE)

  # Divide the patients into two groups.
  moduleData <- moduleDF[,colSums(moduleDF) != 0 ]
  moduleData <- as.matrix(moduleData)
  ConsensusClusterPlus(moduleData, maxK = 3,
                       reps = 1000, pItem = 0.8,
                       pFeature = 1,
                       seed = 10000,
                       clusterAlg = "pam",
                       distance = "pearson",
                       title = "sample.cluster",
                       plot = "png",
                       writeTable=TRUE)

  # Compare the drug sensitivity scores of the two groups of patients.
  drugData <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
  colnames(drugData)[1] <- "Sample"
  cluster <- read.csv("./sample.cluster/sample.cluster.k=2.consensusClass.csv",header = FALSE)
  colnames(cluster) <- c("Sample","Cluster")

  dat <- merge(cluster,drugData,by="Sample")
  drug <- colnames(dat)[3:ncol(dat)]

  pvalue <- c()
  label <- c()
  for(i in 3:ncol(dat)){
    #i=3
    dat.cluster1 <- dat[dat$Cluster== 1,i]
    dat.cluster2 <- dat[dat$Cluster== 2,i]
    wilcox <- wilcox.test(dat.cluster1,dat.cluster2)
    pvalue <- c(pvalue,wilcox$p.value)
    label <- c(label,ifelse(wilcox$p.value < 0.05, "sig","non-sig"))
  }
  drugSig <- data.frame("drug"=drug,"pvalue"=pvalue,"label" = label)

  countSig <- c(sum(drugSig$label=="non-sig"),sum(drugSig$label !="non-sig"))
  sigCount <- data.frame(label = c("non-sig","sig"),count = countSig)
  write.csv(sigCount,"./drug_response_level.csv", row.names = F,quote = F)
  return(sigCount)

}
