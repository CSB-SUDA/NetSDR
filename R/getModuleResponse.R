#' @title Calculate the drug response levels of the modules
#' @description This function predicts the clinical drug response of patients with a specific subtype
#'     based on the expression profile of the module, and calculates the number of drugs whose
#'     sensitivity score is affected by the expression of the module.
#'
#' @param expr_file Expression profile of the module, with rows representing proteins and columns representing samples.
#'
#' @return This function doesn't return anything, but saves the results to the "NetSDR_results/DRN" files.
#' @import utils oncoPredict ConsensusClusterPlus
#' @importFrom stats phyper wilcox.test
#' @export
#'
#' @examples
#' \dontrun{
#'   expr_file <- "inst/extdata/expression_module.txt"
#'   getModuleResponse(expr_file)
#' }
#'

getModuleResponse <- function(expr_file){

  data("GDSC2_Expr")
  data("GDSC2_Res")
  testData <- read.table(expr_file,header=T,sep = "\t",row.names = 1)

  out_dir <- "NetSDR_results/DRN"
  dir.create(out_dir,recursive = T)
  setwd(out_dir)

  # Predict patients' clinical responses based on module expressions.
  testExpr <- as.matrix(testData)
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
  moduleData <- testData[,colSums(testData) != 0 ]
  moduleData <- as.matrix(moduleData)
  ConsensusClusterPlus(moduleData, maxK = 3,
                       reps = 1000, pItem = 0.8,
                       pFeature = 1,
                       seed = 10000,
                       clusterAlg = "pam",
                       distance = "pearson",
                       title = "./sample.cluster",
                       plot = "png",
                       writeTable=TRUE)

  # Compare the drug sensitivity scores of the two groups of patients.
  drugData <- read.csv("./calcPhenotype_Output/DrugPredictions.csv",check.names = F)
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

  message("Clinical response assessment completed!")
  setwd("../..")


}
