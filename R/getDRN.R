#' @title Construct the drug response network
#' @description This function utilizes protein expression within modules and predicted drug sensitivity scores based on
#'    module expression profiles to compare scores between two groups with high and low expression of a particular protein.
#'    If there is a significant difference (p < 0.05), it suggests a potential interaction between the drug and the protein.
#'
#' @param expr_file Expression profile of the module, with rows representing proteins and columns representing samples.
#'
#' @return This function doesn't return anything, but saves the results to the "NetSDR_results/DRN" files.
#' @export
#' @import utils oncoPredict
#' @importFrom stats median wilcox.test
#' @examples
#' \dontrun{
#'   expr_file <- "inst/extdata/expression_module.txt"
#'   getDRN(expr_file)
#' }
#'

getDRN <- function(expr_file){

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
                removeLowVaryingGenes = 0,
                minNumSamples = 10,
                printOutput = TRUE,
                removeLowVaringGenesFrom = 'rawData',
                cc = TRUE)

  moduleData <- data.frame(t(testData))
  moduleLabel <- apply(moduleData, 2, function(x) {ifelse(x > median(x), "high", "low")})

  drugData <- read.csv("./calcPhenotype_Output/DrugPredictions.csv",row.names = 1,check.names = F)
  dat <- merge(moduleLabel,drugData,by="row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]

  protein <- colnames(moduleData)
  drug <- colnames(drugData)

  proteins <- c()
  drugs <- c()
  pvalue <- c()
  for (pro in protein) {
    for (dg in drug) {
      proteins <- c(proteins,pro)
      drugs <- c(drugs,dg)

      temp = dat[,c(pro,dg)]
      dat.high = temp[temp[,1] == "high",2]
      dat.low = temp[temp[,1] != "high",2]
      if(all(is.na(dat.high))|all(is.na(dat.low))){
        pvalue <- c(pvalue,NA)
      }else{
        wilcox <- wilcox.test(dat.high,dat.low)
        pvalue <- c(pvalue,wilcox$p.value)
      }
    }
  }

  DPI_all <- data.frame("protein"=proteins,"drug"=drugs,"pvalue"=pvalue)
  DPI_all$label <- ifelse(DPI_all$pvalue < 0.05,"sig","non-sig")
  DRN <- subset(DPI_all,DPI_all$pvalue < 0.05)
  write.table(DRN,"DRN.txt",sep = "\t",quote = F,row.names = F)
  DRN.info <- data.frame(node = c(protein,drug),type = c(rep("protein",length(protein)),rep("drug",length(drug))))
  write.table(DRN.info,"DRN_info.txt",sep = "\t",quote = F,row.names = F)

  setwd("../..")

}
