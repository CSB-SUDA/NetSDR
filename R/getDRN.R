#' @title Construct the drug response network
#' @description This function utilizes protein expression within modules and predicted drug sensitivity scores based on
#'    module expression profiles to compare scores between two groups with high and low expression of a particular protein.
#'    If there is a significant difference (p < 0.05), it suggests a potential interaction between the drug and the protein.
#'
#' @param moduleDF A data frame storing expression values of the module, with rows representing proteins and columns representing samples.
#' @param drugDF A data frame storing drug sensitivity scores predicted by the module, with rows representing samples and columns representing drugs.
#'
#' @return a list containing the drug response network and its node information
#' @export
#' @import utils
#' @importFrom stats median wilcox.test
#' @examples
#' \dontrun{
#'   drugDF <- read.csv("calcPhenotype_Output/DrugPredictions.csv",row.names = 1,check.names = F)
#'   getDRN(moduleDF,drugDF)
#' }
#'
getDRN <- function(moduleDF,drugDF){

  moduleData <- data.frame(t(moduleDF))
  moduleLabel <- apply(moduleData, 2, function(x) {ifelse(x > median(x), "high", "low")})
  drugData <- drugDF
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
  DRN.info <- data.frame(node = c(protein,drug),type = c(rep("protein",length(protein)),rep("drug",length(drug))))
  DRN_list <- list(DRN=DRN,DRN.info=DRN.info)

  write.table(DRN,"DRN.txt",sep = "\t",quote = F,row.names = F)
  write.table(DRN.info,"DRN_info.txt",sep = "\t",quote = F,row.names = F)

  return(DRN_list)

}
