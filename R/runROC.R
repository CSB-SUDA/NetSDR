#' @title Perform ROC analysis
#' @description This function performs ROC validation on the predicted DPI. Extract DPIs greater than or equal to a certain ps in sequence,
#'     allocate 1 (known) and 0 (unknown) according to the drug-target interactions of TTD and DrugBank, and perform ROC analysis.
#'
#' @param pred_DPI A data frame contains DPIs and their ps, where the first column represents drugs, the second column represents proteins, and the third column represents the ps value.
#'
#' @return A data frame storing AUC from different cut-off.
#' @import utils pROC
#' @export
#'
#' @examples
#' \dontrun{
#'   runROC(pred_DPI)
#' }
#'
runROC <- function(pred_DPI){

  # Load known drug-target interaction information
  data("DTI_TTD")
  data("DTI_DrugBank")
  drug2target <- rbind(TTD,DrugBank)
  drug2target$DrugName <- toupper(drug2target$DrugName)
  drug2target <- unique(drug2target)

  # Load predicted DPI
  pred_dpi <- pred_DPI
  colnames(pred_dpi)[1:3] <- c("Drug.Name","Target.Name","score")
  pred_dpi$Drug.Name <- substring(pred_dpi$Drug.Name,1,nchar(pred_dpi$Drug.Name)-5)
  pred_dpi$Drug.Name <- toupper(pred_dpi$Drug.Name)

  cutoff_ls <- c()
  auc_ls <- c()
  drug_ls <- c()
  dpi_ls <- c()
  seq_ls <- pred_dpi$score

  for(cutoff in seq_ls){
    temp <- pred_dpi[pred_dpi$score>=cutoff,]
    temp$type <- ifelse(temp$Drug.Name%in%drug2target$DrugName & temp$Target.Name%in%drug2target$TargetName,1,0)
    if(length(unique(temp$type))==2){
      cutoff_ls <- c(cutoff_ls,cutoff)
      roc1 <- roc(temp$type, temp$score)
      auc_ls <- c(auc_ls,auc(roc1))
      drug_ls <- c(drug_ls,length(unique(temp$Drug.Name)))
      dpi_ls <- c(dpi_ls,nrow(temp))
    }

  }
  df <- data.frame(cutoff=cutoff_ls,auc=auc_ls,drugCount=drug_ls,dpiCount=dpi_ls)
  return(df)
}
