#' @title Predict the binding affinity between drugs and proteins
#' @description The function uses the DeepPurpose tool to predict the binding affinity between drug-protein pairs,
#'     and the results are saved in the "result/virtual_screening.txt" file.
#'
#' @param smilesDF A data frame with "drug" column as the drug name and "SMILES" column as the simplified molecular input line entry system format of the drug.
#' @param seqDF A data frame with "Protein" column as symbols of proteins and "Sequence" column as the amino acid sequences of the proteins.
#' @param DRN A data frame contains two columns: "drug"(drug name) and "protein"(symbol).
#' @param condaenv The virtual environment path for DeepPurpose installed using conda.
#' @param path The path where the DeepPurpose package is located.
#' @param pretrained_model The DeepPurpose model type used, which defaults to 'MPNN_CNN_BindingDB_IC50'.
#'
#' @return NULL
#' @import utils reticulate
#' @export
#'
#' @examples
#' \dontrun{
#'   data(drug_smiles)
#'   data(protein_sequences)
#'   data(DRN)
#'   getAffinity(smiles,seqs,DRN,condaenv="C:/Users/A/.conda/envs/DeepPurpose",
#'     path="C:/Users/A/DeepPurpose")
#' }
#'
getAffinity <- function(smilesDF,seqDF,DRN,condaenv,path,pretrained_model='MPNN_CNN_BindingDB_IC50'){

  smiles <- smilesDF
  seqs <- seqDF
  DPI <- DRN

  # Integrate SMILES and Sequence to DPI
  DPI$Sequence <- apply(DPI, 1, function(p) seqs[p[1],"Sequence"])
  DPI$SMILES <- apply(DPI, 1, function(p) smiles[p[2],"SMILES"])
  DPI_info <- DPI[,c("drug","protein","Sequence","SMILES")]
  DPI_info <- DPI_info[DPI_info$SMILES !="",]
  colnames(DPI_info) <- c("drug_name","target_name","target_seq","drug_smiles")
  write.table(DPI_info,"DPI_info.txt",sep = "\t",quote = F,row.names = F)

  #-----------
  # Activate the virtual environment for DeepPurpose
  use_condaenv(condaenv)
  py_run_string("import sys")
  py_run_string(paste0("sys.path.append(r'",path,"')"))
  py_run_string("from DeepPurpose import DTI as models")
  py_run_string("import pandas as pd")

  # Import DPI information
  df = py$pd$read_table("./DPI_info.txt", sep="\t")
  target_seq = as.list(df$target_seq)
  drug_smiles =  as.list(df$drug_smiles)
  target_name =  as.list(df$target_name)
  drug_name =  as.list(df$drug_name)
  # Import pre-trained model
  model = py$models$model_pretrained(model = pretrained_model)
  # Predict binding affinity
  my_predict = py$models$virtual_screening(drug_smiles, target_seq, model, drug_name, target_name,verbose=FALSE)

}
