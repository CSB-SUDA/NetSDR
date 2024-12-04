#' @title Predict the binding affinity between drugs and proteins
#' @description The function uses the DeepPurpose tool to predict the binding affinity between drug-protein pairs,
#'     and the results are saved in the "result/virtual_screening.txt" file.
#'
#' @param smiles_file The information of simplified molecular input line entry system format of the drug, where the first column is the name or ID of the drug,
#'     and the second column is its smiles symbol.
#' @param seq_file The amino acid sequences of proteins, with the first column representing the protein symbol and the second column representing the sequence.
#' @param DRN_file Drug-protein interaction data, contains two columns: the first column is "protein" (symbol), and the second column is "drug" (drug name).
#' @param py_env The virtual environment path for python installed using conda.
#' @param pretrained_model The DeepPurpose model type used, which defaults to 'MPNN_CNN_BindingDB_IC50'.
#'
#' @return NULL
#' @import utils reticulate
#' @export
#'
#' @examples
#' \dontrun{
#'   smiles_file <- "inst/extdata/smiles.txt"
#'   seq_file <- "inst/extdata/sequences.txt"
#'   DRN_file <- "inst/extdata/DRN.txt"
#'   py_env <- "D:/Software/anaconda3/envs/py3.9"
#'   getAffinity(smiles_file,seq_file,DRN_file,py_env)
#' }
#'
getAffinity <- function(smiles_file,seq_file,DRN_file,py_env,pretrained_model='MPNN_CNN_BindingDB_IC50'){
  #Drug smiles
  smiles <- read.table(smiles_file,header = T,sep = "\t",comment.char = "")
  colnames(smiles) <- c("drug","SMILES")
  row.names(smiles) <- smiles$drug
  #Protein sequences
  seqs <- read.table(seq_file,header = T,sep = "\t",comment.char = "")
  colnames(seqs) <- c("protein","Sequence")
  row.names(seqs) <- seqs$protein
  #Drug-Protein Interactions
  DPI <- read.table(DRN_file,header = T,sep = "\t")
  colnames(DPI) <- c("protein","drug")

  out_dir <- "NetSDR_results/PS"
  dir.create(out_dir,recursive = T)
  setwd(out_dir)

  # Integrate SMILES and Sequence to DPI
  DPI$Sequence <- apply(DPI, 1, function(p) seqs[p[1],"Sequence"])
  DPI$SMILES <- apply(DPI, 1, function(p) smiles[p[2],"SMILES"])
  DPI_info <- DPI[,c("drug","protein","Sequence","SMILES")]
  DPI_info <- DPI_info[DPI_info$SMILES !="",]
  colnames(DPI_info) <- c("drug_name","target_name","target_seq","drug_smiles")
  write.table(DPI_info,"./DPI_info.txt",sep = "\t",quote = F,row.names = F)

  #-----------
  # Activate the virtual environment for DeepPurpose
  use_condaenv(py_env)
  py_run_string("from DeepPurpose import DTI as models")

  # Import DPI information
  df = read.table("DPI_info.txt",header = T, sep="\t",comment.char = "")
  target_seq = as.list(df$target_seq)
  drug_smiles =  as.list(df$drug_smiles)
  target_name =  as.list(df$target_name)
  drug_name =  as.list(df$drug_name)
  # Import pre-trained model
  model = py$models$model_pretrained(model = pretrained_model)
  # Predict binding affinity
  my_predict = py$models$virtual_screening(drug_smiles, target_seq, model, drug_name, target_name,verbose=FALSE)

  setwd("../..")
}
