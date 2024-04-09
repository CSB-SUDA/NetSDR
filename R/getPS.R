#' @title Calculat drug score(ps) using PRS methods for DRN.
#' @description This function uses the PRS algorithm to calculate the ps of DPIs and saves it in the "prs_dti_score.csv" file.
#'
#' @param BA_file The path to the predicted binding affinity file.
#' @param PPIN_file The path to the PPI network file (or to a certain module).
#' @param condaenv The virtual environment path for DeepPurpose installed using conda.To maintain consistency in Python, we use the Python version in the DeepPurpose environment.
#' @param path The path to the enm package. Its code is downloaded from https://github.com/oacar/enm_package, which can be placed in any path.
#'
#' @return A data frame contains the ps of DPI.
#' @import utils reticulate
#' @export
#'
#' @examples
#' \dontrun{
#'   getPS(BA_file="data/virtual_screening.txt",PPIN_file="data/edge.txt",
#'     condaenv="C:/Users/A/.conda/envs/DeepPurpose",path="D:/enm_package-master")
#' }
#'
getPS <- function(BA_file,PPIN_file,condaenv,path){

  use_condaenv(condaenv)
  py_run_string("import os")
  py_run_string("import numpy as np")
  py_run_string("import networkx as nx")
  py_run_string("import pandas as pd")
  py_run_string("import sys")
  py_run_string(paste0("sys.path.append(r'",path,"')"))
  py_run_string("from enm.Enm import *")

  enm = py$Enm('PPIN')
  enm$read_network(PPIN_file, sep='\t')
  enm$gnm_analysis(normalized=FALSE)
  enm$cluster_matrix(enm$prs_mat)
  write.csv(enm$df,"pcc_df.csv",row.names = T)
  write.table(enm$prs_mat,"prs_df.txt",quote = F)
  write.table(enm$prs_mat_df,"prs_mat_df.txt",quote = F)

  BA <- read.table(BA_file,skip = 3,sep = "|",comment.char ="+")
  BA <- BA[,2:5]
  colnames(BA) <- c("Rank", "Drug.Name", "Target.Name", "Binding.Score")
  BA$Drug.Name <- trimws(BA$Drug.Name)
  BA$Target.Name <- trimws(BA$Target.Name)

  Sens <- read.csv("pcc_df.csv",row.names = 1)
  colnames(Sens)[1] <- "Target.Name"

  ps <- merge(BA,Sens,by="Target.Name")
  ps <- ps[,c('Drug.Name', 'Target.Name', 'Binding.Score', 'sens')]
  ps$score <- ps$Binding.Score*ps$sens
  ps <- ps[order(ps$score,decreasing = T),]
  write.csv(ps,"prs_dti_score.csv",row.names = F)

  return(ps)

}
