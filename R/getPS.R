#' @title Calculat drug score(ps) using PRS methods for DRN.
#' @description This function uses the PRS algorithm to calculate the ps of DPIs and saves it in the "prs_dti_score.csv" file.
#'
#' @param BA_file The path to the predicted binding affinity file.
#' @param PPIN_file The path to the PPI network file (or to a certain module).
#' @param py_env The virtual environment path for python installed using conda.
#'
#' @return A data frame contains the ps of DPI.
#' @import utils reticulate
#' @export
#'
#' @examples
#' \dontrun{
#'   py_env <- "D:/Software/anaconda3/envs/py3.9"
#'   BA_file <- "inst/extdata/virtual_screening.txt"
#'   PPIN_file <- "inst/extdata/edge_module.txt"
#'   getPS(BA_file,PPIN_file,py_env)
#' }
#'
getPS <- function(BA_file,PPIN_file,py_env){

  out_dir <- "NetSDR_results/PS"
  dir.create(out_dir,recursive = T)

  use_condaenv(py_env)
  py_run_string("import os")
  py_run_string("import numpy as np")
  py_run_string("import networkx as nx")
  py_run_string("import pandas as pd")
  py_run_string("import sys")

  enm_path <- paste0(system.file(package = "NetSDR"),"/python")
  py_run_string(paste0("sys.path.append(r'",enm_path,"')"))
  py_run_string("from enm.Enm import *")

  enm = py$Enm('PPIN')
  enm$read_network(PPIN_file, sep='\t')
  enm$gnm_analysis(normalized=FALSE)
  enm$cluster_matrix(enm$prs_mat)
  write.csv(enm$df,paste0(out_dir,"/pcc_df.csv"),row.names = T)
  write.table(enm$prs_mat_df,paste0(out_dir,"/prs_mat_df.txt"),quote = F)

  BA <- read.table(BA_file,skip = 3,sep = "|",comment.char ="+")
  BA <- BA[,2:5]
  colnames(BA) <- c("Rank", "Drug.Name", "Target.Name", "Binding.Score")
  BA$Drug.Name <- trimws(BA$Drug.Name)
  BA$Target.Name <- trimws(BA$Target.Name)

  Sens <- read.csv(paste0(out_dir,"/pcc_df.csv"),row.names = 1)
  colnames(Sens)[1] <- "Target.Name"

  ps <- merge(BA,Sens,by="Target.Name")
  ps <- ps[,c('Drug.Name', 'Target.Name', 'Binding.Score', 'sens')]
  ps$ps <- ps$Binding.Score*ps$sens
  ps <- ps[order(ps$ps,decreasing = T),]
  write.csv(ps,paste0(out_dir,"/prs_dti_score.csv"),row.names = F)

  return(ps)

}
