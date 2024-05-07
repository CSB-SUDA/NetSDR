#' @title Identify subtype-specific network modules
#' @description Calculate subtype-specific signature proteins based on the proteome and map them
#'     onto the human interactome (combined score >0.4) to identify network modules for the subtype.
#'     The signature proteins are saved in the 'signature.csv' file. The 'edges.txt' file stores
#'     the subtype-specific network, while the 'node_Module.txt' and 'edge_Module.txt' files provide
#'     information on the nodes and edges of robust modules, respectively.
#'
#' @param exprDF A data frame storing expression values without log2 transformation, with rows representing proteins and columns representing samples.
#' @param groupDF A data frame storing subtype information, with the first column being samples and the second column being subtype grouping.
#' @param subtype A vector representing the analysis of a specific subtype.
#'
#' @return This function doesn't return anything, but saves the results to the specified files.
#' @import utils igraph
#' @importFrom stats phyper wilcox.test
#' @export
#'
#' @examples
#' \dontrun{
#'   data(Expr_Group)
#'   getModule(exprDF=Expr,groupDF=Group,subtype="G-IV")
#' }
#'
getModule <- function(exprDF,groupDF,subtype){

  expr <- exprDF
  group <- groupDF
  colnames(group) <- c("Sample","Group")
  egSample = group[group$Group == subtype,"Sample"]
  cgSample = group[group$Group != subtype,"Sample"]

  message("Calculate signature proteins...",appendLF = T)
  pvalue = rep(NA,nrow(expr))
  log2FC = rep(NA,nrow(expr))
  Frequency = rep(0,nrow(expr))
  for (i in 1:nrow(expr)) {

    egExpr = as.numeric(expr[i,egSample])
    cgExpr = as.numeric(expr[i,cgSample])

    if(all(is.na(egExpr)) | all(is.na(cgExpr)) | sum(egExpr,na.rm=TRUE)==0 | sum(cgExpr,na.rm=TRUE)==0){
      next
    }else{
      Test = wilcox.test(egExpr,cgExpr)
      pvalue[i] = Test$p.value
      log2FC[i] = log2(mean(egExpr,na.rm=TRUE))-log2(mean(cgExpr,na.rm=TRUE))
    }
    Frequency[i] = sum(egExpr != 0,na.rm =TRUE)
    if(i%%1000 == 0)cat(paste0(i,"/",nrow(expr),"\r"))
  }
  result = data.frame(protein=row.names(expr),log2FoldChange=log2FC,p.value = pvalue,frequency = Frequency)

  freq.cf =  0.1*length(egSample)
  result$label <- rep("non-significant",nrow(result))
  result$label[result$p.value < 0.05 & result$log2FoldChange > 1 & result$frequency > freq.cf] = "up"
  result$label[result$p.value < 0.05 & result$log2FoldChange < -1 & result$frequency > freq.cf] = "down"

  signature <- subset(result,result$label != "non-significant")
  write.csv(signature,"./signature.csv",row.names = F)

  message("Identify network modules...", appendLF = T)
  data("HPPIN")

  # Map signature proteins to human protein-protein interaction network
  network <- HPPIN[HPPIN$node1%in%signature$protein & HPPIN$node2%in%signature$protein,]

  # Create the network
  edges <- network
  g_net <- graph_from_data_frame(edges, directed=FALSE)

  # Perform community detection using GN and LPA
  set.seed(123)
  GN <- cluster_edge_betweenness(g_net,weights=NULL)
  set.seed(123)
  LP <- cluster_label_prop(g_net,weights=NULL)

  # Associate each node with its corresponding module
  GN_label <- data.frame(node = get.vertex.attribute(g_net)[[1]],
                         module = GN$membership)
  LP_label <- data.frame(node = get.vertex.attribute(g_net)[[1]],
                         module = LP$membership)

  # Extract objects in different clusters
  GN_label_list <- split(GN_label$node,GN_label$module)
  LP_label_list <- split(LP_label$node,LP_label$module)

  # Remove outliers that belong to no cluster
  GN_label_list2<- GN_label_list[(lengths(GN_label_list) >1)]
  LP_label_list2<- LP_label_list[(lengths(LP_label_list) >1)]


  GN_LP <- as.data.frame(matrix(nrow=length(LP_label_list2),
                               ncol=length(GN_label_list2)))
  GN_LP_union <- GN_LP
  GN_LP_phyper <- GN_LP

  # Hypergeometric test between each cluster pairs
  for(i in 1:length(GN_label_list2)){
    for(j in 1:length(LP_label_list2)){
      GN_LP[j,i]=length(intersect(GN_label_list2[[i]],
                                  LP_label_list2[[j]]))
      GN_LP_union[j,i]=length(union(GN_label_list2[[i]],
                                    LP_label_list2[[j]]))
    }
  }
  for(i in 1:ncol(GN_LP_phyper)){
    for(j in 1:nrow(GN_LP_phyper)){
      GN_LP_phyper[j,i]=1-phyper(GN_LP[j,i],
                                 length(GN_label_list2[[i]]),GN_LP_union[j,i],
                                 length(LP_label_list2[[j]]))
    }
  }
  colnames(GN_LP_phyper) <- paste("GN",1:length(GN_label_list2),sep="_")
  rownames(GN_LP_phyper) <- paste("LP",1:length(LP_label_list2),sep="_")

  # Extract significantly correlated clusters
  sig <- data.frame(which(GN_LP_phyper < 0.05, arr.ind = TRUE))
  sig_list <- apply(sig,1,function(x) intersect(GN_label_list2[[x[2]]],LP_label_list2[[x[1]]]))

  # Integrate the nodes included in robustness modules
  modules <- data.frame(node = unlist(sig_list),
                        module = rep(1:length(sig_list),times=as.vector(sapply(sig_list, length))))
  row.names(modules) <- modules$node

  # Get module label of each interaction.
  edges <- edges[,1:2]
  edges$module <- apply(edges, 1, function(x) ifelse((!x[1]%in%modules$node)|(!x[2]%in%modules$node),0,
                                                     ifelse(modules[x[1],2]==modules[x[2],2],modules[x[1],2],0)))
  # Extract interactions of robustness modules
  patternM <- edges[!edges$module==0,]

  # Write the modules and edges data frames to separate output files
  write.table(modules,"node_Module.txt",sep="\t",quote=F,row.names=F,col.names=c("node","module"))
  write.table(patternM,"edge_Module.txt",sep="\t",quote=F,row.names=F,col.names=c("node1","node2","module"))
  write.table(edges,"edges.txt",sep="\t",quote=F,row.names=F,col.names=c("node1","node2","module"))
  message("Module identification succeeded!")

}
