#' Explore counts structure
#'
#' Explore counts structure: MDS (edgeR) and clustering
#'
#' @param group factor vector of the condition from which each sample belongs
#' @param gene.selection selection of the features in MDSPlot (\code{"pairwise"} by default)
#' @param col colors used for the PCA/MDS (one per biological condition)
#' @return A list containing the dds object and the results object
#' @author Hugo Varet

exploreCounts <- function(object, group, gene.selection, col, OutDir){

	# Clusterplot
	clusterPlot(group=target[,varInt], OutDir)
	
	# MDS plot
	MDSPlot(group=target[,varInt], gene.selection, col, OutDir) 
  
} 