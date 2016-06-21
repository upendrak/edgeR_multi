#' Summarize edgeR analysis
#'
#' Summarize edgeR analysis: diagnotic plots, dispersions plot, summary of the independent filtering, export results...
#'
#' @param out.edgeR the result of \code{run.edgeR()}
#' @param group factor vector of the condition from which each sample belongs
#' @param counts matrix of raw counts
#' @param alpha significance threshold to apply to the adjusted p-values
#' @param col colors for the plots
#' @return A list containing: (i) a list of \code{data.frames} from \code{exportResults.edgeR()} and (ii) a table summarizing the number of differentially expressed features
#' @author Hugo Varet

summarizeResults.edgeR <- function(out.edgeR, group, counts, alpha, col){  
  
  # boxplots before and after normalisation
  source("/Users/upendrakumardevisetty/Documents/git_repos/edgeR_multifactorial/countsBoxplots.R")
  countsBoxplots(out.edgeR$dge, group, col, OutDir)

  # dispersions
  source("/Users/upendrakumardevisetty/Documents/git_repos/edgeR_multifactorial/BCVPlot.R")
  BCVPlot(dge=out.edgeR$dge, OutDir)
  
  # exporting results of the differential analysis
  source("/Users/upendrakumardevisetty/Documents/git_repos/edgeR_multifactorial/exportResults.edgeR.R")
  complete <- exportResults.edgeR(out.edgeR, group, counts, alpha, OutDir)

  # small table with number of differentially expressed features
  nDiffTotal <- nDiffTotal(complete, alpha)
  cat("Number of features down/up and total:\n")
  print(nDiffTotal, quote=FALSE)
  
  # histograms of raw p-values
  source("/Users/upendrakumardevisetty/Documents/git_repos/edgeR_multifactorial/rawpHist.R")
  rawpHist(complete, OutDir)
  
  # MA-plots
  source("/Users/upendrakumardevisetty/Documents/git_repos/edgeR_multifactorial/MAPlot.R")
  MAPlot(complete, alpha, OutDir)
  
  # Volcano plots
  source("/Users/upendrakumardevisetty/Documents/git_repos/edgeR_multifactorial/volcanoPlot.r")
  volcanoPlot(complete, alpha, OutDir)
  
  return(list(complete, nDiffTotal))
}