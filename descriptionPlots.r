#' Description plots of the counts
#'
#' Description plots of the counts according to the biological condition
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors for the plots (one per biological condition)
#' @return PNG files in the "figures" directory and the matrix of the most expressed sequences
#' @author Hugo Varet

descriptionPlots <- function(counts, n=3, group=target[,varInt], OutDir, col){
  # create the figures directory if does not exist
  dir.create(paste(OutDir,"figures",sep="/"), showWarnings=FALSE)
  
  source("/barplotTotal.R")

  # total number of reads per sample
  barplotTotal(counts, group=target[,varInt], OutDir, col)
  
  source("/barplotNull.R")

  # percentage of null counts per sample
  barplotNull(counts, group=target[,varInt], OutDir, col)
  
  source("/densityPlot.R")

  # distribution of counts per sample
  densityPlot(counts, group=target[,varInt], OutDir, col)
  
  source("/majSequences.R")

  # features which catch the most important number of reads
  majSequences <- majSequences(counts,  n=3, group=target[,varInt], OutDir, col)
  
  source("/pairwiseScatterPlots.R")

  # SERE and pairwise scatter plots
  cat("Matrix of SERE statistics:\n")
  print(tabSERE(counts))
  pairwiseScatterPlots(counts, group=target[,varInt], OutDir)
  
  return(majSequences)
}
