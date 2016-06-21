#' BCV plot (for edgeR dispersions)
#'
#' Biological Coefficient of Variation plot (for edgeR objects)
#'
#' @param dge a \code{DGEList} object
#' @param outfile TRUE to export the figure in a png file
#' @return A file named BCV.png in the figures directory with a BCV plot produced by the \code{plotBCV()} function of the edgeR package
#' @author Marie-Agnes Dillies and Hugo Varet

BCVPlot <- function(dge=out.edgeR$dge, OutDir){	
  nd = paste(OutDir,"figures",sep="/")

  png(filename=paste(nd,"BCV.png",sep="/"), width=1800, height=1800, res=300)
  
  plotBCV(dge, las = 1, main = "BCV plot")
  
  dev.off()	
}
