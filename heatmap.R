Heatmap <- function(OutDir) 
{
		object=out.edgeR$dge

		# normalization
		dge <- calcNormFactors(object, method=normalizationMethod)
		
		# countspermillion 
		countspermi <- cpm(dge, normalized.lib.sizes=TRUE)

		# Now pick the genes with the top variance over all samples:
		rv <- rowVars(countspermi)
		idx <- order(-rv)[1:20]

		nd = paste(OutDir,"figures",sep="/")

		# Plotting
  		png(filename=paste(nd,"heatmap.png",sep="/"),width=1800,height=1800,res=300)
		heatmap(countspermi[idx,], Colv = NA)
		title("Heatmap", line = 3)
		dev.off()
}