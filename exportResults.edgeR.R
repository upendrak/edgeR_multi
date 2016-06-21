#' Export results for edgeR analyses
#'
#' Export counts and edgeR results
#'
#' @param out.edgeR the result of \code{run.edgeR()}
#' @param group factor vector of the condition from which each sample belongs
#' @param counts non-filtered counts (used to keep them in the final table)
#' @param alpha threshold to apply to adjusted p-values
#' @return A list of \code{data.frame} containing counts, pvalues, FDR, log2FC...
#' @details \code{counts} are used as input just in order to export features with null counts too.
#' @author Marie-Agnes Dillies and Hugo Varet

exportResults.edgeR <- function(out.edgeR, group=target[,varInt], counts, alpha=0.05, OutDir){
  
  dge <- out.edgeR$dge
  res <- out.edgeR$results
  
  # comptages bruts, normalisés et baseMean
  tmm <- dge$samples$norm.factors
  N <- colSums(dge$counts)
  f <- tmm * N/mean(tmm * N)
  normCounts <- round(scale(dge$counts, center=FALSE, scale=f))
  base <- data.frame(Id=rownames(counts), counts)
  names(base) <- c("Id", colnames(counts))
  norm.bm <- data.frame(Id=rownames(normCounts),normCounts)
  names(norm.bm) <- c("Id", paste0("norm.",colnames(normCounts)))
  norm.bm$baseMean <- round(apply(scale(dge$counts, center=FALSE, scale=f),1,mean),2)
  for (cond in levels(group)){
    norm.bm[,cond] <- round(apply(as.data.frame(normCounts[,group==cond]),1,mean),0)
  }
  base <- merge(base,norm.bm,by="Id",all=TRUE)

  complete <- list()
  for (name in names(res)){
    complete.name <- base

    # ajout d'elements depuis res
    res.name <- data.frame(Id=rownames(res[[name]]),FC=round(2^(res[[name]][,"logFC"]),3),
                           log2FoldChange=round(res[[name]][,"logFC"],3),pvalue=res[[name]][,"PValue"],
						   padj=res[[name]][,ifelse("FDR" %in% names(res[[name]]), "FDR", "FWER")])
    complete.name <- merge(complete.name, res.name, by="Id", all=TRUE)
    
    # ajout d'elements depuis dge
    dge.add <- data.frame(Id=rownames(dge$counts),tagwise.dispersion=round(dge$tagwise.dispersion,4),
                          trended.dispersion=round(dge$trended.dispersion,4))
    complete.name <- merge(complete.name, dge.add, by="Id", all=TRUE)
    complete[[name]] <- complete.name
	
  # sélection up and down
	up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$log2FoldChange>=0),]
	up.name <- up.name[order(up.name$padj),]
	down.name <- complete.name[which(complete.name$padj <= alpha & complete.name$log2FoldChange<=0),]
	down.name <- down.name[order(down.name$padj),]		

	# exports
	name <- gsub("_","",name)

  dir.create(paste(OutDir,"tables",sep="/"), showWarnings=FALSE)
  nd = paste(OutDir,"tables",sep="/")

  file1 = paste0(name,".complete.txt")
  file2=paste(nd, file1 , sep="/")
  write.table(complete.name, file=file2, sep="\t", row.names=FALSE, dec=".", quote=FALSE)

  file3 = paste0(name,".up.txt")
  file4=paste(nd, file3 , sep="/")
  write.table(up.name, file=file4, row.names=FALSE, sep="\t", dec=".", quote=FALSE)

  file5 = paste0(name,".down.txt")
  file6=paste(nd, file5 , sep="/")
  write.table(down.name, file=file6, row.names=FALSE, sep="\t", dec=".", quote=FALSE)

  }

  return(complete)
}
