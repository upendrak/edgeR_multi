#!/usr/bin/Rscript

################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Upendra Devisetty
### June 16th, 2016
################################################################################

# Load libraries
library(edgeR)
library(genefilter)
library(devtools)
library(SARTools)
library(getopt)

args<-commandArgs(TRUE)

options<-matrix(c('project',  'pn', 1,  "character",    # project name
                  'author', 'au', 1,  "character",      # author of the statistical analysis/report  
                  'Dir',  'r',  1,  "character",        # path to the directory containing raw counts files
                  'OutDir',  'w',  1,  "character",     # path to the output file 
                  'target', 't',  1,  "character",      # path to the design/target file
                  'features', 'fe', 2,  "character",    # names of the features to be removed (specific HTSeq-count information and rRNA for example)
                  'varInt', 'v',  2,  "character",      # factor of interest
                  'condRef',  'c',  2,  "character",    # reference biological condition
                  'batch',  'b',  2,  "character",      # blocking factor: NULL (default) or "batch" for example
                  'cpmCutoff',  'cc', 2,  "integer",    # counts-per-million cut-off to filter low counts
                  'alpha',  'a',  2,  "double",         # threshold of statistical significance
                  'pAdjust',  'p',  2,  "character",    # p-value adjustment method: "BH" (default) or "BY"
                  'genesele',  'gs', 2,  "character",   # selection of the features in MDSPlot
                  'norm',  'n',  2,  "character",       # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"
                  'colors', 'co', 2,  "character",      # vector of colors of each biological condition on the plots
                  'help',   'h',    0,      "logical"),
                            ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

projectName <- ret.opts$project
author  <-  ret.opts$author
targetFile <- ret.opts$target
rawDir <- ret.opts$Dir
outDir <- ret.opts$OutDir
features <- ret.opts$features
varInt  <- ret.opts$varInt
condRef <- ret.opts$condRef
batch <- ret.opts$batch
fitType <- ret.opts$filtType
gene.selection <- ret.opts$genesele
cpmCutoff <- ret.opts$cpmCutoff
normalizationMethod <- ret.opts$norm
alpha <- ret.opts$alpha
pAdjust <- ret.opts$pAdjust
colors <- ret.opts$colors

# checking parameters
checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

# loading target file
# target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
# counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
# majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# edgeR analysis
# out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
#                        batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
#                        pAdjustMethod=pAdjustMethod)

# MDS + clustering
# exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
# summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)


# generating HTML report
# writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
#                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
#                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
#                   condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
#                   gene.selection=gene.selection, normalizationMethod=normalizationMethod)
