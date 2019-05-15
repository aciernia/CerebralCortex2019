#Co-expression similarity and adjacency
#signed, biweight midcorrelation between all genes 

#library("WGCNA")
library(WGCNA)
allowWGCNAThreads(nThreads = 24)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#load data
load("Zhong2018_WGCNA_log2RPKMInputdata.Rdata")

#=====================================================================================
bwnet = blockwiseModules(PFCdatExpr, maxBlockSize = 30000,#by setting max block size this high, forces one block
                         power = 15, TOMType = "signed", 
                         minModuleSize = 20,
                         networkType = "signed",
                         corType = "bicor", #biweight midcorrelation
                         maxPOutliers = 0.05, #forces bicor to never regard more than the specified proportion of samples as outliers.
                         reassignThreshold = 0,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Zhong_TOM",
                         nThreads = 24,
                         verbose = 3)

save(bwnet, file = "NetworkConstruction_Zhong2018.RData")


