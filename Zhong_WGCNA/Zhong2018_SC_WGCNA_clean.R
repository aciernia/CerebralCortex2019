#clean single cell human gene expression matrix

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gplots)
library(biomaRt)
library(nlme)
library(lsmeans)
library(data.table)
library(readxl)

##############################################################################################################
#read in count matrix after UMI filter, etc.
#process to RPKM for each cell 
#average across each cell type (as identified from each paper by clustering)
##############################################################################################################

##############################################################################################################
#Zhong 2018
##############################################################################################################
library(readxl)

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/OverlapLists/HumanSCdata/Zhong2018")
#from GEO and convert to .csv on command line: ssconvert GSE104276_all_pfc_2394_UMI_count_NOERCC.xls GSE104276_all_pfc_2394_UMI_count_NOERCC.csv
#https://stackoverflow.com/questions/10557360/convert-xlsx-to-csv-in-linux-command-line

counts <- read.table("GSE104276_all_pfc_2394_UMI_count_NOERCC_mod.txt")
barcodes <- read.csv("GSE104276_readme_sample_barcode.csv")

counts_df <- counts
genes <- unique(rownames(counts_df))


#count total reads per cell:
Count_Totals <- colSums(counts_df)

#Total Counts/1,000,000
Counts_per_million <- Count_Totals/1000000

#Reads per million: 
#mapply(`/`, data.frame(a), b)

RPM <- mapply(`/`, data.frame(counts_df), Counts_per_million)

RPM <- as.data.frame(RPM)

RPM$gene_ids <- rownames(counts_df)

HumanHg38 <- read.csv("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/8_2018_4CpGDMRs/genelist_lengthcorrOverlaps/Ensemble_hg38genes.csv")

df <- merge(RPM,HumanHg38, by.x = "gene_ids",by.y = "hgnc_symbol")
df$length <- abs(df$end_position - df$start_position)

#divide RPM data by length:
RPKM <- df[,2:(ncol(df)-8)]/df$length

df_RPKM <- cbind(df[,c(1,2399:ncol(df))],RPKM)

df_RPKM2 <- df_RPKM %>% gather(Cells_name, df_RPKM, 7:ncol(df_RPKM))

#add in cell data:
cellinfo <- barcodes[,1:3]


df_RPKM3 <- merge(df_RPKM2,cellinfo, by = "Cells_name", all.x=T)
df_RPKM4 <- df_RPKM3 %>% separate(Cells_name, into = c("GW","Sample","cell"),remove =F)
df_RPKM4$gestationweek_celltype <- paste(df_RPKM4$week, df_RPKM4$cell_types)
df_RPKM4$gestationweek_celltype <- factor(df_RPKM4$gestationweek_celltype)


#average for each gene for each cell type:
RPKM_mean <- df_RPKM4 %>% group_by(ensembl_gene_id,gestationweek_celltype) %>% 
  summarize(meanRPKM = mean(df_RPKM)) 

RPKM_Samplemean <- df_RPKM4 %>% group_by(Sample,ensembl_gene_id,gestationweek_celltype) %>% 
  dplyr::summarize(meanRPKM = mean(df_RPKM)) %>%
 spread(gestationweek_celltype, meanRPKM)


write.csv(RPKM_mean,"Zhong2018_RPKMmeans_byAgeCellType.csv")
write.csv(RPKM_Samplemean,"Zhong2018_RPKM_Samplemeans_byAgeCellType.csv")

#average for each gene for each cell type and log2 tranform
log2RPKM_mean <- df_RPKM3 %>% group_by(ensembl_gene_id,gene_ids, gestationweek_celltype) %>% summarize(meanRPKM = mean(df_RPKM)) %>%
  mutate(log2RPKMmean = log2(meanRPKM+1)) %>% dplyr::select(-meanRPKM) %>%
  spread(gestationweek_celltype, log2RPKMmean)

write.csv(log2RPKM_mean,"Zhong2018_log2meanRPKM+1_byAgeCellType.csv")

save(df_RPKM3,RPKM_mean,log2RPKM_mean, file = "Zhong2018_RPKM_bycelltype.rda")
load("Zhong2018_RPKM_bycelltype.rda")
################################################################################
#load
setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/ZhongSC_WGCNA")
load("Zhong2018_RPKM_bycelltype.rda")


####################################################################################
#WGCNA
####################################################################################

library(WGCNA);
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=12
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#make gene matrix with genes as columns and rows as samples
PFCdatExpr0 <- distinct(RPKM_mean[,c(1,3:44)])
PFCdatExpr0 <- PFCdatExpr0[complete.cases(PFCdatExpr0),]

#duplicated row name: PFCdatExpr0$ensembl_gene_id[(duplicated(PFCdatExpr0$ensembl_gene_id))]
#"ENSG00000230417"
#remove it

PFCdatExpr0 <- PFCdatExpr0[ !grepl("ENSG00000230417",PFCdatExpr0$ensembl_gene_id ) , ]
PFCdatExpr0 <- as.data.frame(PFCdatExpr0)
rownames(PFCdatExpr0) <- PFCdatExpr0$ensembl_gene_id

PFCdatExpr0 <- PFCdatExpr0 [,c(2:ncol(PFCdatExpr0))]

#make gene matrix with genes as columns and rows as samples

PFCdatExpr0T <- t(PFCdatExpr0)
  
#detect genes with missing values or samples with missing values
gsg = goodSamplesGenes(PFCdatExpr0T, verbose = 3);
gsg$allOK # FALSE: Excluding 7311 genes from the calculation due to too many missing samples or zero variance

#remove the genes with too many missing values
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(PFCdatExpr0T)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(PFCdatExpr0T)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  PFCdatExpr0T2 = PFCdatExpr0T[gsg$goodSamples, gsg$goodGenes]
}


#gene with RPKM value of 0.1 or higher in at least one sample:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4851109/
#max(col)>=.1
#https://stackoverflow.com/questions/24212739/how-to-find-the-highest-value-of-a-column-in-a-data-frame-in-r/24212879
colMax <- function(data) sapply(data, max, na.rm = TRUE)

PFCdatExpr0T2 <- as.data.frame(PFCdatExpr0T2)

PFCdatExpr0.1 <- PFCdatExpr0T2[, colMax(PFCdatExpr0T2) >= .001] #leaves 13,157 genes

#log2 transform:
PFCdatExpr <- log2(PFCdatExpr0.1+1)

# median absolute deviation 
#remove if = 0 https://support.bioconductor.org/p/65124/
colMad <- function(data) sapply(data, mad, na.rm = TRUE)
mad.PFCdatExpr = PFCdatExpr[, colMad(PFCdatExpr) == 0]

PFCdatExpr = PFCdatExpr[, colMad(PFCdatExpr) != 0] #leaves 10,895 genes

save(PFCdatExpr,file="Zhong2018_WGCNA_log2RPKMInputdata.Rdata")
#=====================================================================================
#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(PFCdatExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#=====================================================================================
#phenotype data:
traits <- as.data.frame(rownames(PFCdatExpr))
colnames(traits) <- c("Sample")

traits$GW <- substr(traits$Sample,3,4)
traits$GW <- as.numeric(traits$GW)

traits$Sample <- as.character(traits$Sample )

traits$cell_type <- substr(traits$Sample, 6,nchar(traits$Sample))
levels = c("Stem cells","OPC","Neurons","GABAergic neurons","Astrocytes","Microglia")
traits$cell_type  <- factor(traits$cell_type , levels=levels)
traits$cell_type2 <- as.numeric(traits$cell_type)
  
save(PFCdatExpr,traits,file="Zhong2018_WGCNA_log2RPKMInputdata.Rdata")
load("Zhong2018_WGCNA_log2RPKMInputdata.Rdata")
#=====================================================================================


sampleTree2 = hclust(dist(PFCdatExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traits[,c(2,4)], signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traits[,c(2,4)]), 
                    main = "Sample dendrogram and trait heatmap")

ggsave(filename="Sample_clustering_Zhong2018.pdf",width = 8, height = 5, dpi = 300)

#=====================================================================================
#choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency 
# Choose a set of soft-thresholding powers
powers = c(c(1:40), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(PFCdatExpr, powerVector = powers, networkType = "signed",corFnc="bicor" ,verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="black")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="blue")

write.csv(sft, "SoftThresholding_powercalc.csv")
#lowest power for which the scale-free topology fit index of > 0.8

#=====================================================================================
#Run this section on an external server, save the network data and load into R studio
#=====================================================================================
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

#=====================================================================================
#back in R studio: load network data
#=====================================================================================

#bwnet$colors contains the module assignment, and bwnet$MEs contains the module eigengenes of the modules.
load("NetworkConstruction_Zhong2018.RData")


table(bwnet$colors)


# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwnet$colors)


# Plot the dendrogram and the module colors for all blocks
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    main = "Gene dendrogram and module colors in all blocks",
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree, 
#      file = "DLFPC-networkConstruction-blockwise.RData")
#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(PFCdatExpr, colors = moduleColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#=====================================================================================
#merge down to fewer modules based on plot height
MEDissThres = 0.8
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(PFCdatExpr, moduleColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mColors = merge$colors;

#fix colors:
mColors <- gsub("grey60","purple",mColors)
mColors <- gsub("cyan","green",mColors)
mColors <- gsub("darkturquoise","red",mColors)
mColors <- gsub("lightsteelblue1","blue",mColors)

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

colnames(mergedMEs) <- gsub("grey60","purple",colnames(mergedMEs))
colnames(mergedMEs) <- gsub("cyan","green",colnames(mergedMEs))
colnames(mergedMEs) <- gsub("darkturquoise","red",colnames(mergedMEs))
colnames(mergedMEs) <- gsub("lightsteelblue1","blue",colnames(mergedMEs))

#replot
sizeGrWindow(12, 9)
pdf(file = "Zhong_Dendrogram_postMerged.pdf", wi = 9, he = 6)

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(mColors[bwnet$blockGenes[[1]]]),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
# Rename to moduleColors
moduleColors = mColors

modulesize <- table(moduleColors)
library(openxlsx)
write.xlsx(modulesize,"modulesize.xlsx")

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, traits,PFCdatExpr,file = "Zhong-networkConstruction_postmerge.RData")
load("Zhong-networkConstruction_postmerge.RData")
#=====================================================================================
#repeat clustering of new MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
pdf(file = "PFCgeneDendroMerged_postmerge.pdf", wi = 5, he = 5)
plot(METree, main = "Clustering of module eigengenes after merging modules",
     xlab = "", sub = "")

dev.off()

#=====================================================================================
#on server
#=====================================================================================
#Co-expression similarity and adjacency
#signed, biweight midcorrelation between all genes 

library("WGCNA")
allowWGCNAThreads(nThreads = 24)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#=====================================================================================

load("Zhong-networkConstruction_postmerge.RData")

# #loads as TOM
load("Zhong_TOM-block.1.RData")

#=====================================================================================

# 
TOM <- as.matrix(TOM)
colnames(TOM) = colnames(PFCdatExpr)
# 
dissTOM=1-TOM

#MDS plot
cmd1=cmdscale(as.dist(dissTOM),2)
sizeGrWindow(7, 6)

pdf(file = "MDSplot_dissTOM.pdf", wi = 16, he = 11)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(moduleColors), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
dev.off()

save(moduleColors, dissTOM,cmd1,PFCdatExpr,file = "Zhong2018_Tom_postmerge.RData")
#=====================================================================================
#Trait Correlation
#=====================================================================================
#correlate eigengenes with external traits and look for the most signicant associations
traits2 <- traits[,c(2,4)] #numeric
rownames(traits2) <- traits$Sample

# Define numbers of genes and samples
nGenes = ncol(PFCdatExpr);
nSamples = nrow(PFCdatExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(PFCdatExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#FDR correct the pvalues
#sapply(pval,p.adjust,method="fdr") #per column

#FDR correction on entire matrix
FDR <- matrix(p.adjust(as.vector(moduleTraitPvalue), method='fdr'),ncol=ncol(moduleTraitPvalue))
colnames(FDR) <- paste(colnames(moduleTraitPvalue),"FDR pvalue", sep =" ")
rownames(FDR) <- rownames(moduleTraitPvalue)


corout <- cbind(moduleTraitCor,FDR)
write.csv(corout, "PearsonsCorrelations_mod-trait.csv",row.names = T)

#We color code each association by the correlation value:

sizeGrWindow(10,6)

pdf(file = "Module-traitrelationshipsFDRcorrected.pdf", wi = 8.5, he = 11)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\ (", #pearson correlation coefficient, space, (FDR corrected pvalue)
                    signif(FDR, 4), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

#=====================================================================================

#quantify associations of individual genes with our trait of interest (years) by dening Gene Signicance GS as
#(the absolute value of) the correlation between the gene and the trait

#For each module, we also dene a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression prole. This
#allows us to quantify the similarity of all genes on the array to every module.


#remove extra subject column from MEs
MEs <- MEs[,-7]
# names (colors) of the modules
modNames = substring(names(MEs), 3)

##module membership MM 
#correlation between expression of each gene and MEs (uncorrected p values)
geneModuleMembership = as.data.frame(cor(PFCdatExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#or
# calculate the module membership values (aka. module eigengene based
# connectivity kME):
#datKME = signedKME(PFCdatExpr, MEs)

# Define variable age containing the years column of datTrait
GW = as.data.frame(datTraits$GW);
names(GW) = "GW"

#correlate each gene expression with variable of interest
geneTraitSignificanceGW = as.data.frame(cor(PFCdatExpr, GW, use = "p"));
GSPvalueGW = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceGW), nSamples));

names(geneTraitSignificanceGW) = paste("GS.", names(GW), sep="");
names(GSPvalueGW) = paste("p.GS.", names(GW), sep="")

# Define variable age containing the years column of datTrait
cell_type = as.data.frame(datTraits$cell_type2);
names(cell_type) = "cell_type"

#correlate each gene expression with variable of interest
geneTraitSignificancecell_type = as.data.frame(cor(PFCdatExpr, cell_type, use = "p"));
GSPvaluecell_type = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificancecell_type), nSamples));

names(geneTraitSignificancecell_type) = paste("GS.", names(cell_type), sep="");
names(GSPvaluecell_type) = paste("p.GS.", names(cell_type), sep="")

#=====================================================================================
#Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high signicance for weight as well as high module
#membership in interesting modules.


#As an example, we look at the purple module that has the highest association with GW.

########## puprle module and cell type ##########
module = "purple"
column = match(module, modNames);
moduleGenes = moduleColors==module;


pdf(file = "PurpleModule_MMvsGS_celltype.pdf", wi = 5, he = 5,useDingbats = F)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificancecell_type[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Gestational Week",
                   main = paste("Module membership vs. gene significance cell type\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#GS and MM are correlated, illustrating that genes highly significantly associated with a trait are 
#often also the most important (central) elements of modules associated with the trait.
dev.off()

########## red module and cell type ##########
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;


pdf(file = "RedModule_MMvsGS_celltype.pdf", wi = 5, he = 5,useDingbats = F)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificancecell_type[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Cell Type",
                   main = paste("Module membership vs. gene significance cell type\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#GS and MM are correlated, illustrating that genes highly significantly associated with a trait are 
#often also the most important (central) elements of modules associated with the trait.
dev.off()

########## skyblue module and GW ##########
module = "skyblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;


pdf(file = "SkyBlueModule_MMvsGS_GW.pdf", wi = 5, he = 5,useDingbats = F)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificanceGW[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Gestational Week",
                   main = paste("Module membership vs. gene significance GW\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#GS and MM are correlated, illustrating that genes highly significantly associated with a trait are 
#often also the most important (central) elements of modules associated with the trait.
dev.off()

#=====================================================================================

#modules with high association with our trait of interest, and have identified their central players by the Module Membership measure.
#merge this statistical information with gene annotation and write out a file that summarizes the most important results and can be inspected in standard spreadsheet software 
# Create the starting data frame
Ensemble = colnames(PFCdatExpr)

#get matching data for each ensemble id from the row data from Brainspain input
ids <- RPKM_mean[which(RPKM_mean$ensembl_gene_id %in% Ensemble),]

geneInfo0 = data.frame(Ensemble = colnames(PFCdatExpr),
                       geneSymbol = ids$gene_ids,
                       moduleColor = moduleColors,
                       geneTraitSignificanceGW = geneTraitSignificanceGW,
                       GSPvalueGW = GSPvalueGW,
                       geneTraitSignificancecell_type = geneTraitSignificancecell_type,
                       GSPvaluecell_type = GSPvaluecell_type,
                       t(PFCdatExpr))


#write to file:
write.xlsx(geneInfo0, file = "Zhong2018_WGCNA_geneInfo.xlsx")
write.csv(geneInfo0, file = "Zhong2018_WGCNA_geneInfo.csv")

traits3 <-cbind(traits,MEs)
write.xlsx(traits3,file = "Zhong2018_WGCNA_MEs.xlsx")
#geneInfo0 <- read.csv("Zhong2018_WGCNA_geneInfo.csv")

#=====================================================================================
#=====================================================================================
#Extract modules
module_colors = setdiff(unique(moduleColors), "grey")


#Look at expression patterns of these genes, as they are clustered
#heatmap colors:
#install.packages("RColorBrewer")
library("RColorBrewer")
library(gplots)

#gene-module info
GM <- geneInfo0 %>% dplyr::select(Ensemble,moduleColor)
GM <- as.data.frame(GM)
colnames(GM) <- c("EnsembleID","module")

PFCdatExpr2 <- t(PFCdatExpr)

PFCdatExpr2 <- as.data.frame(PFCdatExpr2)

PFCdatExpr2$EnsembleID <- rownames(PFCdatExpr2)

merPFC <- merge(PFCdatExpr2,GM, by = "EnsembleID")
merPFC$module <- as.factor(merPFC$module )
write.csv(merPFC,"Log2RPKM+1_PFC_byModule.csv")

merPFC2 <- merPFC %>% filter(module != "grey") %>% arrange(module)
merPFC2$module <- factor(merPFC2$module )

data <- as.matrix(merPFC2[,2:43])
rownames(data) <- merPFC2$EnsembleID

#reorder by development and then
order <- as.data.frame(colnames(data))
colnames(order) <- c("name")
order$time <- substr(order$name,1,4)
order$time <- as.factor(order$time)

order$name <- as.character(order$name)

order$cell <- substr(order$name, 6,nchar(order$name))
levels = c("Stem cells","OPC","Neurons","GABAergic neurons","Astrocytes","Microglia")
order$cell <- factor(order$cell, levels=levels)

order <- order %>% arrange(cell,time)

#reorder matrix:
data2 <- data[,order$name]

save(data2,order,merPFC2,file="HeatmapMatrix.Rdata")
load("HeatmapMatrix.Rdata")
#=====================================================================================
#Heatmap
#=====================================================================================

#heatmap:
library(pheatmap)
#my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

annotation_col <- data.frame(
  #Sample = factor(d$donor_name),
  GW = order$time,
  Cell_Type = order$cell)

rownames(annotation_col) = colnames(data2)

head(annotation_col)

annotation_row <- data.frame(merPFC2$module)

rownames(annotation_row) = rownames(data2)

colnames(annotation_row) <- c("module")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
colfunc_blue <- colorRampPalette(c("lightblue", "darkblue"))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
#cell type
Var1        <- c("darkorchid3","turquoise2","darkgoldenrod2", "deeppink2","steelblue4","salmon")

names(Var1) <- unique(annotation_col$Cell_Type)

#GW
Var2        <- colfunc_blue(length(unique(annotation_col$GW)))

names(Var2) <- unique(annotation_col$GW)


Var3 <- as.character(unique(merPFC2$module))
names(Var3) <-  unique(merPFC2$module)

anno_colors <- list(Cell_Type = Var1, GW = Var2, Module = Var3)


pdf(file = "Zhongmodule-Heatmap_log2RPKM+1_zscore.pdf", wi = 8, he = 6)
pheatmap::pheatmap(data2, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = F,
                   annotation_names_row=T,
                  color = my_palette, 
                   fontsize = 10,
                   fontsize_row=6, 
                   fontsize_col = 10,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Single Cell Brain Gene Expression Modules")

dev.off()

#=====================================================================================
#=====================================================================================
#line graph of average expression over time
#=====================================================================================
#scale data same as pheatmap: subtract mean and divide by std. dev
#https://www.biostars.org/p/223532/
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}


#scale by row:
merPFCmatrix <- scale_mat(data2, scale = "row")
merPFCmatrix <- as.data.frame(merPFCmatrix)

merPFCmatrix$module <- merPFC2$module

merPFCDF <- merPFCmatrix %>% gather(Sample, scaled_expression, 1:(ncol(merPFCmatrix)-1)) 
# 
merPFCDF <- as.data.frame(merPFCDF)

#add in GW and cell type factors
merPFCDF2 <- merge(merPFCDF,datTraits, by="Sample")

merPFCDF2 <- merPFCDF2 %>% group_by(module,GW,cell_type) %>%
  summarize(meanlog2RPKM = mean(scaled_expression))
# ordercol <- c( "FirstTrimester" ,  "Second_Trimester" ,"Third_Trimester" ,"First_Year" ,"Year_2_4", "Year_8_15", "Year_18_23","Year_30plus")
# 
# merPFCDF$age <- as.factor(merPFCDF$age)
# levels(merPFCDF$age) <- ordercol
# 
# tan <- merPFCDF %>% filter(module == "tan"|module == "blue")

#set colors:
my_palette2 <- unique(merPFCDF2$module)

my_palette2 <- col2hex(my_palette2)

pdf(file = "Zhongmodule_lineplot_log2RPKM+1_scaled.pdf", wi = 11, he = 8.5)

p <- ggplot(merPFCDF2, aes(x=GW, y=meanlog2RPKM, group=module, color=module)) +
  facet_wrap(~cell_type, scales="free")+
  geom_line(size=2) +
  scale_color_manual(values = my_palette2)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))+
  xlab(label = c("Gestational Week")) +
  ylab(label = c("Scaled log2(RPKM+1)"))

p

dev.off()

pdf(file = "Zhongmodule_lineplot_log2RPKM+1_scaledv2.pdf", wi = 11, he = 8.5)

p <- ggplot(merPFCDF2, aes(x=GW, y=meanlog2RPKM, group=cell_type, color=cell_type)) +
  facet_wrap(~module, scales="free")+
  geom_line(size=2) +
  scale_color_manual(values = Var1)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))+
  xlab(label = c("Gestational Week")) +
  ylab(label = c("Scaled log2(RPKM+1)"))

p

dev.off()
#=====================================================================================
#List Enrichment:
#=====================================================================================
###################################################################
#make internal collection for gene lists of interest for ASD and microglia
###################################################################

#install.packages("/Users/annieciernia/Downloads/anRichmentMethods_0.87-1.tar.gz", repos = NULL, type = "source")
#install.packages("/Users/annieciernia/Downloads/anRichment_0.97-1.tar.gz", repos = NULL, type = "source")


library(anRichment)
options(stringsAsFactors = FALSE)

# #need three files: geneSetInfo:
# #geneSetContent
# #groupInfo

###################################################################

#gene lists






###################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)

#get human genes
human = useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

#get all ids
#create vector of chromosomes
my_chr <- c(1:22,'X','Y')

ids <- getBM(
  filters = 'chromosome_name',
  attributes= c("ensembl_gene_id", "hgnc_symbol","entrezgene", "description"),
  values = my_chr,
  mart= human)

write.csv(ids,"HumanHg38_allids.csv")

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/old/8_2018_4CpGDMRs\ copy/CustomCollections/CustomGeneLists")
ids <- read.csv("HumanHg38_allids.csv", header=T)

#function to convert ensemble to entrezgene id
convertEnsbltoEntrez <- function(x){
  
  geneids <- ids[which(ids$ensembl_gene_id %in% x),]
  geneids <- unique(geneids$entrezgene)
  geneids <- geneids[complete.cases(geneids)]
  
  return(geneids)
}

#test2 <- convertEnsbltoEntrez(test)
###################################################################

#load in gene lists
#change 
setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/gene_assignment")
load("Hg38_Conversioninfo.RData")
#ASD DMR genes 
ASD_genes <- read.table("ASD_DMR_EnsemblHG38.txt",header=F)
ASD_genes <- as.character(ASD_genes$V1)
Dup15q_genes <- read.table("Dup15q_DMR_EnsemblHG38.txt",header=F)
Dup15q_genes <- as.character(Dup15q_genes$V1)
RTT_genes <- read.table("Rett_DMR_EnsemblHG38.txt",header=F)
RTT_genes <- as.character(RTT_genes$V1)

ASD_entrez <- convertEnsbltoEntrez(ASD_genes)
Dup15q_entrez <- convertEnsbltoEntrez(Dup15q_genes)
RTT_entrez <- convertEnsbltoEntrez(RTT_genes)

master <- list(ASD_entrez,Dup15q_entrez,RTT_entrez)
names(master) <- c("ASD DMR genes","Dup15q DMR genes","RTT DMR genes")


bbGeneSetASD = newGeneSet(
  geneEntrez = ASD_entrez,
  geneEvidence = "IEP",
  geneSource = paste0("ASD DMR genes"),
  ID = "ASD DMR genes",
  name = "ASD DMR genes",
  description = "ASD DMR genes",
  source = paste0("ASD DMR genes"),
  organism = "human",
  internalClassification = c("PL", "DMR"),
  groups = "PL",
  lastModified = "2018-8-12");


PLgroup = newGroup(name = "PL", description = "PL's experimental group of gene sets",
                   source = "Personal imagination");


bbGeneSetDup15q = newGeneSet(
  geneEntrez = Dup15q_entrez,
  geneEvidence = "IEP",
  geneSource = paste0("Dup15q DMR genes"),
  ID = "Dup15q DMR genes",
  name = "Dup15q DMR genes",
  description = "Dup15q DMR genes",
  source = paste0("Dup15q DMR genes"),
  organism = "human",
  internalClassification = c("PL", "DMR"),
  groups = "PL",
  lastModified = "2018-8-12");

bbGeneSetRTT = newGeneSet(
  geneEntrez = RTT_entrez,
  geneEvidence = "IEP",
  geneSource = paste0("RTT DMR genes"),
  ID = "RTT DMR genes",
  name = "RTT DMR genes",
  description = "RTT DMR genes",
  source = paste0("RTT DMR genes"),
  organism = "human",
  internalClassification = c("PL", "DMR"),
  groups = "PL",
  lastModified = "2018-8-12");

PLcollection = newCollection(dataSets = list(bbGeneSet,bbGeneSetDup15q,bbGeneSetRTT), groups = list(PLgroup));


setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/ZhongSC_WGCNA")
geneInfo0 <- read.csv("Zhong2018_WGCNA_geneInfo.csv")

entrez <- ids[which(ids$ensembl_gene_id %in% geneInfo0$Ensemble),]
merge_entrez <- merge(geneInfo0,entrez, by.x= c("Ensemble"), by.y = c("ensembl_gene_id"))

merge_entrez <- merge_entrez %>% dplyr::select(moduleColor,entrezgene) %>%
  distinct()

Combinedenrichment = enrichmentAnalysis(
  classLabels = merge_entrez$moduleColor, 
  identifiers = merge_entrez$entrezgene,
  refCollection = PLcollection,
  useBackground = "given", 
  threshold = 1.0,
  thresholdType = "FDR",
  getOverlapEntrez = FALSE,
  getOverlapSymbols = TRUE,
  entrySeparator = ",",
  maxReportedOverlapGenes = 1000,
  ignoreLabels = "grey") #ignore genes not assigned to modules (grey)

table <- Combinedenrichment$enrichmentTable

write.xlsx(table,"Zhong_DMR_enrichments.xlsx")

save(Combinedenrichment,table,merge_entrez,PLcollection,file="ASD_celltypeEnrichments,RData")


############################################################################################################
#GO terms:
###########################################################################################################

#calculate the enrichment using the modified:
GOenrichment = enrichmentAnalysis(
  classLabels = geneInfo2$moduleColor, 
  identifiers = geneInfo2$entrezgene,
  refCollection = GOcollection,
  useBackground = "intersection", #intersection of genes in identifiers and the organism database
  threshold = 0.05,
  thresholdType = "FDR",
  getOverlapEntrez = FALSE,
  getOverlapSymbols = TRUE,
  entrySeparator = ",",
  maxReportedOverlapGenes = 500,
  ignoreLabels = "grey")



#Creating enrichment labels for classes
# The function takes as input the enrichment table returned by enrichmentAnalysis in the enrichmentTable component, and lets the user specify the various columns needed to put together the enrichment label

#The function returns a data frame with information about highest enriched terms in each of the groups specified
#in focusOnGroups (“all” is a special keyword meaning all terms irrespective of what group they belong to).

eLabels = enrichmentLabels(
  GOenrichment$enrichmentTable,
  focusOnGroups = c("all", "GO"),
  groupShortNames = c("all", "GO", "CT"),
  minSize = 0.05,
  pValueThreshold = 0.05,
  numericClassLabels = FALSE)

#enrichment results are summarized in the component enrichmentTable:

sigenrichGO <- GOenrichment$enrichmentTable 

sigenrichGO <- sigenrichGO %>% filter(FDR < 0.05)

write.xlsx(sigenrichGO, "SignificantGOtermEnrichment-module.xlsx")

GO <- eLabels %>% dplyr::select(class,highestEnrichedSet.GO)

GObest <- sigenrichGO[which(sigenrichGO$dataSetName %in% GO$highestEnrichedSet.GO ),]

GObestDF <- GObest %>%
  group_by(class) %>%
  top_n(n = -1, wt = rank)  %>% #select lowest rank for each class
  dplyr::select(class,inGroups,dataSetName,FDR,enrichmentRatio)

GOmatrix <- as.matrix(GObestDF$enrichmentRatio)
rownames(GOmatrix) <- GObestDF$class

textMatrix =  paste(GObestDF$dataSetName,"\n", signif(GObestDF$enrichmentRatio, 2), " (", #enrichment value, space, (FDR corrected pvalue)
                    signif(GObestDF$FDR, 2), ")", sep = "")

dim(textMatrix) = dim(GOmatrix)

textMatrix[textMatrix == "NA\n(NA)"] <- c("")


pdf(file = "PFCmodule-GOterm_relationshipsFDRcorrected.pdf", wi = 8.5, he = 11)

# Display the values within a heatmap plot
#values are enrichment ratio (FDR pvalue) and scale is log10(enrichmentratio)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = log10(GOmatrix),
               xLabels = colnames(GOmatrix),
               yLabels = row.names(GOmatrix),
               ySymbols = row.names(GOmatrix),
               colorLabels = FALSE,
               yColorLabels =TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               main = paste("Module-Genelist Enrichments: Log10 Enrichment Ratio")
)

dev.off()


