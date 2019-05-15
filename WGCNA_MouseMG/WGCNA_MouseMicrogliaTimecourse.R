#####################################################################################
# Gene expression from human brain across development: brainspan
####################################################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)

#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gplots)


#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

######################################################################
# Basic function to convert mouse to human gene names
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

convertMouseGeneList <- function(x){
  #for hg38:
  genesV2 = getLDS(attributes = c("mgi_symbol","ensembl_gene_id"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol","ensembl_gene_id","entrezgene"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

######################################################################
#Change to your working directory
setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/MandF2212/githup_code/WGCNA/Mouse_Microglia")
######################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
#Initial data processing, sample cleaning and differential expression analysis 
############################################################################################################################
############################################################################################################################
############################################################################################################################

##############################################################################################
#DE genes, remove sample tube67 for <2,000 reads and tube 71 for <200,000 reads
##############################################################################################
#raw count data per sample from GEO:
bilbo <- read.delim("GSE99622_TranscriptRawReadCounts.txt", header = T,sep="\t", stringsAsFactors=FALSE)
genelength <- read.delim("genelength.txt", header = T,sep="\t", stringsAsFactors=FALSE)
genelength$length <- abs(genelength$mm10.knownGene.txStart - genelength$mm10.knownGene.txEnd)
colnames(genelength)[6] <- c("Transcript")
bilbo3 <- merge(bilbo,genelength, by = "Transcript",all.x=T )

bilbo <- bilbo3 %>% filter(Best == 1)


#remove sample tube67 for <2,000 reads and tube 71 for <200,000 reads
bilbo$X..C9.._tube67_P60_LPS_4.10.14_F.4 <- NULL
bilbo$X..G9.._tube71_P60_LPS_4.10.14_M.3 <- NULL

samples <- colnames(bilbo)[4:63]
samples <- as.data.frame(samples)
library(stringr)
samples2 <- as.data.frame(str_split_fixed(samples$samples, "_", 4))
samples$tube <- samples2$V2
samples$timepoint <- samples2$V3
samples$treatment <- str_split_fixed(samples2$V4, "_", 2)


samples2$sub <- str_sub(samples2$V4, start= -3)
samples2$sub <- sub("4_M", "M.4", samples2$sub) #replace 4_M with 4.M
samples2$sub <- sub("4_F", "F.4", samples2$sub) #replace 4_F with 4.F

samples3 <- substr(samples2$sub, start= 1, stop =1)

samples$sex <- samples3


samples5 <- sub("[0-9]\\.([0-9]|[0-9][0-9])\\.[0-9][0-9]", "none", samples$treatment) #replace date with none
samples5 <- as.data.frame(samples5)
samples$treatment <- samples5$V1
samples$group <- paste(samples$timepoint,samples$sex,samples$treatment, sep="_")

#sub NA for 0
matrix <- as.matrix(bilbo[,4:63])
rownames(matrix) <- bilbo$Gene

matrix2 <- replace(matrix,which(is.na(matrix)),0)

matrix2 <- matrix2[!duplicated(rownames(matrix2)), ]

genes <- rownames(matrix2)
samplenames <- colnames(matrix2)

group <- samples$group

# Make DEGList
d <- DGEList(counts = matrix2, group=group, genes = genes, samples = samplenames)

dim(d)

RGcpm <- d

#reset library sizes
RGcpm$samples$lib.size <- colSums(RGcpm$counts)


#Normalize library
RGnorm <- calcNormFactors(RGcpm)
RGnorm$samples$group

#plot MDS
plotMDS(RGnorm, col = rainbow(length(levels(factor(group))))[factor(group)],cex=1, main="MDS")

#build factors
Timepoint <- as.factor(samples$timepoint)
sex <- as.factor(samples$sex)
treatment <- as.factor(samples$treatment)
treatment <- relevel(treatment, ref="none")


design <- model.matrix(~0+ group,data = samples)

#estimate dispersion
RGnorm <- estimateGLMCommonDisp(RGnorm,design,verbose=TRUE)##RG change
#Disp = 0.0342 , BCV = 0.1849 

RGnorm <- estimateGLMTrendedDisp(RGnorm,design)##RG change
RGnorm <- estimateGLMTagwiseDisp(RGnorm,design)##RG change
plotBCV(RGnorm)##RG change

#png("plotBCV_RGnorm.png")

#DE Expression
fit <- glmFit(RGnorm,design)


#set constrasts
mycontrasts <- makeContrasts(
  E18FvsM = "groupE18_F_none-groupE18_M_none",
  P4FvsM = "groupP4_F_none-groupP4_M_none",
  P14FvsM = "groupP14_F_none-groupP14_M_none",
  P60FSAL_P60MSAL = "groupP60_F_SAL-groupP60_M_SAL",
  P60FLPS_P60MLPS = "groupP60_F_LPS-groupP60_M_LPS",
  P60FSAL_P60FLPS = "groupP60_F_SAL-groupP60_F_LPS",
  P60MSAL_P60MLPS = "groupP60_M_SAL-groupP60_M_LPS",
  F_E18vsP4 = "groupE18_F_none-groupP4_F_none",
  F_E18vsP14 = "groupE18_F_none-groupP14_F_none",
  F_E18vsP60 = "groupE18_F_none-groupP60_F_SAL",
  F_P4vsP14 = "groupP4_F_none-groupP14_F_none",
  F_P4vsP60 = "groupP4_F_none-groupP60_F_SAL",
  F_P14vsP60 = "groupP14_F_none-groupP60_F_SAL",
  M_E18vsP4 = "groupE18_M_none-groupP4_M_none",
  M_E18vsP14 = "groupE18_M_none-groupP14_M_none",
  M_E18vsP60 = "groupE18_M_none-groupP60_M_SAL",
  M_P4vsP14 = "groupP4_M_none-groupP14_M_none",
  M_P4vsP60 = "groupP4_M_none-groupP60_M_SAL",
  M_P14vsP60 = "groupP14_M_none-groupP60_M_SAL",
  levels=design)

##loop for comparisons
comparisons=colnames(mycontrasts)


#run each comparison from contrast through this loop
DEcount2 <-NULL
diff.out2 <- as.data.frame(bilbo$Gene)
colnames(diff.out2)[1] <- c("genes")
diff.out3 <- NULL

for(i in 1:length(comparisons)){
  #comparison name
  print(i)
  comp = comparisons[i]
  #make comparisons 
  diff <- glmLRT(fit,contrast=mycontrasts[,comp])
  o <- order(diff$table$PValue)
  #summary fo diffexp genes
  diffsum <- summary(de <-decideTestsDGE(diff))
  #CPM for diff exp genes per animal
  numberDEgenes <- diffsum[1]+diffsum[3]
  DEUP <- diffsum[1]
  DEDOWN <-diffsum[3]
  DEcount <- cbind(comp,numberDEgenes,DEUP,DEDOWN)
  DEcount2 <- rbind(DEcount2,DEcount)
  
  diff.byanimal <- cpm(RGnorm) [o[1:numberDEgenes],]
  diff.byanimal <- as.data.frame(diff.byanimal)
  diff.byanimal$Geneid <- rownames(diff.byanimal)
  #List of diff exp genes
  #diff.toptags <- topTags(diff)
  #returns all genes - filter by p value
  diff.toptags <- topTags(diff, n=100000, adjust.method="BH")
  diff.toptags <- data.frame(diff.toptags)
  #add Geneid as a column for use in merge step below
  diff.toptags <- cbind(Geneid = rownames(diff.toptags), diff.toptags)
  FDRvalues <- as.data.frame(cbind(diff.toptags$FDR,diff.toptags$genes))
  colnames(FDRvalues)[1] <- paste("FDR")
  colnames(FDRvalues)[2] <- paste("genes")
  FDRvalues$comparison <- paste(comp)
  FDRvalues <- subset(FDRvalues, !duplicated(FDRvalues[,2]))
  
  #diff.out <- left_join(diff.out2, FDRvalues, by = "genes")
  diff.out3 <- rbind(diff.out3,FDRvalues)
  #add in ave CPM per group for the genes in diff.toptags
  #diff.toptags <- merge(x = diff.toptags, y = AveLog2CPM, by = "Geneid", all.x = FALSE)
  
  #Toptags
  write.table(diff.toptags, file = paste(comp, ".toptags.txt", sep="."), sep="\t",row.names = F)
  #Idividual CPM values
  write.table(diff.byanimal, file = paste(comp, ".CPMs.txt", sep="."), sep="\t",row.names = F)
  #Summary of genes up and down
  write.table(diffsum, file = paste(comp, ".diffsum.txt", sep="."), sep="\t")
  
}

#counts
write.table(DEcount2, file = "DEgeneCounts.txt", sep="\t",row.names = F)
############################################################################################################################
#write out FPKM table:
############################################################################################################################

#get normalized count matrix:
nc <- cpm(RGnorm, normalized.lib.sizes=T,log=F, prior.count=1)
nc <-as.data.frame(nc)

#manually change column names
oldnames <- as.data.frame(colnames(nc))
colnames(oldnames)[1] <- c("samples")

samples$name <- paste(samples$group,samples$tube, sep="_")
mergenames <- merge(oldnames,samples, by = "samples")

k <- mergenames[ order(match(mergenames$samples, samples$samples)), ] #reorder so matches dataframe of log10values

#rename columnes for animal and condition
colnames(nc) <- k$name

nc$geneid <- rownames(nc)

write.table(nc,file="Bilbo_normalizedcountmatrix.txt",sep="\t", quote=FALSE,col.names = T,row.names = F)
write.table(samples,file="Bilbo_sampleinfo.txt",sep="\t", quote=FALSE,col.names = T,row.names = F)

############################################################################################################################
############################################################################################################################
############################################################################################################################
#WGCNA
############################################################################################################################
############################################################################################################################
############################################################################################################################
#Biblo Timecourse:no filter
Bilbotimecourse <- read.delim("Bilbo_normalizedcountmatrix.txt", header = T,sep="\t", stringsAsFactors=FALSE)
bilbosamples <- read.delim("Bilbo_sampleinfo.txt", header = T,sep="\t", stringsAsFactors=FALSE)
colnames(Bilbotimecourse)[colnames(Bilbotimecourse)=="geneid"] <- "mgi_symbol"

#convert mgi symbols to human ensemble
geneids <- convertMouseGeneList(Bilbotimecourse$mgi_symbol)

merge <- merge(Bilbotimecourse, geneids, by.x = "mgi_symbol", by.y = "MGI.symbol")

DF <- merge %>% gather(samples,FPKM, 2:61) 

bilbosamples$samples <- paste(bilbosamples$group,bilbosamples$tube,sep ="_")

DF2 <- merge(DF, bilbosamples,by = "samples",all.x=T)

write.csv(DF2,"TransformedDF_BilboMGtimecourse_20001_hg38.csv")

#DF2 <- read.csv("TransformedDF_BilboMGtimecourse_20001_hg38.csv",stringsAsFactors=FALSE)

#make WGCNA table

DFTable <- DF2 %>% dplyr::select(Gene.stable.ID.1,FPKM,samples) %>%
  distinct() 

# #nonduplicates
mydf <- DFTable %>%
   group_by(Gene.stable.ID.1,samples) %>%
   mutate(n = n()) %>%
   filter(n == 1)

#duplicates
mydf2 <- DFTable %>%
  group_by(Gene.stable.ID.1,samples) %>%
  mutate(n = n()) %>%
  filter(n != 1) %>%
  mutate(ident_rnk = min_rank(FPKM)) %>%
  top_n(n = 1,wt = ident_rnk) %>%
  dplyr::select(-ident_rnk)

mydf3 <- rbind(mydf,mydf2)

#make gene matrix with genes as columns and rows as samples
mydf4 <- mydf3 %>% dplyr::select(-n) %>%
  spread(Gene.stable.ID.1, FPKM)

write.csv(mydf4,"Cleaned_DF_BilboMGtimecourse_hg38.csv")
#mydf4 <- read.csv("Cleaned_DF_BilboMGtimecourse_hg38.csv")

#remove LPS samples
mydf4 <- as.data.frame(mydf4)
rownames(mydf4) <- mydf4$samples

DFnoLPS <- mydf4[!grepl("LPS*", rownames(mydf4)), ]

DFnoLPS <- DFnoLPS[,2:ncol(DFnoLPS)]

#remove extra columns:
DFnoLPS <- DFnoLPS %>% dplyr::select(-samples)
DFnoLPS <- as.data.frame(DFnoLPS)
################## process data ########################
#detect genes with missing values or samples with missing values
gsg = goodSamplesGenes(DFnoLPS, verbose = 3);
gsg$allOK # FALSE: Excluding 7311 genes from the calculation due to too many missing samples or zero variance

#remove the genes with too many missing values
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(DFnoLPS)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(DFnoLPS)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  DFnoLPS = DFnoLPS[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(DFnoLPS, verbose = 3);
gsg$allOK # FALSE: Excluding 7311 genes from the calculation due to too many missing samples or zero variance


#gene with RPKM value of 0.25 or higher in at least one sample.
#max(col)>=.25
#https://stackoverflow.com/questions/24212739/how-to-find-the-highest-value-of-a-column-in-a-data-frame-in-r/24212879
colMax <- function(data) sapply(data, max, na.rm = TRUE)

DFnoLPS.1 <- DFnoLPS[, colMax(DFnoLPS) >= .25]

#log2 transform:
Log2DFnoLPS <- log2(DFnoLPS.1+1)

# median absolute deviation 
#remove if = 0 https://support.bioconductor.org/p/65124/
colMad <- function(data) sapply(data, mad, na.rm = TRUE)

Log2DFnoLPS = Log2DFnoLPS[, colMad(Log2DFnoLPS) != 0]

#=====================================================================================
#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(Log2DFnoLPS), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)


pdf(file = "sampleClustering_preCut.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
#remove two samples E18 tube 8 and P60 tube 50
abline(h = 160, col = "red");

dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = Log2DFnoLPS[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#phenotype data:
traits <- bilbosamples
traits$samples2 <- paste(traits$group,traits$tube, sep = "_")

#match sample names and trait data:
Samples = rownames(datExpr)
traitRows = match(Samples, traits$samples2)
datTraits = traits[traitRows, ]

datTraits<- as.data.frame(datTraits)
datTraits$time <- datTraits$timepoint

#change time to numeric
datTraits$time[datTraits$time == "E18"] <- 1
datTraits$time[datTraits$time == "P4"] <- 2
datTraits$time[datTraits$time == "P14"] <- 3
datTraits$time[datTraits$time == "P60"] <- 4

datTraits$sex2 <- as.numeric(as.factor(datTraits$sex)) #M =2 = red,  F= 1 = white

#numeric traits only
numerictraits <- datTraits[,7:8]
rownames(numerictraits) <- datTraits$samples

numerictraits$time <- as.numeric(as.factor(numerictraits$time))

#We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable datTraits. Before we continue with network construction and module detection, we visualize how the clinical traits relate to the sample dendrogram.
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(numerictraits$time, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.

pdf(file = "sampleClustering_withtraits.pdf", width = 12, height = 9)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(numerictraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(numerictraits,datTraits,datExpr, file = "BilboMG-networkConstruction-inputdata.RData")

#load("BilboMG-networkConstruction-inputdata.RData")

write.csv(datExpr,"Bilbo_log2FPKM+1_cleaned.csv")
#=====================================================================================
#choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency 
# Choose a set of soft-thresholding powers
powers = c(c(1:30))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed",corFnc="bicor" ,verbose = 5)
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
#lowest power for which the scale-free topology fit index of > 0.8 ->10

#=====================================================================================
#Run this section on an external server, save the network data and load into R studio
#=====================================================================================
#Co-expression similarity and adjacency
library(WGCNA)
allowWGCNAThreads()

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#load data
load("BilboMG-networkConstruction-inputdata.RData")


Allnet = blockwiseModules(datExpr, maxBlockSize = 50000,
                          power = 10, TOMType = "signed", minModuleSize = 200,
                          networkType = "signed",
                          corType = "bicor", #biweight midcorrelation
                          maxPOutliers = 0.05, #forces bicor to never regard more than the specified proportion of samples as outliers.
                          reassignThreshold = 0,
                          numericLabels = TRUE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = "TOM-BilboMGExpression",
                          verbose = 3)

save(Allnet, file = "NetworkConstruction-auto_BilboMGExpression.RData")
#=====================================================================================
#back in Rstudio:
load("NetworkConstruction-auto_BilboMGExpression.RData")

table(Allnet$colors)

# Convert labels to colors for plotting
AllnetModuleColors = labels2colors(Allnet$colors)


AllnetTree = Allnet$dendrograms[[1]]

#plot the gene dendrogram and the corresponding module colors
sizeGrWindow(8,6);

pdf(file = "AllnetDendrogram_bilboMGtimecourse.pdf", wi = 8, he = 6)
plotDendroAndColors(AllnetTree, AllnetModuleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "All Samples Cluster Dendrogram")
dev.off()


#=====================================================================================

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = AllnetModuleColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

pdf(file = "Pretrim_allnet_dendrogram.pdf", wi = 9, he = 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=0.25, col = "red")

dev.off()



#=====================================================================================
#merge down to fewer modules based on plot height
MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(datExpr, AllnetModuleColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#replot
sizeGrWindow(12, 9)
pdf(file = "Allnet_postmerg.25MEClustering.pdf", wi = 9, he = 6)

plotDendroAndColors(AllnetTree, mColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "All Samples Cluster Dendrogram")
dev.off()

#=====================================================================================
# Rename to moduleColors
moduleColors = mColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, AllnetTree, datTraits, numerictraits,datExpr,Allnet,file = "allsamplesbilboMGtimecourse-networkpostmerge.25.RData")

load("allsamplesbilboMGtimecourse-networkpostmerge.25.RData")
#=====================================================================================
#repeat clustering of new MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
pdf(file = "Pretrim_allnet_dendrogram_postmerge.pdf", wi = 9, he = 6)
plot(METree, main = "Clustering of module eigengenes after merging modules",
     xlab = "", sub = "")

dev.off()

#=====================================================================================
#=====================================================================================
#correlate eigengenes with external traits and look for the most signicant associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, numerictraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#FDR correct the pvalues
#sapply(pval,p.adjust,method="fdr") #per column

#FDR correction on entire matrix
FDR <- matrix(p.adjust(as.vector(moduleTraitPvalue), method='fdr'),ncol=ncol(moduleTraitPvalue))
colnames(FDR) <- colnames(moduleTraitPvalue)
rownames(FDR) <- rownames(moduleTraitPvalue)


corout <- cbind(moduleTraitCor,FDR)
write.csv(corout, "PearsonsCorrelations_mod-trait.csv",row.names = T)

#We color code each association by the correlation value:

sizeGrWindow(10,6)

pdf(file = "module-traitrelationshipsFDRcorrected.pdf", wi = 8.5, he = 11)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\ (", #pearson correlation coefficient, space, (FDR corrected pvalue)
                    signif(FDR, 4), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(numerictraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

#=====================================================================================

#quantify associations of individual genes with our trait of interest (years) by dening Gene Signicance GS as
#(the absolute value of) the correlation between the gene and the trait

#For each module, we also dene a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression prole. This
#allows us to quantify the similarity of all genes on the array to every module.

# Define variable age containing the years column of datTrait
age = as.data.frame(numerictraits$time);
names(age) = "age"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

##module membership MM 
#correlation between expression of each gene and MEs (uncorrected p values)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#or
# calculate the module membership values (aka. module eigengene based
# connectivity kME):
#datKME = signedKME(PFCdatExpr, MEs)

#correlate each gene expression with variable of interest
geneTraitSignificance = as.data.frame(cor(datExpr, age, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(age), sep="");
names(GSPvalue) = paste("p.GS.", names(age), sep="")


#=====================================================================================
#Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high signicance for weight as well as high module
#membership in interesting modules.


#As an example, we look at the brown module that has the highest association with age.
module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for age",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#Clearly, GS and MM are highly correlated, illustrating that genes highly significantly associated with a trait are 
#often also the most important (central) elements of modules associated with the trait.

#=====================================================================================

#We have found modules with high association with our trait of interest, and have identified their central players by the Module Membership measure.
#We now merge this statistical information with gene annotation and write out a file that summarizes the most important results and can be inspected in standard spreadsheet software 
# Create the starting data frame
Ensemble = colnames(datExpr)

df <- DF2 %>% dplyr::select(Gene.stable.ID.1,HGNC.symbol,NCBI.gene.ID) %>% distinct()

#get matching data for each ensemble id from the row data from Brainspain input
ids <- df[which(df$Gene.stable.ID.1 %in% Ensemble),]

#match order

idRows = match(Ensemble, ids$Gene.stable.ID.1)
ids2 = ids[idRows, ]


geneInfo0 = data.frame(Ensemble = colnames(datExpr),
                       geneSymbol = ids2$HGNC.symbol,
                       entrez_id = ids2$NCBI.gene.ID,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue,
                       t(datExpr))
# Order modules by their significance for age
modOrder = order(-abs(cor(MEs, age, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}



write.csv(geneInfo0, file = "geneInfo_BilboMGtimecourse.csv")

#geneInfo0 <- read.csv("geneInfo_BilboMGtimecourse.csv")

#=====================================================================================
#output a character vector of genes, where the genes are the hub gene picked for each module, and the names correspond to the module in which each gene is a hub.

hubs <- chooseTopHubInEachModule(datExpr,moduleColors, type = "signed",power =10) #https://support.bioconductor.org/p/46342/
#power of 2 for unsigned, 4 for signed
write.csv(hubs,"TopHubInEachModule.csv")

#These genes represent the top 10 genes per module based on kME  
topGenesKME = NULL
for (i in 1:length(colnames(geneModuleMembership))){
  KMErank = geneModuleMembership[(order(-geneModuleMembership[,i])),] #order by column, - decreasing
  KMErank$Ensemble <- rownames(KMErank)
  GenesKME = KMErank[c(1:10),c(i,10)] #where column 12 is the ensemble id
  
  #get gene info
  topKMEinfo <- ids2[which(ids2$Gene.stable.ID.1 %in% rownames(GenesKME)),]
  
  #get kME
  merge <- merge(topKMEinfo,GenesKME,by.x="Gene.stable.ID.1", by.y = "Ensemble")
  colnames(merge)[4] <- c("kME")
  
  merge$module <- substr(colnames(geneModuleMembership[i]),3,nchar(colnames(geneModuleMembership[i])))
  
  topGenesKME = rbind(topGenesKME,merge)
}

write.csv(topGenesKME,"Top10HubInEachModule_KME.csv")


#topGenesKME <- read.csv("TopHubInEachModule.csv")

#=====================================================================================
#Extract modules
module_colors = setdiff(unique(moduleColors), "grey")
# for (color in module_colors){
#   module=datExpr[,moduleColors==module_colors]
#   write.table(module, paste("module_",color, "log2RPKM+1expression_DLPFC.txt",sep=""), sep="\t", row.names=T, col.names=T,quote=FALSE)
#   
# }

#Look at expression patterns of these genes, as they are clustered
#heatmap colors:
#install.packages("RColorBrewer")
library("RColorBrewer")
library(gplots)

#gene-module info

heatDF <- geneInfo0 %>% filter(moduleColor != "grey") %>% arrange(moduleColor) %>% dplyr::select(-X)
heatDF$moduleColor <- factor(heatDF$moduleColor)

data <- as.matrix(heatDF[,7:51])
rownames(data) <- heatDF$Ensemble


#data <- t(m)
library(pheatmap)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:
annotation_col <- data.frame(
  #Sample = factor(d$donor_name),
  Age = factor(numerictraits$time, levels = unique(numerictraits$time)),
  Sex = factor(numerictraits$sex2))

rownames(annotation_col) = colnames(data)

head(annotation_col)

annotation_row <- data.frame(heatDF$moduleColor)

rownames(annotation_row) = rownames(data)

colnames(annotation_row) <- c("module")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- sample(col_vector, 4)
names(Var1) <- unique(annotation_col$Age)

Var2        <-  c("cornflowerblue","darksalmon")
names(Var2) <- c("2","1")

Var3 <- as.character(unique(heatDF$moduleColor))
names(Var3) <-  unique(heatDF$moduleColor)

anno_colors <- list(Age = Var1, Sex = Var2, moduleColor = Var3)

pdf(file = "Bilbomodule-Heatmap_log2RPKM+1_zscore.pdf", wi = 20, he = 16)
pheatmap::pheatmap(data, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   color = my_palette, 
                   fontsize = 10,
                   fontsize_row=6, 
                   show_rownames = F,
                   fontsize_col = 10,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "MG Developmental Timecourse Gene Expression Modules")

dev.off()

#=====================================================================================
#average heatmap by age
heatDF2 <- heatDF[,c(1,4,7:51)]

heatDF2 <- heatDF2 %>% gather(samples, log2exp, 3:ncol(heatDF2))
headDF3 <- merge(heatDF2,datTraits, by = "samples", all.x=T)
headDF3$group <- paste(headDF3$time, headDF3$sex)

heatDF4 <- headDF3 %>% group_by(moduleColor,Ensemble,group) %>% summarize(meanlog2exp = mean(log2exp)) %>%
  spread(group,meanlog2exp)%>% arrange(moduleColor)

matrix <- as.matrix(heatDF4[,3:ncol(heatDF4)])
rownames(matrix) <- heatDF4$Ensemble


my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

#define column groups:
age <- substr(colnames(matrix),1,1)
sex <- substr(colnames(matrix),3,3)

annotation_col <- data.frame(
  Age = factor(age),
  Sex = factor(sex))

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

annotation_row <- data.frame(heatDF4$moduleColor)

rownames(annotation_row) = rownames(matrix)

colnames(annotation_row) <- c("moduleColor")

head(annotation_row)

# Specify colors
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))


# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- sample(col_vector, 4)
names(Var1) <- unique(annotation_col$Age)

Var2        <-  c("cornflowerblue","darksalmon")
names(Var2) <- c("M","F")

Var3 <- as.character(unique(heatDF4$moduleColor))
names(Var3) <-  unique(heatDF4$moduleColor)

anno_colors <- list(Age = Var1, Sex = Var2, moduleColor = Var3)

pdf(file = "BilboMGtimecoursemodule-Heatmap_log2RPKM+1_zscore_noReplicates.pdf", wi = 11, he = 8.5)
pheatmap::pheatmap(matrix, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   color = my_palette, 
                   show_rownames = F,
                   fontsize = 12,
                   fontsize_row=6, 
                   fontsize_col = 12,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "MG Developmental Timecourse Gene Expression Modules")

dev.off()


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

#data:
heatDF4 <- headDF3 %>% group_by(moduleColor,Ensemble,group) %>% summarize(meanlog2exp = mean(log2exp)) %>%
  spread(group,meanlog2exp)%>% arrange(moduleColor)

matrix <- as.matrix(heatDF4[,3:ncol(heatDF4)])
rownames(matrix) <- heatDF4$Ensemble

#scale by row:
merPFCmatrix <- scale_mat(matrix, scale = "row")
merPFCmatrix <- as.data.frame(merPFCmatrix)

merPFCmatrix$module <- heatDF4$moduleColor

merPFCDF <- merPFCmatrix %>% gather(age, scaled_expression, 1:(ncol(merPFCmatrix)-1)) %>% group_by(module,age) %>%
  summarize(meanscaled_expression = mean(scaled_expression))

merPFCDF <- as.data.frame(merPFCDF)


#set colors:
my_palette2 <- unique(merPFCDF$module)

my_palette2 <- col2hex(my_palette2)

pdf(file = "MGmodules_lineplot_log2RPKM+1_scaled_noReplicates.pdf", wi = 11, he = 8.5)

p <- ggplot(merPFCDF, aes(x=age, y=meanscaled_expression, group=module,colour=module)) +
  geom_line(size=2) +
  scale_color_manual(values = my_palette2)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))+
  xlab(label = c("Developmental Age")) +
  ylab(label = c("Scaled log2(RPKM+1)"))

p

dev.off()

#facet by module
theme_set(theme_classic(base_size = 24)) 
pdf(file = "MGmodule_lineplot_log2RPKM+1_scaled_noReplicates.pdf", wi = 8, he = 18)

p <- ggplot(merPFCDF, aes(x=age, y=meanscaled_expression, group = module, colour=module)) +
  facet_wrap(~module,ncol = 1) +
  geom_line(size=2) +
  scale_color_manual(values = my_palette2)+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        legend.position="none",
        strip.text.x = element_text(size = 30))+
  xlab(label = c("Developmental Age")) +
  ylab(label = c("Scaled log2(RPKM+1)"))

p

dev.off()

#blue moduel only
blueDF <- merPFCDF %>% filter(module == "blue")
theme_set(theme_classic(base_size = 24)) 
pdf(file = "BlueMGmodule_lineplot_log2RPKM+1_scaled_noReplicates.pdf", wi = 10, he = 6)

p <- ggplot(blueDF, aes(x=age, y=meanscaled_expression, group = module, colour=module)) +
  facet_wrap(~module,ncol = 1) +
  geom_line(size=2) +
  scale_color_manual(values = "blue")+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        legend.position="none",
        strip.text.x = element_text(size = 30))+
  xlab(label = c("Developmental Age")) +
  ylab(label = c("Scaled log2(RPKM+1)"))

p

dev.off()

#yellow module only
yellowDF <- merPFCDF %>% filter(module == "yellow")
theme_set(theme_classic(base_size = 24)) 
pdf(file = "YellowMGmodule_lineplot_log2RPKM+1_scaled_noReplicates.pdf", wi = 10, he = 6)

p <- ggplot(yellowDF, aes(x=age, y=meanscaled_expression, group = module, colour=module)) +
  facet_wrap(~module,ncol = 1) +
  geom_line(size=2) +
  scale_color_manual(values = "yellow")+
  theme(axis.text.x  = element_text(angle=45, vjust=0.5),
        strip.background = element_rect(fill="#B0C4DE"),
        legend.position="none",
        strip.text.x = element_text(size = 30))+
  xlab(label = c("Developmental Age")) +
  ylab(label = c("Scaled log2(RPKM+1)"))

p

dev.off()




#=====================================================================================
#GO term enrichment:https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/
#source("http://bioconductor.org/biocLite.R");
#biocLite(c("impute", "AnnotationDBI", "GO.db", "org.Hs.eg.db", "org.Mm.eg.db"))

#install.packages("~/Box Sync/LaSalle Lab/Experiments/WGCNA/anRichment_0.82-1.tgz", repos = NULL, type = .Platform$pkgType, lib="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
# install.packages("~/Box Sync/LaSalle Lab/Experiments/WGCNA/SCSBrainCellTypeCollection_1.00.tgz", repos = NULL, type = .Platform$pkgType, lib="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

library(anRichment)
options(stringsAsFactors = FALSE)

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/Mouse_MicrogliaWGCNA")
geneInfo0 <- read.csv("geneInfo_BilboMGtimecourse.csv")
geneInfo <- geneInfo0
###################################################################
GOcollection = buildGOcollection(organism = "human")

# evaluates the enrichment of the gene modules in the collection of GO terms

GOenrichment = enrichmentAnalysis(
   classLabels = geneInfo0$moduleColor, identifiers = geneInfo0$entrez_id,
   refCollection = GOcollection,
   useBackground = "intersection", #intersection of genes in identifiers and the organism database
   threshold = 0.05,
   thresholdType = "FDR",
   getOverlapEntrez = TRUE,
   getOverlapSymbols = TRUE,
   maxReportedOverlapGenes = 500,
   ignoreLabels = "grey") #ignore genes not assigned to modules (grey)
# 
 collectGarbage()
# 
 names(GOenrichment)
# 
 #enrichment results are summarized in the component enrichmentTable:
 names(GOenrichment$enrichmentTable)
# 

 write.csv(GOenrichment$enrichmentTable, file = "GOenrichment-enrichmentTable.csv")
# 
# 
 ###################################################################
 #NDD DMR genes 
 ###################################################################
 #make internal collections for DMR associated genes and genetic mutations, etc
 
 #Custom lists of ASD relevant genes
 load("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/CustomCollections/CustomGeneLists/NDD_DMRgenes_collection.Rdata")
 knownGroups(PLcollection, sortBy = "size")

 #calculate the enrichment using the modified collection:
 Combinedenrichment = enrichmentAnalysis(
   classLabels = geneInfo$moduleColor, 
   identifiers = geneInfo$entrez_id,
   refCollection = PLcollection,
   useBackground = "given", 
   threshold = 1.0,
   thresholdType = "FDR",
   getOverlapEntrez = FALSE,
   getOverlapSymbols = TRUE,
   entrySeparator = ",",
   maxReportedOverlapGenes = 1000,
   ignoreLabels = "grey") #ignore genes not assigned to modules (grey)
 
 enrichment_table <- Combinedenrichment$enrichmentTable
# enrichment_table2 <-  enrichment_table %>% filter(class != "grey") %>% filter(class != "red")
# enrichment_table2$FDRpval <- p.adjust(as.numeric(enrichment_table2$pValue),method="fdr")
 write.xlsx(Combinedenrichment$enrichmentTable, file = "NDD_DMRgenesCustomListEnrichment-MGWGCNA.xlsx")
 
###########################################################################################################
 ###################################################################
 #make internal collections for DMR associated genes and genetic mutations, etc
 
 #Custom lists of microglia relevant genes
 load("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/8_2018_4CpGDMRs/CustomCollections/CustomGeneLists/MicrogliaHg38listsCollection.RData")
 
 knownGroups(MG_genelists, sortBy = "size")
 dataSetNames(MG_genelists)

 
 
 ######################################################################################################################################
 
############################################################################################################
MAAcollection = subsetCollection(MG_genelists, tags = c("MAA DE genes","MAA DMR genes","MAA HyperDMR genes","MAA HypoDMR genes"))

dataSetNames(MAAcollection)


#calculate the enrichment using the combined collection:
Combinedenrichment = enrichmentAnalysis(
  classLabels = geneInfo$moduleColor, identifiers = geneInfo$entrez_id,
  refCollection = MAAcollection,
  useBackground = "intersection", #intersection of genes in identifiers and the organism database
  threshold = 1.00,
  thresholdType = "FDR",
  getOverlapEntrez = FALSE,
  getOverlapSymbols = TRUE,
  entrySeparator = ", ",
  maxReportedOverlapGenes = 1000,
  ignoreLabels = "grey") #ignore genes not assigned to modules (grey)



#enrichment results are summarized in the component enrichmentTable:
enrichment_table <- Combinedenrichment$enrichmentTable

write.csv(Combinedenrichment$enrichmentTable, file = "MAAgenelistCustomenrichment_MGWGCNA.csv")


######################################################################################################################################
#custom collections
#Enrichment for Microglial gene lists

#calculate the enrichment using the combined collection:
Combinedenrichment = enrichmentAnalysis(
  classLabels = geneInfo$moduleColor, identifiers = geneInfo$entrez_id,
  refCollection = MG_genelists,
  useBackground = "intersection", #intersection of genes in identifiers and the organism database
  threshold = 0.05,
  thresholdType = "FDR",
  getOverlapEntrez = FALSE,
  getOverlapSymbols = TRUE,
  entrySeparator = ",",
  maxReportedOverlapGenes = 1000,
  ignoreLabels = "grey") #ignore genes not assigned to modules (grey)



#enrichment results are summarized in the component enrichmentTable:
enrichment_table <- Combinedenrichment$enrichmentTable

write.csv(Combinedenrichment$enrichmentTable, file = "MGgenelistCustomenrichment_MGWGCNA.csv")

###########################################################################################################
#GO terms:
#Creating enrichment labels for classes
# The function takes as input the enrichment table returned by enrichmentAnalysis in the enrichmentTable component, and lets the user specify the various columns needed to put together the enrichment label

#The function returns a data frame with information about highest enriched terms in each of the groups specified
#in focusOnGroups (“all” is a special keyword meaning all terms irrespective of what group they belong to).

eLabels = enrichmentLabels(
  GOenrichment$enrichmentTable,
  focusOnGroups = c("GO"),
  groupShortNames = c("GO"),
  minSize = 0.05,
  pValueThreshold = 0.05,
  numericClassLabels = FALSE)


GO <- eLabels %>% dplyr::select(class,highestEnrichedSet.GO)

GObest <- GOenrichment$enrichmentTable[which(GOenrichment$enrichmentTable$dataSetName %in% GO$highestEnrichedSet.GO),]

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



pdf(file = "BiboMGmodule-GOterm_relationshipsFDRcorrected.pdf", wi = 8.5, he = 11)

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




