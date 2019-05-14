#author: Annie Vogel Ciernia
#a.ciernia@gmail.com
#Gene Overlap Enrichment testing by permutation test with correction for gene length
#takes as input: Target gene list for overlapping, DMR genes, and BG DF
#each must have a column for gene length titled "length" and "ensembl_gene_id"
#BG data is first binned into 50 bins by gene length (50 is arbitrary)
#the number of target list genes within each bin is then found
#A random subsample of the BG is then made for each bin so that the gene lengths of the target list is matched
#calculate a t-test between the lengths in susample DF and the actual target gene list
#if the t-test is not significant then it keeps the subsample and performs an intersection with the DMR genes
#saves the counts of the overlaps and the overlapping ids
#repeats this bin-based gene length matched resampling 1000 times
#compares the actual overlap between the target list and DMR genes to the random resampling overlaps = permutation p-value
#outputs a histogram of the background resampling distribution with a line for the actual overlap

##############################################################################################################
library(dplyr)
library(tidyr)
library(cowplot)
library(gplots)
#library(openxlsx)
library(stringr)
library(parallel)
##############################################################################################################
#based on this:
#Gene list overlaps by permutation with correction for gene length
#https://www.biostars.org/p/209717/
# geneListA # vector of genes from organism A
# geneListB # vectof of genes from organism B
# N <- 1000 # simulation iteration
# 
# gA # list of genes of interest from organism A
# gB # list of genes of interest from organism B
# 
# inter.r <- length(intersect(gA,gB)) # number of gene in common 
# 
# getSimIntersect <- function(gA,gB,geneListA,geneListB){
#   gA.sim <- sample(geneListA,length(gA))
#   gB.sim <- sample(geneListB,length(gB))
#   return(length(intersect(gA.sim,gB.sim)))
# }
# 
# 
# #perform purmutations
# inter.sim <- replicate(N,getSimIntersect(gA,gB,geneListA,geneListB))
# 
# #permuted pvalue
# #p value is not = 0 . 
# #If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
# p <- sum(inter.sim>=inter.r)/N

#and this:
#https://stats.stackexchange.com/questions/94716/sampling-data-to-have-specific-mean-and-standard-deviation
# # 
# setwd("/Users/annieciernia/Box Sync/LaSalle Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/genelist_lengthcorrOverlaps")
# # # #################################################################
# # # #load in Datasets
# # # #################################################################
# # #
# # # # vector of genes from all of Hg38 (BioMart), need for getting gene length values and GC values
# load("Hg38_Conversioninfo.RData")
# HumanHg38 <- ids
# # #
# # # # list of genes associated with DMR
# DMRgenes <- read.delim("ASD_DMR_EnsemblHG38.txt",header=F)
# colnames(DMRgenes) <- "ensembl_gene_id"
# DMRlistname <- c("ASD DMRs")
# # #
# # # #background gene list. If want all hg38 then take from HumanHg38
# BG_DF <- read.csv("ConsensusBG_DMR_EnsemblHG38.txt",header=F)
# # # #fix col name
# colnames(BG_DF) <- c("ensembl_gene_id")
# # #
# # # #add gene start,stop and length columns if don't already exist
# #
# setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/GeneOverlapsGClength")
# # # #list of genes to overlap with = target genes (ie. DE expressed, genetic mutation, etc)
# OverlapList <- read.csv("Overlaplist_3_2019.csv")
# #
# OverlapList$length <- abs(OverlapList$end_position - OverlapList$start_position)
# OverlapList <- distinct(OverlapList)
# 
# #combine some datasets in overlap list
# names <- read.csv("OverlaplistNames_3_2019.csv")
# names <- names[,c(1:2)]
# OverlapList2 <- merge(OverlapList, names, by.x="List", by.y="Original.List")
# OverlapList2 <- distinct(OverlapList2)
# OverlapList <- OverlapList2
# #
# # #
# DMRgenes2 <- merge(DMRgenes,HumanHg38, by = c("ensembl_gene_id"),all.x=T)
# DMRgenes2$length <- abs(DMRgenes2$end_position - DMRgenes2$start_position)
# DMRgenes <- DMRgenes2[,c("ensembl_gene_id","length")]
# DMRgenes <- distinct(DMRgenes)
# # #
# BG_DF2 <- merge(BG_DF,HumanHg38, by = c("ensembl_gene_id"),all.x=T)
# BG_DF2$length <- abs(BG_DF2$end_position - BG_DF2$start_position)
# BG_DF2 <- BG_DF2[,c("ensembl_gene_id","length")]
# BG_DF <- distinct(BG_DF2)
# #
# # #add in GC content for DMRs and BG
# load("Ensemble_hg38genes.Rdata")
# #
# GC <- human_genes %>% dplyr::select(ensembl_gene_id,percentage_gene_gc_content)
# DMRmerge <- merge(DMRgenes,GC,by="ensembl_gene_id")
# DMRmerge <- distinct(DMRmerge)
# #
# BGmerge <- merge(BG_DF,GC,by="ensembl_gene_id")
# BGmerge <- distinct(BGmerge)
# #
# Listmerge <- merge(OverlapList,GC,by="ensembl_gene_id")
# Listmerge <- distinct(Listmerge)
# #
# #
# DMRgenes <- DMRmerge
# BG_DF <-BGmerge
# OverlapList <- Listmerge
# # # #make length measure numeric
# BG_DF$length <- as.numeric(BG_DF$length)
# DMRgenes$length <- as.numeric(DMRgenes$length)
# OverlapList$length <- as.numeric(OverlapList$length)
# names(OverlapList)[names(OverlapList) == 'List'] <- 'OriginalList'
# names(OverlapList)[names(OverlapList) == 'Listname'] <- 'List'
# 
# write.csv(OverlapList, "OverlapList2_3_31_19.csv")
# 
# #
# # # #make GC measure numeric
# BG_DF$percentage_gene_gc_content <- as.numeric(BG_DF$percentage_gene_gc_content)
# DMRgenes$percentage_gene_gc_content <- as.numeric(DMRgenes$percentage_gene_gc_content)
# OverlapList$percentage_gene_gc_content <- as.numeric(OverlapList$percentage_gene_gc_content)
# # #
# # # #make ensemble ids characters
# BG_DF$ensembl_gene_id  <- as.character(BG_DF$ensembl_gene_id)
# DMRgenes$ensembl_gene_id  <- as.character(DMRgenes$ensembl_gene_id)
# OverlapList$ensembl_gene_id  <- as.character(OverlapList$ensembl_gene_id)
# # #
# # # #split OverlapList into list of lists
# OverlapList$List <- as.character(OverlapList$List)
# OverlapList <- OverlapList[OverlapList$List !="",]
# OverlapList$List <- as.factor(OverlapList$List)
# #
# List <- split(OverlapList, OverlapList$List)
# #
# save(List, DMRgenes,BG_DF,file= "ASDoverlapInputdata.Rdata")

load("ASDoverlapInputdata.Rdata")
#setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/8_2018_4CpGDMRs/genelist_lengthcorrOverlaps/test")
#################################################################
#modified for resample from background with match of gene length
#################################################################
#input=list of genes targets
#BG = dataframe of background
#DMRgenes = DMR associated genes
#testList <- List[[1]]
#x = testList
#df = BG_DF
#DMRgenes_df=DMRgenes
#n_perm = number of permutations, no commas
#DMRlistname = c("DMR listname"), character value with name of DMR list
Overlap <- function(x,df,DMRgenes_df,n_perm,DMRlistname) {
  listname <- unique(as.character(x$List))
  print(listname)
  # Parameters
  #set overlap output so not zero rows to start
  sim_out <- data.frame(subsample_mean = NULL,
                        subsample_std = NULL,
                        overlap = NULL,
                        subsample_ensemble = NULL)
  
  #remove rows with no length or GC
  x <- x[!is.na(x$length),]
  x <- x[!is.na(x$percentage_gene_gc_content),]
  
  #keep only genes with known length and GC content:
  OverlapList_target_genelength <- nrow(x)

  subsample_size = nrow(x)
  
  #for each gene in target list find another one randomly with ~ the same size and GC:
  
  #assign 50 bins to BG by length
  df$bin <- ntile(df$length, 50)  
  
  #get range of quartiles
  df_range <- df %>% group_by(bin) %>%
    summarize(min = min(length), max= max(length)) %>% na.omit()
  
  #find where gene falls
  list_bins <- NULL
  for (i in x$ensembl_gene_id) {
    tmp <- x %>% filter(ensembl_gene_id == i)
    
    fallsinBin <- findInterval(x = tmp$length,v = df_range$min)
    
    d <- cbind(tmp,fallsinBin)
    
    list_bins <- rbind(d,list_bins)
  }
  
  #how many genes fall into each length bin?
  bin_counts <-  list_bins %>% group_by(fallsinBin) %>% summarize(count= n())
  
  #repeated subsampling from within each bin:
  repeat {
    
    #select that many genes randomly for each bin from the BG dataset
    SubsampleBG <-NULL
    
    for (i in unique(bin_counts$fallsinBin)) {
      
      subsample_size_bin <- bin_counts %>% filter(fallsinBin == i) %>% dplyr::select(count)
      subsample_size_bin <- as.numeric(subsample_size_bin)
      #get BG for the bin
      BGselected_bin <- df %>% filter(bin == i) 
      
      # Pick subsamples from the bin 
      subsample = sample(nrow(BGselected_bin),size=subsample_size_bin) #selects random row numbers
      subsampledf = BGselected_bin[subsample,] #select values based on random row selection
      
      SubsampleBG <- rbind(subsampledf,SubsampleBG)  
    }
    
    #calculate a t-test between the random sample and the target overlap sample for length 
    #If not significantly different, keep the distribution for testing.
    
    ttest <- t.test(x$length, SubsampleBG$length)
    pval <- ttest$p.value
    print(pval)
    
    #and for GC
    #ttestGC <- t.test(x$percentage_gene_gc_content, SubsampleBG$percentage_gene_gc_content)
   # pvalGC <- ttestGC$p.value
    #print(pvalGC)
    
    #if the simulation distribution is not different from the target distribution keep it:
    #if(pval>0.05 & pvalGC>0.05) {
    if(pval>0.05) {
      #overlap target gene ids and subsample gene ids
      simoverlap <- length(intersect(DMRgenes_df$ensembl_gene_id,SubsampleBG$ensembl_gene_id))
      
      
      #output the subsample information and overlap data:
      simdata <- data.frame(subsample_mean = mean(SubsampleBG$length),
                            subsample_std = sd(SubsampleBG$length),
                            overlap = as.numeric(simoverlap),
                            subsample_ensemble = ensembl_gene_id<- str_c(as.character(SubsampleBG$ensembl_gene_id),collapse=','))
      
      sim_out <- rbind(simdata,sim_out)
      print(simoverlap)
      
    }
    # exit if the condition is met: end when the number of kept permutations is 10,000 
    if (nrow(sim_out) >= n_perm) break
  }
  sim_out <- sim_out[complete.cases(sim_out),]
  #return(sim_out)
  
  #how much to they actually overlap?
  inter.r <- length(intersect(DMRgenes$ensembl_gene_id,x$ensembl_gene_id)) # number of gene in common between the two lists
  
  inter.r <- as.numeric(inter.r)
  
  #permuted pvalue
  #p value is not = 0 . 
  #If you do N=10,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
  p <- sum(sim_out$overlap>=inter.r)/length(sim_out$overlap)
  
  overlapids <-str_c(as.character(intersect(DMRgenes$ensembl_gene_id,x$ensembl_gene_id),collapse=','))
  
  DMRlistname <- as.character(DMRlistname)
  DF_pvales <- data.frame(DMRlist = DMRlistname,
                          OverlapList = listname,
                          overlap = inter.r,
                          pvalue = p,
                          number_permutations = n_perm,
                          Number_BG_greater_overlap = sum(sim_out$overlap>=inter.r),
                          Average_ResampleBG_overlap = mean(sim_out$overlap))
  
  print(DF_pvales)
  
  DMRlistname <- gsub(" ","_",DMRlistname)
  listname <- gsub(" ","_",listname)
  
  pdf(file=paste("Histogram",DMRlistname,"_",listname,"_p=",p,".pdf",sep=""),height = 5,width = 5,useDingbats = F)
  hist(sim_out$overlap, 
       main=paste("Histogram for ",DMRlistname," overlap ",listname,"\n p = ",p,sep=""),
       xlab="Overlap Distribution \n from Gene Length & GC Matched Resampling from Background")
  abline(v=inter.r, lwd=2, col="purple")
  dev.off()
  
  write.csv(DF_pvales,file = paste(DMRlistname,"_",listname,".csv",sep=""))
  
  return(DF_pvales)
  
}

#testList <- List[1]
#BGmatched <- lapply(testList,function(x) Overlap(x,df = BG_DF,DMRgenes_df=DMRgenes, n_perm =10, DMRlistname=c("ASD DMRs")))

cores = 8
BGmatched <- mclapply(List,function(x) Overlap(x,df = BG_DF,DMRgenes_df=DMRgenes,n_perm=1000,DMRlistname=c("ASD DMRs")),mc.cores=cores)

save(BGmatched, file= "BGmatched_DF_ASDOverlaps.Rdata")

BGmatched_DF = do.call("rbind", BGmatched)

#BGmatched_DF$FDR <- p.adjust(BGmatched_DF$pvalue,method="fdr")

write.csv(BGmatched_DF,"BGmatched_DF_ASDOverlaps.csv")

#save(BGmatched_DF,BGmatched, file= "BGmatched_DF_ASDOverlaps.Rdata")
