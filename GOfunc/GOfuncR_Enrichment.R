##########################################################################
#AVC a.ciernia@gmail.com
#8/2018
#GOfuncR: Gene Ontology Enrichment Using FUNC
#http://bioconductor.org/packages/devel/bioc/vignettes/GOfuncR/inst/doc/GOfuncR.html#hypergeometric-test-with-genomic-regions-as-input
##########################################################################
library(dplyr)
library(tidyr)
library(xlsx)
library(rtracklayer)
library(ChIPseeker)
library(tools)

#developer version
#install.packages("devtools")
library(sm)
library(devtools)
#install_github("sgrote/GOfuncR")

library(GOfuncR)
#human hg38:
#source("https://bioconductor.org/biocLite.R")
#biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
#biocLite("org.Hs.eg.db")
#source("https://bioconductor.org/biocLite.R")
#biocLite("Homo.sapiens")

library(org.Hs.eg.db)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#########################################################################
#build custom coordinates: gene start and stop + 10kb up and downstream
#########################################################################
library(GenomicRanges)


setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/GOenrichment/gofunc")
#biomart: gene start, gene end, chr and strand
allgenes <- read.delim("Hg38_Genestart_stop.txt",header=T)
allgenes$chr <- paste("chr",allgenes$Chromosome.scaffold.name,sep="")
allgenes$Strand <- gsub("-1","-",allgenes$Strand)
allgenes$Strand <- gsub("1","+",allgenes$Strand)
allgenes <- allgenes[,c(7,4,3,6,1,5)]
names(allgenes) <- c("chr","start","stop","strand","ensembleid","geneid")

genedf <- makeGRangesFromDataFrame(allgenes,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=FALSE,
                                   seqinfo=NULL,
                                   seqnames.field=c("seqnames", "seqname",
                                                    "chromosome", "chrom",
                                                    "chr", "chromosome_name",
                                                    "seqid"),
                                   start.field="start",
                                   end.field=c("end", "stop"),
                                   strand.field="strand",
                                   starts.in.df.are.0based=T)


#add 10000 bp up and downstream of gene start and gene end
#https://support.bioconductor.org/p/78652/
extend <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

newgr <- extend(genedf,upstream=10000,downstream=10000)

head(newgr)

#make to dataframe
df <- data.frame(seqnames=seqnames(newgr),
                 starts=start(newgr)-1,
                 ends=end(newgr),
                 names=c(rep(".", length(newgr))),
                 scores=c(rep(".", length(newgr))),
                 strands=strand(newgr),
                 ensemble=elementMetadata(newgr)$ensembleid,
                 geneid=elementMetadata(newgr)$geneid)

#select: ensembleid, chr,start and end 
Hg38genecord <- df[,c(8,1,2,3)]
colnames(Hg38genecord) <- c("gene","chr","start","end")
Hg38genecord <- Hg38genecord[!duplicated(Hg38genecord$gene),]

Hg38genecord$gene <- as.character(Hg38genecord$gene)
Hg38genecord$chr <- as.character(Hg38genecord$chr)
Hg38genecord$start <- as.integer(Hg38genecord$start)
Hg38genecord$end <- as.integer(Hg38genecord$end)

save(newgr,df,Hg38genecord,file="GRanges_Hg38Genes_10kbextension.RData")

load("GRanges_Hg38Genes_10kbextension.RData")
##########################################################################
#get input data
##########################################################################

#ASD DMRs
ASD_DMRs <- read.delim("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/ASD/ASD_DMRs.bed",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(ASD_DMRs) <- c("chr","start","end")
ASD_DMRs$chr2 <-substr(ASD_DMRs$chr,4,5)
#chr:start-stop, where start always has to be smaller than stop.
ASD_DMRs$cord <- paste(ASD_DMRs$chr2,":",ASD_DMRs$start,"-",ASD_DMRs$end, sep="")

BG_DMRs <- read.delim("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/ConsensusBackgroundRegions_sort.bed",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(BG_DMRs) <- c("chr","start","end")
BG_DMRs$chr2 <-substr(BG_DMRs$chr,4,5)
#chr:start-stop, where start always has to be smaller than stop.
BG_DMRs$cord <- paste(BG_DMRs$chr2,":",BG_DMRs$start,"-",BG_DMRs$end, sep="")

#add in DMRs removed by merging:
ASD_DMRs_notBG <- ASD_DMRs[!(ASD_DMRs$cord %in% BG_DMRs$cord),] 

ASD_DF <- rbind(BG_DMRs,ASD_DMRs_notBG)
ASD_DF$is_candidate <- 0
ASD_DF$is_candidate[which(ASD_DF$cord %in% ASD_DMRs$cord)] <- "1"
table(ASD_DF$is_candidate)

ASD_DF <- ASD_DF[,c(5:6)]

#The input for the hypergeometric test is a dataframe with two columns: 
#(1) a column with gene-symbols and 
#(2) a binary column with 1 for a candidate gene and 0 for a background gene.

names(ASD_DF) <-c("regions","is_candidate")

ASD_DF$is_candidate <- as.numeric(ASD_DF$is_candidate)

#run Enrichment
res_hyper_ASD = go_enrich(ASD_DF, 
                            test='hyper',
                            n_randsets=1000,
                            regions=TRUE,
                            gene_coords = Hg38genecord,
                            circ_chrom=TRUE,
                            orgDb='org.Hs.eg.db',
                            txDb='TxDb.Hsapiens.UCSC.hg38.knownGene')

stats_region = res_hyper_ASD[[1]]
head(stats_region)
stats_region$FDR <- p.adjust(stats_region$raw_p_overrep,method='fdr')

write.csv(stats_region,"ASDBG_GO.csv")
save(res_hyper_ASD,file="ASD_GOdata.RData")



#########################################################################
#Dup15q
#########################################################################
#Dup15q DMRs
Dup15q_DMRs <- read.delim("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/Dup15q/Dup15q_DMRs.bed",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(Dup15q_DMRs) <- c("chr","start","end")
Dup15q_DMRs$chr2 <-substr(Dup15q_DMRs$chr,4,5)
#chr:start-stop, where start always has to be smaller than stop.
Dup15q_DMRs$cord <- paste(Dup15q_DMRs$chr2,":",Dup15q_DMRs$start,"-",Dup15q_DMRs$end, sep="")

BG_DMRs <- read.delim("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/ConsensusBackgroundRegions_sort.bed",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(BG_DMRs) <- c("chr","start","end")
BG_DMRs$chr2 <-substr(BG_DMRs$chr,4,5)
#chr:start-stop, where start always has to be smaller than stop.
BG_DMRs$cord <- paste(BG_DMRs$chr2,":",BG_DMRs$start,"-",BG_DMRs$end, sep="")

#add in DMRs removed by merging:
Dup15q_DMRs_notBG <- Dup15q_DMRs[!(Dup15q_DMRs$cord %in% BG_DMRs$cord),] 

Dup15q_DF <- rbind(BG_DMRs,Dup15q_DMRs_notBG)
Dup15q_DF$is_candidate <- 0
Dup15q_DF$is_candidate[which(Dup15q_DF$cord %in% Dup15q_DMRs$cord)] <- "1"
table(Dup15q_DF$is_candidate)

Dup15q_DF <- Dup15q_DF[,c(5:6)]

#The input for the hypergeometric test is a dataframe with two columns: 
#(1) a column with gene-symbols and 
#(2) a binary column with 1 for a candidate gene and 0 for a background gene.

names(Dup15q_DF) <-c("regions","is_candidate")

Dup15q_DF$is_candidate <- as.numeric(Dup15q_DF$is_candidate)

res_hyper_Dup15q = go_enrich(Dup15q_DF, 
                             test='hyper',
                             n_randsets=1000,
                             regions=TRUE,
                             gene_coords = Hg38genecord,
                             circ_chrom=TRUE,
                             orgDb='org.Hs.eg.db',
                             txDb='TxDb.Hsapiens.UCSC.hg38.knownGene')

stats_region = res_hyper_Dup15q[[1]]
head(stats_region)
stats_region$FDR <- p.adjust(stats_region$raw_p_overrep,method='fdr')

write.csv(stats_region,"Dup15q_GO.csv")
save(res_hyper_Dup15q,file="res_hyper_Dup15q.RData")




#########################################################################
#########################################################################
#RTT
#########################################################################
#RTT DMRs
RTT_DMRs <- read.delim("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/Rett/RTT_DMRs.bed",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(RTT_DMRs) <- c("chr","start","end")
RTT_DMRs$chr2 <-substr(RTT_DMRs$chr,4,5)
#chr:start-stop, where start always has to be smaller than stop.
RTT_DMRs$cord <- paste(RTT_DMRs$chr2,":",RTT_DMRs$start,"-",RTT_DMRs$end, sep="")

BG_DMRs <- read.delim("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/ConsensusBackgroundRegions_sort.bed",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(BG_DMRs) <- c("chr","start","end")
BG_DMRs$chr2 <-substr(BG_DMRs$chr,4,5)
#chr:start-stop, where start always has to be smaller than stop.
BG_DMRs$cord <- paste(BG_DMRs$chr2,":",BG_DMRs$start,"-",BG_DMRs$end, sep="")

#add in DMRs removed by merging:
RTT_DMRs_notBG <- RTT_DMRs[!(RTT_DMRs$cord %in% BG_DMRs$cord),] 

RTT_DF <- rbind(BG_DMRs,RTT_DMRs_notBG)
RTT_DF$is_candidate <- 0
RTT_DF$is_candidate[which(RTT_DF$cord %in% RTT_DMRs$cord)] <- "1"
table(RTT_DF$is_candidate)

RTT_DF <- RTT_DF[,c(5:6)]

#The input for the hypergeometric test is a dataframe with two columns: 
#(1) a column with gene-symbols and 
#(2) a binary column with 1 for a candidate gene and 0 for a background gene.

names(RTT_DF) <-c("regions","is_candidate")

RTT_DF$is_candidate <- as.numeric(RTT_DF$is_candidate)

res_hyper_RTT = go_enrich(RTT_DF, 
                          test='hyper',
                          n_randsets=1000,
                          regions=TRUE,
                          gene_coords = Hg38genecord,
                          circ_chrom=TRUE,
                          orgDb='org.Hs.eg.db',
                          txDb='TxDb.Hsapiens.UCSC.hg38.knownGene')

stats_region = res_hyper_RTT[[1]]
head(stats_region)
stats_region$FDR <- p.adjust(stats_region$raw_p_overrep,method='fdr')

write.csv(stats_region,"RTT_GO.csv")
save(res_hyper_RTT,file="res_hyper_RTT.RData")



