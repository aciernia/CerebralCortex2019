#filter closest gene output from bedtools to limit DMR-gene associations to within genebodies and 5kb up or downstream
#basic plots and stats of gene assignments
#author: Annie Vogel Ciernia
#a.ciernia@gmail.com
#2/2018
##############################################################################################################
setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/gene_assignment")


#read in DMRs with nearest genes added using bedtools
#Dup15q DMRs
DMRgenes <- read.delim("Dup15q_DMRs_closestgenes.txt",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(DMRgenes) <- c("chr","DMRstart","DMRend","chrom","gene_start","gene_end","Ensbl","distance_relhg38")

DMRgenes$distance_relhg38 <- as.numeric(DMRgenes$distance_relhg38)
#filter +/- 5kb from DMR
DMRgenes$label <- c("DMR Genes")

#BG loci DMRs
BGDMRgenes <- read.delim("ConsensusBackgroundRegions_closestgenes.txt",header=F, sep="\t", stringsAsFactors=FALSE)
colnames(BGDMRgenes) <- c("chr","DMRstart","DMRend","chrom","gene_start","gene_end","Ensbl","distance_relhg38")

BGDMRgenes$distance_relhg38 <- as.numeric(BGDMRgenes$distance_relhg38)
BGDMRgenes$label <- c("BG Region Genes")

#all Genes
AllGenes <- rbind(DMRgenes,BGDMRgenes)
AllGenes10kb <- subset(AllGenes, distance_relhg38 <= 10000)
AllGenes10kb <- subset(AllGenes10kb , distance_relhg38 >= -10000)

#histogram of distance
# Basic histogram
library(ggplot2)
plot <- ggplot(AllGenes10kb, aes(x=distance_relhg38,color=label, fill=label)) + 
 # Change the width of bins
  geom_histogram(binwidth=1000) + # Change colors
  #geom_vline(aes(xintercept=0), color="blue", linetype="dashed") +
  scale_color_manual(values=c("#999999", "#E69F00"))+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_y_log10() +
  labs(title="Distance from Nearest Gene",x="Distance (bp)", y = "Log10(count)")

ggsave(plot = plot, filename="Dup15q_Log10CountHistogram_Genes10kb.pdf",width = 8, height = 5, dpi = 300, useDingbats = FALSE)

###################################################################
#annotate
###################################################################

load("Hg38_Conversioninfo.RData")
#ids2 <- subset(ids, ids$ensembl_gene_id %in% AllGenes10kb$Ensbl)

DMRinfomerge <- merge(AllGenes10kb,ids,by.x="Ensbl",by.y="ensembl_gene_id", all.x=T)
DMRinfomerge <- distinct(DMRinfomerge)


###################################################################
#write out files
DMRgenes_filter <- DMRinfomerge %>% filter(label == "DMR Genes") #2794
write.csv(DMRgenes_filter, file="Dup15qDMRgenes_within10kb.csv")

BGgenes_filter <- DMRinfomerge %>% filter(label == "BG Region Genes")#48852
write.csv(BGgenes_filter, file="Dup15q_BGregiongenes_within10kb.csv")

DMRgenes_filter <- DMRinfomerge %>% filter(label == "DMR Genes") %>% 
  dplyr::select("Ensbl") %>% unique() #2367

BGgenes_filter <- DMRinfomerge %>% filter(label == "BG Region Genes") %>% 
  dplyr::select("Ensbl") %>% unique() #18974


write.table(DMRgenes_filter,file='Dup15q_DMR_EnsemblHG38.txt',sep="\t", quote=FALSE,col.names = F,row.names = F)

#write.table(BGgenes_filter,file='Dup15q_BG_DMR_EnsemblHG38.txt',sep="\t", quote=FALSE,col.names = F,row.names = F)

##############################################################################################################
#percent of DMR associated genes per genomic feature (genebody, up or downstream)
##############################################################################################################

#distribution by location up or down 
AllDMRs <- AllGenes10kb
#label each DMR-gene by location
AllDMRs$genomiclocation <- ifelse(AllDMRs$distance_relhg38 == 0, "genebody","not")
AllDMRs$genomiclocation <- ifelse(AllDMRs$distance_relhg38 > 0, "downstream",AllDMRs$genomiclocation)
AllDMRs$genomiclocation <- ifelse(AllDMRs$distance_relhg38 < 0, "upstream",AllDMRs$genomiclocation)

#count table
AllDMRs$DMR <- paste(AllDMRs$chr, AllDMRs$DMRstart,AllDMRs$DMRend)
counts <- table(AllDMRs$label, AllDMRs$genomiclocation)
counts  <- as.data.frame(counts )

#percents of total
counts2 <- counts %>% spread(Var2,Freq)
counts2$total <- rowSums(counts2[,2:4])
counts2$percent_downstream = counts2$downstream/counts2$total *100
counts2$percent_genebody = counts2$genebody/counts2$total *100
counts2$percent_upstream = counts2$upstream/counts2$total *100

DMRcounts <- counts2 %>%
  gather(genomic_location,percent_DMRs, 6:8)

#save to file
write.csv(DMRcounts, "Dup15q_Counts_DMR-gene_perlocation.csv",quote=F, row.names = F)

#order for upstream, genebody,downstream
DMRcounts$genomic_location <- factor(DMRcounts$genomic_location)
levels(DMRcounts$genomic_location) <-  c("percent_upstream","percent_genebody","percent_downstream")

library(ggplot2)
library(cowplot)

ggplot(DMRcounts,aes(x=Var1,y=percent_DMRs,fill=genomic_location)) + 
  geom_bar(stat="identity", position = "dodge") + 
  #ggnice() + 
  labs(title="Genomic Context of DMR Associated Genes") +
  scale_x_discrete(name = "Genomic Feature") +
  ylab("Percentage") +
  scale_y_continuous(limits = c(0, 100))+
  theme(text = element_text(size=20)) 

ggsave(filename="Dup15q_DMR-genes_GeneomicLocation.pdf",width = 8, height = 5, dpi = 300,  useDingbats = FALSE)

##################################################################################

#Fishers exact test enrichment compared to BG loci:
#two sided fisher's exact test
call <- unique(AllDMRs$genomiclocation)
counts <- table(AllDMRs$label, AllDMRs$genomiclocation)
counts  <- as.data.frame(counts )

#DMRs:
out <-NULL
for (i in call) {
  print(i)
  
  DMRs <- filter(counts,Var1=="DMR Genes")
  DMRs_total <- sum(DMRs$Freq)
  tmp <- filter(DMRs,Var2 ==i)
  roiplus <- tmp$Freq
  roineg <- DMRs_total - tmp$Freq
  
  BG <- filter(counts,Var1=="BG Region Genes")
  BG_total <- sum(BG$Freq)
  tmp2 <- filter(BG,call ==i)
  BGplus <- tmp2$Freq
  BGneg <- BG_total - tmp2$Freq
  
  #contingency table
  mytable2 <- rbind(c(BGneg,roineg),c(BGplus,roiplus))
  test <- fisher.test(mytable2,alternative='two.sided')

  #output
  pvalue <- as.numeric(as.character(test$p.value))
  HigherCI <- as.numeric(as.character(test$conf.int[1]))
  LowerCI <- as.numeric(as.character(test$conf.int[2]))
  OddsRatio <- as.numeric(as.character(test$estimate))
  location <- paste(i)
  results <- cbind(location,pvalue,HigherCI,LowerCI,OddsRatio)
  out <- rbind(out,results) 
  
}

out <- as.data.frame(out)
out$pvalue <- as.numeric(as.character(out$pvalue))


#FDR adjustment
out$adj.pvalue <- p.adjust(out$pvalue, method="fdr") #FDR


#write to file
write.csv(out,"Dup15q_TwoTailFishersExact_DMRgenes_GenomicContext.csv",quote=F, row.names = F)
write.csv(counts2,"Dup15q_DMRgenecounts_GenomicContex.csv",quote=F, row.names = F)





