
library(readxl)
library(qdapRegex)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(tools)
library(openxlsx)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- org.Hs.eg.db


# ASD annotation ----------------------------------------------------------

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/ASD")

ASD <- read_excel("DMRs_annotated.xlsx")
ASD <- ASD[,c(1:15)]
ASD$annotation <- rm_between_multiple(ASD$annotation, left="(uc", right=")")
ASD <- as.data.frame(ASD)
ASD$start <- as.numeric(ASD$start)
ASD$end <- as.numeric(ASD$end)
colnames(ASD)[1:3] <- c("chr",  "DMRstart",    "DMRend")
ASD$DMRstart <- ASD$DMRstart+1


#read in gene annotation from bedtools:
ASDgenes <-read.csv("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/gene_assignment/ASDDMRgenes_within10kb.csv")
ASDgenes <- ASDgenes[,c(2:13)]
ASDgenes$description <- rm_between_multiple(ASDgenes$description, left="[", right="]")

ASDgenes$DMRstart <- as.numeric(ASDgenes$DMRstart) #puts DMR back in 0 based coordiante system
ASDgenes$DMRend <- as.numeric(ASDgenes$DMRend)

#merge
ASD_anno <- merge(ASD,ASDgenes,by = c("chr","DMRstart","DMRend") , all.x=T)
ASD_anno <- ASD_anno %>% dplyr::select(-qval,-label)
write.xlsx(ASD_anno,"ASD_DMR_annotations.xlsx")

# Dup15q annotation ----------------------------------------------------------

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/Dup15q")

Dup15q <- read_excel("DMRs_annotated.xlsx")
Dup15q <- Dup15q[,c(1:15)]
Dup15q$annotation <- rm_between_multiple(Dup15q$annotation, left="(uc", right=")")
Dup15q <- as.data.frame(Dup15q)
Dup15q$start <- as.numeric(Dup15q$start)
Dup15q$end <- as.numeric(Dup15q$end)
colnames(Dup15q)[1:3] <- c("chr",  "DMRstart",    "DMRend")
Dup15q$DMRstart <- Dup15q$DMRstart+1


#read in gene annotation from bedtools:
Dup15qgenes <-read.csv("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/gene_assignment/Dup15qDMRgenes_within10kb.csv")
Dup15qgenes <- Dup15qgenes[,c(2:13)]
Dup15qgenes$description <- rm_between_multiple(Dup15qgenes$description, left="[", right="]")

Dup15qgenes$DMRstart <- as.numeric(Dup15qgenes$DMRstart) #puts DMR back in 0 based coordiante system
Dup15qgenes$DMRend <- as.numeric(Dup15qgenes$DMRend)

#merge
Dup15q_anno <- merge(Dup15q,Dup15qgenes,by = c("chr","DMRstart","DMRend") , all.x=T)
Dup15q_anno <- Dup15q_anno %>% dplyr::select(-qval,-label)
write.xlsx(Dup15q_anno,"Dup15q_DMR_annotations.xlsx")

# RTT annotation ----------------------------------------------------------

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/Rett")

RTT <- read_excel("DMRs_annotated.xlsx")
RTT <- RTT[,c(1:15)]
RTT$annotation <- rm_between_multiple(RTT$annotation, left="(uc", right=")")
RTT <- as.data.frame(RTT)
RTT$start <- as.numeric(RTT$start)
RTT$end <- as.numeric(RTT$end)
colnames(RTT)[1:3] <- c("chr",  "DMRstart",    "DMRend")
RTT$DMRstart <- RTT$DMRstart+1


#read in gene annotation from bedtools:
RTTgenes <-read.csv("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/gene_assignment/RettDMRgenes_within10kb.csv")
RTTgenes <- RTTgenes[,c(2:13)]
RTTgenes$description <- rm_between_multiple(RTTgenes$description, left="[", right="]")

RTTgenes$DMRstart <- as.numeric(RTTgenes$DMRstart) #puts DMR back in 0 based coordiante system
RTTgenes$DMRend <- as.numeric(RTTgenes$DMRend)

#merge
RTT_anno <- merge(RTT,RTTgenes,by = c("chr","DMRstart","DMRend"), all.x=T)
RTT_anno <- RTT_anno %>% dplyr::select(-qval,-label)
write.xlsx(RTT_anno,"RTT_DMR_annotations.xlsx")




# BG annotation -----------------------------------------------------------


#annotate consensus background
BG <- readPeakFile("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/ConsensusBackgroundRegions_sort.bed")
BGanno <- annotatePeak(BG, TxDb=txdb, annoDb = "org.Hs.eg.db", overlap = "all")
BGanno <- as.data.frame(BGanno)

BGanno$annotation <- rm_between_multiple(BGanno$annotation, left="(uc", right=")")
BGanno <- as.data.frame(BGanno)
BGanno$start <- as.numeric(BGanno$start)
BGanno$end <- as.numeric(BGanno$end)
colnames(BGanno)[1:3] <- c("chr",  "DMRstart",    "DMRend")
BGanno <- BGanno[,c(1:6)]
BGanno$DMRstart <- BGanno$DMRstart-1

#read in gene annotation from bedtools:
BGgenes <-read.csv("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/gene_assignment/BGregiongenes_within10kb.csv")
BGgenes<- BGgenes[,c(2:13)]
BGgenes$description <- rm_between_multiple(BGgenes$description, left="[", right="]")

BGgenes$DMRstart <- as.numeric(BGgenes$DMRstart) #puts DMR back in 0 based coordiante system
BGgenes$DMRend <- as.numeric(BGgenes$DMRend)

#merge
BG_anno <- merge(BGanno,BGgenes,by = c("chr","DMRstart","DMRend"), all.x=T)
BG_anno <- BG_anno %>% dplyr::select(-label)
write.xlsx(BG_anno,"BG_regions_annotations.xlsx")

# Summaries gene counts---------------------------------------------------------------
ASD_hypergenes <- ASD_anno %>% filter(percentDifference > 0 ) %>%
  dplyr::select(Ensbl) %>%
  distinct() %>%
  summarize(count = n())
ASD_hypergenes$DMR <- c("ASD Hyper DMRs")

ASD_hypogenes <- ASD_anno %>% filter(percentDifference < 0 ) %>%
  dplyr::select(Ensbl) %>%
  distinct() %>%
  summarize(count = n())
ASD_hypogenes$DMR <- c("ASD Hypo DMRs")

Dup15q_hypergenes <- Dup15q_anno %>% filter(percentDifference > 0 ) %>%
  dplyr::select(Ensbl) %>%
  distinct() %>%
  summarize(count = n())
Dup15q_hypergenes$DMR <- c("Dup15q Hyper DMRs")

Dup15q_hypogenes <- Dup15q_anno %>% filter(percentDifference < 0 ) %>%
  dplyr::select(Ensbl) %>%
  distinct() %>%
  summarize(count = n())
Dup15q_hypogenes$DMR <- c("Dup15q Hypo DMRs")

RTT_hypergenes <- RTT_anno %>% filter(percentDifference > 0 ) %>%
  dplyr::select(Ensbl) %>%
  distinct() %>%
  summarize(count = n())
RTT_hypergenes$DMR <- c("RTT Hyper DMRs")

RTT_hypogenes <- RTT_anno %>% filter(percentDifference < 0 ) %>%
  dplyr::select(Ensbl) %>%
  distinct() %>%
  summarize(count = n())
RTT_hypogenes$DMR <- c("RTT Hypo DMRs")

outgenes <- rbind(ASD_hypergenes,ASD_hypogenes,Dup15q_hypergenes,Dup15q_hypogenes,RTT_hypergenes,RTT_hypogenes)
colnames(outgenes)[1] <-c("Gene Counts")

write.xlsx(outgenes,"DMRassociatedGeneCounts.xlsx")

# Summaries DMRs---------------------------------------------------------------

#summarize ASD:
ASD_hyper <- ASD %>% filter(percentDifference > 0 ) %>%
  summarize(median = median(`percentDifference`), 
            mean= mean(`percentDifference`),
            min = min(`percentDifference`),
            max = max(`percentDifference`),
            count = n())

ASD_hyper <- data.frame(ASD_hyper)
ASD_hyper$DMR <- c("ASD Hyper DMRs")

ASD_hypo <-  ASD %>% filter(percentDifference < 0 ) %>%
  summarize(median = median(`percentDifference`), 
            mean= mean(`percentDifference`),
            min = min(`percentDifference`),
            max = max(`percentDifference`),
            count = n())

ASD_hypo <- data.frame(ASD_hypo)
ASD_hypo$DMR <- c("ASD Hypo DMRs")

ASD_all <- rbind(ASD_hyper,ASD_hypo)

ASDqoverall <- ASD %>% mutate(AbsRawDiff = abs(percentDifference)) %>%
  summarize(median = median(AbsRawDiff), 
            mean= mean(AbsRawDiff),
            min = min(AbsRawDiff),
            max = max(AbsRawDiff),
            count = n())
ASDqoverall$DMR <-c("ASD DMRs")


#Dup15q

Dup15q_hyper <- Dup15q %>% filter(percentDifference > 0 ) %>%
  summarize(median = median(`percentDifference`), 
            mean= mean(`percentDifference`),
            min = min(`percentDifference`),
            max = max(`percentDifference`),
            count = n())

Dup15q_hyper <- data.frame(Dup15q_hyper)
Dup15q_hyper$DMR <- c("Dup15q Hyper DMRs")

Dup15q_hypo <-  Dup15q %>% filter(percentDifference < 0 ) %>%
  summarize(median = median(`percentDifference`), 
            mean= mean(`percentDifference`),
            min = min(`percentDifference`),
            max = max(`percentDifference`),
            count = n())

Dup15q_hypo <- data.frame(Dup15q_hypo)
Dup15q_hypo$DMR <- c("Dup15q Hypo DMRs")

Dup15q_all <- rbind(Dup15q_hyper,Dup15q_hypo)

Dup15qqoverall <- Dup15q %>% mutate(AbsRawDiff = abs(percentDifference)) %>%
  summarize(median = median(AbsRawDiff), 
            mean= mean(AbsRawDiff),
            min = min(AbsRawDiff),
            max = max(AbsRawDiff),
            count = n())
Dup15qqoverall$DMR <-c("Dup15q DMRs")


#RTT
RTT_hyper <- RTT %>% filter(percentDifference > 0 ) %>%
  summarize(median = median(`percentDifference`), 
            mean= mean(`percentDifference`),
            min = min(`percentDifference`),
            max = max(`percentDifference`),
            count = n())

RTT_hyper <- data.frame(RTT_hyper)
RTT_hyper$DMR <- c("RTT Hyper DMRs")

RTT_hypo <-  RTT %>% filter(percentDifference < 0 ) %>%
  summarize(median = median(`percentDifference`), 
            mean= mean(`percentDifference`),
            min = min(`percentDifference`),
            max = max(`percentDifference`),
            count = n())

RTT_hypo <- data.frame(RTT_hypo)
RTT_hypo$DMR <- c("RTT Hypo DMRs")

RTT_all <- rbind(RTT_hyper,RTT_hypo)

RTTqoverall <- RTT %>% mutate(AbsRawDiff = abs(percentDifference)) %>%
  summarize(median = median(AbsRawDiff), 
            mean= mean(AbsRawDiff),
            min = min(AbsRawDiff),
            max = max(AbsRawDiff),
            count = n())
RTTqoverall$DMR <-c("RTT DMRs")

#outputs

out <- rbind(ASD_all,Dup15q_all,RTT_all,ASDqoverall,Dup15qqoverall,RTTqoverall)
out <- as.data.frame(out)

write.xlsx(out,"DMRcounts.xlsx")


# Fisher’s exact test for genomic locations -------------------------------

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/GenomicElementErichment")
library(tools)

ASD_anno2 <- ASD_anno %>% dplyr::select(chr,DMRstart,DMRend,annotation) %>% distinct() %>% mutate(DMR = c("ASD DMRs"))
Dup15q_anno2 <- Dup15q_anno %>% dplyr::select(chr,DMRstart,DMRend,annotation) %>% distinct() %>% mutate(DMR = c("Dup15q DMRs"))
RTT_anno2 <- RTT_anno %>% dplyr::select(chr,DMRstart,DMRend,annotation) %>% distinct() %>% mutate(DMR = c("RTT DMRs"))
BG_anno2 <- BG_anno %>% dplyr::select(chr,DMRstart,DMRend,annotation) %>% distinct() %>% mutate(DMR = c("BG Regions"))

AllData<- rbind(ASD_anno2,Dup15q_anno2,RTT_anno2,BG_anno2)


#make df
count_df_allData <- AllData %>% 
  dplyr::group_by(DMR,annotation) %>%
  dplyr::summarise(count = n()) %>%
  mutate(total = sum(count)) %>%
  mutate(percent = count/total*100)

TotalCounts <-as.data.frame(count_df_allData)

#make list by DMR type
TotalCounts$DMR <- factor(TotalCounts$DMR,levels=c("ASD DMRs","Dup15q DMRs","RTT DMRs","BG Regions"))
sumcountsList <- split(TotalCounts, TotalCounts$DMR)
names(sumcountsList)

#write to one file:
library(openxlsx)
write.xlsx(TotalCounts, file = "Annotation_Allfiles_Counts.xlsx")


#two sided fisher's exact test
call <- c("ASD DMRs","Dup15q DMRs" ,"RTT DMRs")
masterfile <- NULL
for (i in call) {
  print(i)
  
  #get background data
  BG <- sumcountsList$`BG Regions`
  BG <- as.data.frame(BG)
  
  #all annotations possible
  annotations <- unique(BG$annotation)
  
  #get dataset
  data <- sumcountsList[[i]]
  data <- as.data.frame(data)
  
  data2 <- merge(BG,data,by ="annotation",all.x=T)
  
  #if missing annotation replace with zero
  data2[is.na(data2)] <- 0
  
  #for each annotation:
  out <- NULL
  for(an in annotations) {
    
    tmp <- filter(data2,annotation == an)
    #for DEpeak list:
    roiplus <- tmp$count.y
    roineg <- tmp$total.y - tmp$count.y
    #for BG peak list:
    BGplus <- tmp$count.x
    BGneg <- tmp$total.x - tmp$count.x
    mytable <- rbind(c(roiplus,BGplus),c(roineg,BGneg))
    test <- fisher.test(mytable,alternative='two.sided')
    
    pvalue <- as.numeric(as.character(test$p.value))
    HigherCI <- as.numeric(as.character(test$conf.int[1]))
    LowerCI <- as.numeric(as.character(test$conf.int[2]))
    OddsRatio <- as.numeric(as.character(test$estimate))
    
    results <- data.frame(annotation=paste(an),pvalue,HigherCI,LowerCI,OddsRatio)
    out <- rbind(out,results) 
  }
  
  #make masteroutfile    
  out <- as.data.frame(out)
  out$comparison <- paste(names(sumcountsList[i]))
  
  masterfile <- rbind(out,masterfile)  
  
}

#masterfile <- filter(masterfile, comparison != "ConsensusBackgroundRegions")

masterfile$FDR <- p.adjust(masterfile$pvalue, method="fdr") #FDR

openxlsx::write.xlsx(masterfile, file = "TwoTailFishersTest_Annotations_ChIPseeker.xlsx")


# #collapse count data down to df
# sumcountsList$ConsensusBackgroundRegions$list <- c("ConsensusBackgroundRegions")
# sumcountsList$`ASD DMRs`$list <- c("ASD DMRs")
# sumcountsList$`Dup15q DMRs`$list <- c("Dup15q DMRs")
# sumcountsList$`Rett DMRs`$list <- c("Rett DMRs")
# countDF <- rbind(sumcountsList$ASDDMRs,sumcountsList$Dup15qDMRs,sumcountsList$RettDMRs,
#                  sumcountsList$ConsensusBackgroundRegions)
# 
# openxlsx::write.xlsx(countDF, file = "Counts_Annotations_ChIPseeker.xlsx")
# 
#graphs
library(ggplot2)
library(cowplot)
countDF <- TotalCounts
countDF$annotation <- as.factor(countDF$annotation)
countDF$DMR <- factor(countDF$DMR, levels= c("BG Regions","RTT DMRs","Dup15q DMRs","ASD DMRs"))
countDF$percent <- as.numeric(as.character(countDF$percent))
# The palette:
#cbPalette <- c("#CC79A7","#E69F00", "#D55E00" ,"#56B4E9", "#009E73")
library(RColorBrewer)
library(randomcoloR)
marker = color = brewer.pal(11, "Paired")
colors()[1:100] 
palette()

marker <- c("darkmagenta","cornflowerblue","pink","darkblue","blanchedalmond","coral","cadetblue","orange","brown3","darkgreen","bisque3")

plot <-ggplot(countDF, aes(fill=annotation, y=percent, x=DMR)) +
  geom_bar( stat="identity") +
  coord_flip()+
  #ggnice() + 
  labs(title="Genomic Context of DMRs")+
  ylab("Percentage") +
  xlab(" ") +
  theme(text = element_text(size=14))+
  scale_fill_manual(values=marker) +
  theme(legend.title = element_text(size =15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=15),
        axis.text.y = element_text(size=15))

plot

pdf("ASD_DMRcontextbar.pdf", height = 4, width = 8, useDingbats=FALSE)    # create PNG for the heat map       
print(plot)
dev.off()



