#Load in gene lists, gene lengths and analyze DMR associated genes vs BG genes 
#author: Annie Vogel Ciernia
#a.ciernia@gmail.com
#3/2019

##############################################################################################################
library(dplyr)
library(tidyr)
library(cowplot)
##############################################################################################################

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/12_6_2018_4CpGDMRs/DMR_analysis.R/gene_assignment")
load("Hg38_Conversioninfo.RData")
#ASD DMR genes 
ASD_genes <- read.csv("ASDDMRgenes_within10kb.csv",header=T)
Dup15q_genes <- read.csv("Dup15qDMRgenes_within10kb.csv",header=T)
RTT_genes <- read.csv("RettDMRgenes_within10kb.csv",header=T)
BG_genes <- read.csv("Consensus_BGregiongenes_within10kb.csv",header=T)


ASD_genes <- ASD_genes %>% mutate(length=gene_end - gene_start) %>%
  dplyr::select(Ensbl,length) %>%
  distinct()

Dup15q_genes <- Dup15q_genes %>% mutate(length=gene_end - gene_start) %>%
  dplyr::select(Ensbl,length) %>%
  distinct()

RTT_genes <- RTT_genes %>% mutate(length=gene_end - gene_start) %>%
  dplyr::select(Ensbl,length) %>%
  distinct()

BG_genes <- BG_genes %>% mutate(length=gene_end - gene_start) %>%
  dplyr::select(Ensbl,length) %>%
  distinct()

ASD_genes$group <- c("ASD DMR genes")
Dup15q_genes$group <- c("Dup15q DMR genes")
RTT_genes$group <- c("RTT DMR genes")
BG_genes$group <- c("BG DMR genes")

master <- rbind(ASD_genes,Dup15q_genes,RTT_genes,BG_genes)

master <- master[!is.na(master$length),]

master$group <- factor(master$group,levels=unique(master$group))

#histogram of length
# Basic histogram
library(ggplot2)
plot <- ggplot(master, aes(x=length,color=group)) + 
 # Change the width of bins
  geom_density() + # Change colors
  #geom_vline(aes(xintercept=0), color="blue", linetype="dashed") +
  #scale_color_manual(values=c("#999999", "#E69F00"))+
 # scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_x_log10() +
  labs(title=" ",x="Log10 Gene Length (bp)", y = "Density")
plot

ggsave(plot = plot, filename="GeneLengthDensityPlot_Genes10kb.pdf",width = 8, height = 5, dpi = 300, useDingbats = FALSE)

###################################################################
load("Hg38_Conversioninfo.RData")
#all human genes in Hg38
HumanHg38 <- ids

#length
HumanHg38$length <- abs(HumanHg38$end_position - HumanHg38$start_position)

#remove rows with no length 
HumanHg38 <- HumanHg38[!is.na(HumanHg38$length),]

HumanHg38 <- distinct(HumanHg38[,c("ensembl_gene_id","length")])

names(HumanHg38) <- c("Ensbl","length")

HumanHg38$group <- c("All Hg38 genes")

master2 <- rbind(master,HumanHg38)
master2$group <- factor(master2$group,levels=unique(master2$group))

colorpal <- c("firebrick1","goldenrod2","blue2","gray48","black")

plot <- ggplot(master2, aes(x=length,color=group)) + 
  # Change the width of bins
  geom_density() + # Change colors
  #geom_vline(aes(xintercept=0), color="blue", linetype="dashed") +
  scale_color_manual(values=colorpal)+
  #scale_fill_manual(values=colorpal)+
  scale_x_log10() +
  labs(title=" ",x="Log10 Gene Length (bp)", y = "Density")
plot

ggsave(plot = plot, filename="GeneLengthDensityPlot_Genes10kb_Hg38.pdf",width = 8, height = 5, dpi = 300, useDingbats = FALSE)

###################################################################
#boxplot
###################################################################
colorpal <- c("firebrick1","goldenrod2","blue2","burlywood4","gray48")
master2$group <- gsub("BG DMR genes", "BG Region genes", master2$group)

master2$group <- factor(master2$group, levels= unique(master2$group))

pdf("Boxplot_DMRassocGenes_log10genelength.pdf", height = 4, width =6)    # create PNG for the heat map       

ggplot(master2, aes(x=group, y=log10(length)),group=group) + 
  #facet_wrap(~DMR_features, scales = "free")+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=group))+ #gives line as median, top of box and bottom of box as 25th and 75th quartile, line = 5th and 9th percentile
  #geom_point(position = position_dodge(width = 0.90),aes(group=group)) + 
  scale_fill_manual(values = colorpal) +
  theme_cowplot(font_size = 15)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "DMR associated genes") + #,limits = order2)+
  scale_y_continuous(name = "Log10 Gene Length (bp)")

dev.off();

###################################################################
#ANOVA: Are the distributions different from each other?
###################################################################
library(lsmeans)
library(nlme)


#means:
meanlength <- master2 %>% group_by(group) %>% summarize(meanlength = mean(length))

#model length by group
tmp <- kruskal.test(length ~ group, data = master2)

kt_DF <- data.frame(method= tmp$method,
                     pval = tmp$p.value,
                     statistic = tmp$statistic,
                    test = tmp$data.name)

#posthoc tests for group at each timepoint
cox <- pairwise.wilcox.test(master2$length, master2$group,
                     p.adjust.method = "BH")

cox_DF <- data.frame(method= cox$method,
                     pval = cox$p.value,
                     p.adjust.method = cox$p.adjust.method
                     )

library(openxlsx)
write.xlsx(kt_DF, file="kruskal.test_genelength_DMRassocgenes.xlsx")
write.xlsx(cox_DF, file="BHpairwise.wilcox.test_genelength_DMRassocgenes.xlsx")
write.xlsx(meanlength, file="Mean_genelength_DMRassocgenes.xlsx")

###################################################################
#Permutation: Are they different from the whole genome?
###################################################################



#function to get mean length distribution randomly from all hg38 genes
getSimDistribution <- function(gA,geneListA){
  gA.sim <- sample(nrow(geneListA),size=length(gA$Ensbl))
  gA.sim.DF <- geneListA[gA.sim,]
  return(median(gA.sim.DF$length))
}

gA <- ASD_genes # list of genes of interest from organism A
geneListA <- HumanHg38 #background genelist
N <- 100000 # simulation iteration


#perform purmutations
inter.sim <- replicate(N,getSimDistribution(gA,geneListA))

#permuted pvalue
#p value is not = 0 .
#If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
inter.r<-median(gA$length)
p_ASD <- sum(inter.sim>=inter.r)/N

pdf(file="Histogram_medianGeneLengthASD_DMRs.pdf",height = 3,width = 5,useDingbats = F)
hist(inter.sim, 
     main="Histogram for ASD DMR Gene Length",
     xlab="Median Gene Length Matched Resampling from all Hg38",
     xlim = c(1000,60000))
abline(v=inter.r, lwd=2, col="purple")
dev.off()

#Dup15q

gA <- Dup15q_genes # list of genes of interest from organism A
geneListA <- HumanHg38 #background genelist
N <- 100000 # simulation iteration


#perform purmutations
inter.sim <- replicate(N,getSimDistribution(gA,geneListA))

#permuted pvalue
#p value is not = 0 .
#If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
inter.r<-median(gA$length)
p_Dup15q <- sum(inter.sim>=inter.r)/N

pdf(file="Histogram_medianGeneLengthDup15q_DMRs.pdf",height = 5,width = 5,useDingbats = F)
hist(inter.sim, 
     main="Histogram for Dup15q DMR Gene Length",
     xlab="Median Gene Length Matched Resampling from all Hg38",
     xlim = c(1000,60000))
abline(v=inter.r, lwd=2, col="purple")
dev.off()

#RTT

gA <- RTT_genes # list of genes of interest from organism A
geneListA <- HumanHg38 #background genelist
N <- 100000 # simulation iteration


#perform purmutations
inter.sim <- replicate(N,getSimDistribution(gA,geneListA))

#permuted pvalue
#p value is not = 0 .
#If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
inter.r<-median(gA$length)
p_RTT <- sum(inter.sim>=inter.r)/N

pdf(file="Histogram_medianGeneLengthRTT_DMRs.pdf",height = 5,width = 5,useDingbats = F)
hist(inter.sim, 
     main="Histogram for RTT DMR Gene Length",
     xlab="Median Gene Length Matched Resampling from all Hg38",
     xlim = c(1000,60000))
abline(v=inter.r, lwd=2, col="purple")
dev.off()

#BG DMRs

gA <- BG_genes # list of genes of interest from organism A
geneListA <- HumanHg38 #background genelist
N <- 100000 # simulation iteration


#perform purmutations
inter.sim <- replicate(N,getSimDistribution(gA,geneListA))

#permuted pvalue
#p value is not = 0 .
#If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
inter.r<-median(gA$length)
p_BG <- sum(inter.sim>=inter.r)/N

pdf(file="Histogram_medianGeneLengthBG_DMRs.pdf",height = 5,width = 5,useDingbats = F)
hist(inter.sim, 
     main="Histogram for BG DMR Gene Length",
     xlab="Median Gene Length Matched Resampling from all Hg38",
     xlim = c(1000,20000))
abline(v=inter.r, lwd=2, col="purple")
dev.off()

#all genes vs Hg38
Allp <- rbind(p_ASD,p_Dup15q,p_RTT,p_BG)
Allp <- as.data.frame(Allp)
names(Allp) <- c("permutation pvalues")
Allp$DMRgenes <- rownames(Allp)
Allp$DMRgenes <- gsub("p_","",Allp$DMRgenes)
#replace p=0 with 1/N
Allp$`permutation pvalues` <- gsub("0",1/N,Allp$`permutation pvalues`)

Allp$backgroundlist <- c("All Hg38 genes")

Allp$`median gene length(bp` <- c(median(ASD_genes$length),
                                  median(Dup15q_genes$length),
                                  median(RTT_genes$length),
                                  median(BG_genes$length))

write.xlsx(Allp,file="PermutationPvalues_vsHg38genes.xlsx")

###################################################################
#Permutation: Are they different from the the BG list?
###################################################################

#function to get mean length distribution randomly from all hg38 genes
getSimDistribution <- function(gA,geneListA){
  gA.sim <- sample(nrow(geneListA),size=length(gA$Ensbl))
  gA.sim.DF <- geneListA[gA.sim,]
  return(median(gA.sim.DF$length))
}

gA <- ASD_genes # list of genes of interest from organism A
geneListA <- BG_genes #background genelist
N <- 100000 # simulation iteration


#perform purmutations
inter.sim <- replicate(N,getSimDistribution(gA,geneListA))

#permuted pvalue
#p value is not = 0 .
#If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
inter.r<-median(gA$length)
p_ASD <- sum(inter.sim>=inter.r)/N

pdf(file="Histogram_medianGeneLengthASD_DMRs_vsBGDMRs.pdf",height = 3,width = 5,useDingbats = F)
hist(inter.sim, 
     main="Histogram for ASD DMR Gene Length",
     xlab="Median Gene Length Matched Resampling from BG regions",
     xlim = c(1000,median(gA$length)+10000))
abline(v=inter.r, lwd=2, col="purple")
dev.off()

#Dup15q

gA <- Dup15q_genes # list of genes of interest from organism A
geneListA <- BG_genes #background genelist
N <- 100000 # simulation iteration


#perform purmutations
inter.sim <- replicate(N,getSimDistribution(gA,geneListA))

#permuted pvalue
#p value is not = 0 .
#If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
inter.r<-median(gA$length)
p_Dup15q <- sum(inter.sim>=inter.r)/N

pdf(file="Histogram_medianGeneLengthDup15q_vsBGDMRs.pdf",height = 5,width = 5,useDingbats = F)
hist(inter.sim, 
     main="Histogram for Dup15q DMR Gene Length",
     xlab="Median Gene Length Matched Resampling from BG regions",
     xlim = c(1000,median(gA$length)+10000))
abline(v=inter.r, lwd=2, col="purple")
dev.off()

#RTT

gA <- RTT_genes # list of genes of interest from organism A
geneListA <- BG_genes #background genelist
N <- 100000 # simulation iteration


#perform purmutations
inter.sim <- replicate(N,getSimDistribution(gA,geneListA))

#permuted pvalue
#p value is not = 0 .
#If you do N=100,000 iterations the only thing you can say is that p < 1e-06 if you do not see any simulated intersection bigger than the real one. Nevertheless you can report the mean simulated intersection also.
inter.r<-median(gA$length)
p_RTT <- sum(inter.sim>=inter.r)/N

pdf(file="Histogram_medianGeneLengthRTT_DMRsvsBGDMRs.pdf",height = 5,width = 5,useDingbats = F)
hist(inter.sim, 
     main="Histogram for RTT DMR Gene Length",
     xlab="Median Gene Length Matched Resampling from BG regions",
     xlim = c(1000,median(gA$length)+10000))
abline(v=inter.r, lwd=2, col="purple")
dev.off()


#all genes vs Hg38
Allp <- rbind(p_ASD,p_Dup15q,p_RTT)
Allp <- as.data.frame(Allp)
names(Allp) <- c("permutation pvalues")
Allp$DMRgenes <- rownames(Allp)
Allp$DMRgenes <- gsub("p_","",Allp$DMRgenes)
#replace p=0 with 1/N
Allp$`permutation pvalues` <- gsub("0",1/N,Allp$`permutation pvalues`)

Allp$backgroundlist <- c("BG region genes")

Allp$`median gene length(bp` <- c(median(ASD_genes$length),
                                  median(Dup15q_genes$length),
                                  median(RTT_genes$length))

write.xlsx(Allp,file="PermutationPvalues_vsBGregions.xlsx")

