#9/5/18
##############################################################################
# MulitQC output from Bismark for alignment data etc
# read in and run ANOVA
#Annie Vogel Ciernia
##############################################################################

#load in and concatinate pm_stats output from wgbs tools
library(ggplot2)
library(lsmeans)
library(nlme)
library(dplyr)
library(tidyr)
library(xlsx)
options("scipen"=100, "digits"=4) #prevent exponents
##############################################################################
#Dup15q samples
path <- "/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/WGBS_QC/Dup15q"
setwd(path)

##############################################################################
#Total % methylation from Bismark
##############################################################################

PerMe <- read.delim("multiqc_bismark_methextract.txt",header=T)

#calcualte percent ch methylation
#combined reads for methylation chh and chg
PerMe$meth_ch_reads <- PerMe$meth_chg+PerMe$meth_chh
PerMe$unmeth_ch_reads <- PerMe$unmeth_chg + PerMe$unmeth_chh
#mCH/CH
PerMe$percent_chh_meth <- (PerMe$meth_ch_reads/(PerMe$unmeth_ch_reads+PerMe$meth_ch_reads))*100

#mCG/CG
PerMe$percentmeth_cpg <- (PerMe$meth_cpg/(PerMe$meth_cpg+PerMe$unmeth_cpg))*100

#add in group info:
#read in subject data
sampleinfo <- read.csv("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/AllSamples.csv",header=T)
sampleinfo <- sampleinfo[which(sampleinfo$filename %in% PerMe$Sample),]

PerMemerge <- merge(PerMe,sampleinfo,by.x="Sample",by.y="filename")


PerMemerge$Group <- factor(PerMemerge$Group)
PerMemerge$Sex <- factor(PerMemerge$Sex)
PerMemerge$Age <- as.numeric(PerMemerge$Age..yrs.)

DFout1 <- PerMemerge[,c(1,6,4,9,13,10,11,12)]

##############################################################################
#Read alignment info from Bismark
##############################################################################

DedupReads <- read.table("multiqc_bismark_dedup.txt",header=T)

Coverage <- read.delim("multiqc_general_stats.txt",header=T)
DF <- merge(DedupReads,Coverage,by="Sample")


DFMemerge <- merge(DF,sampleinfo,by.x="Sample",by.y="filename")

DFMemerge$Group <- factor(DFMemerge$Group)
DFMemerge$Age <- as.numeric(DFMemerge$Age..yrs.)

DFtest <- DFMemerge %>% dplyr::select(Sample,Group,Sex,Age,aligned_reads,Bismark_mqc.generalstats.percent_aligned,
                                      dedup_reads,dup_reads_percent,Bismark_mqc.generalstats.total_c,
                                      Bismark_mqc.generalstats.C_coverage)

colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.percent_aligned"] <- c("percent_aligned")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.total_c"] <- c("TotalCs")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.C_coverage"] <- c("Coverage")

colnames(DFtest)

DFmaster <- merge(PerMemerge,DFtest,by="Sample")

##############################################################################
#ANOVA
##############################################################################

#model with Age as a covariate
tmp <- lme(percentmeth_cpg ~ Group.x+Age.x+Coverage, ~1|Sample, data = DFmaster)
#model indcludes random effects subject ID 
anova <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova$factors <- rownames(anova)
anova <- anova[!grepl('(Intercept)', anova$factors),]
anova$measure <- c("Percent mCG/CG")
anova


#model
tmp2 <- lme(percent_chh_meth ~ Group.x+Age.x+Coverage, ~1|Sample, data = DFmaster)
#model indcludes random effects subject ID 
anova2 <- anova(tmp2, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova2$factors <- rownames(anova2)
anova2 <- anova2[!grepl('(Intercept)', anova2$factors),]
anova2$measure <- c("Percent mCH/CH")
anova2


ANOVAout <- rbind(anova,anova2)

write.xlsx(ANOVAout, file = "Dup15q_RM-ANOVA_PercentMethylation.xlsx")

###############################################################################################
#graphs
###############################################################################################
library(ggplot2)
library(cowplot)

PerMemerge2 <- PerMemerge%>%
  group_by(Sample,Group) %>% gather(Type,PercentMeth,c(13,12))

PerMemerge2$Type <- gsub("percent_chh_meth","Percent mCH/CH",PerMemerge2$Type)
PerMemerge2$Type <- gsub("percentmeth_cpg","Percent mCG/CG",PerMemerge2$Type)
PerMemerge2$Group <- gsub("Control_Dup15q","Control",PerMemerge2$Group)
PerMemerge2$Group <- as.factor(PerMemerge2$Group)
PerMemerge2$Group <- relevel(PerMemerge2$Group,ref="Control")

pdf("Dup15q_GLobaMethylation.pdf", height = 3, width =6)    # create PNG for the heat map       

cbPalette <- c("#0072B2","#D55E00")
ggplot(PerMemerge2, aes(x=Group, y=PercentMeth)) + 
  facet_wrap(~Type, scales = "free")+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Group))+
  geom_point(position = position_dodge(width = 0.90),aes(group=Group)) + 
  scale_fill_manual(values = cbPalette) +
  theme_cowplot(font_size = 15)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  scale_x_discrete(name = "Diagnosis") + #,limits = order2)+
  scale_y_continuous(name = "Global Percent Methylation")+
  ggtitle("Global Percent Methylation") 

dev.off()

##############################################################################
#Read alignment info from Bismark
##############################################################################

DedupReads <- read.table("multiqc_bismark_dedup.txt",header=T)

Coverage <- read.delim("multiqc_general_stats.txt",header=T)
DF <- merge(DedupReads,Coverage,by="Sample")


DFMemerge <- merge(DF,sampleinfo,by.x="Sample",by.y="filename")

DFMemerge$Group <- factor(DFMemerge$Group)
DFMemerge$Age <- as.numeric(DFMemerge$Age..yrs.)

DFtest <- DFMemerge %>% dplyr::select(Sample,Group,Sex,Age,aligned_reads,Bismark_mqc.generalstats.percent_aligned,
                                      dedup_reads,dup_reads_percent,Bismark_mqc.generalstats.total_c,
                                      Bismark_mqc.generalstats.C_coverage)

colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.percent_aligned"] <- c("percent_aligned")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.total_c"] <- c("TotalCs")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.C_coverage"] <- c("Coverage")

colnames(DFtest)

#model t-test
tmp <- t.test(aligned_reads ~ Group, data = DFtest)
test1 <- cbind(statistic=tmp$statistic,
               df=tmp$parameter,
               pvalue=tmp$p.value,
               method=tmp$method,
               test=tmp$alternative,
               measure=c("Bismark Aligned Reads"))

#test2
tmp <- t.test(percent_aligned ~ Group, data = DFtest)
test2 <- cbind(statistic=tmp$statistic,
               df=tmp$parameter,
               pvalue=tmp$p.value,
               method=tmp$method,
               test=tmp$alternative,
               measure=c("Percent Aligned Reads"))


#model Group * sex
tmp <- t.test(dedup_reads ~ Group, data = DFtest)
test3 <- cbind(statistic=tmp$statistic,
               df=tmp$parameter,
               pvalue=tmp$p.value,
               method=tmp$method,
               test=tmp$alternative,
               measure=c("DeDuplicated Reads"))


#model Group * sex
tmp <- t.test(dup_reads_percent ~ Group, data = DFtest)
test4 <- cbind(statistic=tmp$statistic,
               df=tmp$parameter,
               pvalue=tmp$p.value,
               method=tmp$method,
               test=tmp$alternative,
               measure=c("Percent Duplicated Reads"))


#model Group * sex
tmp <- t.test(Coverage ~ Group, data = DFtest)
test5 <- cbind(statistic=tmp$statistic,
               df=tmp$parameter,
               pvalue=tmp$p.value,
               method=tmp$method,
               test=tmp$alternative,
               measure=c("Coverage"))


#model Group * sex
tmp <- t.test(TotalCs ~ Group, data = DFtest)
test6 <- cbind(statistic=tmp$statistic,
               df=tmp$parameter,
               pvalue=tmp$p.value,
               method=tmp$method,
               test=tmp$alternative,
               measure=c("Total C's Analyzed"))

TTESTout <- rbind(test1,test2,test3,test4,test5,test6)
TTESTout <- as.data.frame(TTESTout)
TTESTout$cohort <- c("Dup15q")

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/WGBS_QC/Dup15q")
write.xlsx(TTESTout,file="Dup15q_WGBSQC_T-test.xlsx",row.names = F)

mergeDF <- merge(sampleinfo,DFtest,by.y="Sample",by.x="filename",all.y=T)
mergeDF2 <- merge(mergeDF,DFout1,by.x="filename",by.y="Sample")

write.xlsx(mergeDF2,file="Dup15q_WGBS_sampleinfo.xlsx")

##############################################################################
#PMI and AGE
##############################################################################

DFmaster$Age <- as.numeric(DFmaster$Age..yrs.)
DFmaster$PMI <- as.numeric(DFmaster$PMI)

tmp <- t.test(PMI ~ Group.x, data = DFmaster)
PMI <- cbind(statistic=tmp$statistic,
               df=tmp$parameter,
               pvalue=tmp$p.value,
               method=tmp$method,
               test=tmp$alternative,
               measure=c("PMI"))

tmp <- t.test(Age ~ Group.x, data = DFmaster)
Age <- cbind(statistic=tmp$statistic,
             df=tmp$parameter,
             pvalue=tmp$p.value,
             method=tmp$method,
             test=tmp$alternative,
             measure=c("Age"))


TTESTout <- rbind(PMI,Age)
TTESTout <- as.data.frame(TTESTout)
TTESTout$cohort <- c("Dup15q")

write.xlsx(TTESTout,"Ttest_PMI_Age.xlsx")

#correlatoin global me and pmi
corrMEPMI <- cor.test(DFmaster$PMI,DFmaster$percent_cpg_meth, method="pearson")
coroutput <- data.frame(comparison = c("Pearson Correlation mCG/CG and PMI"),
                        r = corrMEPMI$estimate,
                        pvalue= corrMEPMI$p.value)

write.xlsx(coroutput,"PearsonCorl_PMI_mCG.xlsx")

