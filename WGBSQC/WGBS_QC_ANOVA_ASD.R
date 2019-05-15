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
#ASD samples
path <- "/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/WGBS_QC/ASD"
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
PerMe$Sample <- gsub("_filtered","",PerMe$Sample)
PerMemerge <- merge(PerMe,sampleinfo,by.x="Sample",by.y="filename",all.x=T)


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


DF$Sample <- gsub("_filtered","",DF$Sample)
DFMemerge <- merge(DF,sampleinfo,by.x="Sample",by.y="filename",all.x=T)

DFMemerge$Group <- factor(DFMemerge$Group)
DFMemerge$Sex <- factor(DFMemerge$Sex)
DFMemerge$Age <- as.numeric(DFMemerge$Age..yrs.)

DFtest <- DFMemerge %>% dplyr::select(Sample,Group,Sex,Age,aligned_reads,Bismark_mqc.generalstats.bismark.percent_aligned,
                                      dedup_reads,dup_reads_percent,Bismark_mqc.generalstats.bismark.total_c,
                                      Bismark_mqc.generalstats.bismark.C_coverage)

colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.bismark.percent_aligned"] <- c("percent_aligned")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.bismark.C_coverage"] <- c("Coverage")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.bismark.total_c"] <- c("TotalCs")

colnames(DFtest)

DFmaster <- merge(PerMemerge,DFtest,by="Sample")
##############################################################################
#ANOVA
##############################################################################

#model with sex, Age and coverage as a covariates
tmp <- lme(percentmeth_cpg ~ Group.x+Sex.x+Age.x+Coverage, ~1|Sample, data = DFmaster)
#model indcludes random effects subject ID 
anova <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova$factors <- rownames(anova)
anova <- anova[!grepl('(Intercept)', anova$factors),]
anova$measure <- c("Percent mCG/CG")
anova

# #posthoc tests for group at each timepoint
# refgrid <- ref.grid(tmp)
# tmp2 <- lsmeans(refgrid, ~Group|Sex)
# tmp3 <- summary(pairs(tmp2, adjust = "none"))
# tmp3 <- as.data.frame(tmp3)
# #adjust p value for all comparisons:
# tmp3$p.adjust <- p.adjust(tmp3$p.value,method = "BH") #BH adjusted p values
# tmp3$measure <- c("Percent mCG/CG")
# tmp3
# 
#model
tmp2 <- lme(percent_chh_meth ~ Group.x+Sex.x+Age.x+Coverage, ~1|Sample, data = DFmaster)
#model indcludes random effects subject ID 
anova2 <- anova(tmp2, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova2$factors <- rownames(anova2)
anova2 <- anova2[!grepl('(Intercept)', anova2$factors),]
anova2$measure <- c("Percent mCH/CH")
anova2

# #posthoc tests for group at each timepoint
# refgrid <- ref.grid(tmp2)
# PH2 <- lsmeans(refgrid, ~Group|Sex)
# PH2 <- summary(pairs(PH2, adjust = "none"))
# PH2 <- as.data.frame(PH2)
# #adjust p value for all comparisons:
# PH2$p.adjust <- p.adjust(PH2$p.value,method = "BH") #BH adjusted p values
# PH2$measure <- c("Percent mCH/CH")
# PH2

ANOVAout <- rbind(anova,anova2)
#BHout <- rbind(tmp3,PH2)

#write.xlsx(BHout, file = "ASD_GlobalPercentMethylationPostHocs.xlsx")

write.xlsx(ANOVAout, file = "ASD_RM-ANOVA_PercentMethylation.xlsx")

###############################################################################################
#graphs
###############################################################################################
library(ggplot2)
library(cowplot)

PerMemerge2 <- PerMemerge%>%
  group_by(Sample,Group,Sex) %>% gather(Type,PercentMeth,c(12,11))

PerMemerge2$Type <- gsub("percent_chh_meth","Percent mCH/CH",PerMemerge2$Type)
PerMemerge2$Type <- gsub("percentmeth_cpg","Percent mCG/CG",PerMemerge2$Type)
PerMemerge2$Group <- gsub("Control_ASD","Control",PerMemerge2$Group)
PerMemerge2$Group <- as.factor(PerMemerge2$Group)
PerMemerge2$Group <- relevel(PerMemerge2$Group,ref="Control")

pdf("ASD_GLobaMethylation.pdf", height = 3, width =6)    # create PNG for the heat map       

cbPalette <- c("#0072B2","#D55E00")
ggplot(PerMemerge2, aes(x=Group, y=PercentMeth),group=Sex) + 
  facet_wrap(~Type, scales = "free")+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Group))+
  geom_point(position = position_dodge(width = 0.90),aes(group=Group,color=Sex)) + 
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


DF$Sample <- gsub("_filtered","",DF$Sample)
DFMemerge <- merge(DF,sampleinfo,by.x="Sample",by.y="filename",all.x=T)

DFMemerge$Group <- factor(DFMemerge$Group)
DFMemerge$Sex <- factor(DFMemerge$Sex)
DFMemerge$Age <- as.numeric(DFMemerge$Age..yrs.)

DFtest <- DFMemerge %>% dplyr::select(Sample,Group,Sex,Age,aligned_reads,Bismark_mqc.generalstats.bismark.percent_aligned,
                                      dedup_reads,dup_reads_percent,Bismark_mqc.generalstats.bismark.total_c,
                                      Bismark_mqc.generalstats.bismark.C_coverage)

colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.bismark.percent_aligned"] <- c("percent_aligned")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.bismark.C_coverage"] <- c("Coverage")
colnames(DFtest)[colnames(DFtest)=="Bismark_mqc.generalstats.bismark.total_c"] <- c("TotalCs")

colnames(DFtest)

#model Group + sex
tmp <- lme(aligned_reads ~ Group+Sex, ~1|Sample, data = DFtest)
#model indcludes random effects subject ID 
anova <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova$factors <- rownames(anova)
anova <- anova[!grepl('(Intercept)', anova$factors),]
anova$measure <- c("Bismark Aligned Reads")
anova

#model Group * sex
tmp <- lme(percent_aligned ~ Group+Sex+, ~1|Sample, data = DFtest)
#model indcludes random effects subject ID 
anova1 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova1$factors <- rownames(anova1)
anova1 <- anova1[!grepl('(Intercept)', anova1$factors),]
anova1$measure <- c("Percent Aligned Reads")
anova1

#model Group + sex
tmp <- lme(dedup_reads ~ Group+Sex, ~1|Sample, data = DFtest)
#model indcludes random effects subject ID 
anova2 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova2$factors <- rownames(anova2)
anova2 <- anova2[!grepl('(Intercept)', anova2$factors),]
anova2$measure <- c("DeDuplicated Reads")
anova2

#model Group + sex
tmp <- lme(dup_reads_percent ~ Group+Sex, ~1|Sample, data = DFtest)
#model indcludes random effects subject ID 
anova3 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova3$factors <- rownames(anova3)
anova3 <- anova3[!grepl('(Intercept)', anova3$factors),]
anova3$measure <- c("Percent Duplicated Reads")
anova3

#model Group + sex
tmp <- lme(Coverage ~ Group+Sex, ~1|Sample, data = DFtest)
#model indcludes random effects subject ID 
anova4 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova4$factors <- rownames(anova4)
anova4 <- anova4[!grepl('(Intercept)', anova4$factors),]
anova4$measure <- c("Coverage")
anova4

#model Group + sex
tmp <- lme(TotalCs ~ Group+Sex, ~1|Sample, data = DFtest)
#model indcludes random effects subject ID 
anova5 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova5$factors <- rownames(anova5)
anova5 <- anova5[!grepl('(Intercept)', anova5$factors),]
anova5$measure <- c("Total C's Analyzed")
anova5

ANOVAout <- rbind(anova,anova1,anova2,anova3,anova4,anova5)
ANOVAout$cohort <- c("ASD")

setwd("/Users/annieciernia/Box\ Sync/LaSalle\ Lab/Experiments/humanASD_WGBSDMRs/WGBS_QC/ASD")
write.xlsx(ANOVAout,file="ASD_WGBSQC_ANOVA.xlsx")

mergeDF <- merge(sampleinfo,DFtest,by.y="Sample",by.x="filename",all.y=T)
mergeDF2 <- merge(mergeDF,DFout1,by.x="filename",by.y="Sample")

write.xlsx(mergeDF2,file="ASD_WGBS_sampleinfo.xlsx")

##############################################################################
#Read alignment info from Bismark
##############################################################################

ASDinfo <- filter(sampleinfo) %>% filter(cohort == 1)
ASDinfo$PMI <- as.numeric(ASDinfo$PMI)

DFtest2 <- merge(ASDinfo,DFtest,by.x = "filename",by.y="Sample")
#model Group + sex
tmp <- lme(PMI ~ Group.x+Sex.x, ~1|filename, data = DFtest2)
#model indcludes random effects subject ID 
anova5 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova5$factors <- rownames(anova5)
anova5 <- anova5[!grepl('(Intercept)', anova5$factors),]
anova5$measure <- c("PMI")
anova5

write.xlsx(anova5,"ASD_PMI_ANOVA.xlsx")


#age
DFtest2$Age <- as.numeric(DFtest2$Age..yrs.)

#model Group + sex
tmp <- lme(Age ~ Group.x+Sex.x, ~1|filename, data = DFtest2)
#model indcludes random effects subject ID 
anova5 <- anova(tmp, type = "marginal") #marginal gives Type 3 SS for ANOVA
#write anova to output file
anova5$factors <- rownames(anova5)
anova5 <- anova5[!grepl('(Intercept)', anova5$factors),]
anova5$measure <- c("Age")
anova5

write.xlsx(anova5,"ASD_Age_ANOVA.xlsx")

#correlatoin global me and pmi
DFmaster <- merge(PerMemerge,DFtest,by="Sample")
DFmasterASD <- DFmaster %>% filter(cohort==1)

corrMEPMI <- cor.test(DFmasterASD$PMI,DFmasterASD$percent_cpg_meth, method="pearson")
coroutput <- data.frame(comparison = c("Pearson Correlation mCG/CG and PMI"),
  r = corrMEPMI$estimate,
pvalue= corrMEPMI$p.value)

write.xlsx(coroutput,"PearsonCorl_PMI_mCG.xlsx")

