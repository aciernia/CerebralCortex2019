# CerebralCortex2019
Analysis code associated with Vogel Ciernia et al. Cerebral Cortex 2019.

Code included in this folder was used to generate data for Vogel Ciernia et al. 2019 Cerebral Cortex analysis of Whole Genome Bisulfite Sequencing (WGBS) datasets from human cortex from Autism Spectrum Disorder, Dup15q Syndrome, and Rett Syndrome. Details on the tissue samples, library preparation, and sequencing can be found in the published article. 

Inital processing of fastq files, adapter trimming, QC and Bismark alignment was performed using:https://github.com/ben-laufer/CpG_Me
DMR calling was performed as described here: https://github.com/ben-laufer/DMRichR

Analysis scripts included here are mainly designed to be run interactively in R studio. Some more computationally intensive scripts require either a high performance desktop computer or access to a server. We provide the code and example datasets used in the manuscript. All code is provided "as is", so please use at your own risk. 

Folder contents are described as follows:

WGBSQC: R scripts for QC of WGBS datasets
DMR_Annotation: R scripts describing how genes were assigned to DMRs and assessment of median methylation differences. Data folder included.

GeneLengthAnalysis: R script for permutation based testing of gene lengths of DMR associated genes relative to BG region genes. Data folder included.

GeneOverlaps_lengthcorrected: R script for permutation based testing of gene list overlaps between DMR associated genes and published gene lists. Permutation distribution is taken from the BG regions and matched for gene length to the target (published) gene list. Example data is given for ASD DMR associated genes.

GOfunc: Gene Ontology enrichment analysis using GOFuncR permutation testing. The custom gene-region annotation is included as well as background regions and DMRs in bed format.

methmotifs: R scripts for graphing output data from http://bioinfo-csi.nus.edu.sg/methmotif/ analysis of DMRs. Script for analyzing smoothed methylation values from specific DMRs is also included. Data folder is included. 

ROIsmoothed: R scripts for extracting smoothed methylation values from WGBS data, graphing smoothed methylation values and conducting statistical analysis. Datasets included.

WGCNA_MouseMG: R scripts for cleaning, filtering and conducting differential analysis on RNA-seq data from mouse microglia across development from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5540146/
Raw data from GEO, filtered and normalized data are included. Includes R scripts for WGCNA analysis and enrichment for DMR associated genes. Output networks are also included.

Zhong_WGCNA: Rscripts for WGCNA analysis of single cell RNA-seq data from human brain across development from https://www.ncbi.nlm.nih.gov/pubmed/29539641. Datasets and output networks are included.

