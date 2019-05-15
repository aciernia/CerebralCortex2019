#!/usr/bin/env Rscript

# WGCNAsmooth
# Ben Laufer

# R settings --------------------------------------------------------------

rm(list=ls())
options(scipen=999)

if(length(grep("genomecenter.ucdavis.edu", .libPaths())) > 0){
  .libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
  AnnotationHub::setAnnotationHubOption("CACHE", "/share/lasallelab/programs/DMRichR/R_3.5")
}else{
  sink("DMRichR_log.txt", type = c("output", "message"), append = FALSE, split = TRUE)
}

# Functions ---------------------------------------------------------------

#' packageLoad
#' @description Install and load desired packages
#' @param packages Character string of desired packages
#' @export packageLoad
packageLoad <- function(packages = packages){
  message("\n","Checking for BiocManager and helpers...")
  CRAN <- c("BiocManager", "remotes", "magrittr")
  new.CRAN.packages <- CRAN[!(CRAN %in% installed.packages()[,"Package"])]
  if(length(new.CRAN.packages)>0){
    install.packages(new.CRAN.packages, repos ="https://cloud.r-project.org", quiet = TRUE)
  }
  message("Loading package management...")
  stopifnot(suppressMessages(sapply(CRAN, require, character.only = TRUE)))
  
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    message("\n","Installing missing packages...")
    new.packages <- packages %>%
      gsub("ggbiplot", "vqv/ggbiplot", .) %>% 
      gsub("DMRichR", "ben-laufer/DMRichR", .)
    BiocManager::install(new.packages, ask = FALSE, quiet = TRUE)
  }
  message("Loading packages...")
  stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
  suppressWarnings(BiocManager::valid(fix = TRUE, update = TRUE, ask = FALSE))
}

# Install and update ------------------------------------------------------

cat("\n[DMRichR] Installing and updating packages \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
packageLoad(c("tidyverse", "dmrseq", "annotatr", "rGREAT", "enrichR", "ChIPseeker", "BiocParallel", "ggbiplot",
              "liftOver", "openxlsx", "CMplot", "optparse", "gplots", "RColorBrewer", "broom", "lsmeans", "DMRichR"))

BiocManager::install("kdkorthauer/dmrseq")
library(dmrseq)

BiocManager::install("ben-laufer/DMRichR")
library(DMRichR)

# Consensus ---------------------------------------------------------------
message("Consensus")
setwd("/share/lasallelab/Annie/NDD_DMR_12_6_18/ROIsmoothed")
regions <- readPeakFile("/share/lasallelab/Annie/NDD_DMR_12_6_18/ROIsmoothed/hg38_cpg_shores.bed", as = "data.frame")
regions[2] <- regions[2] - 1
names(regions) <- c("chr", "start", "end")

# ASD ---------------------------------------------------------------------
message("ASD")
setwd("/share/lasallelab/Ben/ASD")
load("bsseq.RData")

getSmooth(bsseq = bs.filtered.bsseq,
          regions = regions,
          out = "/share/lasallelab/Annie/NDD_DMR_12_6_18/ROIsmoothed/ASD_Shores_smoothed_methylation.txt")

rm(bs.filtered.bsseq)

# Rett --------------------------------------------------------------------
message("Rett")
setwd("/share/lasallelab/Ben/Rett/")
load("bsseq.RData")

getSmooth(bsseq = bs.filtered.bsseq,
          regions = regions,
          out = "/share/lasallelab/Annie/NDD_DMR_12_6_18/ROIsmoothed/RTT_Shores_smoothed_methylation.txt")
rm(bs.filtered.bsseq)


# Dup15q ------------------------------------------------------------------
message("Dup15q")
setwd("/share/lasallelab/Ben/Dup15q/")
load("bsseq.RData")

getSmooth(bsseq = bs.filtered.bsseq,
          regions = regions,
          out = "/share/lasallelab/Annie/NDD_DMR_12_6_18/ROIsmoothed/Dup15q_Shores_smoothed_methylation.txt")

rm(bs.filtered.bsseq)

# End ---------------------------------------------------------------------

cat("\n[DMRichR] Finishing \t\t\t\t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")

sessionInfo()
rm(list = ls())
message("\n","Done...")
quit(save = "no", status = 0, runLast = FALSE)

