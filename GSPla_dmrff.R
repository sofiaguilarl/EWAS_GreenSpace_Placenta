##################################################################################################
#DMRs  -dmrff R package
# Sofía Aguilar Lacasaña
##################################################################################################

### Creating a correlation matrix for DMRFF
### This function calculates the correlations between each CpG site and the next 20 CpG sites in the genome.
### The resulting information cannot be used to reveal any information about individual samples in the dataset so it is safe to share with external users.

### Load packages
library(matrixStats)
library(parallel)
#library(meffil) #to install meffil https://github.com/perishky/meffil/wiki/Installation
library(data.table)
library(dmrff)

#### LOAD METHYLATION DATA 
load("BISC_Pla_MethylEPIC_PreprocessedBetas_nooutliers_20231115.RData")
methylation <- betafinal.nooutlier # Name object representing methylation data meth 
str(methylation) 
methylation[1:5,1:5]  #rows are sites and columns are people. Transpose data if needed.rm(betafinal.nooutlier)

#-------------
# HOME
#-------------
# 
NDVI100<-fread("NDVI100_home_all_EPIC_AdjC/NDVI100_home_all_EPIC_AdjC_QCData.txt",data.table=F)
NDVI100<-NDVI100[order(NDVI100$P_VAL),]
head(NDVI100)
ewas_res<-NDVI100
rownames(ewas_res)<-ewas_res$probeID
# 
##We need the CpG site in rownames
## Methylation data
 methylation <- betafinal.nooutlier # Name object representing methylation data meth 
 str(methylation) 
 methylation[1:5,1:5]  #rows are sites and columns are people. Transpose data if needed.rm(betafinal.nooutlier)
# 
## Load annotation data and pull out annotation file.
##For the Illumina 450K microaray, we use the IlluminaHumanMethylation450kanno.ilmn12.hg19 package.
 library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# 
## Create an annotation data frame.
 data(list="IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
 data(Locations)
 data(Other)
 annotation <- cbind(as.data.frame(Locations), as.data.frame(Other))
 dim(annotation)  #485512     17
 annotation$name <- rownames(annotation)  #Rownames to colume to obtain CpG site name as a column
# 
## Match the methylation, ewas result, and annotation data and pull out beta, chr, pos, coef, se, pval.
 methylation<-methylation[na.omit(match(rownames(ewas_res),rownames(methylation))),] 
 ewas_res1 <- ewas_res[na.omit(match(rownames(methylation), rownames(ewas_res))),] 
 ifelse(all(rownames(ewas_res1)==rownames(methylation)), "meth and results data successfully matched :) ","Data not matched :(")
# 
 annotation <- annotation[na.omit(match(rownames(methylation), rownames(annotation))),]
 methylation<-methylation[na.omit(match(rownames(annotation),rownames(methylation))),]
 ifelse(all(rownames(annotation)==rownames(methylation)), "meth and annotation data successfully matched :) ","Data not matched :(")

## Double check of the order
#common.sites <- intersect(rownames(ewas_res), rownames(methylation))
#annotation <- annotation[rownames(annotation) %in% common.sites,]
#annotation <- annotation[order(annotation$chr, annotation$pos),] 
#methylation <- methylation[match(rownames(annotation),rownames(methylation)),]
#ewas_res1 <- ewas_res[match(rownames(annotation),rownames(ewas_res)),]

#
 coef <- ewas_res1$BETA ## coefficients for each CpG site from the EWAS results (note that beta refers to ewas result coefficients)
 se <- ewas_res1$SE ## standard errors for each CpG site from the EWAS results
 pval <- ewas_res1$P_VAL## p-values for each CpG site from the EWAS results
 chr <- strsplit(annotation$chr,"chr") #list with 1 element per CpG, saying the chrm it belongs to 
 chr <- do.call(rbind, chr)
 chr <- as.numeric(chr[,2]) ## chromosome of each cpg site in 'methylation'
 pos <- annotation$pos ## chromosomal position of each cpg site in 'methylation'
# 
# Run DMRs
 dmrs <- dmrff(estimate=coef, ## effect estimate for each CpG site
               se=se,       ## standard error of the estimate for each CpG site
               p.value=pval,  ## p-value
               methylation=methylation, ## methylation matrix
               chr=chr,      ## chromosome of each CpG site
               pos=pos)      ## position of each CpG site
 

 cohort<-"BISC" #define cohort name
 save(dmrs, file=paste0(cohort,"_GSpla_DMR_NDVI100preg_home_results_MAIN", format(Sys.Date(), "%d%m%Y"),".rda")) 
 str(dmrs)
# 
