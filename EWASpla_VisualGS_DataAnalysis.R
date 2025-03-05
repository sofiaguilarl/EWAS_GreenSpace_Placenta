#############################
# Visual access - Categorical
# Sofía Aguilar Lacasaña
# Final 

#Load package
library(PACEanalysis)
library(data.table)
library(dplyr)


phenodataframe<-fread("GSpla_db_final_20240710.txt",data.table=F)

#12w
phenodataframe$GS3_m04_2cat<-as.factor(phenodataframe$GS3_m04_2cat)
levels(phenodataframe$GS3_m04_2cat)
phenodataframe$GS3_m04_2cat<-relevel(phenodataframe$GS3_m04_2cat,ref="less than half")
phenodataframe<-phenodataframe[(!is.na(phenodataframe$GS3_m04_2cat)),]
samples.included<-phenodataframe$Basename


#Covariates
#- Age
class(phenodataframe$Age)
summary(phenodataframe$Age)
#- Sex
phenodataframe$Sex<-as.factor(phenodataframe$Sex)
summary(phenodataframe$Sex)
#- GA
class(phenodataframe$gestage_0y_c_weeks)
summary(phenodataframe$gestage_0y_c_weeks)
#- ethnicity_c_4cat
phenodataframe$ethnicity_c_4cat<-as.factor(phenodataframe$ethnicity_c_4cat)
summary(phenodataframe$ethnicity_c_4cat)
#- matedu_2cat
phenodataframe$matedu_2cat<-as.factor(phenodataframe$matedu_2cat)
summary(phenodataframe$matedu_2cat)
#- ses_income_anhh_2020: proxy of SES NEIGH
class(phenodataframe$ses_income_anhh_2020)
summary(phenodataframe$ses_income_anhh_2020)
#- smoke_any_m
phenodataframe$smoke_any_m<-as.factor(phenodataframe$smoke_any_m)
summary(phenodataframe$smoke_any_m)
#-covid_confinement_m_3cat
phenodataframe$covid_confinement_m_3cat <-as.factor(phenodataframe$covid_confinement_m_3cat)
summary(phenodataframe$covid_confinement_m_3cat)
#- hospital_del_pa_3cat
phenodataframe$hospital_del_pa_3cat<-as.factor(phenodataframe$hospital_del_pa_3cat)
summary(phenodataframe$hospital_del_pa_3cat)
#- Meanlog2oddsContamination
class(phenodataframe$Meanlog2oddsContamination)
summary(phenodataframe$Meanlog2oddsContamination)

phenodataframe$parity_m_2cat<-as.factor(phenodataframe$parity_m_2cat)
summary(phenodataframe$parity_m_2cat)

rownames(phenodataframe)<-phenodataframe$Basename

#Load processed data
load("BISC_Pla_MethylEPIC_Preprocessed_20231115.RData")
load("BISC_Pla_MethylEPIC_PreprocessedBetas_nooutliers_20231115.RData")

#Celltypes
celltypes<-processedOut$Omega
class(celltypes)
head(celltypes)
celltypes<-as.data.frame(celltypes)
head(rownames(celltypes))
#Keep the same number of samples than in the pheno
celltypes.inc<-celltypes[rownames(celltypes) %in% samples.included,]
dim(celltypes.inc)
#Keep the same number of samples than in the pheno
class(betafinal.nooutlier)
betafinal.nooutlier<-as.data.frame(betafinal.nooutlier)
head(colnames(betafinal.nooutlier))
betafinal.nooutlier.inc<-betafinal.nooutlier[,colnames(betafinal.nooutlier) %in% samples.included]
dim(betafinal.nooutlier.inc)

#CHECK ORDER
#####################################################################
#pheno + matrix DNA meth
table(ifelse(phenodataframe$Basename==colnames(betafinal.nooutlier.inc),"Matched","--NOT MATCHED--"))
# Order samples
phenodataframe1<-phenodataframe[order(match(phenodataframe$Basename,colnames(betafinal.nooutlier.inc))),]
table(ifelse(phenodataframe1$Basename==colnames(betafinal.nooutlier.inc),"Matched","--NOT MATCHED--"))

#pheno + celltypes
table(ifelse(phenodataframe1$Basename==rownames(celltypes.inc),"Matched","--NOT MATCHED--"))

celltypes.inc<-as.matrix(celltypes.inc)
dim(celltypes.inc)
head(celltypes.inc)

betafinal.nooutlier.inc<-as.matrix(betafinal.nooutlier.inc)
dim(betafinal.nooutlier.inc)
head(betafinal.nooutlier.inc)

#Run analysis

## if running in parallel, checking the number of available cores 
library(parallel)

#####################
####MAIN
#####################
   

allvarsofinterest=c("GS3_m04_2cat")

for (i in 1:length(allvarsofinterest)){
  
  cat("Exposure:",allvarsofinterest[i],"\n")
  tempresults.MAIN<-dataAnalysis(phenofinal=phenodataframe1,
                  betafinal=betafinal.nooutlier.inc,
                  array="EPIC",
                  maxit=100,
                  robust=TRUE,
                  Omega=celltypes.inc,
                  vartype="ExposureCat",
                  varofinterest=allvarsofinterest[i],
                  Table1vars=c("GS3_32w_m04_2cat_completed_12w","GS3_m04","GS3_m04_2cat","GS3_32w_m04","GS3_32w_m04_2cat","Sex","Age","ethnicity_c_3cat","matedu_2cat","ses_income_anhh_2020","hospital_del_pa_3cat","covid_confinement_m_3cat",
                  "Meanlog2oddsContamination","smoke_any_m","gestage_0y_c_weeks","season_conception_m","pm25_home_hybrid","pm25_total_hybrid","pm25_total_lur","parity_m_2cat", "prematurity_c_2cat","fgr_c_32w","preg_compl","preg_compl_2cat",
                  "preg_compl_pla","preg_compl_pla_2cat","preg_compl_pla_baby","preg_compl_pla_baby_2cat"),
                  StratifyTable1=FALSE,
                  StratifyTable1var=NULL,
                  adjustmentvariables=c("Sex","Age","ethnicity_c_3cat","matedu_2cat","ses_income_anhh_2020","hospital_del_pa_3cat","covid_confinement_m_3cat",
                  "Meanlog2oddsContamination","gestage_0y_c_weeks","smoke_any_m","parity_m_2cat"),
                  RunUnadjusted=TRUE,
                  RunAdjusted=TRUE,
                  RunCellTypeAdjusted=TRUE,
                  RunSexSpecific=TRUE,
                  RunCellTypeInteract=TRUE,
                  RestrictToSubset=FALSE,
                  RestrictionVar=NULL,
                  RestrictToIndicator=NULL,
                  number_cores=4,
                  runparallel=TRUE,
                  destinationfolder="/EWASpla_GreenSpace_SA/results/",
                  savelog=TRUE,
                  cohort="BISC",analysisdate="20240710",
                  analysisname="MAIN")
                  
}

#---------------------------
#---------------------------
#RESTRICTED TO EUROPEANS
#---------------------------
#---------------------------

 allvarsofinterest=c("GS3_m04_2cat")
# 
 for (i in 1:length(allvarsofinterest)){
#   
   cat("Exposure:",allvarsofinterest[i],"\n")
   tempresults.EUR<-dataAnalysis(phenofinal=phenodataframe1,
                   betafinal=betafinal.nooutlier.inc,
                   array="EPIC",
                   maxit=100,
                   robust=TRUE,
                   Omega=celltypes.inc,
                   vartype="ExposureCat",
                   varofinterest=allvarsofinterest[i],
                   Table1vars=c("Visual_access_preg","Sex","Age","ethnicity_c_3cat","matedu_2cat","ses_income_anhh_2020","hospital_del_pa_3cat","covid_confinement_m_3cat",
                   "Meanlog2oddsContamination","smoke_any_m","gestage_0y_c_weeks","season_conception_m","pm25_home_hybrid","pm25_total_hybrid","pm25_total_lur","parity_m_2cat", "prematurity_c_2cat","fgr_c_32w","preg_compl","preg_compl_2cat",
                   "preg_compl_pla","preg_compl_pla_2cat","preg_compl_pla_baby","preg_compl_pla_baby_2cat"),
                   StratifyTable1=FALSE,
                   StratifyTable1var=NULL,
                   adjustmentvariables=c("Sex","Age","matedu_2cat","ses_income_anhh_2020","hospital_del_pa_3cat","covid_confinement_m_3cat",
                   "Meanlog2oddsContamination","gestage_0y_c_weeks","smoke_any_m","parity_m_2cat"),
                   RunUnadjusted=FALSE,
                   RunAdjusted=FALSE,
                   RunCellTypeAdjusted=TRUE,
                   RunSexSpecific=FALSE,
                   RunCellTypeInteract=FALSE,
                   RestrictToSubset=TRUE,
                   RestrictionVar="ethnicity_c_3cat",
                   RestrictToIndicator="eur",
                   number_cores=4,
                   runparallel=TRUE,
                   destinationfolder="/EWASpla_GreenSpace_SA/results/",
                   savelog=TRUE,
                  cohort="BISC",analysisdate="20240710",
                   analysisname="EUR")
                   
 }
