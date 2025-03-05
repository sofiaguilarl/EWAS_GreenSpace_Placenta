#======================================================================
# Update 06.05.2024
# Sofía Aguilar Lacasaña
# Final - RUN EWAS GREEN SPACE EXPOSURE VS PLACENTAL DNA METHYLATION
#======================================================================

######################
# 1) Load packages
######################
library(PACEanalysis)
library(data.table)
library(parallel)

######################################################
# 2) Load phenodata with the exposures and covariates
######################################################
phenodataframe<-fread("EWASpla_GreenSpace_final_with_new_PC_13062024.txt",data.table=F)
samples.included<-phenodataframe$Basename

#-------------------
# Exposures
#-------------------
#"env.tot.ndvi100_preg.iqr"
#"env.tot.ndvi300_preg.iqr"
#"env.tot.ndvi500_preg.iqr"
#"env.home.ndvi100_preg.iqr"
#"env.home.ndvi300_preg.iqr"
#"env.home.ndvi500_preg.iqr"
#"env.tot_mgreen.dist_preg.iqr"
#"env.home_mgreen.dist_preg.iqr"
#"ts_all_GS_preg.iqr"
#"ts_public_GS_preg.iqr"

#--------------
#Covariates
#--------------
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
#Parity 2cat
phenodataframe$parity_m_2cat<-as.factor(phenodataframe$parity_m_2cat)
summary(phenodataframe$parity_m_2cat)
rownames(phenodataframe)<-phenodataframe$Basename

####################################################################
#3) Load methylation data processed data and cell type proportions
####################################################################

#---------------------------------------------------------
#Load the beta matrix after the QC and without outliers
#----------------------------------------------------------
load("BISC_Pla_MethylEPIC_PreprocessedBetas_nooutliers_20231115.RData")
#Keep the same number of samples than in the pheno
class(betafinal.nooutlier)
betafinal.nooutlier<-as.data.frame(betafinal.nooutlier)
head(colnames(betafinal.nooutlier))
betafinal.nooutlier.inc<-betafinal.nooutlier[,colnames(betafinal.nooutlier) %in% samples.included]
dim(betafinal.nooutlier.inc)
#---------------------------------------------------------------------
#Load these preprocessed object which has the celltype proportions
#---------------------------------------------------------------------
load("BISC_Pla_MethylEPIC_Preprocessed_20231115.RData")
celltypes<-processedOut$Omega
class(celltypes)
head(celltypes)
celltypes<-as.data.frame(celltypes)
head(rownames(celltypes))
#Keep the same number of samples than in the pheno
celltypes.inc<-celltypes[rownames(celltypes) %in% samples.included,]
dim(celltypes.inc)


###########################################################
#4) Check order across the samples of the different objects 
###########################################################
#--------------------------
#pheno + matrix DNA meth
#--------------------------
table(ifelse(phenodataframe$Basename==colnames(betafinal.nooutlier.inc),"Matched","--NOT MATCHED--"))
# Order samples
phenodataframe1<-phenodataframe[order(match(phenodataframe$Basename,colnames(betafinal.nooutlier.inc))),]
table(ifelse(phenodataframe1$Basename==colnames(betafinal.nooutlier.inc),"Matched","--NOT MATCHED--"))

#-----------------
#pheno + celltypes
#-----------------
table(ifelse(phenodataframe1$Basename==rownames(celltypes.inc),"Matched","--NOT MATCHED--"))

celltypes.inc<-as.matrix(celltypes.inc)
dim(celltypes.inc)
head(celltypes.inc)

betafinal.nooutlier.inc<-as.matrix(betafinal.nooutlier.inc)
dim(betafinal.nooutlier.inc)
head(betafinal.nooutlier.inc)

#############################
#5 Run analysis
#############################

## if running in parallel, checking the number of available cores 

#-----------------
#-----------------
#MAIN
#-----------------
#-----------------

# Select the exposure of interest
allvarsofinterest=c("env.home.ndvi100_preg.iqr","env.home.ndvi300_preg.iqr","env.home.ndvi500_preg.iqr","env.home_mgreen.dist_preg.iqr","ts_all_GS_preg.iqr")

for (i in 1:length(allvarsofinterest)){
  
  cat("Exposure:",allvarsofinterest[i],"\n")
  tempresults.MAIN<-dataAnalysis(phenofinal=phenodataframe1, #add your final phenodataframe
                  betafinal=betafinal.nooutlier.inc
                  array="EPIC", #specify the array type
                  maxit=100, #default
                  robust=TRUE, # specify "TRUE" if you want to run robust linear regression models
                  Omega=celltypes.inc, #add your final celltype db
                  vartype="ExposureCont", #specify which type are the variables of interest - in this case, we are interested in continuous exposures
                  varofinterest=allvarsofinterest[i],
                  #Add all variables that you want to add in the descriptive table
                  Table1vars=c("env.tot_ndvi.100","env.tot_ndvi.300","env.tot_ndvi.500","env.home_ndvi.100","env.home_ndvi.300","env.home_ndvi.500","env.tot_mgreen.dist","env.home_mgreen.dist","ts_all_GS","ts_all_GS_preg.iqr","ts_public_GS","env.tot.ndvi100_preg.iqr","env.tot.ndvi300_preg.iqr","env.tot.ndvi500_preg.iqr","env.home.ndvi100_preg.iqr","env.home.ndvi300_preg.iqr","env.home.ndvi500_preg.iqr",
                  "env.tot_mgreen.dist_preg.iqr","env.home_mgreen.dist_preg.iqr","ts_all_GS_preg.iqr","ts_public_GS_preg.iqr","Sex","Age","ethnicity_c_3cat","matedu_2cat","ses_income_anhh_2020","hospital_del_pa_3cat","covid_confinement_m_3cat",
                  "Meanlog2oddsContamination","smoke_any_m","gestage_0y_c_weeks","season_conception_m","pm25_home_hybrid","pm25_total_hybrid","parity_m_2cat","Complications"),
                  StratifyTable1=FALSE, #in case you want to create de descriptive table stratifying by an specific variable
                  StratifyTable1var=NULL,
                  #Add the adjusment variables (covariates)
                  adjustmentvariables=c("Sex","Age","ethnicity_c_3cat","matedu_2cat","ses_income_anhh_2020","hospital_del_pa_3cat","covid_confinement_m_3cat",
                  "Meanlog2oddsContamination","gestage_0y_c_weeks","smoke_any_m","parity_m_2cat"),
                  RunUnadjusted=TRUE,#model without adjusting for any covariate
                  RunAdjusted=TRUE,#model without adjusting for celltype proportions
                  RunCellTypeAdjusted=TRUE,#model adjusted for celltype proportions (in my case, the MAIN model)
                  RunSexSpecific=TRUE,# run analyses stratified by sex
                  RunCellTypeInteract=TRUE,#run interaction by celltype proportion
                  RestrictToSubset=FALSE, # If you want to run the analyses restricted to a subset of individuals (e.g. by Europeans) --> see code below
                  RestrictionVar=NULL,
                  RestrictToIndicator=NULL,
                  number_cores=4, #default                  
                  runparallel=TRUE,#default
                  destinationfolder="/EWASpla_GreenSpace_SA/results/",
                  savelog=TRUE,#default
                  cohort="BISC",analysisdate="20240506",
                  analysisname="MAIN")
                  
}


#---------------------------
#---------------------------
#RESTRICTED TO EUROPEANS
#---------------------------
#---------------------------
allvarsofinterest=c("env.home.ndvi100_preg.iqr","env.home.ndvi300_preg.iqr","env.home.ndvi500_preg.iqr","env.home_mgreen.dist_preg.iqr")

for (i in 1:length(allvarsofinterest)){
  
  cat("Exposure:",allvarsofinterest[i],"\n")
  tempresults.EUR<-dataAnalysis(phenofinal=phenodataframe1,
                  betafinal=betafinal.nooutlier.inc,
                  array="EPIC",
                  maxit=100,
                  robust=TRUE,
                  Omega=celltypes.inc,
                  vartype="ExposureCont",
                  varofinterest=allvarsofinterest[i],
                  Table1vars=c("env.tot_ndvi.100","env.tot_ndvi.300","env.tot_ndvi.500","env.home_ndvi.100","env.home_ndvi.300","env.home_ndvi.500","env.tot_mgreen.dist","env.home_mgreen.dist","ts_all_GS","ts_public_GS","env.tot.ndvi100_preg.iqr","env.tot.ndvi300_preg.iqr","env.tot.ndvi500_preg.iqr","env.home.ndvi100_preg.iqr","env.home.ndvi300_preg.iqr","env.home.ndvi500_preg.iqr",
                  "env.tot_mgreen.dist_preg.iqr","env.home_mgreen.dist_preg.iqr","ts_all_GS_preg.iqr","ts_public_GS_preg.iqr","Sex","Age","ethnicity_c_3cat","matedu_2cat","ses_income_anhh_2020","hospital_del_pa_3cat","covid_confinement_m_3cat",
                  "Meanlog2oddsContamination","smoke_any_m","gestage_0y_c_weeks","season_conception_m","pm25_home_hybrid","pm25_total_hybrid","parity_m_2cat"),
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
                  RestrictionVar="ethnicity_c_4cat",
                  RestrictToIndicator="eur",
                  number_cores=4,
                  runparallel=TRUE,
                  destinationfolder="/EWASpla_GreenSpace_SA/results/",
                 savelog=TRUE,
                  cohort="BISC",analysisdate="20240506",
                  analysisname="EUR")
                  
}


