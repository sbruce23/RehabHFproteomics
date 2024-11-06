rm(list=ls())

library(tidyverse)
library(XLConnect)
library(readxl)
library(janitor)
library(gt)
library(gtsummary)
library(reshape2)
library(sjPlot)
library(ggforce)
library(gridExtra)
library(glinternet)
library(bartMachine)
library(optmatch)
library(rpart)
library(rpart.plot)
library(lsmeans)
library(treeClust)
library(pdfCluster)
library(mediation)
library(dplyr)

setwd("~/Library/CloudStorage/OneDrive-SharedLibraries-HarvardUniversity/REHAB-HFpEF Ancillary - General/REHAB-HF proteomics/REHAB-HF proteomics data")

############################################
## 0.1 Data dictionary
#############################################
#age (years)
#sex (1 for female, 0 for male)
#race___4 (1 for white, 0 for non-white)
#hf_cat (1 for hfpef (EF >=45%), 0 for hfref (EF < 45%))
#creatinine_value (mg/dL)
#bl_sppb, fu_sppb (baseline (bl) and follow up (fu) SPPB score)
#bl_smw, fu_smw (baseline (bl) and follow up (fu) six minute walk test distance (m))
#rehosp (rehospitalization)
#death (death)
#event (rehospitalization or death)
#lfu (lost to follow up)
#eGFR (derived from age, sex, and creatinine level 
#(https://www.kidney.org/content/ckd-epi-creatinine-equation-2021))

############################################
## 0.2 Set theme for ggplot2 plots ##
############################################
hw <- theme_gray()+ theme(
  plot.title=element_text(hjust=0.5,size=18,face="bold"),
  plot.subtitle=element_text(hjust=0.5,size=12),
  plot.caption=element_text(hjust=-.5,size=10),
  strip.background=element_rect(fill=rgb(.9,.95,1),
                                colour=gray(.5), linewidth=.2),
  panel.border=element_rect(fill=FALSE,colour=gray(.70)),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.spacing.x = unit(0.2,"cm"),
  panel.spacing.y = unit(0.2,"cm"),
  axis.text=element_text(colour="black",size=10),
  axis.text.y=element_text(margin=ggplot2::margin(0,3,0,3)),
  axis.text.x=element_text(margin=ggplot2::margin(-1,0,3,0)),
  axis.title=element_text(size=16,face="bold"),
  legend.text=element_text(size=14),
  legend.title = element_blank(),
);

############################################
## 1.1 Read in and initial data pre-processing ##
############################################


olink_dta=read_excel("Rehab-HF Olink_deFilippi_CVDIII NPX raw data.xlsx", range = "A8:CO418", col_names = FALSE);
olink_colnames=read_excel("Rehab-HF Olink_deFilippi_CVDIII NPX raw data.xlsx", range = "B4:CO4", col_names = FALSE);
olink_olinkid=read_excel("Rehab-HF Olink_deFilippi_CVDIII NPX raw data.xlsx", range = "B4:CO6", col_names = FALSE);
olink_olinkid = as.data.frame(t(olink_olinkid));
colnames(olink_olinkid) = c('Assay','UniProt','OlinkID')

olink_lod=as.numeric(read_excel("Rehab-HF Olink_deFilippi_CVDIII NPX raw data.xlsx", range = "B420:CO420", col_names = FALSE));
colnames(olink_dta)=c('ID',as.character(olink_colnames));
#15% missingness in SPON1
#4% missingness in CHIT1
#some minor missingness across others

#baseline
olink_dta_bl=olink_dta[grep("BL",olink_dta$ID),]
olink_dta_bl=olink_dta_bl[-grep("BL_2",olink_dta_bl$ID),] #remove the extra BL reading for 20-017
olink_dta_bl$timepoint="Baseline";
olink_dta_bl$subject_id=gsub(" BL","",olink_dta_bl$ID,fixed=TRUE)

#followup
olink_dta_fu=olink_dta[grep("FU",olink_dta$ID),]
olink_dta_fu$timepoint="Follow Up";
olink_dta_fu$subject_id=gsub(" FU","",olink_dta_fu$ID,fixed=TRUE)
olink_dta_fu$subject_id=gsub("FU","",olink_dta_fu$subject_id,fixed=TRUE)

setwd("~/Library/CloudStorage/OneDrive-SharedLibraries-HarvardUniversity/REHAB-HFpEF Ancillary - General/REHAB-HF dataset");

## Read in data
rhfwb <- loadWorkbook("MERGED_FINAL2_dictionary.xlsx", password="RehabHF");
rhfds <- readWorksheet(rhfwb, "MERGED_FINAL2");
rbiowb <- loadWorkbook("hs-cTn Rehab-HF 20210551_deFilippi sample manifest_RESULTS.xls", password="06242022");
rbiods <- readWorksheet(rbiowb, "Sheet1",endRow = 402);

## Clean up column names and data types
rbiods <- clean_names(rbiods);
colnames(rbiods)[3] <- "sample_id";
colnames(rbiods)[4] <- "sample_type";
colnames(rbiods)[5] <- "sample_volume_ml";
colnames(rbiods)[9] <- "hf_category_1_hfpef_0_hfref";
colnames(rbiods)[12] <- "creatinine_mg_dl";
colnames(rbiods)[13] <- "troponin_i_pg_ml";

rhfds$sex <- factor(rhfds$sex,levels=c(0,1));
rhfds$race___4 <- factor(rhfds$race___4,levels=c(0,1));
rhfds$hf_cat <- factor(rhfds$hf_cat,levels=c(0,1));
rhfds$creatinine_value <- as.numeric(rhfds$creatinine_value)
rhfds$bl_sppb <- as.numeric(rhfds$bl_sppb)
rhfds$fu_sppb <- as.numeric(rhfds$fu_sppb)
rhfds$bl_smw <- as.numeric(rhfds$bl_smw)
rhfds$fu_smw <- as.numeric(rhfds$fu_smw)


## Create eGFR
kpa = ifelse(rhfds$sex==1,0.7,0.9);
alpha = ifelse(rhfds$sex==1,-0.241,-0.302)
const = ifelse(rhfds$sex==1,1.012,1)
rhfds$egfr = 142*
  (pmin(as.numeric(rhfds$creatinine_value)/kpa,1)^alpha)*
  (pmax(as.numeric(rhfds$creatinine_value)/kpa,1)^-1.2)*
  (0.9938^rhfds$age)*const

## Clean up entries in rbiods

## Find and remove all (D) indicators
rbiods <- apply(rbiods,2,function(x) gsub("(D)","",x,fixed=TRUE));

## Find and remove all empty spaces in biomarker columns 12-16
rbiods[,12:16] <- apply(rbiods[,12:16],2,function(x) gsub(" ","",x,fixed=TRUE));

## Change QNS to NA and add QNS indicator for all biomarkers (columns 12-16)
qns_ind <- apply(rbiods[,12:16],2,function(x) as.numeric(grepl("QNS",x,fixed=TRUE)));
colnames(qns_ind) <- paste(colnames(qns_ind),"_qns_ind",sep="");
rbiods <- as.data.frame(cbind(rbiods,qns_ind));
rbiods[,12:16] <- apply(rbiods[,12:16],2,function(x) gsub("QNS",NA,x,fixed=TRUE));

## Add censoring indicator columns and change values to NA
cens_ind <- matrix(0,nrow=nrow(rbiods),ncol=5)
cens_ind[apply(rbiods[,12:16],2,function(x) grepl("[><]",x))] <- 
  unlist(apply(rbiods[,12:16],2,function(x) x[grep("[><]",x)]));
colnames(cens_ind) <- paste(colnames(rbiods)[12:16],"_cens_ind",sep="");
rbiods <- as.data.frame(cbind(rbiods,cens_ind));
rbiods[,12:16] <- apply(rbiods[,12:16],2,function(x) gsub("[><]",NA,x));

## Impute censored observations: 
## Halfway between 0 and lower limit of detection for left censored observations
## Upper limit of detection for right censored observations
rbiods$creatinine_mg_dl[rbiods$creatinine_mg_dl_cens_ind != 0] <- 0.15/2;
rbiods$hs_crp_mg_l[rbiods$hs_crp_mg_l_cens_ind != 0] <- 90;
rbiods$nt_pro_bnp[rbiods$nt_pro_bnp_cens_ind != 0] <- 35000;
rbiods$troponin_t_ng_l[rbiods$troponin_t_ng_l_cens_ind != 0] <- 6/2;

## Convert biomarker fields and censoring and QNS indicators to numeric
rbiods[,12:21] <- apply(rbiods[,12:21],2,function(x) as.numeric(x));


#merge (left join olink to rhfds) on subject_id in olink_dta_bl and study_id in rhfds
varlist = c("age","sex","race___4","hf_cat","creatinine_value",
            "bl_sppb","fu_sppb","bl_bal_scr","fu_bal_scr",
            "bl_walk_time","fu_walk_time","rehosp","death","event","lfu","bl_smw","fu_smw","nrehosp_death","egfr")
rprotblds <- merge(x=olink_dta_bl,y=rhfds[,c("study_id",varlist)],
                   by.x='subject_id',by.y='study_id',all.x=TRUE)

rprotfuds <- merge(x=olink_dta_fu,y=rhfds[,c("study_id",varlist)],
                   by.x='subject_id',by.y='study_id',all.x=TRUE)


#add intervention indicator from rbiods
rprotblds <- merge(x=rprotblds,y=rbiods[rbiods$timepoint=='Baseline',
                                        c('subject_id','intervention_1_control_0')],
                   by.x='subject_id',by.y='subject_id',all.x=TRUE)
rprotblds$intervention_1_control_0[rprotblds$subject_id=='10-031'] = 0; #10-031 only had measurement at follow up but is in control
rprotblds$intervention_1_control_0 <- factor(rprotblds$intervention_1_control_0,
                                             levels=c(0,1))

rprotfuds <- merge(x=rprotfuds,y=rbiods[rbiods$timepoint=='Baseline',
                                        c('subject_id','intervention_1_control_0')],
                   by.x='subject_id',by.y='subject_id',all.x=TRUE)
rprotfuds$intervention_1_control_0 <- factor(rprotfuds$intervention_1_control_0,
                                             levels=c(0,1))

# #means and sds of protein expressions
# apply(rprotblds[,3:94],2,mean)
# apply(rprotblds[,3:94],2,sd)

############################################
## 1.2 Baseline biomarker outlier removal and final pre-processing ##
############################################

## Identify and remove outliers in baseline biomarker levels

## How many outliers are there for each biomarker?
## More than 3 standard deviations beyond the mean after log (base 2) transformation
apply(rbiods[rbiods$timepoint=='Baseline',12:16],2,
      function(x) sum((log2(x+0.001)>mean(log2(x+0.001),na.rm=TRUE)+3*sd(log2(x+0.001),na.rm=TRUE))|
                        (log2(x+0.001)<mean(log2(x+0.001),na.rm=TRUE)-3*sd(log2(x+0.001),na.rm=TRUE)),na.rm=TRUE))

## Which rows correspond to outliers for each biomarker?
bmoutliers <- apply(rbiods[rbiods$timepoint=='Baseline',12:16],2,
                    function(x) rbiods[rbiods$timepoint=='Baseline',c("sample_id","timepoint")]
                    [which((log2(x+0.001)>mean(log2(x+0.001),na.rm=TRUE)+
                              3*sd(log2(x+0.001),na.rm=TRUE))|
                             (log2(x+0.001)<mean(log2(x+0.001),na.rm=TRUE)-
                                3*sd(log2(x+0.001),na.rm=TRUE))),])

## Replace outlier values with NA
for (i in 1:length(bmoutliers)){
  bmtmp <- names(bmoutliers)[i]
  if (nrow(bmoutliers[[i]]>0)){
    for (j in 1:nrow(bmoutliers[[i]])){
      print(bmoutliers[[i]][j,])
      idtmp <- bmoutliers[[i]][j,1]
      rbiods[rbiods$timepoint=='Baseline' & rbiods$sample_id == idtmp,bmtmp] <- NA
    }
  }
}

## Merge rhfds and rbiods
rcombds <- merge(x=rbiods,y=rhfds,
                 by.x='subject_id',by.y='study_id',all.x=TRUE);

############################################
## 2.1 Protein EMM variable screening step 1: individual regressions ##
############################################
setwd("/Users/scottbruce/Library/CloudStorage/OneDrive-SharedLibraries-HarvardUniversity/REHAB-HFpEF Ancillary - General/REHAB-HF proteomics/REHAB-HF proteomics analysis")

#sppb
intxn_df = matrix(NA,nrow=92,ncol=2)
colnames(intxn_df) = c("effect","p-value")
for (i in 3:94){
  prot=names(rprotblds)[i]
  formtmp=paste("fu_sppb ~ bl_sppb + age + sex + race___4 + hf_cat + `",prot,
                "`*intervention_1_control_0",sep="")
  tmpmod=summary(lm(data=rprotblds,formula=formtmp))
  intxn_df[i-2,]=tmpmod$coefficients[nrow(tmpmod$coefficients),c(1,4)]
}

intxn_df = data.frame(Protein=names(rprotblds)[3:94],intxn_df,
                      Significance=ifelse(intxn_df[,2]<0.05,"Significant","Not Significant"))
summary(p.adjust(intxn_df$p.value,method='fdr'))
#none are significant after FDR adjustment

volcplot.sppb <- ggplot(data=intxn_df,mapping=aes(x=effect,y=-log10(p.value),color=Significance))+geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dotted")+
  labs(x = "Effect modification", y = "-log10(p-value)") +
  OlinkAnalyze::set_plot_theme() + ggtitle('SPPB Effect Modification (Unadjusted P-values)')+
  ggrepel::geom_label_repel(data = intxn_df[intxn_df$Significance=="Significant",],
                            ggplot2::aes(label = Protein), box.padding = 1, show.legend = FALSE,max.overlaps=Inf) 


pdf("Results/SPPB_volcano_effectmod.pdf",width=7,height=6,)
print(volcplot.sppb)
dev.off()

ggsave(volcplot.sppb,file='Results/SPPB_volcano_effectmod.eps',width=7,height=6,device="eps")



#6mwd
intxn_df = matrix(NA,nrow=92,ncol=2)
colnames(intxn_df) = c("effect","p-value")
for (i in 3:94){
  prot=names(rprotblds)[i]
  formtmp=paste("fu_smw ~ bl_smw + age + sex + race___4 + hf_cat + `",prot,
                "`*intervention_1_control_0",sep="")
  tmpmod=summary(lm(data=rprotblds,formula=formtmp))
  intxn_df[i-2,]=tmpmod$coefficients[nrow(tmpmod$coefficients),c(1,4)]
}

intxn_df = data.frame(Protein=names(rprotblds)[3:94],intxn_df,
                      Significance=ifelse(intxn_df[,2]<0.05,"Significant","Not Significant"))
summary(p.adjust(intxn_df$p.value,method='fdr'))
#none are significant after FDR adjustment

volcplot.6mwd <- ggplot(data=intxn_df,mapping=aes(x=effect,y=-log10(p.value),color=Significance))+geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dotted")+
  labs(x = "Effect modification", y = "-log10(p-value)") +
  OlinkAnalyze::set_plot_theme() + ggtitle('6 Minute Walk Effect Modification (Unadjusted P-values)')+
  ggrepel::geom_label_repel(data = intxn_df[intxn_df$Significance=="Significant",],
                            ggplot2::aes(label = Protein), box.padding = 1, show.legend = FALSE,max.overlaps=Inf) 

pdf("Results/6MWD_volcano_effectmod.pdf",width=7,height=6,)
print(volcplot.6mwd)
dev.off()

ggsave(volcplot.6mwd,file='Results/6MWD_volcano_effectmod.eps',width=7,height=6,device="eps")

############################################
## 2.2 Protein EMM variable screening step 2: LASSO regression with glinternet interaction screening ##
############################################
#sppb
sppb_data <- cbind(rprotblds[,c("fu_sppb","intervention_1_control_0","bl_sppb","age","sex","race___4","hf_cat")],
                   rprotblds[,3:94])

#remove incomplete cases (any rows with missingness, mostly due to lack of follow up)
sppb_data = na.omit(sppb_data);

# ##using glmnet doesn't enforce strong hierarchy
# f <- as.formula(fu_sppb ~ .+intervention_1_control_0*.)
# sppb_y <- sppb_data$fu_sppb
# sppb_x <- model.matrix(f, sppb_data)[, -1]
# 
# fit <- glmnet(sppb_x, sppb_y)
# plot(fit)

##using glinternet, which does enforce strong hierarchy
set.seed(23)
sppb_y <- sppb_data$fu_sppb
#reformat factors as 0/1 vars
sppb_x <- data.frame(intervention_1_control_0=as.numeric(sppb_data$intervention_1_control_0)-1,
                     sppb_data$bl_sppb,sppb_data$age,
                     sex=as.numeric(sppb_data$sex)-1,
                     race___4=as.numeric(sppb_data$race___4)-1,
                     hf_cat=as.numeric(sppb_data$hf_cat)-1,
                     sppb_data[,8:99])
#number of levels for each independent var (1 for continuous vars)
numLevels <- c(2,1,1,2,2,2,rep(1,92))
fit <- glinternet(sppb_x,sppb_y,numLevels=numLevels,
                  interactionCandidates = 1)
coeffs <- coef(fit)

fit <- glinternet.cv(sppb_x,sppb_y,numLevels=numLevels,
                     interactionCandidates = 1)

#looking through results
# plot(fit)
# fit$lambdaHat #0.01029845
# fit$lambdaHat1Std #0.07411649
# fit$activeSet #variables selected using lambdaHat above
# fit$betahat
# fit$family
# coeffs=coef(fit)
# coeffs$mainEffects
# coeffs$interactions
#categorical variables
#1: intervention vs control
#2: sex
#3: race
#4: hf_cat

#continuous variables
#1: bl_sppb
#2: age
#3-94: proteins

#proteins with treatment effect interactions selected
names(sppb_data[,8:99])[coeffs$interactions$catcont[,2]-2]
# [1] "LDL receptor"  "ALCAM"         "BLM hydrolase" "SPON1"        
# [5] "CXCL16"        "GP6"           "CCL16"        

coeffs$interactionsCoef
#give you estimates of the deviation from the main effect slope for control vs treatment

names(sppb_data[,8:99])[coeffs$mainEffects$cont[c(-1,-2)]-2]
unlist(coeffs$mainEffectsCoef$cont)[c(-1,-2)]

# #fit a model
# summary(lm(formula=fu_sppb ~ intervention_1_control_0+hf_cat+
#              bl_sppb+`BLM hydrolase`+`TIMP4`+`TLT-2`+`TR`+SELE+MPO+`IGFBP-1`+
#              PI3+`AP-N`+TNFSF13B+PCSK9+`U-PAR`+`SHPS-1`+CCL15+`CASP-3`+CPB1+
#              CHI3L1+ST2+SCGB3A2+PON3+CTSZ+`LDL receptor`+ALCAM+SPON1+CXCL16+       
#              +GP6+CCL16+
#              `LDL receptor`:intervention_1_control_0+ALCAM:intervention_1_control_0+
#              `BLM hydrolase`:intervention_1_control_0+SPON1:intervention_1_control_0+        
#              `CXCL16`:intervention_1_control_0+GP6:intervention_1_control_0+
#              `CCL16`:intervention_1_control_0        
#            ,data=rprotblds))

#fit models considering lambda values between fit$lambdaHat (0.01029845) and fit$lambdaHat1Std (0.07411649)
#find the model that is a bit more parismonious than the one at lambdaHat by increasing lambda slightly
#to reduce the model size
fit.step.sppb <- glinternet(sppb_x,sppb_y,numLevels=numLevels,lambda=fit$lambda[7:28],
                            interactionCandidates = 1)

lapply(coef(fit.step.sppb),function(x) names(sppb_data[,8:99])[x$interactions$catcont[,2]-2])
lapply(coef(fit.step.sppb),function(x) names(sppb_data[,8:99])[x$interactions$catcont[,2]-2])[[17]]
#[1] "LDL receptor" "ALCAM"        "GP6"          "CCL16"  
fit.step.sppb$lambda[17]
#0.01647607

#lambda slightly higher than lambdaHat, but returns 4 EMMS vs. 7

#6mwd
smw_data <- cbind(rprotblds[,c("fu_smw","intervention_1_control_0","bl_smw","age","sex","race___4","hf_cat")],
                  rprotblds[,3:94])

#remove incomplete cases (any rows with missingness, mostly due to lack of follow up)
smw_data = na.omit(smw_data);

# ##using glmnet, but unfortunately this doesn't enforce strong hierarchy
# f <- as.formula(fu_smw ~ .+intervention_1_control_0*.)
# smw_y <- smw_data$fu_smw
# smw_x <- model.matrix(f, smw_data)[, -1]
# 
# fit <- glmnet(smw_x, smw_y)
# plot(fit)

##using glinternet, which does enforce strong hierarchy
set.seed(34)
smw_y <- smw_data$fu_smw
#reformat factors as 0/1 vars
smw_x <- data.frame(intervention_1_control_0=as.numeric(smw_data$intervention_1_control_0)-1,
                    smw_data$bl_smw,smw_data$age,
                    sex=as.numeric(smw_data$sex)-1,
                    race___4=as.numeric(smw_data$race___4)-1,
                    hf_cat=as.numeric(smw_data$hf_cat)-1,
                    smw_data[,8:99])
#number of levels for each independent var (1 for continuous vars)
numLevels <- c(2,1,1,2,2,2,rep(1,92))
fit <- glinternet(smw_x,smw_y,numLevels=numLevels,
                  interactionCandidates = 1)
coeffs = coef(fit)

fit <- glinternet.cv(smw_x,smw_y,numLevels=numLevels,
                     interactionCandidates = 1)

#looking through results
# plot(fit)
# fit$lambdaHat #0.375042
# fit$lambdaHat1Std #6.288827
# fit$activeSet #variables selected using lambdaHat above
# fit$betahat
# fit$family
# coeffs=coef(fit)
# coeffs$mainEffects
# coeffs$interactions
#categorical variables
#1: intervention vs control
#2: sex
#3: race
#4: hf_cat

#continuous variables
#1: bl_smw
#2: age
#3-94: proteins

#proteins with treatment effect interactions selected
names(smw_data[,8:99])[coeffs$interactions$catcont[,2]-2]
# [1] "LDL receptor"   "ALCAM"          "SELP"           "FABP4"          "CHIT1"         
# [6] "Ep-CAM"         "Gal-4"          "t-PA"           "PDGF subunit A"

coeffs$interactionsCoef
#give you estimates of the deviation from the main effect slope for control vs treatment

names(smw_data[,8:99])[coeffs$mainEffects$cont[c(-1,-2)]-2]
unlist(coeffs$mainEffectsCoef$cont)[c(-1,-2)]

#fit a model
# summary(lm(formula=fu_smw ~ intervention_1_control_0+
#              bl_smw+age+`LDL receptor`+`IL-17RA`+`IL2-RA`+`Gal-3`+`BLM hydrolase`+`Notch 3`+`DLK-1`+
#              CXCL16+`IL-6RA`+GP6+`PSP-D`+`PI3`+`AP-N`+`MMP-2`+`TNFSF13B`+PRTN3+
#              `SHPS-1`+uPA+CHI3L1+SCGB3A2+COL1A1+`vWF`+ALCAM+SELP+FABP4+CHIT1+`Ep-CAM`+`Gal-4`    
#            +`t-PA`+`PDGF subunit A`+
#              `LDL receptor`:intervention_1_control_0+ALCAM:intervention_1_control_0+
#              SELP:intervention_1_control_0+FABP4:intervention_1_control_0+        
#              CHIT1:intervention_1_control_0+`Ep-CAM`:intervention_1_control_0+
#              `Gal-4`:intervention_1_control_0 +`t-PA`:intervention_1_control_0+
#              `PDGF subunit A`:intervention_1_control_0
#            ,data=rprotblds))

#fit models considering lambda values between fit$lambdaHat (0.375042) and fit$lambdaHat1Std (6.288827)
#find the model that is a bit more parismonious than the one at lambdaHat by increasing lambda slightly
#to reduce the model size
fit.step <- glinternet(smw_x,smw_y,numLevels=numLevels,lambda=fit$lambda[2:32],
                       interactionCandidates = 1)

lapply(coef(fit.step),function(x) names(smw_data[,8:99])[x$interactions$catcont[,2]-2])
lapply(coef(fit.step),function(x) names(smw_data[,8:99])[x$interactions$catcont[,2]-2])[[23]]
#[1] "LDL receptor" "ALCAM"        "PLC"          "CHIT1"        "Gal-4" 
fit.step$lambda[23]
#0.7954456

#lambda slightly higher than lambdaHat, but returns 5 EMMS vs. 9

############################################
## 2.3 Protein EMM variable screening step 3: BART with interaction importance ##
############################################

set.seed(45)
options(java.parameters = "-Xmx4000m")
set_bart_machine_num_cores(1)


#model fit
sppb_data <- cbind(rprotblds[,c("fu_sppb","intervention_1_control_0","bl_sppb","age","sex","race___4","hf_cat")],
                   rprotblds[,3:94])

#remove incomplete cases (any rows with missingness, mostly due to lack of follow up)
sppb_data = na.omit(sppb_data);

sppb_y <- sppb_data$fu_sppb
sppb_x <- data.frame(sppb_data$intervention_1_control_0,
                     sppb_data$bl_sppb,sppb_data$age,
                     sppb_data$sex,
                     sppb_data$race___4,
                     sppb_data$hf_cat,
                     sppb_data[,8:99])

bart_machine_sppb <- bartMachine(sppb_x, sppb_y,
                                 num_trees=100,
                                 num_burn_in=500,
                                 num_iterations_after_burn_in=2000,
                                 use_missing_data=TRUE,
                                 mem_cache_for_speed = FALSE,
                                 seed=34)
bart_machine_sppb
pdf("Results/bart_sppb_ConvergenceDiagnostics.pdf", onefile = TRUE)
plot_convergence_diagnostics(bart_machine_sppb)
dev.off()

#if doesn't work, run install.packages('bartMachine') to update package I guess?
reps=100
intxninvest=interaction_investigator(bart_machine_sppb, 
                                     plot=FALSE,
                                     num_var_plot=25,
                                     bottom_margin=20,
                                     num_replicates_for_avg=reps,
                                     num_trees_bottleneck=20)

#average count of interactions between protein and treatment/control
sort(rowSums(intxninvest$interaction_counts_avg[3:94,95:96]),decreasing=TRUE)

#mean and sd of interaction counts
cbind(intxninvest$interaction_counts_avg[3:94,95:96],
      1.96*intxninvest$interaction_counts_sd[3:94,95:96]/sqrt(reps))

#sum of count of interactions
bartintxn_means=apply(intxninvest$interaction_counts[3:94,95:96,],1,sum)/reps
bartintxn_sds=apply(apply(intxninvest$interaction_counts[3:94,95:96,],c(1,3),sum),1,sd)
bartintxn_moe=1.96 * bartintxn_sds / sqrt(reps)
bartintxn_df=data.frame(Protein=rownames(intxninvest$interaction_counts_avg)[3:94],
                        Importance=bartintxn_means,MOE=bartintxn_moe)
bartintxn_df = bartintxn_df[order(bartintxn_df$Importance,decreasing=TRUE),]
#top 20 only
bartintxn_df=bartintxn_df[1:20,]
pdf("Results/InteractionImportance_sppb.pdf", onefile = TRUE)
par(mar = c(7, 6, 3, 0))
bars = barplot(bartintxn_df$Importance, 
               names.arg = bartintxn_df$Protein, 
               las = 2, 
               ylab = "Importance", 
               col = "gray",
               ylim = c(0, max(bartintxn_df$Importance + bartintxn_df$MOE))
)
conf_upper = bartintxn_df$Importance + bartintxn_df$MOE
conf_lower = bartintxn_df$Importance - bartintxn_df$MOE
segments(bars, bartintxn_df$Importance, bars, conf_upper, col = rgb(0.59, 0.39, 0.39), lwd = 3) # Draw error bars
segments(bars, bartintxn_df$Importance, bars, conf_lower, col = rgb(0.59, 0.39, 0.39), lwd = 3)
dev.off()





set.seed(231)
options(java.parameters = "-Xmx4000m")
set_bart_machine_num_cores(1)


#model fit
smw_data <- cbind(rprotblds[,c("fu_smw","intervention_1_control_0","bl_smw","age","sex","race___4","hf_cat")],
                  rprotblds[,3:94])

#remove incomplete cases (any rows with missingness, mostly due to lack of follow up)
smw_data = na.omit(smw_data);

smw_y <- smw_data$fu_smw
smw_x <- data.frame(smw_data$intervention_1_control_0,
                    smw_data$bl_smw,smw_data$age,
                    smw_data$sex,
                    smw_data$race___4,
                    smw_data$hf_cat,
                    smw_data[,8:99])

bart_machine_smw <- bartMachine(smw_x, smw_y,
                                num_trees=100,
                                num_burn_in=2000,
                                num_iterations_after_burn_in=2000,
                                use_missing_data=TRUE,
                                mem_cache_for_speed = FALSE,
                                seed=231)
bart_machine_smw
pdf("Results/bart_smw_ConvergenceDiagnostics.pdf", onefile = TRUE)
plot_convergence_diagnostics(bart_machine_smw)
dev.off()

reps=100
intxninvest=interaction_investigator(bart_machine_smw, 
                                     num_var_plot=25,
                                     bottom_margin=20,
                                     num_replicates_for_avg=reps,
                                     num_trees_bottleneck=20)

#average count of interactions between protein and treatment/control
sort(rowSums(intxninvest$interaction_counts_avg[3:94,95:96]),decreasing=TRUE)

#mean and sd of interaction counts
cbind(intxninvest$interaction_counts_avg[3:94,95:96],
      1.96*intxninvest$interaction_counts_sd[3:94,95:96]/sqrt(reps))

#sum of count of interactions
bartintxn_means=apply(intxninvest$interaction_counts[3:94,95:96,],1,sum)/reps
bartintxn_sds=apply(apply(intxninvest$interaction_counts[3:94,95:96,],c(1,3),sum),1,sd)
bartintxn_moe=1.96 * bartintxn_sds / sqrt(reps)
bartintxn_df=data.frame(Protein=rownames(intxninvest$interaction_counts_avg)[3:94],
                        Importance=bartintxn_means,MOE=bartintxn_moe)
bartintxn_df = bartintxn_df[order(bartintxn_df$Importance,decreasing=TRUE),]
#top 20 only
bartintxn_df=bartintxn_df[1:20,]
pdf("Results/InteractionImportance_6MWD.pdf", onefile = TRUE)
par(mar = c(7, 6, 3, 0))
bars = barplot(bartintxn_df$Importance, 
               names.arg = bartintxn_df$Protein, 
               las = 2, 
               ylab = "Importance", 
               col = "gray",
               ylim = c(0, max(bartintxn_df$Importance + bartintxn_df$MOE))
)
conf_upper = bartintxn_df$Importance + bartintxn_df$MOE
conf_lower = bartintxn_df$Importance - bartintxn_df$MOE
segments(bars, bartintxn_df$Importance, bars, conf_upper, col = rgb(0.59, 0.39, 0.39), lwd = 3) # Draw error bars
segments(bars, bartintxn_df$Importance, bars, conf_lower, col = rgb(0.59, 0.39, 0.39), lwd = 3)
dev.off()

############################################
## 3. Table 1 (Baseline Patient Characteristics) ##
############################################

#LASSO selected vars for 6MWD
# [1] "LDL receptor"   "ALCAM"          "PLC"         "CHIT1"     "Gal-4"         

#BART selected vars for 6MWD: Gal-4 (already selected via LASSO)

#LASSO selected vars for SPPB
# #[1] "LDL receptor" "ALCAM"        "GP6"          "CCL16"      

#BART selected vars for SPPB: ST2, CCL16


protlist = c("LDL receptor","ALCAM","PLC","CHIT1","Gal-4","GP6","CCL16","ST2")


varlist = c(c("intervention_1_control_0","age","sex","race___4",
            "hf_cat","egfr","bl_sppb","bl_smw"),protlist)

tbl_summary(rprotblds[,varlist],
            by=intervention_1_control_0,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)"
            ),label = list(age ~ "Age",
                           sex ~ "Sex (Female=1)",
                           race___4 ~ "Race(White=1)",
                           hf_cat ~ "HF Category (HFpEF=1) (Study)",
                           egfr ~ "EGFR (Study)",
                           bl_sppb ~ "Baseline SPPB",
                           bl_smw ~ "Baseline 6MWD"
            ),
            missing_text = "(Missing)"
) %>% add_p(list(all_continuous() ~ "wilcox.test",all_categorical() ~ "fisher.test")) %>% add_n(
) %>% modify_header(label ~ "**Variable**") %>% 
  modify_caption("**Table 1. Baseline Patient Characteristics**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Control (0) vs. Intervention (1)**") %>%
  as_gt() %>% 
  gtsave(filename = "Results/Table1_Baselinecharacteristics.html")

############################################
## 4. Supplemental Figure 2 - Spearman correlation of baseline proteins ##
############################################
scmat <- round(cor(rprotblds[,which(colnames(rprotblds) %in% protlist)],
    method="spearman",use="pairwise.complete.obs"),2)
scmat[lower.tri(scmat)] <- NA
#rownames(scmat) <- colnames(scmat) <- protlist
scmat_melt <- melt(scmat, na.rm = TRUE)

#compute p-values
for (i in 1:nrow(scmat_melt)){
    b1 <- as.character(scmat_melt$Var1[i])
    b2 <- as.character(scmat_melt$Var2[i])
    scmat_melt$pval[i] <- cor.test(rprotblds[,which(colnames(rprotblds) %in% b1)],
                                   rprotblds[,which(colnames(rprotblds) %in% b2)],
                    method="spearman",use="pairwise.complete.obs")$p.value
}
scmat_melt$sig = scmat_melt$pval<0.01
scmat_melt$text = ifelse(scmat_melt$sig==TRUE & scmat_melt$value != 1,paste0(scmat_melt$value,"*"),scmat_melt$value)
ggheatmap <- ggplot(scmat_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = text), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
# Print the heatmap
pdf("Results/SuppFigure2_SpearmanCorrelation.pdf",width=5,height=5)
print(ggheatmap)
dev.off()

ggsave(ggheatmap,file='Results/SuppFigure2_SpearmanCorrelation.eps',width=5,height=5,device="eps")

############################################
## 5. Table 2 (Linear Regression Results) ##
############################################
#tabmodel won't print out interaction effects for variables with special characters in name, so rename
#columns to remove special characters

tmpdf <- rprotblds
colnames(tmpdf)[which(colnames(rprotblds) %in% protlist)] =  make.names(protlist, unique=TRUE)
protlistclean= make.names(protlist, unique=TRUE)
regresults.sppb=data.frame(Protein=protlist,Unadjusted.Effect=NA,Unadjusted.Pval=NA,
                      Adjusted.Effect=NA,Adjusted.Pval=NA,Interaction.Pval=NA)
regresults.6mwd=regresults.sppb
for (i in 1:length(protlist)){
  
  fname=paste("Results/RegressionAnalyses/SPPBvs",protlist[i],"_regressionanalysis.html",sep="")
  ttl=paste("Follow up SPPB ~ Baseline SPPB + ",protlist[i]," + Intervention + More",sep="")
  mod.sppb.unadj=lm(data=tmpdf,
          formula=paste0("fu_sppb ~ bl_sppb + 
                               ",protlistclean[i]," +
                               intervention_1_control_0"));
  mod.sppb.adj=lm(data=tmpdf,
                  formula=paste0("fu_sppb ~ bl_sppb + 
                               ",protlistclean[i],"*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"));
  print(tab_model(summary(mod.sppb.unadj),
                  summary(mod.sppb.adj),
                  title=ttl,file=fname,
                  digits=3,digits.p=4,show.obs=TRUE))
  
  #store results
  regresults.sppb$Unadjusted.Effect[i]=paste0(round(mod.sppb.unadj$coefficients[protlistclean[i]],4),
                                         "(",
                                         round(confint(mod.sppb.unadj,protlistclean[i])[1],4),
                                         ",",
                                         round(confint(mod.sppb.unadj,protlistclean[i])[2],4),
                                         ")")
  regresults.sppb$Unadjusted.Pval[i]=round(summary(mod.sppb.unadj)$coefficients[protlistclean[i],4],4)
  regresults.sppb$Adjusted.Effect[i]=paste0(round(mod.sppb.adj$coefficients[protlistclean[i]],4),
                                         "(",
                                         round(confint(mod.sppb.adj,protlistclean[i])[1],4),
                                         ",",
                                         round(confint(mod.sppb.adj,protlistclean[i])[2],4),
                                         ")")
  regresults.sppb$Adjusted.Pval[i]=round(summary(mod.sppb.adj)$coefficients[protlistclean[i],4],4)
  regresults.sppb$Interaction.Pval[i]=round(summary(mod.sppb.adj)$coefficients[paste0(protlistclean[i],":intervention_1_control_01"),4],4)
   
  fname=paste("Results/RegressionAnalyses/6MWDvs",protlist[i],"_regressionanalysis.html",sep="")
  ttl=paste("Follow up 6MWD ~ Baseline 6MWD + ",protlist[i]," + Intervention + More",sep="")
  mod.6mwd.unadj=lm(data=tmpdf,
                    formula=paste0("fu_smw ~ bl_smw +  
                               ",protlistclean[i]," +
                               intervention_1_control_0"));
  mod.6mwd.adj=lm(data=tmpdf,
                  formula=paste0("fu_smw ~ bl_smw + 
                               ",protlistclean[i],"*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"));
  print(tab_model(summary(mod.6mwd.unadj),
                  summary(mod.6mwd.adj),
                  title=ttl,file=fname,
                  digits=3,digits.p=4,show.obs=TRUE))

  #store results
  regresults.6mwd$Unadjusted.Effect[i]=paste0(round(mod.6mwd.unadj$coefficients[protlistclean[i]],4),
                                              "(",
                                              round(confint(mod.6mwd.unadj,protlistclean[i])[1],4),
                                              ",",
                                              round(confint(mod.6mwd.unadj,protlistclean[i])[2],4),
                                              ")")
  regresults.6mwd$Unadjusted.Pval[i]=round(summary(mod.6mwd.unadj)$coefficients[protlistclean[i],4],4)
  regresults.6mwd$Adjusted.Effect[i]=paste0(round(mod.6mwd.adj$coefficients[protlistclean[i]],4),
                                            "(",
                                            round(confint(mod.6mwd.adj,protlistclean[i])[1],4),
                                            ",",
                                            round(confint(mod.6mwd.adj,protlistclean[i])[2],4),
                                            ")")
  regresults.6mwd$Adjusted.Pval[i]=round(summary(mod.6mwd.adj)$coefficients[protlistclean[i],4],4)
  regresults.6mwd$Interaction.Pval[i]=round(summary(mod.6mwd.adj)$coefficients[paste0(protlistclean[i],":intervention_1_control_01"),4],4)
  
}  
write.csv(regresults.sppb,"Results/regsummary_sppb.csv")
write.csv(regresults.6mwd,"Results/regsummary_6mwd.csv")
############################################
## 6. Figures 1 and 2 (Linear Regression Results) ##
############################################

#create long form dataset
rprotblds_long=pivot_longer(data=rprotblds[,
                                       c("subject_id","timepoint","intervention_1_control_0",
                                          "bl_smw","fu_smw","bl_sppb","fu_sppb",protlist)],
                             cols=protlist,
                          names_to='Protein',values_to='ExpressionLevel')
#rprotblds_long$intervention_1_control_0 <- factor(rcombds_long_bl$intervention_1_control_0,levels=c(0,1))
rprotblds_long$sppb_chg = rprotblds_long$fu_sppb - rprotblds_long$bl_sppb
rprotblds_long$smw_chg = rprotblds_long$fu_smw - rprotblds_long$bl_smw
colnames(rprotblds_long)[3] <- "Intervention"
rprotblds_long$Protein = factor(rprotblds_long$Protein,levels=protlist)

protlist.6mwd <- c("LDL receptor","ALCAM","PLC","CHIT1","Gal-4")
protlist.sppb <- c("LDL receptor","ALCAM","GP6","CCL16","ST2")

#SPPB

fig1 <- ggplot(data=rprotblds_long[rprotblds_long$Protein %in% protlist.sppb,],
               mapping=aes(x=ExpressionLevel,y=sppb_chg,color=Intervention))+
  geom_smooth(method='lm')+
  geom_hline(yintercept=0,color='black',linetype='dashed')+
  facet_wrap('Protein',scales='free',ncol=3)+
  labs(x="Expression Level",y="Change in SPPB Score")+
  scale_color_discrete(breaks=c(1,0),
                       labels=c("Rehabilitation\nIntervention","Attention\nControl"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85,0.15),
        legend.box.background = element_rect(color = "black"))

#6MWD
fig2 <- ggplot(data=rprotblds_long[rprotblds_long$Protein %in% protlist.6mwd,],
               mapping=aes(x=ExpressionLevel,y=smw_chg,color=Intervention))+
  geom_smooth(method='lm')+
  geom_hline(yintercept=0,color='black',linetype='dashed')+
  facet_wrap('Protein',scales='free',ncol=3)+
  labs(x="Expression Level",y="Change in 6MWD (meters)")+
  scale_color_discrete(breaks=c(1,0),
                       labels=c("Rehabilitation\nIntervention","Attention\nControl"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85,0.15),
        legend.box.background = element_rect(color = "black"))


# Print figures
pdf("Results/Figure1_ChangeinSPPB.pdf",width=10,height=6)
print(fig1)
dev.off()
ggsave(fig1,file='Results/Figure1_ChangeinSPPB.eps',width=6,height=3,device="eps")

pdf("Results/Figure2_Changein6MWD.pdf",width=10,height=6)
print(fig2)
dev.off()
ggsave(fig2,file='Results/Figure2_Changein6MWD.eps',width=6,height=3,device="eps")

############################################
## 7. Supplemental Table 1 compare those with biomarkers vs. unmeasured ## - DO WE NEED THIS BLOCK?  I DON'T THINK SO
############################################
# #supplemental table 1: comparing those with baseline biomarkers measured vs. no baseline biomarkers
# rhfds$bm_measured = (rhfds$study_id %in% 
#                        rcombds$subject_id[rcombds$timepoint=='Baseline'])
# rhfds$bm_measured = factor(rhfds$bm_measured,levels=c("TRUE","FALSE"))
# tbl_summary(rhfds[,c("bm_measured","age","sex","race___4","hf_cat","egfr","bl_sppb","bl_smw")],
#             by=bm_measured,
#             type = all_continuous() ~ "continuous2",
#             statistic = list(all_continuous() ~ c("{mean} ({sd})",
#                                                   "{median} ({p25}, {p75})", 
#                                                   "{min}, {max}"),
#                              all_categorical() ~ "{n} / {N} ({p}%)"
#             ),label = list(age ~ "Age",
#                            sex ~ "Sex (Female=1)",
#                            race___4 ~ "Race(White=1)",
#                            hf_cat ~ "HF Category (HFpEF=1) (Study)",
#                            egfr ~ "EGFR (Study)",
#                            bl_sppb ~ "SPPB Score (Baseline)",
#                            bl_smw ~ "6MWD (meters) (Baseline)"
#             ),
#             missing_text = "(Missing)"
# ) %>% add_p(list(all_continuous() ~ "wilcox.test",all_categorical() ~ "fisher.test")) %>% add_n(
# ) %>% modify_header(label ~ "**Variable**") %>% 
#   modify_caption("**Supplemental Table 1. Baseline Patient Characteristics**") %>%
#   modify_spanning_header(c("stat_1", "stat_2") ~ "**Measured vs. Unmeasured**") %>%
#   as_gt() %>% 
#   gtsave(filename = "SuppTable1_measuredvsunmeasured.html")

############################################
## 8. Supplemental Table 3 - secondary outcomes  ##
############################################

tbl_summary(rprotblds[,
                    c("intervention_1_control_0","event","rehosp","death","lfu")],
            by=intervention_1_control_0,
            type = list(all_continuous() ~ "continuous2"
            ),
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)"
            ),label = list(rehosp ~ "Rehospitalization",
                           death ~ "Death",
                           event ~ "Rehospitalization or Death",
                           lfu ~ "Lost to follow up"
            ),
            missing_text = "(Missing)"
) %>% add_p(list(all_continuous() ~ "wilcox.test",all_categorical() ~ "fisher.test")) %>% add_n(
) %>% modify_header(label ~ "**Variable**") %>% 
  modify_caption("**Supplemental Table 2**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Control (0) vs. Intervention (1)**") %>%
  as_gt() %>% 
  gtsave(filename = "Results/SuppTable3_Secondaryoutcomes.html")

############################################
## 9. Supplemental Table 4 (Regression Results - Secondary Outcomes) ##
############################################
#tabmodel won't print out interaction effects for variables with special characters in name, so rename
#columns to remove special characters

tmpdf <- rprotblds
colnames(tmpdf)[which(colnames(rprotblds) %in% protlist)] =  make.names(protlist, unique=TRUE)
protlistclean= make.names(protlist, unique=TRUE)
regresults.death=data.frame(Protein=protlist,Unadjusted.Effect=NA,Unadjusted.Pval=NA,
                           Adjusted.Effect=NA,Adjusted.Pval=NA,Interaction.Pval=NA)
regresults.rehosp=regresults.death
regresults.both=regresults.death

for (i in 1:length(protlist)){
  
  fname=paste("Results/RegressionAnalyses/DeathorRehospvs",protlist[i],"_regressionanalysis.html",sep="")
  ttl=paste("logit(Death or Rehosp) ~ ",protlist[i]," + Intervention + More",sep="")
  mod.both.unadj=glm(data=tmpdf,
                     formula=paste0("event ~  
                               ",protlistclean[i]," +
                               intervention_1_control_0"), family=binomial(link="logit"))
  mod.both.adj=glm(data=tmpdf,
                   formula=paste0("event ~  
                               ",protlistclean[i],"*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"), family=binomial(link="logit"))
  print(tab_model(mod.both.unadj,
                  mod.both.adj,
                  title=ttl,
                  file=fname,
                  digits=4,digits.p=4,show.obs=TRUE))
  
  #store results
  regresults.both$Unadjusted.Effect[i]=paste0(round(exp(mod.both.unadj$coefficients[protlistclean[i]]),4),
                                              "(",
                                              round(exp(confint(mod.both.unadj,protlistclean[i])[1]),4),
                                              ",",
                                              round(exp(confint(mod.both.unadj,protlistclean[i])[2]),4),
                                              ")")
  regresults.both$Unadjusted.Pval[i]=round(summary(mod.both.unadj)$coefficients[protlistclean[i],4],4)
  regresults.both$Adjusted.Effect[i]=paste0(round(exp(mod.both.adj$coefficients[protlistclean[i]]),4),
                                            "(",
                                            round(exp(confint(mod.both.adj,protlistclean[i])[1]),4),
                                            ",",
                                            round(exp(confint(mod.both.adj,protlistclean[i])[2]),4),
                                            ")")
  regresults.both$Adjusted.Pval[i]=round(summary(mod.both.adj)$coefficients[protlistclean[i],4],4)
  regresults.both$Interaction.Pval[i]=round(summary(mod.both.adj)$coefficients[paste0(protlistclean[i],":intervention_1_control_01"),4],4)
  
  
  fname=paste("Results/RegressionAnalyses/Deathvs",protlist[i],"_regressionanalysis.html",sep="")
  ttl=paste("logit(Death) ~ ",protlist[i]," + Intervention + More",sep="")

  mod.death.unadj=glm(data=tmpdf,
                      formula=paste0("death ~  
                               ",protlistclean[i]," +
                               intervention_1_control_0"), family=binomial(link="logit"))
  mod.death.adj=glm(data=tmpdf,
                    formula=paste0("death ~  
                               ",protlistclean[i],"*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"), family=binomial(link="logit"))
  print(tab_model(mod.death.unadj,
                  mod.death.adj,
                  title=ttl,
                  file=fname,
                  digits=4,digits.p=4,show.obs=TRUE))
  
  #store results
  regresults.death$Unadjusted.Effect[i]=paste0(round(exp(mod.death.unadj$coefficients[protlistclean[i]]),4),
                                              "(",
                                              round(exp(confint(mod.death.unadj,protlistclean[i])[1]),4),
                                              ",",
                                              round(exp(confint(mod.death.unadj,protlistclean[i])[2]),4),
                                              ")")
  regresults.death$Unadjusted.Pval[i]=round(summary(mod.death.unadj)$coefficients[protlistclean[i],4],4)
  regresults.death$Adjusted.Effect[i]=paste0(round(exp(mod.death.adj$coefficients[protlistclean[i]]),4),
                                            "(",
                                            round(exp(confint(mod.death.adj,protlistclean[i])[1]),4),
                                            ",",
                                            round(exp(confint(mod.death.adj,protlistclean[i])[2]),4),
                                            ")")
  regresults.death$Adjusted.Pval[i]=round(summary(mod.death.adj)$coefficients[protlistclean[i],4],4)
  regresults.death$Interaction.Pval[i]=round(summary(mod.death.adj)$coefficients[paste0(protlistclean[i],":intervention_1_control_01"),4],4)
  
  
  fname=paste("Results/RegressionAnalyses/Rehospvs",protlist[i],"_regressionanalysis.html",sep="")
  ttl=paste("logit(Rehosp) ~ ",protlist[i]," + Intervention + More",sep="")
  mod.rehosp.unadj=glm(data=tmpdf,
                       formula=paste0("rehosp ~  
                               ",protlistclean[i]," +
                               intervention_1_control_0"), family=binomial(link="logit"))
  mod.rehosp.adj=glm(data=tmpdf,
                     formula=paste0("rehosp ~  
                               ",protlistclean[i],"*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"), family=binomial(link="logit"))
  print(tab_model(mod.rehosp.unadj,
                  mod.rehosp.adj,
                  title=ttl,
                  file=fname,
                  digits=4,digits.p=4,show.obs=TRUE))
  
  #store results
  regresults.rehosp$Unadjusted.Effect[i]=paste0(round(exp(mod.rehosp.unadj$coefficients[protlistclean[i]]),4),
                                               "(",
                                               round(exp(confint(mod.rehosp.unadj,protlistclean[i])[1]),4),
                                               ",",
                                               round(exp(confint(mod.rehosp.unadj,protlistclean[i])[2]),4),
                                               ")")
  regresults.rehosp$Unadjusted.Pval[i]=round(summary(mod.rehosp.unadj)$coefficients[protlistclean[i],4],4)
  regresults.rehosp$Adjusted.Effect[i]=paste0(round(exp(mod.rehosp.adj$coefficients[protlistclean[i]]),4),
                                             "(",
                                             round(exp(confint(mod.rehosp.adj,protlistclean[i])[1]),4),
                                             ",",
                                             round(exp(confint(mod.rehosp.adj,protlistclean[i])[2]),4),
                                             ")")
  regresults.rehosp$Adjusted.Pval[i]=round(summary(mod.rehosp.adj)$coefficients[protlistclean[i],4],4)
  regresults.rehosp$Interaction.Pval[i]=round(summary(mod.rehosp.adj)$coefficients[paste0(protlistclean[i],":intervention_1_control_01"),4],4)
  
  
}  
write.csv(regresults.death,"Results/regsummary_death.csv")
write.csv(regresults.rehosp,"Results/regsummary_rehosp.csv")
write.csv(regresults.both,"Results/regsummary_deathrehosp.csv")

############################################
## 10. Propensity matching  ##
############################################

ps.df=rprotblds
ps.df$trt=as.numeric(ps.df$intervention_1_control_0)-1
ps.df=ps.df[complete.cases(ps.df[,c("LDL receptor","ALCAM","PLC","CHIT1","Gal-4",
                                    "GP6","CCL16","ST2",
                                    'bl_smw','fu_smw')]),]

#propensity score for treatment assignment
ppty=glm(trt~age+sex+race___4+hf_cat
         ,family=binomial(link="logit"),data=ps.df)
prop.score=predict.glm(ppty,type="response",na.action=na.exclude)
ps.df=as.data.frame(cbind(ps.df,prop.score))
ppty.distance=match_on(ppty)

######################
## 11. Matching Tree Functions
######################

#tree pruning and selection functions from Matching Tree paper (Zhang et al., 2021)
#https://doi.org/10.1016/j.csda.2021.107188

#MT (m)------
#prune tree function
parent <- function(x) {
  if (x[1] != 1)
    c(Recall(if (x %% 2 == 0L) x / 2 else (x - 1) / 2), x) else x
}

jesse.tree.select=function(input,set1.caliper.tree){
  fit.overlap=TRUE
  ptree1=input
  while (fit.overlap==TRUE){
    if (dim(ptree1$frame)[1]==1) break
    rownames=as.numeric(rownames(ptree1$frame))
    size=dim(ptree1$frame)[1]
    rownamesmax=rownames[which.max(rownames)]
    rownamessecond=rownames[which(rownames == 
                                    sort(unique(rownames),partial=size-1)[size-1])]
    
    subset1=set1.caliper.tree[ptree1$where==which(rownames==rownamesmax),]
    subset2=set1.caliper.tree[ptree1$where==which(rownames==rownamessecond),]
    dummy=c(rep("1",dim(subset1)[1]),rep("2",dim(subset2)[1]))
    subsetf=cbind(rbind(subset1,subset2),dummy)
    
    subsetf$dummy=as.factor(subsetf$dummy)
    subsetf$trt=as.factor(subsetf$trt)
    ls.lm=lm(outcome~dummy*trt,data=subsetf)
    ls.list=lsmeans(ls.lm,list(pairwise~dummy|trt,pairwise~trt|dummy))
    out1=as.matrix(summary(ls.list[[4]]))
    
    q.value=1.96 #could change q.value for different methods
    CI1=c(as.numeric(out1[1,3])-q.value*as.numeric(out1[1,4]),
          as.numeric(out1[1,3])+q.value*as.numeric(out1[1,4]))
    CI2=c(as.numeric(out1[2,3])-q.value*as.numeric(out1[2,4]),
          as.numeric(out1[2,3])+q.value*as.numeric(out1[2,4]))  
    fit.overlap=max(CI1[1],CI2[1])<min(CI1[2],CI2[2])
    fit.overlap[is.na(fit.overlap)] <- FALSE
    if (fit.overlap==FALSE) break
    ptree1=snip.rpart(ptree1,
                      toss=tail(head(parent(rownamesmax), -1),n=1))}
  if (dim(ptree1$frame)[1]>1){
    return(ptree1)
  } else {
    return (input)
  }
}

######################
## 12. 6MWD Matching Tree
######################

#match on 6MWD vars
# [1] "LDL receptor"   "ALCAM"          "PLC"         "CHIT1"     "Gal-4"         

set.trt=ps.df[ps.df[, "trt"] == 1,]
set.cont=ps.df[ps.df[,"trt"]==0,]
match.set=c()
calp.sd=1
for (i in 1:dim(set.trt)[1]){
  pair.set1=c()
  var.select1=ifelse(abs(set.cont$`LDL receptor`-
                           set.trt$`LDL receptor`[i])<
                       calp.sd*sd(set.trt$`LDL receptor`),1,0)
  var.select2=ifelse(abs(set.cont$ALCAM-
                           set.trt$ALCAM[i])<
                       calp.sd*sd(set.trt$ALCAM),1,0)
  var.select3=ifelse(abs(set.cont$PLC-
                           set.trt$PLC[i])<
                       calp.sd*sd(set.trt$PLC),1,0)
  var.select4=ifelse(abs(set.cont$CHIT1-
                           set.trt$CHIT1[i])<
                       calp.sd*sd(set.trt$CHIT1),1,0)
  var.select5=ifelse(abs(set.cont$`Gal-4`-
                           set.trt$`Gal-4`[i])<
                       calp.sd*sd(set.trt$`Gal-4`),1,0)
  var.selectf=ifelse(var.select1+var.select2+var.select3+
                       var.select4+var.select5==5,1,0)
  id=c(1:length(var.select1))
  lenset.cont2=as.data.frame(cbind(
    id,set.cont,var.selectf))
  lenset.cont3=lenset.cont2[lenset.cont2[,"var.selectf"]==1,]
  if (dim(lenset.cont3)[1]>0){
    cont.index=which.min(abs(lenset.cont3$prop.score-set.trt$prop.score[i]))
    pair.set=rbind(set.trt[i,],
                   lenset.cont3[cont.index,c(-1,-dim(lenset.cont3)[2])])
    pairin=rep(i,2)
    pair.set1=cbind(pair.set,pairin)
    set.cont = set.cont[!set.cont$prop.score == 
                          lenset.cont3$prop.score[cont.index],]
  } else {
    pair.set1=c()
    set.cont = set.cont
  }
  match.set=rbind(match.set,pair.set1)
}

#calculate the difference between treated and control
treat.set=match.set[match.set[,"trt"]==1,]
control.set=match.set[match.set[,"trt"]==0,]
outcome.diff=(treat.set$fu_smw-treat.set$bl_smw)-
  (control.set$fu_smw-control.set$bl_smw)
set.total=as.data.frame(rbind(cbind(treat.set,outcome.diff),cbind(control.set,outcome.diff)))
matchedset.6mwd=set.total #for later use
#build matching tree using all data
tree.rpart=rpart(outcome.diff~`LDL receptor`+ALCAM+PLC+CHIT1+`Gal-4`
                 # +age+sex+race___4+hf_cat
                 ,method="anova",data=set.total,
                 control=rpart.control(minbucket = 20))
#this is the result after prune
set.total$outcome = set.total$fu_smw - set.total$bl_smw
tree.all=jesse.tree.select(tree.rpart,set.total)

pdf("Results/Figure3b_Matchingtree_6MWD.pdf")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()
setEPS()
postscript("Results/Figure3b_Matchingtree_6MWD.eps")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()

#get confidence intervals
set.total$smw_node=rpart.predict.leaves(tree.all,set.total,type="where")
modscore.smw=lm(data=set.total,formula= outcome.diff~0+as.factor(smw_node))
confint(modscore.smw)
# 2.5 %   97.5 %
#   as.factor(smw_node)3 -53.26891 13.49727
# as.factor(smw_node)4 -11.20765 51.84635
# as.factor(smw_node)5  48.75768 99.11928

summary(modscore.smw)
# Call:
#   lm(formula = outcome.diff ~ 0 + as.factor(smw_node), data = set.total)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -322.03  -53.61   14.66   59.63  223.86 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# as.factor(smw_node)3   -19.89      16.87  -1.179    0.241    
# as.factor(smw_node)4    20.32      15.93   1.276    0.204    
# as.factor(smw_node)5    73.94      12.72   5.811 4.85e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 96.9 on 125 degrees of freedom
# Multiple R-squared:  0.2274,	Adjusted R-squared:  0.2088 
# F-statistic: 12.26 on 3 and 125 DF,  p-value: 4.368e-07

#Leave one out stability analysis
set.seed(324)
folds <- cut(seq(1,nrow(treat.set)),breaks=nrow(treat.set),labels=FALSE)
#shuffle folds
folds<- folds[sample(nrow(treat.set))]

#store trees and subpops
tree.2=list()
tree.plot=list()
tree.leaves=list()
tree.mod.train=list()
tree.mod.test=list()
for (i in 1:max(folds)){
  
  #identify pairs in the train
  ids <- treat.set$subject_id[folds != i]
  pairs.train <- set.total$pairin[set.total$subject_id %in% ids]
  
  #this is the result from CART (update with different rpart.control params)
  tmp=rpart(outcome.diff~`LDL receptor`+ALCAM+PLC+CHIT1+`Gal-4`
            # +age+sex+race___4+hf_cat
            ,method="anova",data=set.total[set.total$pairin %in% pairs.train,],
            control=rpart.control(minbucket = 20))
  
  #this is the result after prune
  set.total$outcome = set.total$fu_smw - set.total$bl_smw
  tree.2[[i]]=jesse.tree.select(tmp,set.total)
  tree.plot[[i]]=rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
  
  #numbers are average outcome diff
  #get leaves for all observations
  tree.leaves[[i]]=rpart.predict.leaves(tree.2[[i]],set.total,type="where")
  
  df=cbind(set.total,tree.leaves[[i]])
  colnames(df)
  tree.mod.train[[i]]=summary(lm(outcome.diff ~ 0 + as.factor(`tree.leaves[[i]]`), data=df[df$pairin %in% pairs.train,]))
  
}

pdf("Results/Matchingtree_6MWD_LOO.pdf", onefile = TRUE)
for (i in 1:length(tree.2)) {
  rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
}
dev.off()

######################
## 13. SPPB Matching Tree
######################

#match on sppb vars
# [1] "LDL receptor" "ALCAM"        "GP6"          "CCL16"   "ST2"     

set.trt=ps.df[ps.df[, "trt"] == 1,]
set.cont=ps.df[ps.df[,"trt"]==0,]
match.set=c()
calp.sd=1
for (i in 1:dim(set.trt)[1]){
  pair.set1=c()
  var.select1=ifelse(abs(set.cont$`LDL receptor`-
                           set.trt$`LDL receptor`[i])<
                       calp.sd*sd(set.trt$`LDL receptor`),1,0)
  var.select2=ifelse(abs(set.cont$ALCAM-
                           set.trt$ALCAM[i])<
                       calp.sd*sd(set.trt$ALCAM),1,0)
  var.select3=ifelse(abs(set.cont$GP6-
                           set.trt$GP6[i])<
                       calp.sd*sd(set.trt$GP6),1,0)
  var.select4=ifelse(abs(set.cont$CCL16-
                           set.trt$CCL16[i])<
                       calp.sd*sd(set.trt$CCL16),1,0)
  var.select5=ifelse(abs(set.cont$`ST2`-
                           set.trt$`ST2`[i])<
                       calp.sd*sd(set.trt$`ST2`),1,0)
  var.selectf=ifelse(var.select1+var.select2+var.select3+
                       var.select4+var.select5==5,1,0)
  id=c(1:length(var.select1))
  lenset.cont2=as.data.frame(cbind(
    id,set.cont,var.selectf))
  lenset.cont3=lenset.cont2[lenset.cont2[,"var.selectf"]==1,]
  if (dim(lenset.cont3)[1]>0){
    cont.index=which.min(abs(lenset.cont3$prop.score-set.trt$prop.score[i]))
    pair.set=rbind(set.trt[i,],
                   lenset.cont3[cont.index,c(-1,-dim(lenset.cont3)[2])])
    pairin=rep(i,2)
    pair.set1=cbind(pair.set,pairin)
    set.cont = set.cont[!set.cont$prop.score == 
                          lenset.cont3$prop.score[cont.index],]
  } else {
    pair.set1=c()
    set.cont = set.cont
  }
  match.set=rbind(match.set,pair.set1)
}

#calculate the difference between treated and control
treat.set=match.set[match.set[,"trt"]==1,]
control.set=match.set[match.set[,"trt"]==0,]
outcome.diff=(treat.set$fu_sppb-treat.set$bl_sppb)-
  (control.set$fu_sppb-control.set$bl_sppb)
set.total=as.data.frame(rbind(cbind(treat.set,outcome.diff),cbind(control.set,outcome.diff)))
matchedset.sppb=set.total #for later use
#build matching tree using all data
tree.rpart=rpart(outcome.diff~`LDL receptor`+ALCAM+GP6+CCL16+ST2
                 # +age+sex+race___4+hf_cat
                 ,method="anova",data=set.total,
                 control=rpart.control(minbucket = 20))
#this is the result after prune
set.total$outcome = set.total$fu_sppb - set.total$bl_sppb
tree.all=jesse.tree.select(tree.rpart,set.total)

pdf("Results/Figure3a_Matchingtree_SPPB.pdf")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()
setEPS()
postscript("Results/Figure3a_Matchingtree_SPPB.eps")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()

#get confidence intervals
set.total$sppb_node=rpart.predict.leaves(tree.all,set.total,type="where")
modscore.sppb=lm(data=set.total,formula= outcome.diff~0+as.factor(sppb_node))
confint(modscore.sppb)
# 2.5 %    97.5 %
#   as.factor(sppb_node)2 -1.475802 0.5424683
# as.factor(sppb_node)3  1.572931 2.7381799

summary(modscore.sppb)
# Call:
#   lm(formula = outcome.diff ~ 0 + as.factor(sppb_node), data = set.total)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.1556 -2.1556 -0.1556  1.5611  6.8444 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# as.factor(sppb_node)2  -0.4667     0.5096  -0.916    0.362    
# as.factor(sppb_node)3   2.1556     0.2942   7.326 3.18e-11 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.791 on 118 degrees of freedom
# Multiple R-squared:  0.316,	Adjusted R-squared:  0.3044 
# F-statistic: 27.26 on 2 and 118 DF,  p-value: 1.854e-10

#Leave one out stability analysis
set.seed(67)
folds <- cut(seq(1,nrow(treat.set)),breaks=nrow(treat.set),labels=FALSE)
#shuffle folds
folds<- folds[sample(nrow(treat.set))]

#store trees and subpops
tree.2=list()
tree.plot=list()
tree.leaves=list()
tree.mod.train=list()
tree.mod.test=list()
for (i in 1:max(folds)){
  
  #identify pairs in the train
  ids <- treat.set$subject_id[folds != i]
  pairs.train <- set.total$pairin[set.total$subject_id %in% ids]
  
  #this is the result from CART (update with different rpart.control params)
  tmp=rpart(outcome.diff~`LDL receptor`+ALCAM+GP6+CCL16+ST2
            # +age+sex+race___4+hf_cat
                   ,method="anova",data=set.total[set.total$pairin %in% pairs.train,],
                   control=rpart.control(minbucket = 20))

  #this is the result after prune
  set.total$outcome = set.total$fu_sppb - set.total$bl_sppb
  tree.2[[i]]=jesse.tree.select(tmp,set.total)
  tree.plot[[i]]=rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
  
  #numbers are average outcome diff
  #get leaves for all observations
  tree.leaves[[i]]=rpart.predict.leaves(tree.2[[i]],set.total,type="where")
  
  df=cbind(set.total,tree.leaves[[i]])
  colnames(df)
  tree.mod.train[[i]]=summary(lm(outcome.diff ~ 0 + as.factor(`tree.leaves[[i]]`), data=df[df$pairin %in% pairs.train,]))

}

pdf("Results/Matchingtree_SPPB_LOO.pdf", onefile = TRUE)
for (i in 1:length(tree.2)) {
  rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
}
dev.off()

######################
## 15. Supplemental table 2 matched vs. unmatched
######################
varlist=c("intervention_1_control_0","age","sex","race___4",
          "hf_cat","egfr","troponin_t_ng_l","troponin_i_pg_ml","nt_pro_bnp","hs_crp_mg_l",
          "creatinine_mg_dl")
tmpdf = rbind(data.frame(set.total[,varlist],dataset="matched"),
              data.frame(setdiff(ps.df[,varlist],set.total[,varlist]),dataset="unmatched"))
tmpdf$dataset=factor(tmpdf$dataset,levels=c("unmatched","matched"))
tmpdf$grp=str_c(tmpdf$intervention_1_control_0,tmpdf$dataset)
tmpdf$grp=as.factor(tmpdf$grp)
tmpdf$grp=fct_collapse(tmpdf$grp, Unmatched = c("0unmatched","1unmatched"), MatchedControl = "0matched",MatchedIntervention = "1matched")
tmpdf$grp = factor(tmpdf$grp,levels = c("MatchedControl","MatchedIntervention","Unmatched"))
tmpdf=tmpdf[,-which(names(tmpdf) %in% c("intervention_1_control_0","dataset"))]

tbl_summary(tmpdf,
            by=grp,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)"
            ),label = list(age ~ "Age",
                           sex ~ "Sex (Female=1)",
                           race___4 ~ "Race(White=1)",
                           hf_cat ~ "HF Category (HFpEF=1) (Study)",
                           egfr ~ "EGFR (Study)",
                           troponin_t_ng_l ~ "Hs-cTnT (Inova)",
                           troponin_i_pg_ml ~ "Hs-cTnI (Inova)",
                           nt_pro_bnp ~ "NT-proBNP (Inova)",
                           hs_crp_mg_l ~ "Hs-CRP (Inova)",
                           creatinine_mg_dl ~ "Creatinine (Inova)"
            ),
            missing_text = "(Missing)"
)  %>% modify_header(label ~ "**Variable**") %>% 
  modify_caption("**Supplemental Table 2. Baseline Measurements (Unmatched vs. Matched)**") %>%
  as_gt() %>% 
  gtsave(filename = "SuppTable2_matchedvsunmatched.html")

######################
## 16. Supplemental figure 3 - protein change (baseline to follow up) by treatment
######################
#create long form dataset
rprotfuds_long=pivot_longer(data=rprotfuds[,
                                           c("subject_id","timepoint","intervention_1_control_0",
                                             "bl_smw","fu_smw","bl_sppb","fu_sppb",protlist)],
                            cols=protlist,
                            names_to='Protein',values_to='ExpressionLevel')
rprotfuds_long$sppb_chg = rprotfuds_long$fu_sppb - rprotfuds_long$bl_sppb
rprotfuds_long$smw_chg = rprotfuds_long$fu_smw - rprotfuds_long$bl_smw
colnames(rprotfuds_long)[3] <- "Intervention"
rprotfuds_long$Protein = factor(rprotfuds_long$Protein,levels=protlist)


rprotds_wide=merge(rprotblds_long,
                   rprotfuds_long,by=c('subject_id','Protein','Intervention',"bl_smw","fu_smw","bl_sppb","fu_sppb","sppb_chg","smw_chg"),
                   all.x=TRUE)
rprotds_wide$Protein_change = (rprotds_wide$ExpressionLevel.y) - (rprotds_wide$ExpressionLevel.x)

supfig3 <- ggplot(data=rprotds_wide,
             mapping=aes(x=Intervention, y=ExpressionLevel.y/ExpressionLevel.x, fill=Intervention))+
        geom_boxplot()+facet_wrap_paginate('Protein',scales='free',ncol=3,
                                          page=1)+hw+
        labs(y="Ratio of change from baseline to 12-weeks")+
        geom_hline(yintercept=1,linetype=2)+
        scale_y_continuous(trans='log2',labels = function(x) sprintf("%.3f", x),
                           breaks=seq(0.5,1.5,length.out=6))+
        scale_x_discrete(breaks=c(1,0),
                         labels=c("RI","AC"))+
        scale_fill_discrete(breaks=c(0,1),
                            labels=c("AC","RI"))+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(face="bold"),
              legend.title = element_blank(),
              legend.position = c(0.85,0.15),
              legend.box.background = element_rect(color = "black"))


pdf("Results/SuppFigure3_Proteinchange.pdf",height=9,width=7,onefile = TRUE)
supfig3
dev.off()
ggsave(supfig3,file='Results/SuppFigure3_Proteinchange.eps',height=9,width=7,device="eps")


##################### protein change from baseline to follow up overall then by treatment group
dftmp=rprotds_wide[!is.na(rprotds_wide$timepoint.y),]
for (i in 1:length(protlist)){
  print(protlist[i])
  print(paste(round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                       probs=c(0.5)),3),
              "(",
              round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                       probs=c(0.25)),3),
              ",",
              round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                       probs=c(0.75)),3),
              ")",sep=""))
  
  print(wilcox.test(x=dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                    dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"])$p.value)
}

#attention control
dftmp=rprotds_wide[!is.na(rprotds_wide$timepoint.y) & rprotds_wide$Intervention==0,]
for (i in 1:length(protlist)){
  print(protlist[i])
  print(paste(round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                             probs=c(0.5)),3),
              "(",
              round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                             probs=c(0.25)),3),
              ",",
              round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                             probs=c(0.75)),3),
              ")",sep=""))
  
  print(wilcox.test(x=dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                      dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"])$p.value)
}


#rehab intervention
dftmp=rprotds_wide[!is.na(rprotds_wide$timepoint.y) & rprotds_wide$Intervention==1,]
for (i in 1:length(protlist)){
  print(protlist[i])
  print(paste(round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                             probs=c(0.5)),3),
              "(",
              round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                             probs=c(0.25)),3),
              ",",
              round(quantile(dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                               dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"],
                             probs=c(0.75)),3),
              ")",sep=""))
  
  print(wilcox.test(x=dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.y"]-
                      dftmp[dftmp$Protein==protlist[i],"ExpressionLevel.x"])$p.value)
}

#t-test to compare difference across groups
dftmp=rprotds_wide[!is.na(rprotds_wide$timepoint.y),]
for (i in 1:length(protlist)){
  print(protlist[i])
  print(t.test(Protein_change ~ Intervention, data=dftmp[dftmp$Protein==protlist[i],],var.equal=TRUE)$p.value)
}  


#mediation analysis
#6mwd
rprotblds_long2=pivot_longer(data=rprotblds[,
                                           c("subject_id","timepoint","intervention_1_control_0",
                                             "bl_smw","fu_smw","bl_sppb","fu_sppb",
                                             "age","sex","race___4","hf_cat",protlist)],
                            cols=protlist,
                            names_to='Protein',values_to='ExpressionLevel')
rprotfuds_long2=pivot_longer(data=rprotfuds[,
                                            c("subject_id","timepoint","intervention_1_control_0",
                                              "bl_smw","fu_smw","bl_sppb","fu_sppb",
                                              "age","sex","race___4","hf_cat",protlist)],
                             cols=protlist,
                             names_to='Protein',values_to='ExpressionLevel')
rprotblds_long2$sppb_chg = rprotblds_long2$fu_sppb - rprotblds_long2$bl_sppb
rprotblds_long2$sppb_logchg = log2(rprotblds_long2$fu_sppb) - log2(rprotblds_long2$bl_sppb)
rprotblds_long2$smw_chg = rprotblds_long2$fu_smw - rprotblds_long2$bl_smw
rprotblds_long2$smw_logchg = log2(rprotblds_long2$fu_smw) - log2(rprotblds_long2$bl_smw)
rprotfuds_long2$sppb_chg = rprotfuds_long2$fu_sppb - rprotfuds_long2$bl_sppb
rprotfuds_long2$sppb_logchg = log2(rprotfuds_long2$fu_sppb) - log2(rprotfuds_long2$bl_sppb)
rprotfuds_long2$smw_chg = rprotfuds_long2$fu_smw - rprotfuds_long2$bl_smw
rprotfuds_long2$smw_logchg = log2(rprotfuds_long2$fu_smw) - log2(rprotfuds_long2$bl_smw)
colnames(rprotblds_long2)[3] <- "Intervention"
colnames(rprotfuds_long2)[3] <- "Intervention"
rprotblds_long2$Protein = factor(rprotblds_long2$Protein,levels=protlist)
rprotfuds_long2$Protein = factor(rprotfuds_long2$Protein,levels=protlist)

dftmp.all=merge(rprotblds_long2,
                   rprotfuds_long2,by=c('subject_id','Protein','Intervention',
                                       "bl_smw","fu_smw","bl_sppb","fu_sppb","sppb_chg","smw_chg",
                                       "sppb_logchg","smw_logchg",
                                       "age","sex","race___4","hf_cat"),
                   all.x=TRUE)
dftmp.all$Protein_change = (dftmp.all$ExpressionLevel.y) - (dftmp.all$ExpressionLevel.x)

pdf("Results/mediationanalysis_6mwd.pdf",onefile=TRUE)
for (i in 1:length(protlist)){
  #limited to propensity matched participants
  dftmp = dftmp.all[!is.na(dftmp.all$smw_chg) & 
                      !is.na(dftmp.all$Protein_change) & 
                      dftmp.all$subject_id %in% matchedset.6mwd$subject_id &
                      dftmp.all$Protein==protlist[i],]
  
  #model the effect of the intervention on the outcome
  total_model <- lm(smw_chg ~  Intervention+ExpressionLevel.x + age + sex + race___4 + hf_cat,
                       data = dftmp)

  #model the effect of the intervention on the mediator (change in  level)
  mediator_model <- lm(Protein_change ~  Intervention+ExpressionLevel.x + age + sex + race___4 + hf_cat,
                        data = dftmp)

  #model the outcome with both the intervention and mediator (protein) included
  outcome_model <- lm(smw_chg ~  Intervention+Protein_change +ExpressionLevel.x+age + sex + race___4 + hf_cat, 
                       data = dftmp)

  # Estimate the direct and indirect effects using the mediation package
  med_results <- mediate(mediator_model, outcome_model, 
                         treat = "Intervention", mediator = "Protein_change",
                         boot = TRUE, sims = 1000)  # Bootstrap for confidence intervals
  
  # Print summary of mediation results
  print(summary(med_results))
  plot(med_results,main=protlist[i])
}
dev.off()
#sppb
pdf("Results/mediationanalysis_sppb.pdf",onefile=TRUE)
for (i in 1:length(protlist)){
  #limited to propensity matched participants
  dftmp = dftmp.all[!is.na(dftmp.all$sppb_chg) & 
                      !is.na(dftmp.all$Protein_change) & 
                      dftmp.all$subject_id %in% matchedset.sppb$subject_id &
                      dftmp.all$Protein==protlist[i],]
  
  #model the effect of the intervention on the outcome
  total_model <- lm(sppb_chg ~  Intervention+ExpressionLevel.x + age + sex + race___4 + hf_cat,
                    data = dftmp)
  
  #model the effect of the intervention on the mediator (change in  level)
  mediator_model <- lm(Protein_change ~  Intervention+ExpressionLevel.x + age + sex + race___4 + hf_cat,
                       data = dftmp)
  
  #model the outcome with both the intervention and mediator (protein) included
  outcome_model <- lm(sppb_chg ~  Intervention+Protein_change +ExpressionLevel.x+age + sex + race___4 + hf_cat, 
                      data = dftmp)
  
  # Estimate the direct and indirect effects using the mediation package
  med_results <- mediate(mediator_model, outcome_model, 
                         treat = "Intervention", mediator = "Protein_change",
                         boot = TRUE, sims = 1000)  # Bootstrap for confidence intervals
  
  # Print summary of mediation results
  print(summary(med_results))
  plot(med_results,main=protlist[i])
}
dev.off()
