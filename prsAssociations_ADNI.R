# Load packages ####
library(data.table)
library(ADNIMERGE) # Must be downloaded from ADNI database
library(jtools)
library(dplyr)
library(nonnest2)
library(ggplot2)
library(broom)
library(sjlabelled)
library(tableone)
library(clusterProfiler)
library(ggpubr)
library(corrplot)
source("/home/jelman/Github/misc/Meff_Based_FDR.R")

# Create mapping of new pathway names with PRSice outputs
pathwayLookup = c(Base = "Global PRS",
                  tau_protein_localization = "protein localization",
                  cholesterol_efflux_transport = "cholesterol transport",
                  amyloid_protein_process = "amyloid protein processing",
                  immune_inhibiting_signal = "immune signalling",
                  leukocyte_activation_inflammatory = "inflammatory response",
                  fibril_negative_endocytosis = "endocytosis and fibril regulation",
                  humoral_mediated_circulating = "humoral immune response",
                  receptor_metabolic_process = "receptor metabolic process",
                  response_misfolded_protein = "response misfolded protein",
                  phototransduction = "phototransduction",
                  regulation_cell_junction = "regulation cell junction",
                  regulation_protein_tyrosine = "regulation protein tyrosine",
                  mitophagy = "mitophagy")

# Define functions ####
# Function to create common ordering of dataframes and rename pathways
renamePathways <- function(df, lookup_vector) {
  df  = df %>% mutate(term = factor(term),
                      term = recode_factor(term, !!!lookup_vector))
  df
}

# Plot coefficients from regression models
plotCoefficients <- function(df, plotTitle) {
  ggplot(df, aes(x=estimate, y=term, color=forcats::fct_rev(Model), alpha=FDRsig, xmin=conf.low, xmax=conf.high)) +
    geom_vline(xintercept = 0, color = "dark gray") +
    facet_grid(. ~ Model, labeller = labeller(Model = label_wrap_gen(12))) +
    geom_pointrange(size=ifelse(df$FDRsig==TRUE, 1, .5), position=position_dodge(width=1)) + 
    scale_y_discrete(limits=rev) +
    scale_x_continuous(labels = scales::label_number(accuracy = .1)) + 
    theme_nice(base_family = "Arial") + 
    ggsci::scale_color_jama(palette="default", limit=rev, guide="none") +
    scale_alpha_discrete(range = c(0.4, 1), guide="none") +
    labs(x="Estimate", y="", alpha="FDR sig.") + ggtitle(plotTitle) +
    theme(text=element_text(size=24, face="plain"), 
          axis.text.y = element_text(face="plain", hjust=0),
          axis.text.x = element_text(size=14, family = "Arial Narrow", face="plain"),
          axis.title.x = element_text(face="plain"),
          strip.text.x = element_text(face="plain"),
          plot.title = element_text(size=32, face="plain", hjust = 0.5))
}

# Calculate biomarker positivity using standard thresholds
applyCutpoints <- function(dat){
  require(dplyr)
  ### UPENN BIOMARKER 9
  # ABETA < 977
  # PTAU > 21.8
  
  ### PET
  # PIB > 1.44 (Landau, 2013) 
  # AV45 > 1.11 (Royse, 2021)
  # FBB > 1.08 (Royse, 2021)
  # TAU >=1.78 (Weigand, 2022)
  
  # Convert columns to numeric
  numVars = c("PTAU", "TAU", "ABETA", "PIB", "AV45", "FBB", "META_TEMPORAL_SUVR")
  dat[numVars] = sapply(dat[numVars],as.numeric)
  # Create variables for biomarker positivity
  dat$ABETA.pos = ifelse(dat$ABETA<977, 1, 0)
  dat$PTAU.pos = ifelse(dat$PTAU>21.8, 1, 0)
  dat$AV45.pos = ifelse(dat$AV45>1.11, 1, 0)
  dat$PIB.pos = ifelse(dat$PIB>1.44, 1, 0)  
  dat$FBB.pos = ifelse(dat$FBB>1.08, 1, 0)  
  
  # PVC-corrected tau PET
  dat$META_TEMPORAL_SUVR.pos = ifelse(dat$META_TEMPORAL_SUVR>=1.78, 1, 0) 
  
  # return a dataframe of variables indicating positivity on each biomarker
  biocuts = dat %>% 
    mutate(EXAMDATE = as.character(EXAMDATE)) %>%
    dplyr::select(PTID, EXAMDATE, ends_with(".pos")) %>%
    sjlabelled::remove_all_labels()
  biocuts
}


### Load data ####
prs = fread("/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data/PRSADNI_Kunkle_ClusteredPathways_q25_o5_p1_NoAPOE.best", colClasses=c(IID="character",FID="character"), header = TRUE)
prs = prs %>% 
  select(PTID=IID, everything(), -FID) %>%
  filter(!grepl("HG|NA", PTID)) # Remove 1KG subjects if included

# Clean up PRSice output
prs = prs %>% select(-In_Regression)
names(prs) = gsub("_1","",names(prs))
prs_names = names(prs)[2:ncol(prs)]
# Expected variables names: 
# [1] "Base"                              "protein_tyrosine_kinase"           "negative_precursor_catabolic"     
# [4] "regulation_component_organization" "protein_lipid_complex"             "microglial_cell_activation"       
# [7] "negative_regulation_fibril"        "lipoprotein_particle_receptor"     "humoral_immune_response"  

# Load genetic principal components
pcs = read.csv("/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data/ADNI_genetic_data/ADNI_1_GO2_3_PC1-10.csv")
pcs = pcs %>% rename(PTID=IID)

# Z-score PCs
prs = prs %>% left_join(pcs, by="PTID") %>%
  mutate_if(is.numeric, scale)
str(prs)


## Load ADNI data ####
demodf = adnimerge %>% 
  mutate(AGE = AGE + Years.bl) %>% 
  filter(!(RID==1195 & VISCODE=="m42"), !(RID==4960 & VISCODE=="m48")) %>%
  select(PTID, RID, VISCODE, EXAMDATE, ORIGPROT, COLPROT, Years.bl,
         DX, DX.bl, APOE4, PTGENDER, AGE, PTEDUCAT,
         PIB, AV45, ABETA, PTAU, Hippocampus, ICV, FLDSTRENGB) %>%
  mutate(RID = as.integer(RID),
         PTID = as.character(PTID),
         EXAMDATE = as.character(EXAMDATE),
         APOE4pos = factor(ifelse(APOE4>0, 1, 0)),
         APOE4 = as.integer(APOE4),
         FLDSTRENG = factor(FLDSTRENG),
         PTGENDER = factor(PTGENDER),
         DX = factor(DX),
         DX.bl = factor(DX.bl),
         ABETA = as.numeric(ABETA),
         PTAU = as.numeric(PTAU)) %>%
  sjlabelled::remove_all_labels()

str(demodf)


# Associations with diagnosis ####

## Create dx dataset ----

# Merge ADNI data with PRS
prs.demo.df = prs %>%
  left_join(demodf, by="PTID") %>%
  sjlabelled::remove_all_labels() %>%
  mutate(EXAMDATE = as.character(EXAMDATE),
         RID = as.character(RID))

# Select timepoint of interest
subdf = prs.demo.df %>%
  group_by(PTID) %>%
  arrange(EXAMDATE) %>%
  filter(!is.na(DX)) %>%
  filter(row_number()==1)%>%
  ungroup() 

# Filter for CN vs AD (exclude MCI)
dxdf = subdf %>% 
  filter(DX!="MCI") %>%
  mutate(DX = factor(DX))


## Run associations with dx ----

# Get full results of Base PRS model
dx_base = glm(DX ~ Base + APOE4 + PTGENDER + AGE + PC1 + PC2 + PC3, family="binomial", data=dxdf)
dx_base_summ = summ(dx_base, confint=TRUE, scale=TRUE)
broom::tidy(dx_base_summ)

## Pathway-specific PRS analyses: Controlling for APOE
# Get results for all pathway scores
dx_res_apoe = lapply(prs_names, function(prsScore){
  fmla = paste0("DX ~ ", prsScore," + APOE4 + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=dxdf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
dx_res_apoe = bind_rows(dx_res_apoe) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results and format results
dx_meff_apoe = Meff_based_FDR(dxdf %>% select(one_of(prs_names)), dx_res_apoe %>% select(term, p.value))
dx_res_apoe = dx_res_apoe %>% inner_join(select(dx_meff_apoe, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-negative only
# Get results for all pathway scores
dx_res_apoe4neg = lapply(prs_names, function(prsScore){
  fmla = paste0("DX ~ ", prsScore," + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=subset(dxdf, APOE4==0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
dx_res_apoe4neg = bind_rows(dx_res_apoe4neg) %>% mutate(Model="APOE-e4 neg. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
dx_meff_apoe4neg = Meff_based_FDR(dxdf %>% select(one_of(prs_names)), dx_res_apoe4neg %>% select(term, p.value))
dx_res_apoe4neg = dx_res_apoe4neg %>% inner_join(select(dx_meff_apoe4neg, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-positive only
# Get results for all pathway scores
dx_res_apoe4pos = lapply(prs_names, function(prsScore){
  fmla = paste0("DX ~ ", prsScore," + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=subset(dxdf, APOE4>0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
dx_res_apoe4pos = bind_rows(dx_res_apoe4pos) %>% mutate(Model="APOE-e4 pos. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
dx_meff_apoe4pos = Meff_based_FDR(dxdf %>% select(one_of(prs_names)), dx_res_apoe4pos %>% select(term, p.value))
dx_res_apoe4pos = dx_res_apoe4pos %>% inner_join(select(dx_meff_apoe4pos, term, FDRsig=`Significance based on FDR_E`), by="term")

## Bind results from all APOE, APOE4+ only, and APOE4- only
dx_res_all = bind_rows(dx_res_apoe, dx_res_apoe4neg, dx_res_apoe4pos) %>%
  mutate(Model = factor(Model, levels = c("Controlling for APOE-e4", "APOE-e4 neg. only", "APOE-e4 pos. only")))
dx_res_all = renamePathways(dx_res_all, pathwayLookup)

# Save out analysis of diagnosis
write.csv(dx_res_all, row.names = FALSE, file="/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/results/dx_coef_table.csv")

## Compare model with Base PRS against model with all pathway PRS using Vuong test
fmla = paste0("DX ~ ", paste(prs_names[2:length(prs_names)], collapse=" + "), " + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3")
dx_paths = glm(fmla, family="binomial", data=dxdf)
vuongtest(dx_base, dx_paths)

## APOE x PRS interaction
# Get results for all pathway scores
dx_res_int = lapply(prs_names, function(prsScore){
  fmla = paste0("DX ~ ", prsScore,"*APOE4pos + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=dxdf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[9,]
})

# Bind results from all models
dx_res_int = bind_rows(dx_res_int) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
dx_meff_int = Meff_based_FDR(dxdf %>% select(one_of(prs_names)), dx_res_int %>% select(term, p.value))
dx_res_int = dx_res_int %>% inner_join(select(dx_meff_int, term, FDRsig=`Significance based on FDR_E`), by="term")



#------------------------------------------#
#   START CREATING BIOMARKER DATASET    ####
#------------------------------------------#

# Creates variables ending in ".pos" for biomarker measurements to indicate 
# whether the measure is positive or not based on standard cutpoints

# Function to calculate SUVR
calcSUVR <- function(x, ref = INFERIORCEREBELLUM_SUVR) (x / ref)

# Calculate SUVR for AV1451 (tau PET)
av1451 = ucberkeleyav1451_pvc %>%
  select(RID, EXAMDATE, META_TEMPORAL_SUVR, INFERIOR_CEREBGM_SUVR) %>% 
  mutate_at(vars(META_TEMPORAL_SUVR), ~calcSUVR(., ref=INFERIOR_CEREBGM_SUVR)) %>%
  mutate(EXAMDATE = as.character(EXAMDATE),
         RID = as.character(RID))

# Merge tau PET data with ADNI data
tmpdf = adnimerge %>% 
  mutate(RID = as.character(RID),
         EXAMDATE = as.character(EXAMDATE)) %>%
  sjlabelled::remove_all_labels() %>% 
  left_join(av1451, by=c("RID","EXAMDATE"))

# Calculate biomarker positivity across all measures
biocuts = applyCutpoints(tmpdf)
# Calculate amyloid-positivity across PET and CSF
biocuts$anyAB.pos = biocuts %>% select(ABETA.pos, PIB.pos, AV45.pos, FBB.pos) %>% apply(., 1, max, na.rm=T) 
biocuts$anyAB.pos = replace(biocuts$anyAB.pos, is.infinite(biocuts$anyAB.pos), NA)
# Calculate tau-positivity across PET and CSF
biocuts$anyPTAU.pos = biocuts %>% select(PTAU.pos, META_TEMPORAL_SUVR.pos) %>% apply(., 1, max, na.rm=T)
biocuts$anyPTAU.pos = replace(biocuts$anyPTAU.pos, is.infinite(biocuts$anyPTAU.pos), NA)

# Merge cutpoint data with adnimerge dataset
df = prs.demo.df %>% 
  filter(!(RID==1195 & VISCODE=="m42")) %>%
  left_join(biocuts, by=c("PTID","EXAMDATE")) %>%
  left_join(av1451, by=c("RID","EXAMDATE"))


#----------------------------------#
# Associations with biomarkers ####
#----------------------------------#

## Create dataset for biomarker analyses ----
# Get first observation where subject has both amyloid, tau and hippo volume
biodf = df %>%
  group_by(RID) %>%
  arrange(RID, EXAMDATE) %>%
  filter(!is.na(PTAU.pos) & !is.na(anyAB.pos) & !is.na(Hippocampus)) %>%
  filter(row_number()==1)%>%
  mutate(FLDSTRENG = ifelse(is.na(FLDSTRENG) & COLPROT!="ADNI1", 2, FLDSTRENG),
         Hippo_ICV = (Hippocampus / ICV) * 100) %>%
  ungroup() 


## ANALYSIS OF AMYLOID POSITIVITY ####

## Get full results of Base PRS model
ab_base = glm(anyAB.pos ~ Base + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3, family="binomial", data=biodf)
ab_base_summ = summ(ab_base, confint=TRUE, scale=TRUE) 
broom::tidy(ab_base_summ)

## Pathway-specific PRS analyses: Controlling for APOE
# Get results for all pathway scores
ab_res_apoe = lapply(prs_names, function(prsScore){
  fmla = paste0("anyAB.pos ~ ", prsScore," + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=biodf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
ab_res_apoe = bind_rows(ab_res_apoe) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ab_meff_apoe = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ab_res_apoe %>% select(term, p.value))
ab_res_apoe = ab_res_apoe %>% inner_join(select(ab_meff_apoe, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-negative only
# Get results for all pathway scores
ab_res_apoe4neg = lapply(prs_names, function(prsScore){
  fmla = paste0("anyAB.pos ~ ", prsScore," + PTGENDER  + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=subset(biodf, APOE4==0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
ab_res_apoe4neg = bind_rows(ab_res_apoe4neg) %>% mutate(Model="APOE-e4 neg. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ab_meff_apoe4neg = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ab_res_apoe4neg %>% select(term, p.value))
ab_res_apoe4neg = ab_res_apoe4neg %>% inner_join(select(ab_meff_apoe4neg, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-positive only
# Get results for all pathway scores
ab_res_apoe4pos = lapply(prs_names, function(prsScore){
  fmla = paste0("anyAB.pos ~ ", prsScore," + PTGENDER  + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=subset(biodf, APOE4>0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
ab_res_apoe4pos = bind_rows(ab_res_apoe4pos) %>% mutate(Model="APOE-e4 pos. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ab_meff_apoe4pos = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ab_res_apoe4pos %>% select(term, p.value))
ab_res_apoe4pos = ab_res_apoe4pos %>% inner_join(select(ab_meff_apoe4pos, term, FDRsig=`Significance based on FDR_E`), by="term")

## Bind results from all APOE, APOE4+ only, and APOE4- only
ab_res_all = bind_rows(ab_res_apoe, ab_res_apoe4neg, ab_res_apoe4pos) %>%
  mutate(Model = factor(Model, levels = c("Controlling for APOE-e4", "APOE-e4 neg. only", "APOE-e4 pos. only")))
ab_res_all = renamePathways(ab_res_all, pathwayLookup)

# Save out results of amyloid analyses
write.csv(ab_res_all, row.names = FALSE, file="/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/results/ab_coef_table.csv")

## Compare model with Base PRS against model with all pathway PRS using Vuong test
fmla = paste0("anyAB.pos ~ ", paste(prs_names[2:length(prs_names)], collapse=" + "), " + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3")
ab_paths = glm(fmla, family="binomial", data=biodf)
vuongtest(ab_base, ab_paths)

## APOE x PRS interaction
# Get results for all pathway scores
ab_res_int = lapply(prs_names, function(prsScore){
  fmla = paste0("anyAB.pos ~ ", prsScore,"*APOE4pos + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=biodf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[9,]
})
ab_res_int = bind_rows(ab_res_int) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ab_meff_int = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ab_res_int %>% select(term, p.value))
ab_meff_int
ab_res_int = ab_res_int %>% inner_join(select(ab_meff_int, term, FDRsig=`Significance based on FDR_E`), by="term")


## ANALYSIS OF PTAU POSITIVITY ####

## Get full results of Base PRS model
ptau_base = glm(PTAU.pos ~ Base + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3, family="binomial", data=biodf)
ptau_base_summ = summ(ptau_base, confint=TRUE, scale=TRUE) 
broom::tidy(ptau_base_summ)

## Pathway-specific PRS analyses: Controlling for APOE
# Get results for all pathway scores
ptau_res_apoe = lapply(prs_names, function(prsScore){
  fmla = paste0("PTAU.pos ~ ", prsScore," + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=biodf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
ptau_res_apoe = bind_rows(ptau_res_apoe) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ptau_meff_apoe = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ptau_res_apoe %>% select(term, p.value))
ptau_res_apoe = ptau_res_apoe %>% inner_join(select(ptau_meff_apoe, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-negative only
# Get results for all pathway scores
ptau_res_apoe4neg = lapply(prs_names, function(prsScore){
  fmla = paste0("anyPTAU.pos ~ ", prsScore," + PTGENDER  + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=subset(biodf, APOE4==0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
ptau_res_apoe4neg = bind_rows(ptau_res_apoe4neg) %>% mutate(Model="APOE-e4 neg. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ptau_meff_apoe4neg = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ptau_res_apoe4neg %>% select(term, p.value))
ptau_res_apoe4neg = ptau_res_apoe4neg %>% inner_join(select(ptau_meff_apoe4neg, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-positive only
# Get results for all pathway scores
ptau_res_apoe4pos = lapply(prs_names, function(prsScore){
  fmla = paste0("PTAU.pos ~ ", prsScore," + PTGENDER  + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=subset(biodf, APOE4>0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
ptau_res_apoe4pos = bind_rows(ptau_res_apoe4pos) %>% mutate(Model="APOE-e4 pos. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ptau_meff_apoe4pos = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ptau_res_apoe4pos %>% select(term, p.value))
ptau_res_apoe4pos = ptau_res_apoe4pos %>% inner_join(select(ptau_meff_apoe4pos, term, FDRsig=`Significance based on FDR_E`), by="term")

## Bind results from all APOE, APOE4+ only, and APOE4- only
ptau_res_all = bind_rows(ptau_res_apoe, ptau_res_apoe4neg, ptau_res_apoe4pos) %>%
  mutate(Model = factor(Model, levels = c("Controlling for APOE-e4", "APOE-e4 neg. only", "APOE-e4 pos. only")))
ptau_res_all = renamePathways(ptau_res_all, pathwayLookup)

# Save out analyses of ptau positivity
write.csv(ptau_res_all, row.names = FALSE, file="/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/results/ptau_coef_table.csv")

## Compare model with Base PRS against model with all pathway PRS using Vuong test
fmla = paste0("anyPTAU.pos ~ ", paste(prs_names[2:length(prs_names)], collapse=" + "), " + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3")
ptau_paths = glm(fmla, family="binomial", data=biodf)
vuongtest(ptau_base, ptau_paths)

## APOE x PRS interaction
ptau_res_int = lapply(prs_names, function(prsScore){
  fmla = paste0("anyPTAU.pos ~ ", prsScore,"*APOE4pos + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(glm(fmla, family="binomial", data=biodf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[9,]
})
ptau_res_int = bind_rows(ptau_res_int) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
ptau_meff_int = Meff_based_FDR(biodf %>% select(one_of(prs_names)), ptau_res_int %>% select(term, p.value))
ptau_meff_int
ptau_res_int = ptau_res_int %>% inner_join(select(ptau_meff_int, term, FDRsig=`Significance based on FDR_E`), by="term")


## ANALYSIS OF HIPPOCAMPAL VOLUME ####

## Get full results of Base PRS model
hippo_base = lm(scale(Hippocampus) ~ Base + APOE4 + PTGENDER  + AGE + PC1 + PC2 + PC3, data=biodf)
hippo_base_summ = summ(hippo_base, confint=TRUE, scale=TRUE) 
broom::tidy(hippo_base_summ)

## Pathway-specific PRS analyses: Controlling for APOE
# Get results for all pathway scores
hippo_res_apoe = lapply(prs_names, function(prsScore){
  fmla = paste0("scale(Hippocampus) ~ ", prsScore," + APOE4 + FLDSTRENG + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(lm(fmla, data=biodf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
hippo_res_apoe = bind_rows(hippo_res_apoe) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
hippo_meff_apoe = Meff_based_FDR(biodf %>% select(one_of(prs_names)), hippo_res_apoe %>% select(term, p.value))
hippo_res_apoe = hippo_res_apoe %>% inner_join(select(hippo_meff_apoe, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-negative only
# Get results for all pathway scores
hippo_res_apoe4neg = lapply(prs_names, function(prsScore){
  fmla = paste0("scale(Hippocampus) ~ ", prsScore," + FLDSTRENG + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(lm(fmla, data=subset(biodf, APOE4==0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
hippo_res_apoe4neg = bind_rows(hippo_res_apoe4neg) %>% mutate(Model="APOE-e4 neg. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
hippo_meff_apoe4neg = Meff_based_FDR(biodf %>% select(one_of(prs_names)), hippo_res_apoe4neg %>% select(term, p.value))
hippo_res_apoe4neg = hippo_res_apoe4neg %>% inner_join(select(hippo_meff_apoe4neg, term, FDRsig=`Significance based on FDR_E`), by="term")

## Pathway-specific PRS analyses: APOE4-positive only
# Get results for all pathway scores
hippo_res_apoe4pos = lapply(prs_names, function(prsScore){
  fmla = paste0("scale(Hippocampus) ~ ", prsScore," + FLDSTRENG + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(lm(fmla, data=subset(biodf, APOE4>0)), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[2,]
})
hippo_res_apoe4pos = bind_rows(hippo_res_apoe4pos) %>% mutate(Model="APOE-e4 pos. only") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
hippo_meff_apoe4pos = Meff_based_FDR(biodf %>% select(one_of(prs_names)), hippo_res_apoe4pos %>% select(term, p.value))
hippo_res_apoe4pos = hippo_res_apoe4pos %>% inner_join(select(hippo_meff_apoe4pos, term, FDRsig=`Significance based on FDR_E`), by="term")

## Bind all together
hippo_res_all = bind_rows(hippo_res_apoe, hippo_res_apoe4neg, hippo_res_apoe4pos) %>%
  mutate(Model = factor(Model, levels = c("Controlling for APOE-e4", "APOE-e4 neg. only", "APOE-e4 pos. only")))
hippo_res_all = renamePathways(hippo_res_all, pathwayLookup)

# Save out analysis of hippocampal volume
write.csv(hippo_res_all, row.names = FALSE, file="/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/results/hippo_coef_table.csv")

## Compare model with Base PRS against model with all pathway PRS using Vuong test
fmla = paste0("scale(Hippocampus) ~ ", paste(prs_names[2:length(prs_names)], collapse=" + "), " + APOE4 + FLDSTRENG + PTGENDER  + AGE + PC1 + PC2 + PC3")
hippo_paths = lm(fmla, data=biodf)
vuongtest(hippo_base, hippo_paths)

## APOE x PRS interaction
hippo_res_int = lapply(prs_names, function(prsScore){
  fmla = paste0("scale(Hippocampus) ~ ", prsScore,"*APOE4pos + FLDSTRENG + PTGENDER + AGE + PC1 + PC2 + PC3")
  summ = summ(lm(fmla, data=biodf), confint=TRUE, scale=TRUE)
  summ = broom::tidy(summ)
  summ[10,]
})
hippo_res_int = bind_rows(hippo_res_int) %>% mutate(Model="Controlling for APOE-e4") %>% as.data.frame()
# Apply Li and Ji FDR correction and format results
hippo_meff_int = Meff_based_FDR(biodf %>% select(one_of(prs_names)), hippo_res_int %>% select(term, p.value))
hippo_meff_int
hippo_res_int = hippo_res_int %>% inner_join(select(hippo_meff_int, term, FDRsig=`Significance based on FDR_E`), by="term")

#---------------------#
# Plotting results ####
#---------------------#

# Rename pathway PRS variables 
dx_plot_data = renamePathways(dx_res_all, pathwayLookup)
ptau_plot_data = renamePathways(ptau_res_all, pathwayLookup)  
ab_plot_data = renamePathways(ab_res_all, pathwayLookup)
hippo_plot_data = renamePathways(hippo_res_all, pathwayLookup)

# Plot coefficients for dx, amyloid, ptau, hippocampal volume
plotCoefficients(dx_plot_data, "Alzheimer's diagnosis")
plotCoefficients(ab_plot_data, "Amyloid positivity")
plotCoefficients(ptau_plot_data, "Ptau positivity")
plotCoefficients(hippo_plot_data, "Hippocampal volume")

#-----------------------------------------------#
# Demographic tables and additional analyses ####
#-----------------------------------------------#

# Create demographic table for diagnosis analysis
CreateTableOne(vars=c("PTGENDER", "AGE", "PTEDUCAT", "APOE4pos"), 
               strata="DX", test=TRUE, data=dxdf)

# Create demographic tables for biomarker analyses
CreateTableOne(vars=c("PTGENDER", "AGE", "PTEDUCAT", "APOE4pos", "DX", "anyAB.pos", "anyPTAU.pos", "Hippo_ICV"), 
               factorVars = c("anyAB.pos", "anyPTAU.pos"),
               data=biodf)

# List how many SNPs were used to calculate each PRS
prsSnps = fread("/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data/PRS/ADNI/2022-02-01/ADNI_Kunkle_ClusteredPathways_q25_o5_p1_NoAPOE.snp")
prsSnps %>% select(one_of(prs_names)) %>% colSums()

# Calculate how many pathways each SNP contributes to and plot results. 
prsSnps = prsSnps %>% 
  rowwise(SNP) %>% 
  mutate(nPathways = sum(c_across(prs_names[!prs_names=="Base"]))) %>%
  ungroup()
snpPaths = prsSnps %>% 
  count(nPathways) %>% 
  filter(nPathways>0) %>%
  mutate(nPathways = factor(nPathways))
snpOverlap = ggbarplot(snpPaths, x="nPathways", y="n", 
                       fill="gray60", color="gray60", 
                       label=TRUE, title = "SNPs",
                       xlab= "n Pathways", ylab = "count") +
  font("title", size = 18)

# Calculate how many pathways each gene contributes to and plot results.
genesets = clusterProfiler::read.gmt("/home/jelman/netshare/BETELSHARE/PROJ/PathwayPRS_ATN/data/cytoscape/Kunkle_Magma_BaderLab_Trimmed_10-1000_q25_o5_ClusteredPathwayGeneSets_Symbols.gmt")
genePaths = genesets %>% count(gene, name="nPathways") %>%
  count(nPathways) %>%
  filter(nPathways>0) %>%
  mutate(nPathways = factor(nPathways))
geneOverlap = ggbarplot(genePaths, x="nPathways", y="n", 
                        fill="gray60", color="gray60", 
                        label=TRUE, title = "Genes",
                        xlab= "n Pathways", ylab = "count") +
  font("title", size = 18)  

# Calculate correlation between all PRS. Rename columns first
prsCorrMat = prs %>% 
  select(one_of(prs_names)) %>% 
  rename_with(~ pathwayLookup[.x], everything()) %>% GGally::ggcorr(.)
cor(., use="pairwise.complete") 

prsCorrPlot = corrplot(prsCorrMat, method = "color",
                       type = "upper",
                       addCoef.col = "black", # Add coefficient of correlation
                       tl.col = "black", tl.srt = 45, # Text label color and rotation
                       diag = FALSE) # hide correlation coefficient on the principal diagonal
