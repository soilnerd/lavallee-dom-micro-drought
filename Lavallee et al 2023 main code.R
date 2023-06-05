###     Code for Lavallee et al. (2023) Land management shapes drought    ###
###     responses of dominant soil microbial taxa across grasslands       ###
setwd("...")

# Load libraries ----
library(tidyverse)
library(dplyr)
library(sjmisc)
library(vegan) 
library(lavaan)
library(nlme)
# Read data ----
# OTU tables
load("OTU_tab_transposed_bact.RData")
load("OTU_tab_transposed_fung.RData")

# OTU representative sequences and taxonomy
load("OTU_tab_rep_seqs_tax_bact.RData")
load("OTU_tab_rep_seqs_tax_fung.RData")

# Env data
load("Lavallee_et_al_2023_env_data.RData")


# Define dominant taxa ----
top_rel_abun <- 0.1 # top proportion of OTUs when ranked by mean relative abundance across all samples
pres_level <- "site" # OTU must always be present at X level to be defined as dominant
sample_pres_cutoff <- 1 # OTU must be present in at least X sample per pres_level


OTUs_dom_bact <- OTU_tab_bact %>%
  summarize(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "OTU") %>% 
  slice_max(V1, prop = top_rel_abun) %>% # filter top X% of OTUs by mean abundance across all samples
  filter(OTU %in% (OTU_tab_bact %>%
                     mutate(across(where(is.numeric), ~1 * (. > 0))) %>% # convert to presence/absence (0/1)
                     group_by(across(all_of(pres_level))) %>%
                     summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% # calculate proportion of samples within a site where OTU is present (1)
                     select(-where(~is.numeric(.x) && any(.x < sample_pres_cutoff))) %>% # filter any OTUs not present at pres_level
                     names())) %>%
  pull(OTU)


OTUs_dom_fung <- OTU_tab_fung %>%
  summarize(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column(var ="OTU") %>%
  slice_max(V1, prop = top_rel_abun) %>%
  filter(OTU %in% (OTU_tab_fung %>%
                     mutate(across(where(is.numeric), ~1 * (. > 0))) %>%
                     group_by(across(all_of(pres_level))) %>%
                     summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% # / n() *100)) %>% # calculate percent of samples within a site where OTU is present (1)
                     select(-where(~is.numeric(.x) && any(.x < sample_pres_cutoff))) %>%
                     names())) %>%
  pull(OTU)

OTU_tab_dom_bact <- OTU_tab_bact %>%
  select(c(sample, region, site, mgmt, field, trt, day), 
         all_of(OTUs_dom_bact))

OTU_tab_dom_fung <- OTU_tab_fung %>%
  select(c(sample, region, site, mgmt, field, trt, day), 
         all_of(OTUs_dom_fung))

# Assign drought response and resilience categories  ----
# Can check source code below, but it may be slow and may not run properly on
# different machines or different package versions (see note in code). Instead, 
# you can load the resulting file which includes the model run for each OTU.
# source("Lavallee et al 2023 Resis Resil code.R")
load("OTU_drought_classifications_and_abundances.RData")

# Calculate indices for each drought response category ----
# Standardize the relative abundance data so that all taxa to contribute equally
# to each drought response category index. Otherwise, the most abundant taxa 
# will have more weight in each index than the less abundant ones.

# calculate min and max counts across all samples, and standardize for each 
# OTU in each sample:
# (#counts - min(all counts of OTU)) / 
# ((max(all counts of OTU) - min(all counts of OTU))

index_dom_bact <-  OTU_tab_dom_bact %>%
  select(-c(region:day)) %>%
  rotate_df(rn = "OTU", cn = TRUE) %>%
  rowwise() %>%
  mutate(maxCount = max(across(where(is.numeric))),
         minCount = min(across(where(is.numeric)))) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), 
                ~(. - minCount)/ (maxCount-minCount))) %>%
  select(-c(maxCount, minCount)) %>%
  data.frame() %>%
  rotate_df(rn = "sample", cn = TRUE) %>%
  # calculate group-wise indices using categories from other data frame
  rowwise() %>%
  mutate(index_bact_opp = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of opportunistic OTUs
                                                 filter(Group == "16S" & Resistance_class == "Opportunistic") %>%
                                                 select(OTU_ID) %>%
                                                 pull()))),
         index_bact_sens = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of sensitive OTUs
                                                  filter(Group == "16S" & Resistance_class == "Sensitive") %>%
                                                  select(OTU_ID) %>%
                                                  pull()))),
         index_bact_resis = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of resistant OTUs
                                                 filter(Group == "16S" & Resistance_class == "Resistant") %>%
                                                 select(OTU_ID) %>%
                                                 pull()))),
         index_bact_resil = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of resilient OTUs
                                                   filter(Group == "16S" & Resilience_class == "Resilient") %>%
                                                   select(OTU_ID) %>%
                                                   pull()))),
         index_bact_not_resil = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of non-resilient OTUs
                                                       filter(Group == "16S" & Resilience_class == "Not Resilient") %>%
                                                       select(OTU_ID) %>%
                                                       pull()))),
         opp_sens_ratio_bact = index_bact_opp/ index_bact_sens,
         resil_not_resil_ratio_bact = index_bact_resil / index_bact_not_resil
  ) %>%
  ungroup() %>%
  select(sample, index_bact_opp, index_bact_sens, index_bact_resis, 
         index_bact_resil, opp_sens_ratio_bact, resil_not_resil_ratio_bact)


index_dom_fung <- OTU_tab_dom_fung %>%
  select(-c(region:day)) %>%
  rotate_df(rn = "OTU", cn = TRUE) %>%
  rowwise() %>%
  mutate(maxCount = max(across(where(is.numeric))),
         minCount = min(across(where(is.numeric)))) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), 
                ~(. - minCount)/ (maxCount-minCount))) %>%
  select(-c(maxCount, minCount)) %>%
  data.frame() %>%
  rotate_df(rn = "sample", cn = TRUE) %>%
  # calculate group-wise indices using categories from other data frame
  rowwise() %>%
  mutate(index_fung_opp = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of opportunistic OTUs
                                                 filter(Group == "ITS" & Resistance_class == "Opportunistic") %>%
                                                 select(OTU_ID) %>%
                                                 pull()))),
         index_fung_sens = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of sensitive OTUs
                                                  filter(Group == "ITS" & Resistance_class == "Sensitive") %>%
                                                  select(OTU_ID) %>%
                                                  pull()))),
         index_fung_resis = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of resistant OTUs
                                                 filter(Group == "ITS" & Resistance_class == "Resistant") %>%
                                                 select(OTU_ID) %>%
                                                 pull()))),
         index_fung_resil = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of resilient OTUs
                                                   filter(Group == "ITS" & Resilience_class == "Resilient") %>%
                                                   select(OTU_ID) %>%
                                                   pull()))),
         index_fung_not_resil = mean(c_across(all_of(drought_responses_and_abundances %>% # make vector of non-resilient OTUs
                                                       filter(Group == "ITS" & Resilience_class == "Not Resilient") %>%
                                                       select(OTU_ID) %>%
                                                       pull()))),
         opp_sens_ratio_fung = index_fung_opp / index_fung_sens,
         resil_not_resil_ratio_fung = index_fung_resil / index_fung_not_resil
  ) %>%
  ungroup() %>%
  select(sample, index_fung_opp, index_fung_sens, index_fung_resis, 
         index_fung_resil, opp_sens_ratio_fung, resil_not_resil_ratio_fung)


# Join into one dataset
dat_indices <- dat_env %>%
  left_join(index_dom_bact) %>%
  left_join(index_dom_fung) %>%
  mutate(region = fct_relevel(region, c("De", "Yo", "Sc")))

# Create dataset for SEM with soil variable ----
# Aggregate data at treatment level (average three replicates at each sampling)
dat_sem <- dat_indices %>%
  filter(day == "0" | day == "60") %>%
  select(-abv_plant_biomass) %>% # only have data for day 0
  filter(complete.cases(.)) %>%
  group_by(region, site, management, treatment, day) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
  ungroup()

# Create matrix for soil variable, with lat/long
soil_var_matrix <- dat_sem %>% 
  select(sand_pct, soil_tot_C, soil_tot_N, soil_temp)

coords <- dat_sem %>%
  select(lat, lon)

#library(vegan)
dist <- vegdist(coords, method = "euclidean")

# Define  PCNM vectors for soil variable
pcnm_vectors <- pcnm(dist)
PCNMs <- as.data.frame(pcnm_vectors$vectors)
soil_var_matrix_STD <- decostand(soil_var_matrix, method = "standardize")
resid_detrend <- residuals(rda(soil_var_matrix_STD ~ coords$lat + coords$lon))

# Fit PCNM vectors to Soil Data Matrix after detrending
mod0_space.capAgg <- capscale(resid_detrend ~ 1, PCNMs)  # Model with intercept only
mod1_space.capAgg <- capscale(resid_detrend ~ ., PCNMs) 

step.res <- ordiR2step(mod0_space.capAgg, mod1_space.capAgg, perm.max = 200)
step.res$anova  # Summary table

# Need PCNM 5, 6, and 4 to model spatial patterns in environmental variables 
PCNM_5 <- pcnm_vectors$vec[,5]
PCNM_6 <- pcnm_vectors$vec[,6]
PCNM_4 <- pcnm_vectors$vec[,4]

# Final spatial multivariate linear model for Soil Variables
mod_space <- capscale(soil_var_matrix_STD ~ coords$lat + coords$lon + 
                       PCNM_5 + PCNM_6 + PCNM_4)
resid_space <- residuals(mod_space)

summary(mod_space) # 44% of variance is spatially structured
anova(mod_space, by = "terms")

mod_space_noXY <- rda(soil_var_matrix_STD ~ PCNM_5 + PCNM_6 + PCNM_4) 
summary(mod_space_noXY) ## 15% accounted by autocorrelation

# Create a soil variable that accounts for spatial autocorrelation 
# (but not latitudinal effects to be fitted by SEM)
PCoA_ord_soil <- capscale(residuals(mod_space_noXY)~1) 
summary(PCoA_ord_soil)$cont 
PCoA_axes <- scores(PCoA_ord_soil, choices = c(1,2,3,4))
ord_axes_matrix <- as.data.frame(PCoA_axes$sites)

# Add PCoA ordination soil variable MDS1 to field data
dat_sem <- dat_sem %>%
  bind_cols(ord_axes_matrix %>%
              rownames_to_column("rowname") %>%
              select(MDS1) %>%
              rename(soil_var = MDS1))


# Test for drought effect on full microbial communities immediately following drought (day 0) ----
#
# Create Bray-Curtis distance data sets and associated labels
# Subset day 0 and remove labels 
OTU_tab_bact_D0_nolab <- OTU_tab_bact %>%
  filter(day == "0") %>%
  select(-c(sample, region, site, mgmt, field, trt, day)) %>%
  remove_rownames() %>%
  column_to_rownames("id")

OTU_tab_fung_D0_nolab <- OTU_tab_fung %>%
  filter(day == "0") %>%
  select(-c(sample, region, site, mgmt, field, trt, day)) %>%
  remove_rownames() %>%
  column_to_rownames("id")

### Pull labels separately
OTU_tab_bact_D0_lab <- OTU_tab_bact %>%
  filter(day == "0") %>%
  select(c(sample, region, site, mgmt, field, trt, day)) 

OTU_tab_fung_D0_lab <- OTU_tab_fung %>%
  filter(day == "0") %>%
  select(c(sample, region, site, mgmt, field, trt, day))

## Create distance matrix
#library(vegan)
bact_D0_dist = as.matrix((vegdist(OTU_tab_bact_D0_nolab, "bray")))
fung_D0_dist = as.matrix((vegdist(OTU_tab_fung_D0_nolab, "bray")))

# ADONIS tests
bact_Adonis <- adonis(bact_D0_dist ~  region * mgmt * trt + field, 
                      strata = OTU_tab_bact_D0_lab$field, data = OTU_tab_bact_D0_lab) 
bact_Adonis

fungAdonis <- adonis(fung_D0_dist ~  region * mgmt * trt + field, 
                     strata = OTU_tab_fung_D0_lab$field, data = OTU_tab_fung_D0_lab) 
fungAdonis

# SEMs ----
# library(lavaan)
## BACTERIA DAY 0 ----
model_sem_B <- '
   # soil properties and pH separately
   soil_var ~ lat + management
   pH ~ lat + management
   
   # soil moisture
   soil_moisture ~ treatment + management  + lat + soil_var
   
   # bacterial groups
   index_bact_opp ~ treatment + management  + pH + lat + soil_moisture + soil_var
   index_bact_resis ~ treatment + management + pH + lat + soil_moisture + soil_var
   index_bact_sens ~ treatment + management + pH + lat + soil_moisture + soil_var
   
   # extrinsic covariances
   management ~~ 0*treatment
   pH ~~ soil_var
   pH ~~ soil_moisture
   '

fit_model_sem_B0 <- sem(model_sem_B,
                        data = dat_sem %>%
                          filter(day == "0"),
                        std.lv = TRUE, std.ov = TRUE,
                        fixed.x = FALSE)
summary(fit_model_sem_B0, fit.measures = TRUE, standardized = T, rsq = T)  

# bootstrap for bootstrap test statistics
B0_boot <- bootstrapLavaan(fit_model_sem_B0, R = 1000L, type = "bollen.stine", verbose= FALSE,
                           FUN = fitMeasures, fit.measures = c("chisq"))
# compute a bootstrap based p-value
B0_fit_orig <- fitMeasures(fit_model_sem_B0, "chisq")
B0_boot_pvalue <- length(which(B0_boot > B0_fit_orig))/length(B0_boot)



## BACTERIA DAY 60 ----
fit_model_sem_B60 <- sem(model_sem_B,
                         data = dat_sem %>%
                           filter(day == "60") #%>%
                         #select(-c(index_fung_opp, index_fung_sens, index_fung_tol))
                         ,
                         std.lv = TRUE, std.ov = TRUE,
                         fixed.x = FALSE)
varTable(fit_model_sem_B60)  
summary(fit_model_sem_B60, fit.measures = TRUE, standardized = T, rsq = T)  

# bootstrap for bootstrap test statistics
B60_boot <- bootstrapLavaan(fit_model_sem_B60, R = 1000L, type = "bollen.stine", verbose= FALSE,
                            FUN = fitMeasures, fit.measures = c("chisq"))
# compute a bootstrap based p-value
B60_fit_orig <- fitMeasures(fit_model_sem_B60, "chisq")
B60_boot_pvalue <- length(which(B60_boot > B60_fit_orig))/length(B60_boot)

## FUNGI DAY 0 ----
model_sem_F <- '
   # soil properties and pH separately
   soil_var ~ lat + management
   pH ~ lat + management
   
   # soil moisture
   soil_moisture ~ treatment + management  + lat + soil_var
   
   # fungal groups
   index_fung_opp ~ treatment + management  + pH + lat + soil_moisture + soil_var
   index_fung_resis ~ treatment + management + pH + lat + soil_moisture + soil_var
   index_fung_sens ~ treatment + management + pH + lat + soil_moisture + soil_var
   
   # extrinsic covariances
   management ~~ 0*treatment
   pH ~~ soil_var
   pH ~~ soil_moisture
   '

fit_model_sem_F0 <- sem(model_sem_F,
                        data = dat_sem %>%
                          filter(day == "0"),
                        std.lv = TRUE, std.ov = TRUE,
                        fixed.x = FALSE)
summary(fit_model_sem_F0, fit.measures = TRUE, standardized = T, rsq = T)  


# bootstrap for bootstrap test statistics
F0_boot <- bootstrapLavaan(fit_model_sem_F0, R = 1000L, type = "bollen.stine", verbose= FALSE,
                           FUN = fitMeasures, fit.measures = c("chisq"))
# compute a bootstrap based p-value
F0_fit_orig <- fitMeasures(fit_model_sem_F0, "chisq")
F0_boot_pvalue <- length(which(F0_boot > F0_fit_orig))/length(F0_boot)


## FUNGI DAY 60 ----
fit_model_sem_F60 <- sem(model_sem_F,
                         data = dat_sem %>%
                           filter(day == "60"),
                         std.lv = TRUE, std.ov = TRUE,
                         fixed.x = FALSE)
summary(fit_model_sem_F60, fit.measures = TRUE, standardized = T, rsq = T)  

# bootstrap for bootstrap test statistics
F60_boot <- bootstrapLavaan(fit_model_sem_F60, R = 1000L, type = "bollen.stine", verbose= FALSE,
                            FUN = fitMeasures, fit.measures = c("chisq"))
# compute a bootstrap based p-value
F60_fit_orig <- fitMeasures(fit_model_sem_F60, "chisq")
F60_boot_pvalue <- length(which(F60_boot > F60_fit_orig))/length(F60_boot)


# Mixed effects models of management effects on indices ----
#library(nlme)
vfIdent <- varIdent(form = ~1|fieldID)
# BACTERIA, ratio of sensitive to opportunistic, both timepoints, control plots
mod_lme_bact_ratio <- lme(log(opp_sens_ratio_bact) ~   management,
                          random = ~1|region/site/management, 
                          weights = vfIdent,
                          method = "REML",  # must use REML for final estimation
                          data = dat_indices %>% # aggregate the two timepoints
                            filter(treatment == "C") %>% # only in control plots
                            group_by(region, site, fieldID, management, 
                                     treatment, plot) %>%
                            summarise(across(where(is.numeric), 
                                             ~mean(.x, na.rm = TRUE))) %>%
                            ungroup(), 
                          na.action = na.omit)
plot(mod_lme_bact_ratio) 
summary(mod_lme_bact_ratio) # effect is 0.28, p = 0.0001

# BACTERIA, ratio of resilient to not resilient, both timepoints, control plots
mod_lme_bact_ratio_resil <- lme(log(resil_not_resil_ratio_bact) ~   management,
                                random = ~1|region/site/management, 
                                weights = vfIdent,
                                method = "REML",  # must use REML for final estimation
                                data = dat_indices %>% # aggregate the two timepoints
                                  filter(treatment == "C") %>% # only in control plots
                                  group_by(region, site, fieldID, management, 
                                           treatment, plot) %>%
                                  summarise(across(where(is.numeric), 
                                                   ~mean(.x, na.rm = TRUE))) %>%
                                  ungroup(), 
                                na.action = na.omit)
plot(mod_lme_bact_ratio_resil) 
summary(mod_lme_bact_ratio_resil) # effect is 0.13, p = 0.0037


# FUNGI, ratio of sensitive to opportunistic, both timepoints, control plots
mod_lme_fung_ratio <- lme(log(opp_sens_ratio_fung) ~   management,
                          random = ~1|region/site/management, 
                          weights = vfIdent,
                          method = "REML",  # must use REML for final estimation
                          data = dat_indices %>% # aggregate the two timepoints
                            filter(treatment == "C") %>% # only in control plots
                            group_by(region, site, fieldID, management, 
                                     treatment, plot) %>%
                            summarise(across(where(is.numeric), 
                                             ~mean(.x, na.rm = TRUE))) %>%
                            ungroup(), 
                          na.action = na.omit)
plot(mod_lme_fung_ratio)
summary(mod_lme_fung_ratio) # effect is -0.41, p = 0.0003

# FUNGI, ratio of resilient to not resilient
mod_lme_fung_ratio_resil <- lme(log(resil_not_resil_ratio_fung) ~   management,
                                random = ~1|region/site/management, 
                                weights = varIdent(form = ~1|region*site),
                                method = "REML",  # must use REML for final estimation
                                data = dat_indices %>% # aggregate the two timepoints
                                  filter(treatment == "C") %>% # only in control plots
                                  group_by(region, site, fieldID, management, 
                                           treatment, plot) %>%
                                  summarise(across(where(is.numeric),
                                                   ~mean(.x, na.rm = TRUE))) %>%
                                  ungroup(), 
                                na.action = na.omit)
plot(mod_lme_fung_ratio_resil)
summary(mod_lme_fung_ratio_resil) # effect is -0.377, p = 0.0451


