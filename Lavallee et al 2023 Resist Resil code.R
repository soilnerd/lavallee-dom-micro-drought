###     Code for Lavallee et al. (2023) Land management shapes drought    ###
###     responses of dominant soil microbial taxa across grasslands       ###
###     Resistance and Resilience drought response models                 ###


# Load libraries ----
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(purrr)
library(magrittr)

# Split OTU tables by sampling day ----
OTU_tab_dom_bact_d0 <- OTU_tab_dom_bact %>%
  filter(day == "0")

OTU_tab_dom_bact_d60 <- OTU_tab_dom_bact %>%
  filter(day == "60")

OTU_tab_dom_fung_d0 <- OTU_tab_dom_fung %>%
  filter(day == "0")

OTU_tab_dom_fung_d60 <- OTU_tab_dom_fung %>%
  filter(day == "60")

# BACTERIA, RESISTANCE ----
## Define model functions ---- 
apply_glmmTMB_poiss_no_disp <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_bact_d0,
          family = poisson)
}

apply_glmmTMB_nbinom_no_disp <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_bact_d0,
          family = nbinom1)
}

apply_glmmTMB_nbinom_disp <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          dispformula = ~ site, 
          data = OTU_tab_dom_bact_d0,
          family = nbinom1)
}

apply_glmmTMB_binomial <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          #dispformula = ~ site, 
          data = OTU_tab_dom_bact_d0 %>% 
            mutate(across(where(is.numeric), ~1 * (. > 0))),
          family = binomial)
}


## POISSON MODELS, NO DISPERSION PARAMETER ----
### First pass to see which OTUs throw warnings
### ! NOTE ! THIS MAY NOT WORK PROPERLY IF WARNINGS OTHER THAN MODEL-CONVERGENCE 
### WARNINGS ARE GIVEN! #Change 'warning' to conditionMessage(w) and inspect to 
### make sure warnings are convergence-related and not something else!
models_poisson_no_disp <-
  OTU_tab_dom_bact_d0 %>%
  select(-c(sample:day)) %>%
  names() %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_poiss_no_disp(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') }) #Change 'warning' to conditionMessage(w) and inspect to make sure warnings are convergence-related and not something else!
    })

### Pull OTUs with warnings and no warnings
OTUs_dom_bact_poisson_no_disp_warn <- models_poisson_no_disp %>%
  extract(. == 'warning') %>% 
  names()

OTUs_dom_bact_poisson_no_disp_no_warn <-  OTU_tab_dom_bact_d0 %>%
  select(-c(sample:day)) %>%
  names() %>%
  extract(!(. %in% OTUs_dom_bact_poisson_no_disp_warn))

# Simulate residuals for models where there weren't warnings 
models_resids_poisson_no_disp <-  models_poisson_no_disp[OTUs_dom_bact_poisson_no_disp_no_warn] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests, get p values
plot.new()
model_ps_poisson_no_disp <-
  models_resids_poisson_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_dom_bact_poisson_no_disp_violate <- model_ps_poisson_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05)) 

OTUs_dom_bact_poisson_no_disp_good <- model_ps_poisson_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) > 0.05)) %>%
  select(id) %>%
  pull()

# Produce the summary stats for models that pass
models_summaries_poisson_no_disp <-  
  models_poisson_no_disp %>%
  extract(all_of(OTUs_dom_bact_poisson_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_poisson_no_disp <- models_summaries_poisson_no_disp %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")),
         model = "Poisson with no dispersion")


## NEGATIVE BINOMIAL WITH NO DISPERSION PARAMETER, for all that didn't meet assumptions or gave warnings ----   
# Run models and check models for warnings
models_nbinom_no_disp <-
  c(OTUs_dom_bact_poisson_no_disp_violate %>% 
      pull(id),
    OTUs_dom_bact_poisson_no_disp_warn) %>% 
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_nbinom_no_disp(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_dom_bact_nbinom_no_disp_warn <- models_nbinom_no_disp %>%
  extract(. == 'warning') %>% 
  names()

OTUs_dom_bact_nbinom_no_disp_no_warn <- c(OTUs_dom_bact_poisson_no_disp_violate %>%
                                            pull(id),
                                          OTUs_dom_bact_poisson_no_disp_warn) %>%
  extract(!(. %in% OTUs_dom_bact_nbinom_no_disp_warn))

# Simulate residuals for models where there weren't warnings
models_resids_nbinom_no_disp <-  models_nbinom_no_disp[all_of(OTUs_dom_bact_nbinom_no_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_nbinom_no_disp <-
  models_resids_nbinom_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_dom_bact_nbinom_no_disp_violate <- model_ps_nbinom_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))    

OTUs_dom_bact_nbinom_no_disp_good <- model_ps_nbinom_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  select(id) %>%
  pull()


# Produce the summary stats for models that pass
models_summaries_nbinom_no_disp <-  
  models_nbinom_no_disp %>%
  extract(all_of(OTUs_dom_bact_nbinom_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_nbinom_no_disp <- models_summaries_nbinom_no_disp %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")),
         model = "Neg binomial with no dispersion")


## NEGATIVE BINOMIAL WITH DISPERSION PARAMETER ----    
models_nbinom_disp <-
  c(OTUs_dom_bact_nbinom_no_disp_violate %>%
      filter(!(id %in% c("Bac_OTU3166", "Bac_OTU12563"))) %>%
      pull(id),
    OTUs_dom_bact_nbinom_no_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_nbinom_disp(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_dom_bact_nbinom_disp_warn <- models_nbinom_disp %>%
  extract(. == 'warning') %>% 
  names()


OTUs_dom_bact_nbinom_disp_no_warn <- c(OTUs_dom_bact_nbinom_no_disp_violate %>%
                                         filter(!(id %in% c("Bac_OTU3166", "Bac_OTU12563"))) %>%
                                         pull(id),
                                       OTUs_dom_bact_nbinom_no_disp_warn) %>%
  extract(!(. %in% OTUs_dom_bact_nbinom_disp_warn))



# Simulate residuals for models where there weren't warnings 
models_resids_nbinom_disp <-  models_nbinom_disp[all_of(OTUs_dom_bact_nbinom_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_nbinom_disp <-
  models_resids_nbinom_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_dom_bact_nbinom_disp_violate <- model_ps_nbinom_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))    

OTUs_dom_bact_nbinom_disp_good <- model_ps_nbinom_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  pull(id)


# Produce summary stats for models that pass
models_summaries_nbinom_disp <-  
  models_nbinom_disp %>%
  extract(all_of(OTUs_dom_bact_nbinom_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_nbinom_disp <- models_summaries_nbinom_disp %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")),
         model = "Neg binomial with dispersion")

## BINOMIAL MODELS ----
models_binomial <-
  c(OTUs_dom_bact_nbinom_disp_violate %>%
      pull(id),
    OTUs_dom_bact_nbinom_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_binomial(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_dom_bact_binomial_warn <- models_binomial %>%
  extract(. == 'warning') %>% 
  names()


OTUs_dom_bact_binomial_no_warn <-c(OTUs_dom_bact_nbinom_disp_violate %>%
                                     pull(id),
                                   OTUs_dom_bact_nbinom_disp_warn) %>%
  extract(!(. %in% OTUs_dom_bact_binomial_warn))

# Simulate residuals for models where there weren't warnings
models_resids_binomial <-  models_binomial[all_of(OTUs_dom_bact_binomial_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_binomial <-
  models_resids_binomial %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d0$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_dom_bact_binomial_violate <- model_ps_binomial %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))   

OTUs_dom_bact_binomial_good <- model_ps_binomial %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  pull(id)


# Produce the summary stats for models that pass
models_summaries_binomial <-  
  models_binomial %>%
  extract(all_of(OTUs_dom_bact_binomial_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_binomial <- models_summaries_binomial %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")),
         model = "Neg binomial with no dispersion")

## Custom models for last OTUs ----
# Make a data frame to keep track of models and stats
OTUs_drought_class_custom <- data.frame(matrix(ncol = 6, nrow=0, 
                                               dimnames = list(NULL, c("id", "model", "Estimate", "Std.Error", "z value", "Pr(>|z|)"
                                               ))))

# One that didn't run before is Bac_OTU3166
mod_Bac_OTU3166 <- glmmTMB(Bac_OTU3166 ~ trt*region + (1|site/field),
                           #dispformula = ~ region,
                           ziformula = ~ region, # tried with site and region, but only region works to avoid violating levene test.
                           data =
                             OTU_tab_dom_bact_d0,
                           family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU3166)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(data.frame(id = "Bac_OTU3166", 
                                                        model = "~ trt*region + (1|site/field), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                             summary(mod_Bac_OTU3166)$coefficients$cond[c(2,8,14,20)] %>%
                                               data.frame() %>%
                                               rename(value = names(.)[1]) %>%
                                               mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                               pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))
# One that didn't run before is Bac_OTU12563
mod_Bac_OTU12563 <- glmmTMB(Bac_OTU12563 ~ trt*region + (1|site/field),
                           #dispformula = ~ region,
                           ziformula = ~ region, # tried with site and region, but only region works to avoid violating levene test.
                           data =
                             OTU_tab_dom_bact_d0,
                           family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU12563)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(data.frame(id = "Bac_OTU12563", 
                                                        model = "~ trt*region + (1|site/field), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                             summary(mod_Bac_OTU3166)$coefficients$cond[c(2,8,14,20)] %>%
                                               data.frame() %>%
                                               rename(value = names(.)[1]) %>%
                                               mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                               pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


# Look at those that violated assumptions
OTUs_dom_bact_binomial_violate%>% pull(id)

mod_Bac_OTU3156 <- glmmTMB(Bac_OTU3156 ~ trt*region + (1|site/field),
                           ziformula = ~ region, 
                           data =   
                             OTU_tab_dom_bact_d0,
                           family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU3156)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(data.frame(id = "Bac_OTU3156", 
                                                        model = "~ trt*region + (1|site/field), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                             summary(mod_Bac_OTU3156)$coefficients$cond[c(2,8,14,20)] %>%
                                               data.frame() %>%
                                               rename(value = names(.)[1]) %>%
                                               mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                               pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#Bac_OTU186
mod_Bac_OTU186 <- glmmTMB(Bac_OTU186 ~ trt*region + (1|site/field),
                           ziformula = ~ region, 
                           data =   
                             OTU_tab_dom_bact_d0,
                           family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU186)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(data.frame(id = "Bac_OTU186", 
                                                        model = "~ trt*region + (1|site/field), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                             summary(mod_Bac_OTU186)$coefficients$cond[c(2,8,14,20)] %>%
                                               data.frame() %>%
                                               rename(value = names(.)[1]) %>%
                                               mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                               pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU1965"
mod_Bac_OTU1965 <- glmmTMB(Bac_OTU1965 ~ trt*region + (1|site/field),
                          ziformula = ~ region, # still violating levene test but not terribly
                          data =   
                            OTU_tab_dom_bact_d0,
                          family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU1965)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(data.frame(id = "Bac_OTU1965", 
                                                        model = "~ trt*region + (1|site/field), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                             summary(mod_Bac_OTU1965)$coefficients$cond[c(2,8,14,20)] %>%
                                               data.frame() %>%
                                               rename(value = names(.)[1]) %>%
                                               mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                               pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU317" 
mod_Bac_OTU317 <- glmmTMB(Bac_OTU317 ~ trt*region + (1|field),
                           dispformula = ~ region,
                           data =   
                             OTU_tab_dom_bact_d0,
                           family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU317)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(data.frame(id = "Bac_OTU317", 
                                                        model = "~ trt*region + (1|field), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                             summary(mod_Bac_OTU317)$coefficients$cond[c(2,8,14,20)] %>%
                                               data.frame() %>%
                                               rename(value = names(.)[1]) %>%
                                               mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                               pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


## Go through those that had warnings
OTUs_dom_bact_binomial_warn

mod_Bac_OTU15057 <-  glmmTMB(Bac_OTU15057 ~ trt*region + (1|site/field),
                             ziformula = ~ site, # tried with site and region, but only site works without warning
                             data =   
                               OTU_tab_dom_bact_d0,
                             family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU15057)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU15057", 
                                                model = "~ trt*region + (1|site/field), ziformula = ~ site, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU15057)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))
#Bac_OTU15598
mod_Bac_OTU15598 <-  glmmTMB(Bac_OTU15598 ~ trt + (1|site/field),
                             dispformula = ~ site,
                             ziformula = ~ site,
                             data =   
                               OTU_tab_dom_bact_d0,
                             family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU15598)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU15598", 
                                                model = "~ trt + (1|site/field), dispformula = ~ site, ziformula = ~ site, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU15598)$coefficients$cond[c(2,4,6,8)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#Bac_OTU5895
mod_Bac_OTU5895 <-  glmmTMB(Bac_OTU5895 ~ trt*region + (1|site/field),
                             dispformula = ~ region,
                             data =   
                               OTU_tab_dom_bact_d0,
                             family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU5895)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU5895", 
                                                model = "~ trt + (1|site/field), dispformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU5895)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))



# Bac_OTU3160
mod_Bac_OTU3160 <-  glmmTMB(Bac_OTU3160 ~ trt*region + (1|site),
                            ziformula = ~ region,
                            data = OTU_tab_dom_bact_d0,
                            family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU3160)
plot(res)


# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU3160", 
                                                model = "~ trt*region + (1|site), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU3160)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


# Bac_OTU57
mod_Bac_OTU57 <-  glmmTMB(Bac_OTU57 ~ trt*region + (1|site),
                          ziformula = ~ region,
                          data = OTU_tab_dom_bact_d0,
                          family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU57)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU57", 
                                                model = "~ trt*region + (1|site), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU57)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))




# Bac_OTU323
mod_Bac_OTU323 <-  glmmTMB(Bac_OTU323 ~ trt*region + (1|site/field),
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU323)
plot(res)


# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU323", 
                                                model = "~ trt*region + (1|site/field), data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU323)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

# Bac_OTU9974
mod_Bac_OTU9974 <- glmmTMB(Bac_OTU9974 ~ trt*region + (1|site),
                           dispformula = ~ site,
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU9974)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU9974", 
                                                model = "~ trt*region + (1|site), dispformula = ~ site, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU9974)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


# Bac_OTU183
mod_Bac_OTU183 <- glmmTMB(Bac_OTU183 ~ trt*region + (1|site),
                          ziformula = ~ region,
                          data = OTU_tab_dom_bact_d0,
                          family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU183)
plot(res)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU183", 
                                                model = "~ trt*region + (1|site), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU183)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))



# Bac_OTU17814 -- this one has an issue with residuals. It's not responding either way so leave it here.
mod_Bac_OTU17814 <- glmmTMB(Bac_OTU17814 ~ trt*region + (1|site/field),
                            dispformula = ~ region,
                            data = OTU_tab_dom_bact_d0,
                            family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU17814)
plot(res)
testZeroInflation(mod_Bac_OTU17814)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU17814", 
                                                model = "~ trt*region + (1|site/field), dispformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU17814)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


# "Bac_OTU1872"  
mod_Bac_OTU1872 <- glmmTMB(Bac_OTU1872 ~ trt*region + (1|site/field),
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU1872)
plot(res)


# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU1872", 
                                                model = "~ trt*region + (1|site/field), data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU1872)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU7383"  
mod_Bac_OTU7383 <- glmmTMB(Bac_OTU7383 ~ trt*region + (1|site/field),
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU7383)
plot(res)


# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU7383", 
                                                model = "~ trt*region + (1|site/field), data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU7383)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU7388"  
mod_Bac_OTU7388 <- glmmTMB(Bac_OTU7388 ~ trt*region + (1|site),
                           dispformula = ~ region,
                           ziformula = ~ region, 
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU7388)
plot(res)
testZeroInflation(mod_Bac_OTU7388)
testDispersion(mod_Bac_OTU7388)

# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU7388", 
                                                model = "~ trt*region + (1|site), dispformula = ~ region, ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU7388)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU7397"  
mod_Bac_OTU7397 <- glmmTMB(Bac_OTU7397 ~ trt*region + (1|site),
                           ziformula  = ~ region, 
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU7397)
plot(res) 

# Add it to the data frame
OTUs_drought_class_custom <- bind_rows(OTUs_drought_class_custom,
                                       bind_cols(
                                         data.frame(id = "Bac_OTU7397", 
                                                    model = "~ trt*region + (1|site), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                         summary(mod_Bac_OTU7397)$coefficients$cond[c(2,8,14,20)] %>%
                                           data.frame() %>%
                                           rename(value = names(.)[1]) %>%
                                           mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                           pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU793"   
mod_Bac_OTU793 <- glmmTMB(Bac_OTU793 ~ trt*region + (1|site),
                          dispformula  = ~ region, 
                          data = OTU_tab_dom_bact_d0,
                          family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU793)
plot(res) 


# Add it to the data frame
OTUs_drought_class_custom <- bind_rows(OTUs_drought_class_custom,
                                       bind_cols(
                                         data.frame(id = "Bac_OTU793", 
                                                    model = "~ trt*region + (1|site), dispformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                         summary(mod_Bac_OTU793)$coefficients$cond[c(2,8,14,20)] %>%
                                           data.frame() %>%
                                           rename(value = names(.)[1]) %>%
                                           mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                           pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU11536" 
mod_Bac_OTU11536 <- glmmTMB(Bac_OTU11536 ~ trt*region + (1|site),
                            #dispformula  = ~ region, 
                            ziformula = ~ region,
                            data = OTU_tab_dom_bact_d0,
                            family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU11536)
plot(res) 

# Add it to the data frame
OTUs_drought_class_custom <- bind_rows(OTUs_drought_class_custom,
                                       bind_cols(
                                         data.frame(id = "Bac_OTU11536", 
                                                    model = "~ trt*region + (1|site), ziformula = ~ region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                         summary(mod_Bac_OTU11536)$coefficients$cond[c(2,8,14,20)] %>%
                                           data.frame() %>%
                                           rename(value = names(.)[1]) %>%
                                           mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                           pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


#"Bac_OTU7382" 
mod_Bac_OTU7382 <- glmmTMB(Bac_OTU7382 ~ trt*region + (1|site),
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU7382)
plot(res) 


# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU7382", 
                                                model = "~ trt*region + (1|site), data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU7382)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU9971" 
mod_Bac_OTU9971 <- glmmTMB(Bac_OTU9971 ~ trt*region + (1|site/field),
                           dispformula = ~region,
                           data = OTU_tab_dom_bact_d0,
                           family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU9971)
plot(res) 
testZeroInflation(mod_Bac_OTU9971)


# Add it to the data frame
OTUs_drought_class_custom <- rbind(OTUs_drought_class_custom,
                                   bind_cols(
                                     data.frame(id = "Bac_OTU9971", 
                                                model = "~ trt*region + (1|site/field), dispformula = ~region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                     summary(mod_Bac_OTU9971)$coefficients$cond[c(2,8,14,20)] %>%
                                       data.frame() %>%
                                       rename(value = names(.)[1]) %>%
                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU17732" 
mod_Bac_OTU17732 <- glmmTMB(Bac_OTU17732 ~ trt*region + (1|site),
                            dispformula = ~region,
                            data = OTU_tab_dom_bact_d0,
                            family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU17732)
plot(res) 
testZeroInflation(mod_Bac_OTU17732)
testDispersion(mod_Bac_OTU17732)


# Add it to the data frame
OTUs_drought_class_custom <- bind_rows(OTUs_drought_class_custom,
                                       bind_cols(
                                         data.frame(id = "Bac_OTU17732", 
                                                    model = "~ trt*region + (1|site), dispformula = ~region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                         summary(mod_Bac_OTU17732)$coefficients$cond[c(2,8,14,20)] %>%
                                           data.frame() %>%
                                           rename(value = names(.)[1]) %>%
                                           mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                           pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


#"Bac_OTU18641" 
mod_Bac_OTU18641 <- glmmTMB(Bac_OTU18641 ~ trt*region + (1|site/field),
                            dispformula = ~region,
                            data = OTU_tab_dom_bact_d0,
                            family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU18641)
plot(res) 


# Add it to the data frame
OTUs_drought_class_custom <- bind_rows(OTUs_drought_class_custom,
                                       bind_cols(
                                         data.frame(id = "Bac_OTU18641", 
                                                    model = "~ trt*region + (1|site/field), dispformula = ~region, data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                         summary(mod_Bac_OTU18641)$coefficients$cond[c(2,8,14,20)] %>%
                                           data.frame() %>%
                                           rename(value = names(.)[1]) %>%
                                           mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                           pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


#"Bac_OTU14995"
mod_Bac_OTU14995 <- glmmTMB(Bac_OTU14995 ~ trt*region + (1|site/field),
                            data = OTU_tab_dom_bact_d0,
                            family = nbinom1)
res  <- simulateResiduals(mod_Bac_OTU14995)
plot(res) 

# Add it to the data frame
OTUs_drought_class_custom <- bind_rows(OTUs_drought_class_custom,
                                       bind_cols(
                                         data.frame(id = "Bac_OTU14995", 
                                                    model = "~ trt*region + (1|site/field), data = OTU_tab_dom_bact_d0, family = nbinom1"),
                                         summary(mod_Bac_OTU14995)$coefficients$cond[c(2,8,14,20)] %>%
                                           data.frame() %>%
                                           rename(value = names(.)[1]) %>%
                                           mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                           pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))



## Put results together (Resistance) ----
OTUs_drought_class_custom <- OTUs_drought_class_custom %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")))
OTUs_drought_class_all <- bind_rows(OTUs_drought_class_poisson_no_disp,
                                    OTUs_drought_class_nbinom_no_disp,
                                    OTUs_drought_class_nbinom_disp,
                                    OTUs_drought_class_binomial,
                                    OTUs_drought_class_custom)


sum_drought_class <- OTUs_drought_class_all %>%
  group_by(Resistance_class) %>%
  summarize(count = n())

# The following should not produce any results. If it does, there are duplicates
# (implies that models ran/passed this time that didn't run previously, likely 
# due to differences in package versions)
OTUs_drought_class_all %>%
  group_by(id) %>%
  summarise(n=n()) %>%
  filter(n > 1) %>%
  pull(id)

 
# BACTERIA, RESILIENCE ----
## model functions ---- 
apply_glmmTMB_poiss_no_disp_60 <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_bact_d60,
          family = poisson)
}

apply_glmmTMB_nbinom_no_disp_60 <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_bact_d60,
          family = nbinom1)
}

apply_glmmTMB_nbinom_disp_60 <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          dispformula = ~ site, 
          data = OTU_tab_dom_bact_d60,
          family = nbinom1)
}

apply_glmmTMB_binomial_60 <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_bact_d60 %>% 
            mutate(across(where(is.numeric), ~1 * (. > 0))),
          family = binomial)
}

## POISSON MODELS WITH NO DISPERSION ----
## First pass to see which OTUs throw warnings
models_60_poisson_no_disp <-
  OTUs_drought_class_all %>%
  filter(Resistance_class != "Resistant") %>%
  pull(id) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_poiss_no_disp_60(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

### Pull OTUs with warnings and no warnings
OTUs_bact_60_poisson_no_disp_warn <- models_60_poisson_no_disp %>%
  extract(. == 'warning') %>% 
  names()

OTUs_bact_60_poisson_no_disp_no_warn <-     OTUs_drought_class_all %>%
  filter(Resistance_class != "Resistant") %>%
  pull(id)  %>%
  extract(!(. %in% OTUs_bact_60_poisson_no_disp_warn))

# Simulate residuals for models where there weren't warnings
models_resids_60_poisson_no_disp <-  models_60_poisson_no_disp[OTUs_bact_60_poisson_no_disp_no_warn] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_60_poisson_no_disp <-
  models_resids_60_poisson_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_bact_60_poisson_no_disp_violate <- model_ps_60_poisson_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05)) 

OTUs_bact_60_poisson_no_disp_good <- model_ps_60_poisson_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  select(id) %>%
  pull()

# Produce the summary stats for models that pass
models_summaries_60_poisson_no_disp <-  
  models_60_poisson_no_disp %>%
  extract(all_of(OTUs_bact_60_poisson_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_60_poisson_no_disp <- models_summaries_60_poisson_no_disp %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient",
                                    "Not Resilient"),
         model = "Poisson with no dispersion")

## NEGATIVE BINOMIAL WITH NO DISP OR ZI ----
# Run models and check models for warnings
models_60_nbinom_no_disp <-
  c(OTUs_bact_60_poisson_no_disp_violate %>%
      pull(id),
    OTUs_bact_60_poisson_no_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_nbinom_no_disp_60(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_bact_60_nbinom_no_disp_warn <- models_60_nbinom_no_disp %>%
  extract(. == 'warning') %>% 
  names()


OTUs_bact_60_nbinom_no_disp_no_warn <- c(OTUs_bact_60_poisson_no_disp_violate %>%
                                           pull(id),
                                         OTUs_bact_60_poisson_no_disp_warn) %>%
  extract(!(. %in% OTUs_bact_60_nbinom_no_disp_warn))

# Simulate residuals for models where there weren't warnings
models_resids_60_nbinom_no_disp <-  models_60_nbinom_no_disp[all_of(OTUs_bact_60_nbinom_no_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_60_nbinom_no_disp <-
  models_resids_60_nbinom_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter any OTUs where p values are <0.05
OTUs_bact_60_nbinom_no_disp_violate <- model_ps_60_nbinom_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))    

OTUs_bact_60_nbinom_no_disp_good <- model_ps_60_nbinom_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  select(id) %>%
  pull()


# Produce the summary stats for models that pass
models_summaries_60_nbinom_no_disp <-  
  models_60_nbinom_no_disp %>%
  extract(all_of(OTUs_bact_60_nbinom_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_60_nbinom_no_disp <- models_summaries_60_nbinom_no_disp %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient",
                                    "Not Resilient"),
         model = "Neg binomial with no dispersion")

## NEGATIVE BINOMIAL WITH DISPERSION PARAMETER ----
models_60_nbinom_disp <-
  c(OTUs_bact_60_nbinom_no_disp_violate %>%
      pull(id), 
    OTUs_bact_60_nbinom_no_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_nbinom_disp_60(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_bact_60_nbinom_disp_warn <- models_60_nbinom_disp %>%
  extract(. == 'warning') %>% 
  names()


OTUs_bact_60_nbinom_disp_no_warn <- c(OTUs_bact_60_nbinom_no_disp_violate %>% 
                                        pull(id), 
                                      OTUs_bact_60_nbinom_no_disp_warn) %>%
  extract(!(. %in% OTUs_bact_60_nbinom_disp_warn))


# Simulate residuals for models where there weren't warnings 
models_resids_60_nbinom_disp <-  models_60_nbinom_disp[all_of(OTUs_bact_60_nbinom_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_60_nbinom_disp <-
  models_resids_60_nbinom_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_bact_60_nbinom_disp_violate <- model_ps_60_nbinom_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))    

OTUs_bact_60_nbinom_disp_good <- model_ps_60_nbinom_disp %>% # only 4 that are good here.
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  pull(id)


# Produce the summary stats for models that pass
models_summaries_60_nbinom_disp <-  
  models_60_nbinom_disp %>%
  extract(all_of(OTUs_bact_60_nbinom_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")


OTUs_drought_class_60_nbinom_disp <- models_summaries_60_nbinom_disp %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient",
                                    "Not Resilient"),
         model = "Neg binomial with dispersion")


## BINOMIAL MODELS ----
models_60_binomial <-
  c(OTUs_bact_60_nbinom_disp_violate %>%
      pull(id),
    OTUs_bact_60_nbinom_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_binomial_60(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_bact_60_binomial_warn <- models_60_binomial %>%
  extract(. == 'warning') %>% 
  names()

OTUs_bact_60_binomial_no_warn <-c(OTUs_bact_60_nbinom_disp_violate %>%
                                    pull(id),
                                  OTUs_bact_60_nbinom_disp_warn) %>%
  extract(!(. %in% OTUs_bact_60_binomial_warn))

# Simulate residuals for models that pass
models_resids_60_binomial <-  models_60_binomial[all_of(OTUs_bact_60_binomial_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_60_binomial <-
  models_resids_60_binomial %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_bact_d60$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_bact_60_binomial_violate <- model_ps_60_binomial %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))   

OTUs_bact_60_binomial_good <- model_ps_60_binomial %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  pull(id)


# Produce the summary stats for models that pass
models_summaries_60_binomial <-  
  models_60_binomial %>%
  extract(all_of(OTUs_bact_60_binomial_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")


OTUs_drought_class_60_binomial <- models_summaries_60_binomial %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient",
                                    "Not Resilient"),
         model = "Binomial no disp no zi")

## Custom models for last 5 OTUs ----
OTUs_bact_60_binomial_violate %>% pull(id)
OTUs_bact_60_binomial_warn
#"Bac_OTU19912" "Bac_OTU15598" "Bac_OTU9776"  "Bac_OTU14995"
# Make data frame to keep track of models and stats
OTUs_drought_class_60_custom <- data.frame(matrix(ncol = 6, nrow=0, 
                                                  dimnames = list(NULL, c("id", "model", "Estimate", "Std.Error", "z value", "Pr(>|z|)"
                                                  ))))

# "Bac_OTU620" violates assumptions
mod_Bac_OTU620 <- glmmTMB(Bac_OTU620 ~ trt*region + (1|site),
                          #dispformula = ~ region,
                          ziformula = ~ region, # tried with site and region, but only region works to avoid violating levene test.
                          data =   
                            OTU_tab_dom_bact_d60,
                          family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU620)
plot(res)

# Add it to the data frame
OTUs_drought_class_60_custom <- rbind(OTUs_drought_class_60_custom,
                                      bind_cols(data.frame(id = "Bac_OTU620", 
                                                           model = "~ trt*region + (1|site), ziformula = ~ region, data = OTU_tab_dom_bact_d60, family = nbinom1"),
                                                summary(mod_Bac_OTU620)$coefficients$cond[c(2,8,14,20)] %>%
                                                  data.frame() %>%
                                                  rename(value = names(.)[1]) %>%
                                                  mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                  pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU19912"
mod_Bac_OTU19912 <- glmmTMB(Bac_OTU19912 ~ trt*region + (1|site/field),
                            data =
                              OTU_tab_dom_bact_d60,
                            family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU19912)
plot(res)

# Add it to the data frame
OTUs_drought_class_60_custom <- rbind(OTUs_drought_class_60_custom,
                                      bind_cols(data.frame(id = "Bac_OTU19912", 
                                                           model = "~ trt*region + (1|site/field), data = OTU_tab_dom_bact_d60, family = nbinom1"),
                                                summary(mod_Bac_OTU19912)$coefficients$cond[c(2,8,14,20)] %>%
                                                  data.frame() %>%
                                                  rename(value = names(.)[1]) %>%
                                                  mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                  pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

# "Bac_OTU15598" 
mod_Bac_OTU15598 <- glmmTMB(Bac_OTU15598 ~ trt*region + (1|site/field),
                            data =
                              OTU_tab_dom_bact_d60,
                            family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU15598)
plot(res)

# Add it to the data frame
OTUs_drought_class_60_custom <- rbind(OTUs_drought_class_60_custom,
                                      bind_cols(data.frame(id = "Bac_OTU15598", 
                                                           model = "~ trt*region + (1|site/field), data = OTU_tab_dom_bact_d60, family = nbinom1"),
                                                summary(mod_Bac_OTU15598)$coefficients$cond[c(2,8,14,20)] %>%
                                                  data.frame() %>%
                                                  rename(value = names(.)[1]) %>%
                                                  mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                  pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU9776"
mod_Bac_OTU9776 <- glmmTMB(Bac_OTU9776 ~ trt + (1|site),
                           dispformula = ~ region,
                           data =
                             OTU_tab_dom_bact_d60,
                           family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU9776)
plot(res)


# Add it to the data frame
OTUs_drought_class_60_custom <- rbind(OTUs_drought_class_60_custom,
                                      bind_cols(data.frame(id = "Bac_OTU9776", 
                                                           model = "~ trt*region + (1|site), dispformula = ~ region, data = OTU_tab_dom_bact_d60, family = nbinom1"),
                                                summary(mod_Bac_OTU9776)$coefficients$cond[c(2,4,6,8)] %>%
                                                  data.frame() %>%
                                                  rename(value = names(.)[1]) %>%
                                                  mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                  pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Bac_OTU14995"
mod_Bac_OTU14995 <- glmmTMB(Bac_OTU14995 ~ trt*region + (1|site/field),
                            data =
                              OTU_tab_dom_bact_d60,
                            family = nbinom1)
res <- simulateResiduals(mod_Bac_OTU14995)
plot(res)

# Add it to the data frame
OTUs_drought_class_60_custom <- rbind(OTUs_drought_class_60_custom,
                                      bind_cols(data.frame(id = "Bac_OTU14995", 
                                                           model = "~ trt*region + (1|site/field), data = OTU_tab_dom_bact_d60, family = nbinom1"),
                                                summary(mod_Bac_OTU14995)$coefficients$cond[c(2,8,14,20)] %>%
                                                  data.frame() %>%
                                                  rename(value = names(.)[1]) %>%
                                                  mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                  pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

OTUs_drought_class_60_custom <- OTUs_drought_class_60_custom %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient","Not Resilient"),
         model = "custom")



## Put results together (Resilience) ----
OTUs_drought_class_60_all <- bind_rows(OTUs_drought_class_60_poisson_no_disp,
                                       OTUs_drought_class_60_nbinom_no_disp,
                                       OTUs_drought_class_60_nbinom_disp,
                                       OTUs_drought_class_60_binomial,
                                       OTUs_drought_class_60_custom) %>%
  rename_with(~paste0(., "_d60"), c("trtD_Estimate":"trtD_Pr(>|z|)", "model"))



sum_drought_class_60 <- OTUs_drought_class_60_all %>%
  group_by(Resilience_class) %>%
  summarize(count = n())



# FUNGI, RESISTANCE ----
## Define model structures ----
apply_glmmTMB_poiss_fung <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_fung_d0,
          family = poisson)
}

apply_glmmTMB_nbinom_no_disp_fung <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_fung_d0,
          family = nbinom1)
}

apply_glmmTMB_binomial_fung <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_fung_d0 %>% 
            mutate(across(where(is.numeric), ~1 * (. > 0))),
          family = binomial)
}

## POISSON MODELS, NO DISPERSION PARAMETER ----
### First pass to see which OTUs throw warnings
models_fung_poisson_no_disp <-
  OTU_tab_dom_fung_d0 %>%
  select(-c(sample:day)) %>%
  names() %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_poiss_fung(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

### Pull OTUs with warnings and no warnings
OTUs_dom_fung_poisson_no_disp_warn <- models_fung_poisson_no_disp %>%
  extract(. == 'warning') %>% 
  names()

OTUs_dom_fung_poisson_no_disp_no_warn <-  OTU_tab_dom_fung_d0 %>%
  select(-c(sample:day)) %>%
  names()%>%
  extract(!(. %in% OTUs_dom_fung_poisson_no_disp_warn))

# Simulate residuals for models where there weren't warnings
models_resids_fung_poisson_no_disp <-  models_fung_poisson_no_disp[all_of(OTUs_dom_fung_poisson_no_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_fung_poisson_no_disp <-
  models_resids_fung_poisson_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_fung_d0$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_fung_d0$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_dom_fung_poisson_no_disp_violate <- model_ps_fung_poisson_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05)) 

OTUs_dom_fung_poisson_no_disp_good <- model_ps_fung_poisson_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) > 0.05)) %>%
  select(id) %>%
  pull()

# Produce the summary stats for models that pass
models_summaries_fung_poisson_no_disp <-  
  models_fung_poisson_no_disp %>%
  extract(all_of(OTUs_dom_fung_poisson_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_fung_poisson_no_disp <- models_summaries_fung_poisson_no_disp %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")),
         model = "Poisson with no dispersion")

## NEGATIVE BINOMIAL WITH NO DISPERSION PARAMETER, for all that didn't meet assumptions or gave warnings ----   
# Run models and check models for warnings
models_fung_nbinom_no_disp <-
  c(OTUs_dom_fung_poisson_no_disp_violate %>%
      pull(id),
    OTUs_dom_fung_poisson_no_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_nbinom_no_disp_fung(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_dom_fung_nbinom_no_disp_warn <- models_fung_nbinom_no_disp %>%
  extract(. == 'warning') %>% 
  names()

OTUs_dom_fung_nbinom_no_disp_no_warn <-   c(OTUs_dom_fung_poisson_no_disp_violate %>%
                                              pull(id),
                                            OTUs_dom_fung_poisson_no_disp_warn) %>%
  extract(!(. %in% OTUs_dom_fung_nbinom_no_disp_warn))

# Simulate residuals for models where there weren't warnings
models_resids_fung_nbinom_no_disp <-  models_fung_nbinom_no_disp[all_of(OTUs_dom_fung_nbinom_no_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests 
model_ps_fung_nbinom_no_disp <-
  models_resids_fung_nbinom_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_fung_d0$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_fung_d0$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_dom_fung_nbinom_no_disp_violate <- model_ps_fung_nbinom_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))    

OTUs_dom_fung_nbinom_no_disp_good <- model_ps_fung_nbinom_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  select(id) %>%
  pull()


# Produce the summary stats for models that pass
models_summaries_fung_nbinom_no_disp <-  
  models_fung_nbinom_no_disp %>%
  extract(all_of(OTUs_dom_fung_nbinom_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_fung_nbinom_no_disp <- models_summaries_fung_nbinom_no_disp %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")),
         model = "Neg binomial with no dispersion")

## BINOMIAL MODELS ----
models_fung_binomial <-
  c(OTUs_dom_fung_nbinom_no_disp_violate %>%
      pull(id),
    OTUs_dom_fung_nbinom_no_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_binomial_fung(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_dom_fung_binomial_warn <- models_fung_binomial %>%
  extract(. == 'warning') %>% 
  names()

OTUs_dom_fung_binomial_no_warn <- c(OTUs_dom_fung_nbinom_no_disp_violate %>% # 23 total
                                      pull(id),
                                    OTUs_dom_fung_nbinom_no_disp_warn) %>%
  extract(!(. %in% OTUs_dom_fung_binomial_warn))

# Simulate residuals for models where there weren't warnings
models_resids_fung_binomial <-  models_fung_binomial[all_of(OTUs_dom_fung_binomial_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_fung_binomial <-
  models_resids_fung_binomial %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_fung_d0$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_fung_d0$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_dom_fung_binomial_violate <- model_ps_fung_binomial %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))   

OTUs_dom_fung_binomial_good <- model_ps_fung_binomial %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  pull(id)


# NProduce the summary stats for models that pass
models_summaries_fung_binomial <-  
  models_fung_binomial %>%
  extract(all_of(OTUs_dom_fung_binomial_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_fung_binomial <- models_summaries_fung_binomial %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")),
         model = "Neg binomial with no dispersion")

## Custom models for last OTUs ----
OTUs_dom_fung_binomial_warn 
OTUs_dom_fung_binomial_violate%>% pull(id)

# Make data frame to keep track of models and stats
OTUs_drought_class_fung_custom <- data.frame(matrix(ncol = 6, nrow=0, 
                                                    dimnames = list(NULL, c("id", "model", "Estimate", "Std.Error", "z value", "Pr(>|z|)"
                                                    ))))

# One left with warning  is Fun_OTU12108
mod_Fun_OTU12108 <- glmmTMB(Fun_OTU12108 ~ trt*region + (1|site/field),
                            dispformula = ~ trt*region, 
                            data =   
                              OTU_tab_dom_fung_d0,
                            family = nbinom1)
res <- simulateResiduals(mod_Fun_OTU12108)
plot(res)
testZeroInflation(res)

# Add it to the data frame
OTUs_drought_class_fung_custom <- rbind(OTUs_drought_class_fung_custom,
                                        bind_cols(data.frame(id = "Fun_OTU12108", 
                                                             model = "~ trt*region + (1|site/field), dispformula = ~ trt*region, data = OTU_tab_dom_fung_d0, family = nbinom1"),
                                                  summary(mod_Fun_OTU12108)$coefficients$cond[c(2,8,14,20)] %>%
                                                    data.frame() %>%
                                                    rename(value = names(.)[1]) %>%
                                                    mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                    pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


# One violating assumptions is "Fun_OTU11684" 
mod_Fun_OTU11684 <- glmmTMB(Fun_OTU11684 ~ trt + (1|site/field),
                            data =   
                              OTU_tab_dom_fung_d0 %>%
                              mutate(across(where(is.numeric), ~1 * (. > 0))),
                            family = binomial)
res <- simulateResiduals(mod_Fun_OTU11684)
plot(res)
testZeroInflation(res)
testDispersion(res)

# Add it to the data frame
OTUs_drought_class_fung_custom <- rbind(OTUs_drought_class_fung_custom,
                                        bind_cols(data.frame(id = "Fun_OTU11684", 
                                                             model = "~ trt + (1|site/field), data = OTU_tab_dom_fung_d0 %>% mutate(across(where(is.numeric), ~1 * (. > 0))), family = binomial"),
                                                  summary(mod_Fun_OTU11684)$coefficients$cond[c(2,4,6,8)] %>%
                                                    data.frame() %>%
                                                    rename(value = names(.)[1]) %>%
                                                    mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                    pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


# "Fun_OTU4042"
mod_Fun_OTU4042 <-  glmmTMB(Fun_OTU4042 ~ trt*region + (1|site/field),
                            dispformula = ~ region, # tried with site and region, but only region works without warning
                            data =   
                              OTU_tab_dom_fung_d0,
                            family = nbinom1)
res  <- simulateResiduals(mod_Fun_OTU4042)
plot(res)
testZeroInflation(res)

# Add it to the data frame
OTUs_drought_class_fung_custom <- rbind(OTUs_drought_class_fung_custom,
                                        bind_cols(
                                          data.frame(id = "Fun_OTU4042", 
                                                     model = "~ trt*region + (1|site/field), dispformula = ~ region, data = OTU_tab_dom_fung_d0, family = nbinom1"),
                                          summary(mod_Fun_OTU4042)$coefficients$cond[c(2,8,14,20)] %>%
                                            data.frame() %>%
                                            rename(value = names(.)[1]) %>%
                                            mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                            pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


OTUs_drought_class_fung_custom <- OTUs_drought_class_fung_custom %>%
  mutate(Resistance_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resistant",
                                    if_else(trtD_Estimate > 0,
                                            "Opportunistic",
                                            "Sensitive")))


## Put results together (Resistance) ----
OTUs_drought_class_fung_all <- bind_rows(OTUs_drought_class_fung_poisson_no_disp,
                                         OTUs_drought_class_fung_nbinom_no_disp,
                                         OTUs_drought_class_fung_binomial,
                                         OTUs_drought_class_fung_custom)


sum_drought_class_fung <- OTUs_drought_class_fung_all %>%
  group_by(Resistance_class) %>%
  summarize(count = n())

# FUNGI, RESILIENCE ----
## model functions ---- 
apply_glmmTMB_poiss_no_disp_fung_60 <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_fung_d60,
          family = poisson)
}

apply_glmmTMB_nbinom_no_disp_fung_60 <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_fung_d60,
          family = nbinom1)
}

apply_glmmTMB_binomial_fung_60 <- function(i) {
  glmmTMB(as.formula(paste0(i, "~ trt + (1|region/site/field)")), 
          data = OTU_tab_dom_fung_d60 %>% 
            mutate(across(where(is.numeric), ~1 * (. > 0))),
          family = binomial)
}

## POISSON MODELS WITH NO DISPERSION ----
## First pass to see which OTUs throw warnings
models_fung_60_poisson_no_disp <-
  OTUs_drought_class_fung_all %>%
  filter(Resistance_class != "Resistant") %>%
  pull(id) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_poiss_no_disp_fung_60(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_fung_60_poisson_no_disp_warn <- models_fung_60_poisson_no_disp %>%
  extract(. == 'warning') %>% 
  names()

OTUs_fung_60_poisson_no_disp_no_warn <-     OTUs_drought_class_fung_all %>%
  filter(Resistance_class != "Resistant") %>%
  pull(id)  %>%
  extract(!(. %in% OTUs_fung_60_poisson_no_disp_warn))

# Simulate residuals for models where there weren't warnings
models_resids_fung_60_poisson_no_disp <-  models_fung_60_poisson_no_disp[all_of(OTUs_fung_60_poisson_no_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_fung_60_poisson_no_disp <-
  models_resids_fung_60_poisson_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_fung_d60$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_fung_d60$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_fung_60_poisson_no_disp_violate <- model_ps_fung_60_poisson_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05)) 

OTUs_fung_60_poisson_no_disp_good <- model_ps_fung_60_poisson_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  select(id) %>%
  pull()

# Produce the summary stats for models that pass
models_summaries_fung_60_poisson_no_disp <-  
  models_fung_60_poisson_no_disp %>%
  extract(all_of(OTUs_fung_60_poisson_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_fung_60_poisson_no_disp <- models_summaries_fung_60_poisson_no_disp %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient",
                                    "Not Resilient"),
         model = "Poisson with no dispersion")
## NEGATIVE BINOMIAL WITH NO DISP OR ZI ----
# Run models and check models for warning
models_fung_60_nbinom_no_disp <-
  c(OTUs_fung_60_poisson_no_disp_violate %>%
      pull(id),
    OTUs_fung_60_poisson_no_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_nbinom_no_disp_fung_60(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_fung_60_nbinom_no_disp_warn <- models_fung_60_nbinom_no_disp %>%
  extract(. == 'warning') %>% 
  names()

OTUs_fung_60_nbinom_no_disp_no_warn <- c(OTUs_fung_60_poisson_no_disp_violate %>%
                                           pull(id),
                                         OTUs_fung_60_poisson_no_disp_warn) %>%
  extract(!(. %in% OTUs_fung_60_nbinom_no_disp_warn))

# Simulate residuals for models where there weren't warnings
models_resids_fung_60_nbinom_no_disp <-  models_fung_60_nbinom_no_disp[all_of(OTUs_fung_60_nbinom_no_disp_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_fung_60_nbinom_no_disp <-
  models_resids_fung_60_nbinom_no_disp %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_fung_d60$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_fung_d60$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_fung_60_nbinom_no_disp_violate <- model_ps_fung_60_nbinom_no_disp %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))    

OTUs_fung_60_nbinom_no_disp_good <- model_ps_fung_60_nbinom_no_disp %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  select(id) %>%
  pull()


# Produce the summary stats for models that pass
models_summaries_fung_60_nbinom_no_disp <-  
  models_fung_60_nbinom_no_disp %>%
  extract(all_of(OTUs_fung_60_nbinom_no_disp_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD row.
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")

OTUs_drought_class_fung_60_nbinom_no_disp <- models_summaries_fung_60_nbinom_no_disp %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient",
                                    "Not Resilient"),
         model = "Neg binomial with no dispersion")


## BINOMIAL MODELS ----
models_fung_60_binomial <-
  c(OTUs_fung_60_nbinom_no_disp_violate %>% # 30 total
      pull(id),
    OTUs_fung_60_nbinom_no_disp_warn) %>%
  purrr::set_names() %>% 
  map(
    function(i) {
      tryCatch(apply_glmmTMB_binomial_fung_60(i),
               warning = function(w) {
                 if(is.character(w$message))
                   return('warning') })
    })

# Pull OTUs with warnings and no warnings
OTUs_fung_60_binomial_warn <- models_fung_60_binomial %>%
  extract(. == 'warning') %>% 
  names()

OTUs_fung_60_binomial_no_warn <-c(OTUs_fung_60_nbinom_no_disp_violate %>%
                                    pull(id),
                                  OTUs_fung_60_nbinom_no_disp_warn) %>%
  extract(!(. %in% OTUs_fung_60_binomial_warn))

# Simulate residuals for models where there weren't warnings
models_resids_fung_60_binomial <-  models_fung_60_binomial[all_of(OTUs_fung_60_binomial_no_warn)] %>%
  map(
    function(m) {
      simulateResiduals(m)
    })

# Run diagnostic tests
model_ps_fung_60_binomial <-
  models_resids_fung_60_binomial %>%
  map_dfr(
    function(.x) {
      data.frame(
        # Dispersion test
        disp_stat = testDispersion(.x, plot = FALSE)$statistic,
        disp_p = testDispersion(.x, plot = FALSE)$p.value,
        # Zero inflation test
        zi_p = testZeroInflation(.x, plot = FALSE)$p.value,
        # Uniformity test (K-S test)
        ks_p = testUniformity(.x, plot = FALSE)$p.value, 
        # Levene test result for homogeneity of variance
        lev_p = testCategorical(.x, catPred = OTU_tab_dom_fung_d60$trt, plot = FALSE)$homogeneity$`Pr(>F)`[1], # The second value is NA 
        bind_cols(
          # within-group K-S test for uniformity
          testCategorical(.x, catPred = OTU_tab_dom_fung_d60$trt, plot = FALSE)$uniformity$p.value.cor %>%
            data.frame() %>%
            rename(ks_group_p = names(.)[1]) %>%
            mutate(trt = c("C", "D")) %>%
            pivot_wider(values_from = ks_group_p, names_from = trt, names_prefix = "ks_group_p_"))
      )},
    .id = "id")

# filter OTUs where p values are <0.05
OTUs_fung_60_binomial_violate <- model_ps_fung_60_binomial %>%
  rowwise %>% 
  filter(any(c_across(where(is.numeric)) < 0.05))   

OTUs_fung_60_binomial_good <- model_ps_fung_60_binomial %>%
  rowwise %>% 
  filter(all(c_across(where(is.numeric)) >= 0.05)) %>%
  pull(id)


# Produce the summary stats for models that pass
models_summaries_fung_60_binomial <-  
  models_fung_60_binomial %>%
  extract(all_of(OTUs_fung_60_binomial_good)) %>%
  map_dfr( function(m) {
    ## Get Estimate, Std.Error, z value, Pr(>|z|) for trtD
    summary(m)$coefficients$cond[c(2,4,6,8)] %>%
      data.frame() %>%
      rename(value = names(.)[1]) %>%
      mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
      pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")
  },
  .id = "id")


OTUs_drought_class_fung_60_binomial <- models_summaries_fung_60_binomial %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient",
                                    "Not Resilient"),
         model = "Binomial no disp no zi")

## Custom models for last 5 OTUs ----
OTUs_fung_60_binomial_violate %>% pull(id) #"Fun_OTU1792"  "Fun_OTU10217" "Fun_OTU11669"
OTUs_fung_60_binomial_warn #"Fun_OTU12108" "Fun_OTU11913"

# Make data frame to keep track of models and stats
OTUs_drought_class_fung_60_custom <- data.frame(matrix(ncol = 6, nrow=0, 
                                                       dimnames = list(NULL, c("id", "model", "Estimate", "Std.Error", "z value", "Pr(>|z|)"
                                                       ))))

# "Fun_OTU1792" violates assumptions
mod_Fun_OTU1792 <- glmmTMB(Fun_OTU1792 ~ trt*region + (1|site/field),
                           data =   
                             OTU_tab_dom_fung_d60,
                           family = nbinom1)
res <- simulateResiduals(mod_Fun_OTU1792)
plot(res)

# Add it to the data frame
OTUs_drought_class_fung_60_custom <- rbind(OTUs_drought_class_fung_60_custom,
                                           bind_cols(data.frame(id = "Fun_OTU1792", 
                                                                model = "~ trt*region + (1|site/field), data = OTU_tab_dom_fung_d60, family = nbinom1"),
                                                     summary(mod_Fun_OTU1792)$coefficients$cond[c(2,8,14,20)] %>%
                                                       data.frame() %>%
                                                       rename(value = names(.)[1]) %>%
                                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

# "Fun_OTU10217"  
mod_Fun_OTU10217 <- glmmTMB(Fun_OTU10217 ~ trt*region + (1|site/field),
                            data =   
                              OTU_tab_dom_fung_d60,
                            family = nbinom1)
res <- simulateResiduals(mod_Fun_OTU10217)
plot(res)

# Add it to the data frame
OTUs_drought_class_fung_60_custom <- rbind(OTUs_drought_class_fung_60_custom,
                                           bind_cols(data.frame(id = "Fun_OTU10217", 
                                                                model = "~ trt*region + (1|site/field), data = OTU_tab_dom_fung_d60, family = nbinom1"),
                                                     summary(mod_Fun_OTU10217)$coefficients$cond[c(2,8,14,20)] %>%
                                                       data.frame() %>%
                                                       rename(value = names(.)[1]) %>%
                                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Fun_OTU11669"
mod_Fun_OTU11669 <- glmmTMB(Fun_OTU11669 ~ trt*region + (1|site/field),
                            dispformula = ~ site,
                            data =
                              OTU_tab_dom_fung_d60,
                            family = nbinom1)
res <- simulateResiduals(mod_Fun_OTU11669)
plot(res)

# Add it to the data frame
OTUs_drought_class_fung_60_custom <- rbind(OTUs_drought_class_fung_60_custom,
                                           bind_cols(data.frame(id = "Fun_OTU11669", 
                                                                model = "~ trt*region + (1|site/field), dispformula = ~ site, data = OTU_tab_dom_fung_d60, family = nbinom1"),
                                                     summary(mod_Fun_OTU11669)$coefficients$cond[c(2,8,14,20)] %>%
                                                       data.frame() %>%
                                                       rename(value = names(.)[1]) %>%
                                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

# "Fun_OTU12108" 
mod_Fun_OTU12108 <- glmmTMB(Fun_OTU12108 ~ trt + (1|site/field),
                            ziformula = ~ 1, # tried with site and region, but only region works to avoid violating levene test.
                            data =
                              OTU_tab_dom_fung_d60 %>%
                              mutate(across(where(is.numeric), ~1 * (. > 0))),
                            family = binomial)
res <- simulateResiduals(mod_Fun_OTU12108)
plot(res)

# Add it to the data frame
OTUs_drought_class_fung_60_custom <- rbind(OTUs_drought_class_fung_60_custom,
                                           bind_cols(data.frame(id = "Fun_OTU12108", 
                                                                model = "~ trt + (1|site/field), data = OTU_tab_dom_fung_d60 %>% mutate(across(where(is.numeric), ~1 * (. > 0))), family = binomial"),
                                                     summary(mod_Fun_OTU12108)$coefficients$cond[c(2,4,6,8)] %>%
                                                       data.frame() %>%
                                                       rename(value = names(.)[1]) %>%
                                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))

#"Fun_OTU11913"
mod_Fun_OTU11913 <- glmmTMB(Fun_OTU11913 ~ trt*region + (1|site/field),
                            data =
                              OTU_tab_dom_fung_d60,
                            family = nbinom1)
res <- simulateResiduals(mod_Fun_OTU11913)
plot(res)

# Add it to the data frame
OTUs_drought_class_fung_60_custom <- rbind(OTUs_drought_class_fung_60_custom,
                                           bind_cols(data.frame(id = "Fun_OTU11913", 
                                                                model = "~ trt*region + (1|site/field), data = OTU_tab_dom_fung_d60, family = nbinom1"),
                                                     summary(mod_Fun_OTU11913)$coefficients$cond[c(2,8,14,20)] %>%
                                                       data.frame() %>%
                                                       rename(value = names(.)[1]) %>%
                                                       mutate(stat = c("Estimate", "Std.Error", "z value", "Pr(>|z|)")) %>%
                                                       pivot_wider(values_from = value, names_from = stat, names_prefix = "trtD_")))


OTUs_drought_class_fung_60_custom <- OTUs_drought_class_fung_60_custom %>%
  mutate(Resilience_class = if_else(`trtD_Pr(>|z|)` > 0.05,
                                    "Resilient","Not Resilient"),
         model = "custom")


## Put results together (Resilience) ----
OTUs_drought_class_fung_60_all <- bind_rows(OTUs_drought_class_fung_60_poisson_no_disp,
                                            OTUs_drought_class_fung_60_nbinom_no_disp,
                                            OTUs_drought_class_fung_60_binomial,
                                            OTUs_drought_class_fung_60_custom) %>%
  rename_with(~paste0(., "_d60"), c("trtD_Estimate":"trtD_Pr(>|z|)", "model"))



sum_drought_class_fung_60 <- OTUs_drought_class_fung_60_all %>%
  group_by(Resilience_class) %>%
  summarize(count = n())


## Combine everything and save data files ----

response_results_all <- OTUs_drought_class_all %>%
  mutate(Group = "16S") %>%
  rename("Model_resistance" = "model") %>%
  left_join(OTUs_drought_class_60_all %>%
              rename("Model_resilience" = "model_d60"),
            by = "id") %>%
  bind_rows(OTUs_drought_class_fung_all %>%
              mutate(Group = "ITS") %>%
              rename("Model_resistance" = "model") %>%
              left_join(OTUs_drought_class_fung_60_all  %>%
                          rename("Model_resilience" = "model_d60"))) %>%
  select(Group, OTU_ID = id, Resistance_class, Resilience_class, 
         "Main_drought_effect_day_0" = "trtD_Estimate",
         "Drought_p_val_day_0" = "trtD_Pr(>|z|)",
         "Main_drought_effect_day_60" = "trtD_Estimate_d60", 
         "Drought_p_val_day_60" = "trtD_Pr(>|z|)_d60",
         "Model_resistance", "Model_resilience") 


