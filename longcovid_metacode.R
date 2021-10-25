#------ Generalised Linear Mixed Models for long COVID meta-analyses ------ #

# code produced by C. Watson and K. Jansen
# last updated on 25/10/21 by K. Jansen (introductory comments added)

# this code contains all meta-analytic calculations for the paper:
# "Persistent neuropsychiatric symptoms after COVID-19: a systematic review and meta-analysis"

# The code is divided into four major sections: 

# PREPARATION: Here, the packages necessary for the analysis are loaded
# and the data frame is prepared for the analysis.

# MAIN ANALYSIS: In this section, the main analyses of the different symptoms are conducted
# using generalized linear mixed models. Forest plots for all symptoms are produced and saved
# in the working directory.

# SENSITIVITY ANALYSIS: In this section, the main analyses are repeated using an inverse-variance model
# with the Freeman-Tukey double-arcsine transformation. This is done as a sensitivity analysis.

# SUBGROUP ANALYSIS: This section contains the subgroup analyses of 1.) hospitalisation status,
# 2.) symptom severity, 3.) time since discharge and 4.) time since symptom onset.
# z-values and p-values are computed for each subgroup comparison. 
# finally, a forest plot of each subgroup analysis is produced and saved in the working directory.


## PREPARATION -------------------------------------------------------------------------------------------------------------------

# Load packages ------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(metafor)
library(readr)
library(ggplot2)
library(tidyverse)


# set working directory ----------------------------------------------------------------------------------------------------------
# set working directory manually

# setwd("")

# Load and prepare data ----------------------------------------------------------------------------------------------------------
csv_name <- "longcovid_analysis.csv"
all_meta_data <- read.csv(csv_name, sep = ";", na.strings = "")

# only keep rows in which there are data points
all_meta_data <- all_meta_data %>% filter(!is.na(Reference))

# recode NAs
all_meta_data[all_meta_data == "<NA>"] <- NA

# change time variable labels
all_meta_data$long_covid_time_edited <- ifelse(!is.na(all_meta_data$long_covid_time_edited),
                                               ifelse(all_meta_data$long_covid_time_edited == "<12weeks", "<12 weeks", "12+ weeks"), NA)

all_meta_data$subgroup_meta = factor(all_meta_data$subgroup_meta) # make factor from subgroup_meta variable


## ANALYSES ---------------------------------------------------------------------------------------------------------------------

# Prepare data frames for full meta-analysis and subgroup meta-analyses ---------------------------------------------------------

# ----------- full meta-analysis of symptom prevalence 
meta_data <- all_meta_data %>% filter(full_meta == "Y") # this was changed from excluding those in the severity subgroups 

# ----------- subgroup analyses:
# subgroup 1: hospitalised vs non-hospitalised. This will include all studies with column D marked '1'. Column E 'hospitalised' and 'non hospitalised' defines the two groups within this analysis
meta_sample_subs <- all_meta_data %>% filter(subgroup_meta == "1")

# subgroup 2: disease severity (WHO criteria) - defined in 3 categories (A) ITU/WHO critical (B) WHO severe (C) non-ITU/non-WHO-critical
meta_sample_sev <- all_meta_data %>% filter(subgroup_meta == "2") # This will include all studies with column D marked '2'. Column E 'ITU/critical/severe' and 'non-ITU/non-severe' defines the two groups in this analysis.

# subgroup 3: sex 
meta_sample_sex <- all_meta_data %>% filter(subgroup_meta == "3") # subgroup 3: male vs female. This will include all studies with column D marked '3'. Column E 'male' and 'female' defines the two groups within this analysis

# ----------- time analyses:

# time 1: post-discharge <12 vs >12 weeks. This will include all studies with column F as 'post discharge'.  Column G '<12' and '>12' defines the two groups within this analysis.
meta_sample_time1 <- all_meta_data %>% filter(long_covid_duration_analysis == "post discharge") 

# time 2: post-symptom <12wk vs >12wk. This will include all studies with column F as 'post symptom'.  Column G '<12' and '>12' defines the two groups within this analysis.
meta_sample_time2 <- all_meta_data %>% filter(long_covid_duration_analysis == "post symptom or PCR") 

## MAIN ANALYSIS -----------------------------------------------------------------------------------------------

# 1: Cognitive dysfunction objective - Full meta-analysis - 6 studies ------------------------------------------

# prepare data frame for meta-analysis
cog_df <- meta_data %>% 
  filter(!is.na(n_reported_cognitive_dsyfunction_1)) %>% 
  mutate(cog_prop = n_reported_cognitive_dsyfunction_1/n) %>% 
  arrange(cog_prop)


# estimate glmm for full meta-analysis
cog_glmm <- rma.glmm(xi = n_reported_cognitive_dsyfunction_1, ni = n, data = cog_df,
                      slab = paste(Reference),
                      measure = "PLO")
cog_glmm

# get prevalence estimates
predict(cog_glmm, transf = transf.ilogit) 

# get sample size
sum(cog_df$n)

# draw forest plot for full objective cognitive dysfunction meta-analysis (and save it)
pdf(file="obj_cog_full.pdf", width = 11)

forest(cog_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit) 

title(main = "Proportion of patients identified with objective cognitive dysfunction")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                             .(formatC(cog_glmm$I2, digits=1, format="f")), "%)")))

dev.off()

# 2: Cognitive dysfunction subjective - Full meta-analysis -----------------------------------------------------------------------

# prepare data frame for meta-analysis
cog2_df <- meta_data %>% 
  filter(!is.na(n_reported_cognitive_dsyfunction_2)) %>% 
  mutate(cog2_prop = n_reported_cognitive_dsyfunction_2/n) %>% 
  arrange(cog2_prop)


# estimate glmm for full meta-analysis
cog2_glmm <- rma.glmm(xi = n_reported_cognitive_dsyfunction_2, ni = n, data = cog2_df,
                      slab = paste(Reference),
                      measure = "PLO")
cog2_glmm

# get prevalence estimates
predict(cog2_glmm, transf = transf.ilogit)

# get sample size
sum(cog2_df$n)

# draw forest plot for full meta-analysis of subjective cognitive dysfunction (and save it)
pdf(file="subj_cog_full.pdf", width = 11)

forest(cog2_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit)

title(main = "Proportion of patients subjectively reporting cognitive dysfunction")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(cog2_glmm$I2, digits=1, format="f")), "%)")))

dev.off()


# 3. Sensorimotor issues - Full meta-analysis ------------------------------------------------------------------------------------

# prepare data for meta-analysis
sm_df <- meta_data %>% 
  filter(!is.na(n_sensmotor1)) %>% 
  mutate(sm_prop = n_sensmotor1/n) %>% 
  arrange(sm_prop)

# estimate glmm
sm_glmm <- rma.glmm(xi = n_sensmotor1, ni = n, data = sm_df,
                      slab = paste(Reference),
                      measure = "PLO")
sm_glmm

# get prevalence
predict(sm_glmm, transf = transf.ilogit)

# get sample size
sum(sm_df$n)

# draw and save forest plot
pdf(file="sensorimotor_full.pdf", width = 11)

forest(sm_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit)

title(main = "Proportion of patients reporting sensorimotor dysfunction")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(sm_glmm$I2, digits=1, format="f")), "%)")))

dev.off()

# 4. Dizziness/ vertigo - Full meta-analysis ------------------------------------------------------------------------------------ 

# prepare data for meta-analysis
diz_df <- meta_data %>% 
  filter(!is.na(n_dizzivertigo))  %>% 
  mutate(diz_prop = n_dizzivertigo/n) %>% 
  arrange(diz_prop)

# estimate glmm
dizglmm <- rma.glmm(xi = n_dizzivertigo, ni = n, data = diz_df,
                    slab = paste(Reference),
                    measure = "PLO")
dizglmm

# get prevalence estimate
predict(dizglmm, transf = transf.ilogit) 

# get sample size
sum(diz_df$n)


# draw and save forest plot
pdf(file="dizzy_full.pdf", width = 11)

forest(dizglmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit)

title(main = "Proportion of patients reporting dizziness or vertigo")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(dizglmm$I2, digits=1, format="f")), "%)")))
dev.off()

# 5. Sleep problems - full meta-analysis ----------------------------------------------------------------------------------------

# prepare data for meta-analysis
sleep_df <- meta_data %>% 
  filter(!is.na(n_sleepprob)) %>% 
  mutate(sleep_prop = n_sleepprob/n) %>% 
  arrange(sleep_prop)

# estimate glmm
sleep_glmm <- rma.glmm(xi = n_sleepprob, ni = n, data = sleep_df,
                    slab = paste(Reference),
                    measure = "PLO")
sleep_glmm

# get prevalence estimate
predict(sleep_glmm, transf = transf.ilogit)

# get sample size
sum(sleep_df$n)


# draw and save forest plot
pdf(file="sleep_full.pdf", width = 11)

forest(sleep_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit)

title(main = "Proportion of patients reporting sleep problems")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(sleep_glmm$I2, digits=1, format="f")), "%)")))
dev.off()


# 6. Depression - full meta-analysis --------------------------------------------------------------------------------------------
#meta_data$n_depression <- as.numeric(as.character(meta_data$n_depression))

# prepare data for meta-analysis
dep_df <- meta_data %>% 
  filter(!is.na(n_depression)) %>% 
  mutate(dep_prop = n_depression/n) %>% 
  arrange(dep_prop)

# estimate glmm
dep_glmm <- rma.glmm(xi = n_depression, ni = n, data = dep_df,
                    slab = paste(Reference),
                    measure = "PLO")
dep_glmm

# get prevalence estimate
predict(dep_glmm, transf = transf.ilogit) 

# get sample size
sum(dep_df$n)

# draw and save forest plot
pdf(file="depression_full.pdf", width = 11)

forest(dep_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit) 

title(main = "Proportion of patients with depression")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(dep_glmm$I2, digits=1, format="f")), "%)")))
dev.off()

# 7. Anxiety - full meta-analysis -----------------------------------------------------------------------------------------------

# prepare data for meta-analysis
anx_df <- meta_data %>% 
  filter(!is.na(n_anxiety)) %>% 
  mutate(anx_prop = n_anxiety/n) %>% 
  arrange(anx_prop)

# estimate glmm
anx_glmm <- rma.glmm(xi = n_anxiety, ni = n, data = anx_df,
                    slab = paste(Reference),
                    measure = "PLO")
anx_glmm

# get prevalence estimate
predict(anx_glmm, transf = transf.ilogit) 

# get sample size
sum(anx_df$n)

# draw and save forest plot
pdf(file="anxiety_full.pdf", width = 12)

forest(anx_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit)

title(main = "Proportion of patients with anxiety")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(anx_glmm$I2, digits=1, format="f")), "%)")))
dev.off()

# 8. PTSD: full meta-analysis ---------------------------------------------------------------------------------------------------
#meta_data$n_PTSD <- as.numeric(as.character(meta_data$n_PTSD))

# prepare data for meta-analysis
pts_df <- meta_data %>% 
  filter(!is.na(n_PTSD)) %>% 
  mutate(pts_prop = n_PTSD/n) %>% 
  arrange(pts_prop)

# estimate glmm
pts_glmm <- rma.glmm(xi = n_PTSD, ni = n, data = pts_df,
                     slab = paste(Reference),
                     measure = "PLO")
pts_glmm

# get prevalence estimate
predict(pts_glmm, transf = transf.ilogit) 

# get sample size
sum(pts_df$n)

# draw and save forest plot
pdf(file="pts_full.pdf", width = 11)

forest(pts_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit) 

title(main = "Proportion of patients reporting PTSD/PTS in long covid")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(pts_glmm$I2, digits=1, format="f")), "%)")))
dev.off()

# 9. Altered taste - full meta-analysis ------------------------------------------------------------------------------------------

# prepare data for meta-analysis
taste_df <- meta_data %>% 
  filter(!is.na(n_dysguesia)) %>% 
  mutate(taste_prop = n_dysguesia/n) %>% 
  arrange(taste_prop)

# estimate glmm
taste_glmm <- rma.glmm(xi = n_dysguesia, ni = n, data = taste_df,
                    slab = paste(Reference),
                    measure = "PLO")
taste_glmm

# get prevalence estimate
predict(taste_glmm, transf = transf.ilogit) 

# get sample size
sum(taste_df$n)

# draw and save forest plot
pdf(file="taste_full.pdf", width = 11)

forest(taste_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit) 

title(main = "Proportion of patients reporting dysgeusia")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(taste_glmm$I2, digits=1, format="f")), "%)")))
dev.off()

# 10. Altered smell: full meta-analysis ----------------------------------------------------------------------------------------

# prepare data for meta-analysis
smell_df <- meta_data %>% 
  filter(!is.na(n_dysnosmia)) %>% 
  mutate(smell_prop = n_dysnosmia/n) %>% 
  arrange(smell_prop)

# estimate glmm
smell_glmm <- rma.glmm(xi = n_dysnosmia, ni = n, data = smell_df,
                    slab = paste(Reference),
                    measure = "PLO")
smell_glmm

# get prevalence estimate
predict(smell_glmm, transf = transf.ilogit) 

# get sample size
sum(smell_df$n)

# draw and save forest plot
pdf(file="smell_full.pdf", width = 12)

forest(smell_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit) 

title(main = "Proportion of patients reporting dysosmia")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(smell_glmm$I2, digits=1, format="f")), "%)")))
dev.off()

# 11. Fatigue - full meta-analysis ----------------------------------------------------------------------------------------------

# prepare data for meta-analysis
fat_df <- meta_data %>% 
  filter(!is.na(n_fatigue)) %>% 
  mutate(fat_prop = n_fatigue/n) %>% 
  arrange(fat_prop)

# estimate glmm
fat_glmm <- rma.glmm(xi = n_fatigue, ni = n, data = fat_df,
                    slab = paste(Reference),
                    measure = "PLO")
fat_glmm

# get prevalence estimate
predict(fat_glmm, transf = transf.ilogit) 

# get sample size
sum(fat_df$n)


# draw and save forest plot
pdf(file="fatigue_full.pdf", width = 11)

forest(fat_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit)

title(main = "Proportion of patients reporting fatigue in long covid")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(fat_glmm$I2, digits=1, format="f")), "%)")))
dev.off()

# 12 Headache - full meta-analysis ----------------------------------------------------------------------------------------------

# prepare data for meta-analysis
head_df <- meta_data %>% 
  filter(!is.na(n_headache)) %>% 
  mutate(head_prop = n_headache/n) %>% 
  arrange(head_prop)

# estimate glmm
head_glmm <- rma.glmm(xi = n_headache, ni = n, data = head_df,
                     slab = paste(Reference),
                     measure = "PLO")
head_glmm

# get prevalence estimate
predict(head_glmm, transf = transf.ilogit) 

# get sample size
sum(head_df$n)

# draw and save forest plot
pdf(file="head_full.pdf", width = 12)

forest(head_glmm, 
       xlim=c(-0.25,1),
       transf = transf.ilogit) 

title(main = "Proportion of patients reporting headache")

text(-0.15, -1, pos=4, cex=1, bquote(paste("(", I^2, " = ",
                                           .(formatC(head_glmm$I2, digits=1, format="f")), "%)")))
dev.off()

## SENSITIVITY ANALYSIS ----------------------------------------------------------------------------------------

# objective cognitive dysfunction
sensitivity_cog1 <- rma.uni(xi = n_reported_cognitive_dsyfunction_1, ni = n, data = cog_df,
                            slab = paste(Reference), 
                            measure = "PFT",
                            method = "REML")
sensitivity_cog1

predict(sensitivity_cog1, transf=transf.ipft.hm, targs=list(ni=sensitivity_cog1$ni))


# subjective cognitive dysfunction
sensitivity_cog2 <- rma.uni(xi = n_reported_cognitive_dsyfunction_2, ni = n, data = cog2_df,
                            slab = paste(Reference), 
                            measure = "PFT",
                            method = "REML")
sensitivity_cog2

predict(sensitivity_cog2, transf=transf.ipft.hm, targs=list(ni=sensitivity_cog2$ni))

# sensorimotor issues
sensitivity_sm <- rma.uni(xi = n_sensmotor1, ni = n, data = sm_df,
                            slab = paste(Reference), 
                            measure = "PFT",
                            method = "REML")
sensitivity_sm

predict(sensitivity_sm, transf=transf.ipft.hm, targs=list(ni=sensitivity_sm$ni))

# dizziness/vertigo
sensitivity_dizzy <- rma.uni(xi = n_dizzivertigo, ni = n, data = diz_df,
                          slab = paste(Reference), 
                          measure = "PFT",
                          method = "REML")
sensitivity_dizzy

predict(sensitivity_dizzy, transf=transf.ipft.hm, targs=list(ni=sensitivity_dizzy$ni))


# sleep
sensitivity_sleep <- rma.uni(xi = n_sleepprob, ni = n, data = sleep_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_sleep

predict(sensitivity_sleep, transf=transf.ipft.hm, targs=list(ni=sensitivity_sleep$ni))

# depression
sensitivity_depr <- rma.uni(xi = n_depression, ni = n, data = dep_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_depr

predict(sensitivity_depr, transf=transf.ipft.hm, targs=list(ni=sensitivity_depr$ni))

# anxiety
sensitivity_anx <- rma.uni(xi = n_anxiety, ni = n, data = anx_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_anx

predict(sensitivity_anx, transf=transf.ipft.hm, targs=list(ni=sensitivity_anx$ni))

# ptsd
sensitivity_pts <- rma.uni(xi = n_PTSD, ni = n, data = pts_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_pts

predict(sensitivity_pts, transf=transf.ipft.hm, targs=list(ni=sensitivity_pts$ni))

# dysgeusia
sensitivity_taste <- rma.uni(xi = n_dysguesia, ni = n, data = taste_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_taste

predict(sensitivity_taste, transf=transf.ipft.hm, targs=list(ni=sensitivity_taste$ni))

# dysosmia
sensitivity_smell <- rma.uni(xi = n_dysnosmia, ni = n, data = smell_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_smell

predict(sensitivity_smell, transf=transf.ipft.hm, targs=list(ni=sensitivity_smell$ni))

# Fatigue
sensitivity_fat <- rma.uni(xi = n_fatigue, ni = n, data = fat_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_fat

predict(sensitivity_fat, transf=transf.ipft.hm, targs=list(ni=sensitivity_fat$ni))


# Headache
sensitivity_head <- rma.uni(xi = n_headache, ni = n, data = head_df,
                             slab = paste(Reference), 
                             measure = "PFT",
                             method = "REML")
sensitivity_head

predict(sensitivity_head, transf=transf.ipft.hm, targs=list(ni=sensitivity_head$ni))


## SUBGROUP ANALYSIS -------------------------------------------------------------------------------------------


# 13.1 Subgroup analysis: hospitalisation status --------------------------------------------------------------

# subjective cognitive dysfunction
hosp_subjCog <- rma.glmm(xi = n_reported_cognitive_dsyfunction_2, ni = n, data = meta_sample_subs,
                         subset = (subgroup_descrip == "hospitalised"),
                         slab = paste(Reference),
                         measure = "PLO")
nonhosp_subjCog <- rma.glmm(xi = n_reported_cognitive_dsyfunction_2, ni = n, data = meta_sample_subs,
                            subset = (subgroup_descrip == "non-hospitalised"),
                            slab = paste(Reference),
                            measure = "PLO")

# Dizziness
hosp_dizzy <- rma.glmm(xi = n_dizzivertigo, ni = n, data = meta_sample_subs,
                       subset = (subgroup_descrip == "hospitalised"),
                       slab = paste(Reference),
                       measure = "PLO")
nonhosp_dizzy <- rma.glmm(xi = n_dizzivertigo, ni = n, data = meta_sample_subs,
                          subset = (subgroup_descrip == "non-hospitalised"),
                          slab = paste(Reference),
                          measure = "PLO" )

# Sleep
hosp_sleep <- rma.glmm(xi = n_sleepprob, ni = n, data = meta_sample_subs,
                       subset = (subgroup_descrip == "hospitalised"),
                       slab = paste(Reference),
                       measure = "PLO")
nonhosp_sleep <- rma.glmm(xi = n_sleepprob, ni = n, data = meta_sample_subs,
                          subset = (subgroup_descrip == "non-hospitalised"),
                          slab = paste(Reference),
                          measure = "PLO")

# Depression
hosp_depr <- rma.glmm(xi = n_depression, ni = n, data = meta_sample_subs,
                      subset = (subgroup_descrip == "hospitalised"),
                      slab = paste(Reference),
                      measure = "PLO")
nonhosp_depr <- rma.glmm(xi = n_depression, ni = n, data = meta_sample_subs,
                         subset = (subgroup_descrip == "non-hospitalised"),
                         slab = paste(Reference),
                         measure = "PLO")

# Anxiety
hosp_anx <- rma.glmm(xi = n_anxiety, ni = n, data = meta_sample_subs,
                     subset = (subgroup_descrip == "hospitalised"),
                     slab = paste(Reference),
                     measure = "PLO")
nonhosp_anx <- rma.glmm(xi = n_anxiety, ni = n, data = meta_sample_subs,
                        subset = (subgroup_descrip == "non-hospitalised"),
                        slab = paste(Reference),
                        measure = "PLO")

# PTSD
hosp_ptsd <- rma.glmm(xi = n_PTSD, ni = n, data = meta_sample_subs,
                      subset = (subgroup_descrip == "hospitalised"),
                      slab = paste(Reference),
                      measure = "PLO")
nonhosp_ptsd <- rma.glmm(xi = n_PTSD, ni = n, data = meta_sample_subs,
                         subset = (subgroup_descrip == "non-hospitalised"),
                         slab = paste(Reference),
                         measure = "PLO")

# Fatigue
hosp_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_subs,
                     subset = (subgroup_descrip == "hospitalised"),
                     slab = paste(Reference),
                     measure = "PLO")
nonhosp_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_subs,
                        subset = (subgroup_descrip == "non-hospitalised"),
                        slab = paste(Reference),
                        measure = "PLO")

# Headache
hosp_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_subs,
                      subset = (subgroup_descrip == "hospitalised"),
                      slab = paste(Reference),
                      measure = "PLO")
nonhosp_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_subs,
                         subset = (subgroup_descrip == "non-hospitalised"),
                         slab = paste(Reference),
                         measure = "PLO")
# subjCog
zval_subjCog = (coef(hosp_subjCog)-coef(nonhosp_subjCog))/sqrt(hosp_subjCog$se^2+nonhosp_subjCog$se^2)
pval_subjCog = 2*pnorm(abs(zval_subjCog), lower.tail = FALSE)

# dizzy
zval_dizzy = (coef(hosp_dizzy)-coef(nonhosp_dizzy))/sqrt(hosp_dizzy$se^2+nonhosp_dizzy$se^2)
pval_dizzy = 2*pnorm(abs(zval_dizzy), lower.tail = FALSE)

# sleep
zval_sleep = (coef(hosp_sleep)-coef(nonhosp_sleep))/sqrt(hosp_sleep$se^2+nonhosp_sleep$se^2)
pval_sleep = 2*pnorm(abs(zval_sleep), lower.tail = FALSE)

# depr
zval_depr= (coef(hosp_depr)-coef(nonhosp_depr))/sqrt(hosp_depr$se^2+nonhosp_depr$se^2)
pval_depr = 2*pnorm(abs(zval_depr), lower.tail = FALSE)

# anx
zval_anx = (coef(hosp_anx)-coef(nonhosp_anx))/sqrt(hosp_anx$se^2+nonhosp_anx$se^2)
pval_anx = 2*pnorm(abs(zval_anx), lower.tail = FALSE)

# ptsd
zval_ptsd = (coef(hosp_ptsd)-coef(nonhosp_ptsd))/sqrt(hosp_ptsd$se^2+nonhosp_ptsd$se^2)
pval_ptsd = 2*pnorm(abs(zval_ptsd), lower.tail = FALSE)

# fat
zval_fat = (coef(hosp_fat)-coef(nonhosp_fat))/sqrt(hosp_fat$se^2+nonhosp_fat$se^2)
pval_fat = 2*pnorm(abs(zval_fat), lower.tail = FALSE)

# head
zval_head = (coef(hosp_head)-coef(nonhosp_head))/sqrt(hosp_head$se^2+nonhosp_head$se^2)
pval_head = 2*pnorm(abs(zval_head), lower.tail = FALSE)


data_hospitalisation <- bind_rows(
  as.data.frame(predict(hosp_anx, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_anx, transf = transf.ilogit)),
  as.data.frame(predict(hosp_depr, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_depr, transf = transf.ilogit)),
  as.data.frame(predict(hosp_sleep, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_sleep, transf = transf.ilogit)),
  as.data.frame(predict(hosp_subjCog, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_subjCog, transf = transf.ilogit)),
  as.data.frame(predict(hosp_dizzy, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_dizzy, transf = transf.ilogit)),
  as.data.frame(predict(hosp_head, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_head, transf = transf.ilogit)),
  as.data.frame(predict(hosp_fat, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_fat, transf = transf.ilogit)),
  as.data.frame(predict(hosp_ptsd, transf = transf.ilogit)),
  as.data.frame(predict(nonhosp_ptsd, transf = transf.ilogit))
)

data_hospitalisation$symptom = rep(c("anxiety", "depression", "sleep problems", "cog dysf subj.",
                                     "dizziness", "headache", "fatigue", "PTSD/PTSS"), each = 2)
data_hospitalisation$subgroup = rep(c("hospitalised", "non-hospitalised"), 8)

data_hospitalisation$I2 = c(hosp_anx$I2, nonhosp_anx$I2, hosp_depr$I2, nonhosp_depr$I2, hosp_sleep$I2, nonhosp_sleep$I2,
                            hosp_subjCog$I2, nonhosp_subjCog$I2, hosp_dizzy$I2, nonhosp_dizzy$I2, hosp_head$I2, nonhosp_head$I2,
                            hosp_fat$I2, nonhosp_fat$I2, hosp_ptsd$I2, nonhosp_ptsd$I2)

data_hospitalisation$zval = round(c(NA, zval_anx, NA, zval_depr, NA, zval_sleep, 
                                    NA, zval_subjCog, NA, zval_dizzy, NA, zval_head,
                                    NA, zval_fat, NA, zval_ptsd),2)

data_hospitalisation$pval = round(c(NA, pval_anx, NA, pval_depr, NA, pval_sleep, 
                                    NA, pval_subjCog, NA, pval_dizzy, NA, pval_head,
                                    NA, pval_fat, NA, pval_ptsd),3)

data_hospitalisation$annotation = ifelse(is.na(data_hospitalisation$zval), "", 
                                         paste0("z = ", format(data_hospitalisation$zval, nsmall = 2), ", p = ",
                                                format(data_hospitalisation$pval, nsmall = 3)))
data_hospitalisation$annotation = ifelse(data_hospitalisation$pval == 0.000,
                                         paste0("z = ", data_hospitalisation$zval, ", p < 0.001"),
                                         data_hospitalisation$annotation)

data_hospitalisation$id = 1:16

superplot_hosp <- ggplot(data_hospitalisation, aes(x = pred, xmin = ci.lb, xmax = ci.ub, y = subgroup))+
  geom_point(color = "black")+
  geom_errorbarh(height = .1)+
  geom_text(aes(label = annotation, x = 0.69, y = 1.5), size = 3, na.rm = TRUE)+
  scale_x_continuous(limits = c(0, 0.75), name = "Prevalence")+ 
  facet_grid(symptom  ~. )+
  ggtitle("Subgroup analysis: Hospitalisation status") +
  ylab("")+
  theme_bw()+
  theme( plot.title = element_text(size = 8),
         axis.title = element_text(size = 7),
         axis.text = element_text(size = 6),
         strip.text = element_text(size = 6))

superplot_hosp

plot_name = "superplot_hosp.pdf"
ggsave(plot_name, device = "pdf", height = 22, width = 19, units = "cm")


# 13.2 Subgroup analysis: severity ----------------------------------------------------------------------------

# subjective cognitive dysfunction:
crit_subjCog <- rma.glmm(xi = n_reported_cognitive_dsyfunction_2, ni = n, data = meta_sample_sev,
                         subset = (subgroup_descrip == "ITU or critical or severe"),
                         slab = paste(Reference),
                         measure = "PLO")
noncrit_subjCog <- rma.glmm(xi = n_reported_cognitive_dsyfunction_2, ni = n, data = meta_sample_sev,
                            subset = (subgroup_descrip == "non-ITU or non-severe"),
                            slab = paste(Reference),
                            measure = "PLO")

# Sleep
crit_sleep <- rma.glmm(xi = n_sleepprob, ni = n, data = meta_sample_sev,
                       subset = (subgroup_descrip == "ITU or critical or severe"),
                       slab = paste(Reference),
                       measure = "PLO")
noncrit_sleep <- rma.glmm(xi = n_sleepprob, ni = n, data = meta_sample_sev,
                          subset = (subgroup_descrip == "non-ITU or non-severe"),
                          slab = paste(Reference),
                          measure = "PLO")

# Anxiety
crit_anx <- rma.glmm(xi = n_anxiety, ni = n, data = meta_sample_sev,
                     subset = (subgroup_descrip == "ITU or critical or severe"),
                     slab = paste(Reference),
                     measure = "PLO")
noncrit_anx <- rma.glmm(xi = n_anxiety, ni = n, data = meta_sample_sev,
                        subset = (subgroup_descrip == "non-ITU or non-severe"),
                        slab = paste(Reference),
                        measure = "PLO")

# Dysguesia
crit_taste <- rma.glmm(xi = n_dysguesia, ni = n, data = meta_sample_sev,
                       subset = (subgroup_descrip == "ITU or critical or severe"),
                       slab = paste(Reference),
                       measure = "PLO")
noncrit_taste <- rma.glmm(xi = n_dysguesia, ni = n, data = meta_sample_sev,
                          subset = (subgroup_descrip == "non-ITU or non-severe"),
                          slab = paste(Reference),
                          measure = "PLO")

# Dysosmia
crit_smell <- rma.glmm(xi = n_dysnosmia, ni = n, data = meta_sample_sev,
                       subset = (subgroup_descrip == "ITU or critical or severe"),
                       slab = paste(Reference),
                       measure = "PLO")
noncrit_smell <- rma.glmm(xi = n_dysnosmia, ni = n, data = meta_sample_sev,
                          subset = (subgroup_descrip == "non-ITU or non-severe"),
                          slab = paste(Reference),
                          measure = "PLO")

# Fatigue
crit_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_sev,
                     subset = (subgroup_descrip == "ITU or critical or severe"),
                     slab = paste(Reference),
                     measure = "PLO")
noncrit_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_sev,
                        subset = (subgroup_descrip == "non-ITU or non-severe"),
                        slab = paste(Reference),
                        measure = "PLO")

# Headache
crit_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_sev,
                      subset = (subgroup_descrip == "ITU or critical or severe"),
                      slab = paste(Reference),
                      measure = "PLO")
noncrit_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_sev,
                         subset = (subgroup_descrip == "non-ITU or non-severe"),
                         slab = paste(Reference),
                         measure = "PLO")

# subjCog
zval_subjCog = (coef(crit_subjCog)-coef(noncrit_subjCog))/sqrt(crit_subjCog$se^2+noncrit_subjCog$se^2)
pval_subjCog = 2*pnorm(abs(zval_subjCog), lower.tail = FALSE)

# sleep
zval_sleep = (coef(crit_sleep)-coef(noncrit_sleep))/sqrt(crit_sleep$se^2+noncrit_sleep$se^2)
pval_sleep = 2*pnorm(abs(zval_sleep), lower.tail = FALSE)

# anx
zval_anx = (coef(crit_anx)-coef(noncrit_anx))/sqrt(crit_anx$se^2+noncrit_anx$se^2)
pval_anx = 2*pnorm(abs(zval_anx), lower.tail = FALSE)

# dysgeusia
zval_taste = (coef(crit_taste)-coef(noncrit_taste))/sqrt(crit_taste$se^2+noncrit_taste$se^2)
pval_taste = 2*pnorm(abs(zval_taste), lower.tail = FALSE)

# dysosmia
zval_smell = (coef(crit_smell)-coef(noncrit_smell))/sqrt(crit_smell$se^2+noncrit_smell$se^2)
pval_smell = 2*pnorm(abs(zval_smell), lower.tail = FALSE)

# fat
zval_fat = (coef(crit_fat)-coef(noncrit_fat))/sqrt(crit_fat$se^2+noncrit_fat$se^2)
pval_fat = 2*pnorm(abs(zval_fat), lower.tail = FALSE)

# head
zval_head = (coef(crit_head)-coef(noncrit_head))/sqrt(crit_head$se^2+noncrit_head$se^2)
pval_head = 2*pnorm(abs(zval_head), lower.tail = FALSE)


data_severity <- bind_rows(
  as.data.frame(predict(crit_anx, transf = transf.ilogit)),
  as.data.frame(predict(noncrit_anx, transf = transf.ilogit)),
  as.data.frame(predict(crit_sleep, transf = transf.ilogit)),
  as.data.frame(predict(noncrit_sleep, transf = transf.ilogit)),
  as.data.frame(predict(crit_subjCog, transf = transf.ilogit)),
  as.data.frame(predict(noncrit_subjCog, transf = transf.ilogit)),
  as.data.frame(predict(crit_head, transf = transf.ilogit)),
  as.data.frame(predict(noncrit_head, transf = transf.ilogit)),
  as.data.frame(predict(crit_fat, transf = transf.ilogit)),
  as.data.frame(predict(noncrit_fat, transf = transf.ilogit)),
  as.data.frame(predict(crit_taste, transf = transf.ilogit)),
  as.data.frame(predict(noncrit_taste, transf = transf.ilogit)),
  as.data.frame(predict(crit_smell, transf = transf.ilogit)),
  as.data.frame(predict(noncrit_smell, transf = transf.ilogit)),
)

data_severity$symptom = rep(c("anxiety", "sleep problems", "cog dysf subj.",
                              "headache", "fatigue", "dysgeusia", "dysosmia"), each = 2)
data_severity$subgroup = rep(c("ITU or critical or severe", "non-ITU or non-severe"), 7)

data_severity$I2 = c(crit_anx$I2, noncrit_anx$I2, crit_sleep$I2, noncrit_sleep$I2, crit_subjCog$I2, noncrit_subjCog$I2,
                     crit_head$I2, noncrit_head$I2, crit_fat$I2, noncrit_fat$I2, crit_taste$I2, noncrit_taste$I2,
                     crit_smell$I2, noncrit_smell$I2)

data_severity$zval = round(c(NA, zval_anx, NA, zval_sleep, NA, zval_subjCog, 
                             NA, zval_head, NA, zval_fat, NA, zval_taste,
                             NA, zval_smell),2)

data_severity$pval = round(c(NA, pval_anx, NA, pval_sleep, NA, pval_subjCog, 
                             NA, pval_head, NA, pval_fat, NA, pval_taste,
                             NA, pval_smell),3)

data_severity$annotation = ifelse(is.na(data_severity$zval), "", 
                                  paste0("z = ", format(data_severity$zval, nsmall = 2), ", p = ",
                                         format(data_severity$pval, nsmall = 3)))
data_severity$annotation = ifelse(data_severity$pval == 0.000,
                                  paste0("z = ", data_severity$zval, ", p < 0.001"),
                                  data_severity$annotation)

data_severity$id = 1:14


superplot_sev <- ggplot(data_severity, aes(x = pred, xmin = ci.lb, xmax = ci.ub, y = subgroup))+
  geom_point(color = "black")+
  geom_errorbarh(height = .1)+
  geom_text(aes(label = annotation, x = 0.69, y = 1.5), size = 3, na.rm = TRUE)+
  scale_x_continuous(limits = c(0, 0.75), name = "Prevalence")+ 
  facet_grid(symptom  ~. )+
  ggtitle("Subgroup analysis: Severity") +
  ylab("")+
  theme_bw()+
  theme( plot.title = element_text(size = 8),
         axis.title = element_text(size = 7),
         axis.text = element_text(size = 6),
         strip.text = element_text(size = 6))

superplot_sev

plot_name = "superplot_sev.pdf"
ggsave(plot_name, device = "pdf", height = 19.25, width = 19, units = "cm")

# 13.3 Subgroup analysis: time 1 ------------------------------------------------------------------------------

# Sleep
time1more_sleep <- rma.glmm(xi = n_sleepprob, ni = n, data = meta_sample_time1,
                            subset = (long_covid_time_edited == "12+ weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
time1less_sleep <- rma.glmm(xi = n_sleepprob, ni = n, data = meta_sample_time1,
                            subset = (long_covid_time_edited == "<12 weeks"),
                            slab = paste(Reference),
                            measure = "PLO")

# Depression
time1more_depr <- rma.glmm(xi = n_depression, ni = n, data = meta_sample_time1,
                           subset = (long_covid_time_edited == "12+ weeks"),
                           slab = paste(Reference),
                           measure = "PLO")
time1less_depr <- rma.glmm(xi = n_depression, ni = n, data = meta_sample_time1,
                           subset = (long_covid_time_edited == "<12 weeks"),
                           slab = paste(Reference),
                           measure = "PLO")
# Anxiety
time1more_anx <- rma.glmm(xi = n_anxiety, ni = n, data = meta_sample_time1,
                          subset = (long_covid_time_edited == "12+ weeks"),
                          slab = paste(Reference),
                          measure = "PLO")
time1less_anx <- rma.glmm(xi = n_anxiety, ni = n, data = meta_sample_time1,
                          subset = (long_covid_time_edited == "<12 weeks"),
                          slab = paste(Reference),
                          measure = "PLO")
# PTSD
time1more_ptsd <- rma.glmm(xi = n_PTSD, ni = n, data = meta_sample_time1,
                           subset = (long_covid_time_edited == "12+ weeks"),
                           slab = paste(Reference),
                           measure = "PLO")
time1less_ptsd <- rma.glmm(xi = n_PTSD, ni = n, data = meta_sample_time1,
                           subset = (long_covid_time_edited == "<12 weeks"),
                           slab = paste(Reference),
                           measure = "PLO")
# Dysguesia
time1more_taste <- rma.glmm(xi = n_dysguesia, ni = n, data = meta_sample_time1,
                            subset = (long_covid_time_edited == "12+ weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
time1less_taste <- rma.glmm(xi = n_dysguesia, ni = n, data = meta_sample_time1,
                            subset = (long_covid_time_edited == "<12 weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
# Dysosmia
time1more_smell <- rma.glmm(xi = n_dysnosmia, ni = n, data = meta_sample_time1,
                            subset = (long_covid_time_edited == "12+ weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
time1less_smell <- rma.glmm(xi = n_dysnosmia, ni = n, data = meta_sample_time1,
                            subset = (long_covid_time_edited == "<12 weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
# Fatigue
time1more_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_time1,
                          subset = (long_covid_time_edited == "12+ weeks"),
                          slab = paste(Reference),
                          measure = "PLO")
time1less_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_time1,
                          subset = (long_covid_time_edited == "<12 weeks"),
                          slab = paste(Reference),
                          measure = "PLO")
# Headache
time1more_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_time1,
                           subset = (long_covid_time_edited == "12+ weeks"),
                           slab = paste(Reference),
                           measure = "PLO")
time1less_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_time1,
                           subset = (long_covid_time_edited == "<12 weeks"),
                           slab = paste(Reference),
                           measure = "PLO")


# sleep
zval_sleep = (coef(time1more_sleep)-coef(time1less_sleep))/sqrt(time1more_sleep$se^2+time1less_sleep$se^2)
pval_sleep = 2*pnorm(abs(zval_sleep), lower.tail = FALSE)

# depr
zval_depr= (coef(time1more_depr)-coef(time1less_depr))/sqrt(time1more_depr$se^2+time1less_depr$se^2)
pval_depr = 2*pnorm(abs(zval_depr), lower.tail = FALSE)

# anx
zval_anx = (coef(time1more_anx)-coef(time1less_anx))/sqrt(time1more_anx$se^2+time1less_anx$se^2)
pval_anx = 2*pnorm(abs(zval_anx), lower.tail = FALSE)

# ptsd
zval_ptsd = (coef(time1more_ptsd)-coef(time1less_ptsd))/sqrt(time1more_ptsd$se^2+time1less_ptsd$se^2)
pval_ptsd = 2*pnorm(abs(zval_ptsd), lower.tail = FALSE)

# taste
zval_taste = (coef(time1more_taste)-coef(time1less_taste))/sqrt(time1more_taste$se^2+time1less_taste$se^2)
pval_taste = 2*pnorm(abs(zval_taste), lower.tail = FALSE)

# smell
zval_smell = (coef(time1more_smell)-coef(time1less_smell))/sqrt(time1more_smell$se^2+time1less_smell$se^2)
pval_smell = 2*pnorm(abs(zval_smell), lower.tail = FALSE)

# fat
zval_fat = (coef(time1more_fat)-coef(time1less_fat))/sqrt(time1more_fat$se^2+time1less_fat$se^2)
pval_fat = 2*pnorm(abs(zval_fat), lower.tail = FALSE)

# head
zval_head = (coef(time1more_head)-coef(time1less_head))/sqrt(time1more_head$se^2+time1less_head$se^2)
pval_head = 2*pnorm(abs(zval_head), lower.tail = FALSE)


data_time1 <- bind_rows(
  as.data.frame(predict(time1less_anx, transf = transf.ilogit)),
  as.data.frame(predict(time1more_anx, transf = transf.ilogit)),
  as.data.frame(predict(time1less_depr, transf = transf.ilogit)),
  as.data.frame(predict(time1more_depr, transf = transf.ilogit)),
  as.data.frame(predict(time1less_sleep, transf = transf.ilogit)),
  as.data.frame(predict(time1more_sleep, transf = transf.ilogit)),
  as.data.frame(predict(time1less_head, transf = transf.ilogit)),
  as.data.frame(predict(time1more_head, transf = transf.ilogit)),
  as.data.frame(predict(time1less_fat, transf = transf.ilogit)),
  as.data.frame(predict(time1more_fat, transf = transf.ilogit)),
  as.data.frame(predict(time1less_taste, transf = transf.ilogit)),
  as.data.frame(predict(time1more_taste, transf = transf.ilogit)),
  as.data.frame(predict(time1less_smell, transf = transf.ilogit)),
  as.data.frame(predict(time1more_smell, transf = transf.ilogit)),
  as.data.frame(predict(time1less_ptsd, transf = transf.ilogit)),
  as.data.frame(predict(time1more_ptsd, transf = transf.ilogit))
)

data_time1$symptom = rep(c("anxiety", "depression", "sleep problems",
                           "headache", "fatigue", "dysgeusia", "dysosmia", "PTSD/PTSS"), each = 2)
data_time1$subgroup = rep(c("<12 weeks", "12+weeks"), 8)

data_time1$zval = round(c(NA, zval_anx, NA, zval_depr, NA, zval_sleep, 
                          NA, zval_head, NA, zval_fat, NA, zval_taste,
                          NA, zval_smell, NA, zval_ptsd),2)

data_time1$pval = round(c(NA, pval_anx, NA, pval_depr, NA, pval_sleep, 
                          NA, pval_head, NA, pval_fat, NA, pval_taste,
                          NA, pval_smell, NA, pval_ptsd),3)

data_time1$annotation = ifelse(is.na(data_time1$zval), "", 
                               paste0("z = ", format(data_time1$zval, nsmall = 2), ", p = ",
                                      format(data_time1$pval, nsmall = 3)))
data_time1$annotation = ifelse(data_time1$pval == 0.000,
                               paste0("z = ", data_time1$zval, ", p < 0.001"),
                               data_time1$annotation)

data_time1$id = 1:16


superplot_time1 <- ggplot(data_time1, aes(x = pred, xmin = ci.lb, xmax = ci.ub, y = subgroup))+
  geom_point(color = "black")+
  geom_errorbarh(height = .1)+
  geom_text(aes(label = annotation, x = 0.69, y = 1.5), size = 3, na.rm = TRUE)+
  scale_x_continuous(limits = c(0, 0.75), name = "Prevalence")+ 
  facet_grid(symptom  ~. )+
  ggtitle("Subgroup analysis: Time since discharge") +
  ylab("")+
  theme_bw()+
  theme( plot.title = element_text(size = 8),
         axis.title = element_text(size = 7),
         axis.text = element_text(size = 6),
         strip.text = element_text(size = 6))

superplot_time1

plot_name = "superplot_time1.pdf"
ggsave(plot_name, device = "pdf", height = 22, width = 19, units = "cm")

# 13.4 Subgroup analysis: time 2 ------------------------------------------------------------------------------

# Dysguesia
time2more_taste <- rma.glmm(xi = n_dysguesia, ni = n, data = meta_sample_time2,
                            subset = (long_covid_time_edited == "12+ weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
time2less_taste <- rma.glmm(xi = n_dysguesia, ni = n, data = meta_sample_time2,
                            subset = (long_covid_time_edited == "<12 weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
# Dysosmia
time2more_smell <- rma.glmm(xi = n_dysnosmia, ni = n, data = meta_sample_time2,
                            subset = (long_covid_time_edited == "12+ weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
time2less_smell <- rma.glmm(xi = n_dysnosmia, ni = n, data = meta_sample_time2,
                            subset = (long_covid_time_edited == "<12 weeks"),
                            slab = paste(Reference),
                            measure = "PLO")
# Fatigue
time2more_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_time2,
                          subset = (long_covid_time_edited == "12+ weeks"),
                          slab = paste(Reference),
                          measure = "PLO")
time2less_fat <- rma.glmm(xi = n_fatigue, ni = n, data = meta_sample_time2,
                          subset = (long_covid_time_edited == "<12 weeks"),
                          slab = paste(Reference),
                          measure = "PLO")
# Headache
time2more_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_time2,
                           subset = (long_covid_time_edited == "12+ weeks"),
                           slab = paste(Reference),
                           measure = "PLO")
time2less_head <- rma.glmm(xi = n_headache, ni = n, data = meta_sample_time2,
                           subset = (long_covid_time_edited == "<12 weeks"),
                           slab = paste(Reference),
                           measure = "PLO")

# taste
zval_taste = (coef(time2more_taste)-coef(time2less_taste))/sqrt(time2more_taste$se^2+time2less_taste$se^2)
pval_taste = 2*pnorm(abs(zval_taste), lower.tail = FALSE)

# smell
zval_smell = (coef(time2more_smell)-coef(time2less_smell))/sqrt(time2more_smell$se^2+time2less_smell$se^2)
pval_smell = 2*pnorm(abs(zval_smell), lower.tail = FALSE)

# fat
zval_fat = (coef(time2more_fat)-coef(time2less_fat))/sqrt(time2more_fat$se^2+time2less_fat$se^2)
pval_fat = 2*pnorm(abs(zval_fat), lower.tail = FALSE)

# head
zval_head = (coef(time2more_head)-coef(time2less_head))/sqrt(time2more_head$se^2+time2less_head$se^2)
pval_head = 2*pnorm(abs(zval_head), lower.tail = FALSE)


data_time2 <- bind_rows(
  as.data.frame(predict(time2less_head, transf = transf.ilogit)),
  as.data.frame(predict(time2more_head, transf = transf.ilogit)),
  as.data.frame(predict(time2less_fat, transf = transf.ilogit)),
  as.data.frame(predict(time2more_fat, transf = transf.ilogit)),
  as.data.frame(predict(time2less_taste, transf = transf.ilogit)),
  as.data.frame(predict(time2more_taste, transf = transf.ilogit)),
  as.data.frame(predict(time2less_smell, transf = transf.ilogit)),
  as.data.frame(predict(time2more_smell, transf = transf.ilogit)),
  
)

data_time2$symptom = rep(c("headache", "fatigue", "dysgeusia", "dysosmia"), each = 2)
data_time2$subgroup = rep(c("<12 weeks", "12+weeks"), 4)


data_time2$zval = round(c( NA, zval_head, NA, zval_fat, NA, zval_taste,
                           NA, zval_smell),2)

data_time2$pval = round(c(NA, pval_head, NA, pval_fat, NA, pval_taste,
                          NA, pval_smell),3)

data_time2$annotation = ifelse(is.na(data_time2$zval), "", 
                               paste0("z = ", format(data_time2$zval, nsmall = 2), ", p = ",
                                      format(data_time2$pval, nsmall = 3)))
data_time2$annotation = ifelse(data_time2$pval == 0.000,
                               paste0("z = ", data_time2$zval, ", p < 0.001"),
                               data_time2$annotation)

data_time2$id = 1:8


superplot_time2 <- ggplot(data_time2, aes(x = pred, xmin = ci.lb, xmax = ci.ub, y = subgroup))+
  geom_point(color = "black")+
  geom_errorbarh(height = .1)+
  geom_text(aes(label = annotation, x = 0.69, y = 1.5), size = 3, na.rm = TRUE)+
  scale_x_continuous(limits = c(0, 0.75), name = "Prevalence")+ 
  facet_grid(symptom  ~. )+
  ggtitle("Subgroup analysis: Time since symptom onset") +
  ylab("")+
  theme_bw()+
  theme( plot.title = element_text(size = 8),
         axis.title = element_text(size = 7),
         axis.text = element_text(size = 6),
         strip.text = element_text(size = 6))

superplot_time2

plot_name = "superplot_time2.pdf"
ggsave(plot_name, device = "pdf", height = 11, width = 19, units = "cm")

