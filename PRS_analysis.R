# ============================================================
# Polygenic Risk Score (PRS) analysis
# Extracted from the original analysis notebook and prepared
# for public sharing. Personal author metadata has been removed.
#
# NOTE: This script assumes analysis datasets have already been
# cleaned and loaded into the R environment. No individual-level
# data or local file paths are included here.
# ============================================================


# ------------------------------------------------------------
# Overview
# ------------------------------------------------------------

# Packages and Data Loading
library(dplyr)
library(ggplot2)
library(tidyverse)
library(psych)
library(Rcpp)
library(MASS)
library(lavaan)
library(nnet)
library(glmnet)
library(boot)
library(broom)
library(survival)
library(caret)

# ------------------------------------------------------------
# Creating a Multivariable Predictor — Regression Models for Univariate Analysis
# ------------------------------------------------------------

PISAallmultinoCov <- lm(formula = fi_new_standardised ~ F1PRS_standardised + F2PRS_standardised + 
                           F3PRS_standardised + F4PRS_standardised + F5PRS_standardised + 
                           F6PRS_standardised + BFPRS_standardised, 
                        data = PISAall)
summary(PISAallmultinoCov)

PISAall$multi <- predict(PISAallmultinoCov, newdata = PISAall)
PISAall$multiscale <- scale(PISAall$multi)


# ------------------------------------------------------------
# Regression Models with Covariates — Regression Models for Univariate Analysis
# ------------------------------------------------------------

modelmulticov <- lm(formula = fi_new_standardised ~ age + sex + PC1 + PC2 + PC3 + 
                       PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                    data = PISAall)
summary(modelmulticov)

PISAllmodelmulti <- lm(formula = fi_new_standardised ~ multiscale + age + sex + PC1 + 
                          PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                       data = PISAall)
summary(PISAllmodelmulti)


# ------------------------------------------------------------
# Regression Models with Individual PRS Predictors
# ------------------------------------------------------------

predictors <- c("F6_PRS_standardised", "F5_PRS_standardised", "F4_PRS_standardised", 
                "F3_PRS_standardised", "F2_PRS_standardised", "F1_PRS_standardised", 
                "BF_PRS_standardised", "FP_PRS_standardised", "FI_PRS_standardised", "multiscale")

fit_model_extract_results <- function(predictor, cov_model) {
  formula <- as.formula(paste("Frailtyindex_standardised ~", predictor, 
                              "+ Sex + Age_years + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                              PC7 + PC8 + PC9 + PC10"))
  model <- lm(formula, data = PISAall)
  r_squared_diff <- summary(model)$r.squared - summary(cov_model)$r.squared

  model_results <- tidy(model) %>%
    filter(term == predictor) %>%
    mutate(r_squared_diff = r_squared_diff)

  conf_ints <- confint(model)[predictor, ]
  model_results <- model_results %>%
    mutate(`CI Lower` = conf_ints[1], `CI Upper` = conf_ints[2])
  return(model_results)
}

PISAcov <- lm(Frailtyindex_standardised ~ Sex + Age_years + PC1 + PC2 + PC3 + PC4 + 
                 PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PISAall)

results_list <- lapply(predictors, fit_model_extract_results, cov_model = PISAcov)
combined_results <- bind_rows(results_list) %>%
  mutate(p_value_FDR = p.adjust(p.value, method = "fdr"))

write.csv(combined_results, file = "regression_results.csv")


# ------------------------------------------------------------
# Preparing Data — Elastic Net Regression
# ------------------------------------------------------------

elasticnetPISAall <- subset(PISAall, select = c(Frailtyindex_standardised, Age_years, 
                                                 Sex, F1_PRS_standardised, F2_PRS_standardised, 
                                                 F3_PRS_standardised, F4_PRS_standardised, 
                                                 F5_PRS_standardised, F6_PRS_standardised, 
                                                 BF_PRS_standardised, PC1, PC2, PC3, PC4, 
                                                 PC5, PC6, PC7, PC8, PC9, PC10))
elasticnetPISAall <- na.omit(elasticnetPISAall)


# ------------------------------------------------------------
# Elastic Net Model Training — Elastic Net Regression
# ------------------------------------------------------------

train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 10, 
                              summaryFunction = defaultSummary)

elastic_net_cv <- train(Frailtyindex_standardised ~ ., data = elasticnetPISAall, 
                        method = "glmnet", trControl = train_control, tuneLength = 10)
best_lambda <- elastic_net_cv$bestTune$lambda

boot_results <- boot(data = elasticnetPISAall, 
                     statistic = function(data, indices) {
                       d <- data[indices, ]
                       predictors <- as.matrix(d[, -1])
                       fit <- glmnet(predictors, d$Frailtyindex_standardised, 
                                     alpha = elastic_net_cv$bestTune$alpha, 
                                     lambda = best_lambda)
                       coef(fit, s = best_lambda)[-1]
                     }, 
                     R = 1000)


# ------------------------------------------------------------
# Creating Quintiles — Quintile Analysis
# ------------------------------------------------------------

prs_vars <- c("F1_PRS_standardised", "F2_PRS_standardised", "F3_PRS_standardised", 
              "F4_PRS_standardised", "F5_PRS_standardised", "F6_PRS_standardised", 
              "FI_PRS_standardised", "FP_PRS_standardised", "multiscale")

for (prs in prs_vars) {
  quintile_col <- paste0(prs, "_quintile")
  PISAall[[quintile_col]] <- cut(
    PISAall[[prs]], 
    breaks = quantile(PISAall[[prs]], probs = seq(0, 1, 0.2), na.rm = TRUE),
    include.lowest = TRUE, 
    labels = 1:5
  )
}


# ------------------------------------------------------------
# Creating Quintiles — Quintile Analysis
# ------------------------------------------------------------

prs_vars <- c("F1_quintile", "F2_quintile", "F3_quintile", "F4_quintile", 
              "F5_quintile", "F6_quintile", "FI_quintile", "FP_quintile", "multi_quintile")

chisq_results_list <- list()

for (prs in prs_vars) {
    if (prs %in% colnames(PISAall_clean)) {
        quintile_col <- prs  

     
        contingency_table <- table(PISAall_clean[[quintile_col]], PISAall_clean$Sex)
        
    
        print(paste("Contingency table for", quintile_col))
        print(contingency_table)
        
 
        if (min(contingency_table) > 0) {
            chisq_test <- chisq.test(contingency_table)
            chisq_results_list[[prs]] <- chisq_test$p.value
        } else {
            chisq_results_list[[prs]] <- "Insufficient data for Chi-squared test"
        }
    } else {
        print(paste(prs, "not found in the dataframe."))
    }
}


print(chisq_results_list)


# ------------------------------------------------------------
# Creating Quintiles — Quintile Analysis
# ------------------------------------------------------------

prs_vars <- c("F1_quintile", "F2_quintile", "F3_quintile", "F4_quintile", 
              "F5_quintile", "F6_quintile", "FI_quintile", "FP_quintile", "multi_quintile")

Age_years_results_list <- list()

for (prs in prs_vars) {

    if (prs %in% colnames(PISAall_clean)) {

        anova_result <- aov(Age_years ~ PISAall_clean[[prs]], data = PISAall_clean)
        

        Age_years_results_list[[prs]] <- summary(anova_result)
    } else {
        print(paste(prs, "not found in the dataframe."))
        Age_years_results_list[[prs]] <- "Column not found"
    }
}


Age_years_results_list


# ------------------------------------------------------------
# Creating Quintiles — Quintile Analysis
# ------------------------------------------------------------

results_list <- list()

for (quintile_var in quintile_vars) {
  formula_str <- paste("Frailtyindex_standardised ~", 
                       paste(quintile_var, "+ Age_years + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10", sep = " "))
  
  model <- lm(as.formula(formula_str), data = PISAall)
  
  summary_model <- tidy(model)
  
  quintile_effects <- summary_model %>% 
    filter(grepl(quintile_var, term) | term %in% c("Age_years", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC10")) %>%
    mutate(quintile_var = quintile_var) %>%  #
    select(term, estimate, std.error, statistic, p.value, quintile_var)
  
  results_list[[quintile_var]] <- quintile_effects
}

results_summary <- bind_rows(results_list)

results_summary <- results_summary %>%
  mutate(fdr_p.value = p.adjust(p.value, method = "fdr"))

print(results_summary)


# Set an output filename before running:
# write.csv(results_summary, file = "quintile_regression_results.csv")


# ------------------------------------------------------------
# Creating Quintiles — Quintile Analysis
# ------------------------------------------------------------

FactorNames = names(Filter(is.factor, df))
df[FactorNames] = lapply(df[FactorNames], as.character)
df[FactorNames] = lapply(df[FactorNames], as.numeric)


# ------------------------------------------------------------
# Creating Quintiles — Quintile Analysis
# ------------------------------------------------------------

library('dplyr')
library('lavaan')

df_mod <- mutate(df,
                   blkdes_w1 = blkdes_w1/2,
                   blkdes_w2 = blkdes_w2/2,
                   blkdes_w3 = blkdes_w3/2,
                   blkdes_w4 = blkdes_w4/2,
                   blkdes_w5 = blkdes_w5/2,
                   vftot_w1 = vftot_w1/2,
                   vftot_w2 = vftot_w2/2,
                   vftot_w3 = vftot_w3/2,
                   vftot_w4 = vftot_w4/2,
                   vftot_w5 = vftot_w5/2,
                   lmtotal_w1 = lmtotal_w1/3,
                   lmtotal_w2 = lmtotal_w2/3,
                   lmtotal_w3 = lmtotal_w3/3,
                   lmtotal_w4 = lmtotal_w4/3,
                   lmtotal_w5 = lmtotal_w5/3,
                   digback_w1 = 3*digback_w1,
                   digback_w2 = 3*digback_w2,
                   digback_w3 = 3*digback_w3,
                   digback_w4 = 3*digback_w4,
                   digback_w5 = 3*digback_w5,
                   digsym_w1 = digsym_w1/2,
                   digsym_w2 = digsym_w2/2,
                   digsym_w3 = digsym_w3/2,
                   digsym_w4 = digsym_w4/2,
                   digsym_w5 = digsym_w5/2,
                   ittotal_w1 = ittotal_w1/2,
                   ittotal_w2 = ittotal_w2/2,
                   ittotal_w3 = ittotal_w3/2,
                   ittotal_w4 = ittotal_w4/2,
                   ittotal_w5 = ittotal_w5/2,
                   crtmean_w1 = -50 * crtmean_w1,
                   crtmean_w2 = -50 * crtmean_w2,
                   crtmean_w3 = -50 * crtmean_w3,
                   crtmean_w4 = -50 * crtmean_w4,
                   crtmean_w5 = -50 * crtmean_w5)

#==================================================================
# LBC1936 factor of curves model
#==================================================================
model <- '
# test growth curves
Imatreas =~ 1*matreas_w1 + 1*matreas_w2 + 1*matreas_w3 + 1*matreas_w4 + 1*matreas_w5
Smatreas =~ 0*matreas_w1 + 2.98*matreas_w2 + 6.75*matreas_w3 + 9.82*matreas_w4 + 12.54*matreas_w5

Iblkdes =~ 1*blkdes_w1 + 1*blkdes_w2 + 1*blkdes_w3 + 1*blkdes_w4 + 1*blkdes_w5
Sblkdes=~ 0*blkdes_w1 + 2.98*blkdes_w2 + 6.75*blkdes_w3 + 9.82*blkdes_w4 + 12.54*blkdes_w5

Ispantot =~ 1*spantot_w1 + 1*spantot_w2 + 1*spantot_w3 + 1*spantot_w4 + 1*spantot_w5
Sspantot=~ 0*spantot_w1 + 2.98*spantot_w2 + 6.75*spantot_w3 + 9.82*spantot_w4 + 12.54*spantot_w5

Inart =~ 1*nart_w1 + 1*nart_w2 + 1*nart_total_w3 + 1*nart_total_w4 + 1*nart_total_w5
Snart =~ 0*nart_w1 + 2.98*nart_w2 + 6.75*nart_total_w3 + 9.82*nart_total_w4 + 12.54*nart_total_w5

Iwtar =~ 1*wtar_w1 + 1*wtar_w2 + 1*wtar_total_w3 + 1*wtar_total_w4 + 1*wtar_total_w5
Swtar =~ 0*wtar_w1 + 2.98*wtar_w2 + 6.75*wtar_total_w3 + 9.82*wtar_total_w4 + 12.54*wtar_total_w5

Ivftot =~ 1*vftot_w1 + 1*vftot_w2 + 1*vftot_w3 + 1*vftot_w4 + 1*vftot_w5
Svftot =~ 0*vftot_w1 + 2.98*vftot_w2 + 6.75*vftot_w3 + 9.82*vftot_w4 + 12.54*vftot_w5

Ivpatotal =~ 1*vpatotal_w1 + 1*vpatotal_w2 + 1*vpatotal_w3 + 1*vpatotal_w4 + 1*vpa_total_w5
Svpatotal =~ 0*vpatotal_w1 + 2.98*vpatotal_w2 + 6.75*vpatotal_w3 + 9.82*vpatotal_w4 + 12.54*vpa_total_w5

Ilmtotal =~ 1*lmtotal_w1 + 1*lmtotal_w2 + 1*lmtotal_w3 + 1*lmtotal_w4 + 1*lmtotal_w5
Slmtotal =~ 0*lmtotal_w1 + 2.98*lmtotal_w2 + 6.75*lmtotal_w3 + 9.82*lmtotal_w4 + 12.54*lmtotal_w5

Idigback =~ 1*digback_w1 + 1*digback_w2 + 1*digback_w3 + 1*digback_w4 + 1*digback_w5
Sdigback =~ 0*digback_w1 + 2.98*digback_w2 + 6.75*digback_w3 + 9.82*digback_w4 + 12.54*digback_w5

Isymsear =~ 1*symsear_w1 + 1*symsear_w2 + 1*symsear_w3 + 1*symsear_w4 + 1*symsear_w5
Ssymsear =~ 0*symsear_w1 + 2.98*symsear_w2 + 6.75*symsear_w3 + 9.82*symsear_w4 + 12.54*symsear_w5

Idigsym =~ 1*digsym_w1 + 1*digsym_w2 + 1*digsym_w3 + 1*digsym_w4 + 1*digsym_w5
Sdigsym =~ 0*digsym_w1 + 2.98*digsym_w2 + 6.75*digsym_w3 + 9.82*digsym_w4 + 12.54*digsym_w5

Iittotal =~ 1*ittotal_w1 + 1*ittotal_w2 + 1*ittotal_w3 + 1*ittotal_w4 + 1*ittotal_w5
Sittotal =~ 0*ittotal_w1 + 2.98*ittotal_w2 + 6.75*ittotal_w3 + 9.82*ittotal_w4 + 12.54*ittotal_w5

Icrtmean =~ 1*crtmean_w1 + 1*crtmean_w2 + 1*crtmean_w3 + 1*crtmean_w4 + 1*crtmean_w5
Scrtmean =~ 0*crtmean_w1 + 2.98*crtmean_w2 + 6.75*crtmean_w3 + 9.82*crtmean_w4 + 12.54*crtmean_w5

# latent g intercept and slope 
Ig =~  Iblkdes + Imatreas  + Ispantot + Ivftot + Ivpatotal + Ilmtotal +
  Idigback + Isymsear + Idigsym + Iittotal + Icrtmean + Inart + Iwtar

Sg =~ Sblkdes + Smatreas + Sspantot + Svftot + Svpatotal + Slmtotal +
  Sdigback + Ssymsear + Sdigsym + Sittotal + Scrtmean + Snart + Swtar

#indicator as scaling reference: loading=1, int=0
Iblkdes ~ 0*1
Sblkdes ~ 0*1 

# within-wave covariances between nart and wtar
nart_w1 ~~ wtar_w1
nart_w2 ~~ wtar_w2
nart_total_w3 ~~ wtar_total_w3
nart_total_w4 ~~ wtar_total_w4
nart_total_w5 ~~ wtar_total_w5

# within-test intercept-slope covariances
Imatreas ~~ Smatreas
Iblkdes ~~ Sblkdes
#Ispantot ~~Sspantot
Inart ~~ Snart
Iwtar ~~ Swtar
Ivftot ~~ Svftot
Ivpatotal ~~ Svpatotal
Ilmtotal ~~ Slmtotal
Idigback ~~ Sdigback
Isymsear ~~ Ssymsear
Idigsym ~~ Sdigsym
Iittotal ~~ Sittotal
Icrtmean ~~ Scrtmean


# within-domain intercept-intercept and slope-slope covariances
Iblkdes ~~ Imatreas # Visuospatial domain
Iblkdes ~~ Ispantot
Imatreas ~~ Ispantot
Sblkdes ~~ Smatreas 
#Sblkdes ~~ Sspantot
#Smatreas ~~ Sspantot

Inart ~~ Ivftot #Crystalized domain
Iwtar ~~ Ivftot
Iwtar ~~ Inart
Snart ~~ Svftot
Swtar ~~ Svftot
Swtar ~~ Snart

Ilmtotal ~~ Ivpatotal # Verbal memory domain
Ilmtotal ~~ Idigback
Ivpatotal ~~ Idigback
Slmtotal ~~ Svpatotal
Slmtotal ~~ Sdigback
Svpatotal ~~ Sdigback

Iittotal ~~ Idigsym #Processing speed domain
Iittotal ~~ Isymsear
Iittotal ~~ Icrtmean
Idigsym ~~ Isymsear
Idigsym ~~ Icrtmean
Isymsear ~~ Icrtmean
Sittotal ~~ Sdigsym 
Sittotal ~~ Ssymsear
Sittotal ~~ Scrtmean
Sdigsym ~~ Ssymsear
Sdigsym ~~ Scrtmean
Ssymsear ~~ Scrtmean

#fixed negative residual variance to 0 
Sspantot ~~ 0*Sspantot
'

fit <- growth(model = model, df_mod,  missing = "ml.x")
summary(fit, fit.measures = T, standardized = T)
print(model)


# ------------------------------------------------------------
# Creating Quintiles — Quintile Analysis
# ------------------------------------------------------------

cogscores <- df_mod %>% select(contains("matreas") | contains("blkdes") | contains("spantot") | contains("nart") | contains("vftot") | contains("vpa") | contains("lmtotal") | contains("digback") | contains("symsear") | contains("digsym") | contains("ittotal") | contains("crtmean") )
wavesindex <- rep(c(1,2,3,4,5), 12) # index the waves
 
# find out which rows have data for only one and for two+ waves
w1 <- matrix(NA, dim(cogscores)[1], 1)
w2 <- matrix(NA, dim(cogscores)[1], 1)
w3 <- matrix(NA, dim(cogscores)[1], 1)
w4 <- matrix(NA, dim(cogscores)[1], 1)
w5 <- matrix(NA, dim(cogscores)[1], 1)
for (i in 1:dim(cogscores)[1]) {
if (!all(is.na(cogscores[i, which(wavesindex == 1)]))) {
w1[i,1] <- 1 }
if (!all(is.na(cogscores[i, which(wavesindex == 2)]))) {
w2[i,1] <- 1 }
if (!all(is.na(cogscores[i, which(wavesindex == 3)]))) {
w3[i,1] <- 1 }
if (!all(is.na(cogscores[i, which(wavesindex == 4)]))) {
w4[i,1] <- 1 }
if (!all(is.na(cogscores[i, which(wavesindex == 5)]))) {
w5[i,1] <- 1 }
}
 
wavesample <- cbind(w1, w2, w3, w4, w5)
rowSums(wavesample, na.rm = T)
onewave <- which(rowSums(wavesample, na.rm = T) < 2) # these people were only tested at one wave
 
# extract factor scores
factorScores <- data.frame(lavPredict(fit, df))
factorScores$Sg[onewave] <- NA # set slopes for participants who were only tested at one wave to NA
factorScores$Ig[which(is.na(w1))] <- NA # set intercepts for participants who were not tested at w1 to NA
LBC1936FactorScores <- data.frame(ID = df[1], Ig = factorScores$Ig, Sg = factorScores$Sg)
 
cbind(c("N intercepts included", "N intercepts excluded", "N slopes included", "N slopes excluded"), as.numeric(c(length(which(is.na(LBC1936FactorScores$Ig) == F)), length(which(is.na(LBC1936FactorScores$Ig) == T)), length(which(is.na(LBC1936FactorScores$Sg) == F)), length(which(is.na(LBC1936FactorScores$Sg) == T)))))


# ------------------------------------------------------------
# motoric cognitive risk — Quintile Analysis
# ------------------------------------------------------------

# Illustrative Cox regression template from the original notebook.
# Replace the placeholder covariates before running.
# for (pred in predictors) {
#   formula <- as.formula(paste(
#     "Surv(event_time, event_status) ~", pred,
#     "+ Age_years + Sex + PC1 + [additional covariates]"
#   ))
#   model <- coxph(formula, data = subset_data)
# }
