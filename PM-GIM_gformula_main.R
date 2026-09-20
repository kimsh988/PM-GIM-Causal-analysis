
################################################################################
# Packages
################################################################################
library(dplyr)
library(haven)
library(tidyverse)
library(gfoRmula)
library(data.table)

################################################################################
# Load data
################################################################################
data1 <- readRDS("data/data1")
data2 <- readRDS("data/data2")

################################################################################
# Data Preprocessing 1
################################################################################
data <- data1 %>% 
  left_join(
    data2 %>% 
      select(-all_of(setdiff(intersect(names(data1), names(data2)), c("patient_id", "dov")))),
    by = c("patient_id", "dov")
  )

data$hpylori_med <- ifelse(data$hpylori_med == "1", 1, 0)

data <- data %>% 
  mutate(
    education = replace(education, is.na(education), 9),
    smk = replace(smk, is.na(smk), 9),
    ipaq_cat = replace(ipaq_cat, is.na(ipaq_cat), 9),
    alc_cat = case_when(
      is.na(alc_amount_grams) ~ 9,
      alc_amount_grams == 0 ~ 0,
      alc_amount_grams > 0 & alc_amount_grams <= 10 ~ 1,
      alc_amount_grams > 10 & alc_amount_grams <= 20 ~ 2,
      alc_amount_grams > 20 ~ 3),
    diabetes_gluc = case_when(
      (!is.na(gluc) & gluc >= 126) | history_diabetes == 1 | med_diabetes == 1 ~ 1,
      !is.na(gluc) ~ 0,
      TRUE ~ NA_real_
    ),
    diabetes_hba1c = case_when(
      (!is.na(hba1c) & hba1c >= 6.5) | history_diabetes == 1 | med_diabetes == 1 ~ 1,
      !is.na(hba1c) ~ 0,
      TRUE ~ NA_real_
    ),
    diabetes = ifelse(
      is.na(diabetes_gluc) & is.na(diabetes_hba1c), NA,
      ifelse(
        (!is.na(diabetes_gluc) & diabetes_gluc == 1) |
          (!is.na(diabetes_hba1c) & diabetes_hba1c == 1), 1, 0)
    ),
    rename(
      hpylori_positive = fdi9_2,
      hpylori_med = fdi9_3
    )
  )

# Analysis cohort
data <- data %>% 
  filter(dov >= as.Date("2011-01-01")) %>%    # after 2011 
  filter(dov < as.Date("2020-01-01"))         # until 2019

# Exclusion rule - patients (cases)
# Missing endoscopy data at baseline
data <- data %>% 
  mutate(base_endoscopy_na = as.integer(
    patient_id %in%
      (data %>% 
         arrange(patient_id, dov) %>% 
         group_by(patient_id) %>% 
         slice(1) %>% 
         filter(rowSums(is.na(across(endoscopy1:endoscopy5))) == 5 |
                  egd_incomplete == 1 |
                  egd_missing == 1) %>% 
         pull(patient_id))
  ))

length(unique(data$patient_id[data$base_endoscopy_na == 1])) # patients number
sum(data$base_endoscopy_na == 1, na.rm = TRUE) # case number

# History of stomach cancer before the visit date
data <- data %>% 
  mutate(base_stomach_ca = as.integer(
    patient_id %in% 
      (data %>% 
         arrange(patient_id, dov) %>% 
         group_by(patient_id) %>% 
         slice(1) %>% 
         filter(history_cancer_stomach == 1) %>% 
         pull(patient_id))
  ))

length(unique(data$patient_id[data$base_stomach_ca == 1])) 
sum(data$base_stomach_ca == 1, na.rm = TRUE) 

# History of cancer at baseline
data <- data %>% 
  mutate(base_history_cancer = as.integer(
    patient_id %in% 
      (data %>% 
         arrange(patient_id, dov) %>% 
         group_by(patient_id) %>% 
         slice(1) %>% 
         filter(history_cancer == 1) %>% 
         pull(patient_id))
  ))

length(unique(data$patient_id[data$base_history_cancer == 1])) 
sum(data$base_history_cancer == 1, na.rm = TRUE) 

# Missing BMI data at baseline
data <- data %>% 
  mutate(base_bmi_na = as.integer(
    patient_id %in% 
      (data %>% 
         arrange(patient_id, dov) %>% 
         group_by(patient_id) %>% 
         slice(1) %>% 
         filter(is.na(bmi)) %>% 
         pull(patient_id))
  ))

length(unique(data$patient_id[data$base_bmi_na == 1])) 
sum(data$base_bmi_na == 1, na.rm = TRUE) 

# Data exclusion
data3 <- data %>% 
  filter(base_endoscopy_na == 0) %>% 
  filter(base_stomach_ca == 0) %>% 
  filter(base_history_cancer == 0) %>% 
  filter(base_bmi_na == 0) %>% 
  filter(base_ag_im == 0)

# Outcome variable
# Atrophic gastritis (AG
data3$AG <- ifelse(
  data3$egd_st_gastritis_ch_atrophic == 1, 1, 0)

# Intestinal metaplasia (IM)
data3$IM <- ifelse(
  rowSums(data3[,paste0("endobiopsy", 1:6)] == 42, na.rm = TRUE) > 0 | 
    data3$egd_st_gastritis_ch_metapl == 1, 1, 0)

# Save data
saveRDS(data3, file = "data3.rds")

################################################################################
# Data Preprocessing 2
################################################################################

# Select variables for main analysis
cov <- c("age", "sex", "education", "center", "yov", "smk", "alc_cat", "ipaq_cat",
         "bmi", "hypertension", "diabetes", "family_cancer_stomach", "hpylori_med",
         "med_hyperlipidemia")
air <- c("TEMP_365d", "RH_365d")
exposure <- c("pm2.5_365d", "pm10_365d")
outcome <- c("IM")

# Data frame
df <- data3[,c("patient_id", "dov", exposure, outcome, cov, air)]

# Delete NA, make survival data
df <- df %>% 
  filter(if_all(everything(), ~!is.na(.))) %>% 
  arrange(patient_id, dov) %>% 
  group_by(patient_id) %>% 
  mutate(visit_order = row_number() - 1) %>% 
  mutate(first_event_order = ifelse(any(IM == 1),
                                    min(visit_order[IM == 1]),
                                    Inf)) %>% 
  filter(visit_order <= first_event_order) %>% 
  select(-first_event_order) %>% 
  ungroup()

df <- df[, c("patient_id", "dov", "visit_order", exposure, outcome, cov, air)]

# Base covariates
df_v2 <- df %>% 
  group_by(patient_id) %>% 
  arrange(dov, .by_group = TRUE) %>% 
  mutate(family_cancer_stomach = {
    y <- family_cancer_stomach
    if (all(is.na(y))) NA else as.integer(any(y == 1, na.rm = TRUE))},
    base_age = first(age),
    base_yov = first(yov),
    education = first(education),
    center = first(center)
    ) %>% 
  ungroup()

# Retain only the earliest visit within the same year
df_v2 <- df_v2 %>% 
  group_by(patient_id, yov) %>% 
  arrange(dov, .by_group = TRUE) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(patient_id) %>% 
  arrange(dov, .by_group = TRUE) %>% 
  mutate(visit_order = row_number() - 1) %>% 
  ungroup()

# Save survival data
saveRDS(df_v2, file = "df_v2.rds")

################################################################################
# Data Preprocessing 3
################################################################################
df_v3 <- df_v2 %>% 
  group_by(patient_id) %>% 
  group_modify(~{dat <- .x %>% arrange(yov)
  
  start_yov <- min(dat$yov, na.rm = TRUE)
  
  # 1. Year of the first occurrence of the event
  stop_im <- if (any(dat$IM == 1, na.rm = TRUE)) {
    min(dat$yov[dat$IM == 1], na.rm = TRUE)
  } else {Inf}
  
  # 2. If the gap between consecutive observed visits exceeds 3 years,
  #    allow follow-up only through 2 years after the previous visit
  gap_idx <- which(diff(dat$yov) > 3)
  
  stop_gap <- if (length(gap_idx) > 0) {
    dat$yov[gap_idx[1]] + 2
  } else {Inf}
  
  # 3. If no event occurs at the last observed visit,
  #    allow follow-up only through 2 years after the last visit
  last_yov <- max(dat$yov, na.rm = TRUE)
  last_event <- dat$IM[which.max(dat$yov)]
  
  stop_yov <- if (!is.na(last_event) && last_event == 0) {
    last_yov + 2
  } else {Inf}
  
  # Determine the final year of follow-up
  end_yov <- min(2019, stop_im, stop_gap, stop_yov)
  
  dat %>% 
    filter(yov <= end_yov) %>% 
    complete(yov == seq(start_yov, end_yov, by = 1)) %>% 
    arrange(yov) %>% 
    mutate(visit = if_else(is.na(dov), 0L, 1L),
           t_yov = yov - start_yov) %>% 
    fill(-c(yov, dov, visit, t_yov), .direction = "down")
  }) %>% 
  ungroup() %>% 
  arrange(patient_id, yov)

# Save the processed data
saveRDS(df_v3, file = "df_v3.rds")

################################################################################
# Main analysis - Parametric gformula
################################################################################

# Setting
setDT(df_v3)
df_v3[, sex := as.integer(haven::zap_labels(sex))]
df_v3[, education := factor(haven::zap_labels(education))]
df_v3[, center := as.integer(haven::zap_labels(center))]
df_v3[, smk := factor(haven::zap_labels(smk))]
df_v3[, yov := factor(haven::zap_labels(yov))]
df_v3[, alc_cat := factor(haven::zap_labels(alc_cat))]
df_v3[, ipaq_cat := factor(haven::zap_labels(ipaq_cat))]
df_v3[, hypertension := as.integer(haven::zap_labels(hypertension))]
df_v3[, diabetes := as.integer(haven::zap_labels(diabetes))]
df_v3[, med_hyperlipidemia := as.integer(haven::zap_labels(med_hyperlipidemia))]
sapply(df_v3, class)

# Parameter
id <- "patient_id"
time_points <- 9
time_name <- "t_yov"
basecovs <- c("sex", "education", "center", "base_age", "base_yov", "family_cancer_stomach")
covnames <- c("visit",
              "TEMP_365d", "RH_365d", "pm2.5_365d",
              "smk", "alc_cat", "ipaq_cat",
              "bmi", "hypertension", "diabetes", "med_hyperlipidemia",
              "hpylori_med")
covtypes <- c("binary",
              "normal", "bounded normal", "normal",
              "categorical", "categorical", "categorical",
              "normal", "binary", "binary", "binary",
              "binary")
visitprocess <- list(c('visit', "TEMP_365d",3),
                     c('visit', "RH_365d",3),
                     c('visit', "pm2.5_365d",3),
                     c('visit', "smk",3),
                     c('visit', "alc_cat",3),
                     c('visit', "ipaq_cat",3),
                     c('visit', "bmi",3),
                     c('visit', "hypertension",3),
                     c('visit', "diabetes",3),
                     c('visit', "hpylori_med",3),
                     c('visit', "med_hyperlipidemia",3))
outcome_name <- "IM"
outcome_type <- "survival"
histories <- c(lagged)
histvars <- list(covnames)

# Time-varying covariates regression models
covparams <- list(covmodels = c(
  # visit
  visit ~ lag1_visit + factor(t_yov),
  
  # TEMP
  TEMP_365d ~ lag1_TEMP_365d + RH_365d + lag1_visit + factor(t_yov),
  
  # Relative humidity
  RH_365d = RH_365d ~ lag1_RH_365d + TEMP_365d + lag1_visit + factor(t_yov),
  
  # PM2.5
  pm2.5_365d ~ lag1_pm2.5_365d + TEMP_365d + lag1_visit + factor(t_yov),
  
  # Smoking
  smk ~ lag1_smk +
    lag1_diabetes + lag1_hypertension + lag1_med_hyperlipidemia +
    alc_cat + ipaq_cat + bmi +
    base_age + sex + lag1_visit + factor(t_yov),
  
  # Alcohol
  alc_cat ~ lag1_alc_cat +
    lag1_diabetes + lag1_hypertension + lag1_med_hyperlipidemia +
    smk + ipaq_cat + bmi +
    base_age + sex + lag1_visit + factor(t_yov),
  
  # Physical activity
  ipaq_cat ~ lag1_ipaq_cat +
    lag1_diabetes + lag1_hypertension + lag1_med_hyperlipidemia +
    smk + alc_cat + bmi +
    base_age + sex + lag1_visit + factor(t_yov),
  
  # BMI
  bmi ~ lag1_bmi +
    lag1_diabetes + lag1_hypertension + lag1_med_hyperlipidemia +
    smk + alc_cat + ipaq_cat +
    base_age + sex + lag1_visit + factor(t_yov),
  
  # Hypertension
  hypertension ~ lag1_hypertension + 
    lag1_diabetes + lag1_med_hyperlipidemia +
    smk + alc_cat + ipaq_cat + bmi +
    base_age + sex + lag1_visit + factor(t_yov),
  
  # Diabetes
  diabetes ~ lag1_diabetes +
    lag1_hypertension + lag1_med_hyperlipidemia +
    smk + alc_cat + ipaq_cat + bmi +
    base_age + sex + lag1_visit + factor(t_yov),
  
  # Medication of hyperlipidemia
  med_hyperlipidemia ~ lag1_med_hyperlipidemia +
    lag1_diabetes + lag1_hypertension +
    smk + alc_cat + ipaq_cat + bmi +
    base_age + sex + lag1_visit + factor(t_yov),
  
  # Medication of helicobacter pylori
  hpylori_med ~ lag1_hpylori_med +
    smk + alc_cat + ipaq_cat + bmi +
    base_age + sex + lag1_visit + factor(t_yov)
    
))

# Outcome regression model
ymodel <- IM ~ lag1_pm2.5_365d + lag1_smk + lag1_alc_cat + lag1_ipaq_cat + lag1_bmi +
  lag1_hypertension + lag1_diabetes + lag1_hpylori_med + lag1_med_hyperlipidemia +
  lag1_TEMP_365d + lag1_RH_365d + lag1_visit +
  pm2.5_365d + base_age + base_yov + education + sex + center + smk + alc_cat +
  ipaq_cat + bmi + hypertension + diabetes + hpylori_med + med_hyperlipidemia +
  family_cancer_stomach + TEMP_365d + RH_365d + factor(t_yov) + visit

# Intervention
intvars <- list(c("pm2.5_365d"), c("pm2.5_365d"))
interventions <- list(list(c(threshold, 0, 15)),
                      list(c(threshold, 15, Inf)))
int_descript <- c("PM2.5 cap at 15", "PM2.5 above 15")

nsimul <- 500000
ncores <- 10
print(Sys.time())
gform <- gformula(obs_data = df_v3,
                  id = id,
                  time_points = time_points,
                  time_name = time_name,
                  covnames = covnames,
                  outcome_name = outcome_name,
                  outcome_type = outcome_type,
                  covtypes = covtypes,
                  visitprocess = visitprocess,
                  covparams = covparams,
                  ymodel = ymodel,
                  intvars = intvars,
                  interventions = interventions,
                  int_descript = int_descript,
                  histories = histories,
                  histvars = histvars,
                  basecovs = basecovs,
                  seed = 1234,
                  parallel = FALSE,
                  nsamples = 200,
                  nsimul = nsimul,
                  boot_diag = TRUE,
                  show_progress = TRUE)
print(Sys.time())
print(gform)

saveRDS(gform, file = "gform.rds")
