
library(tidyverse)
library(ggpubr)
library(survival)
library(survminer)

############################################################
## Input data and selecting relevant columns
############################################################

# Read tab-delimited file with ages and phenotype diagnosis info
df_age_diag <- read.delim('all_data_with_diag_age_stringent_permissive_wgs_only.tsv', stringsAsFactors = FALSE)

df_plot_1 <- df_age_diag %>%
  select(all_of(c('eid', 'recruit_age', 'dob', 'stringent_carrier', 
                  'stringent_final', 'stringent_agediag',
                  'permissive_final', 'permissive_agediag')))

#Reference year is the year at the last HES/cancer registry update occured
# Reference date representing “last known follow-up / registry update” for HES/cancer records.
# Used as a proxy for the censoring date.
reference_year <- ymd('2023-07-01')
df_age_diag$dob <- as.Date(df_age_diag$dob) # Ensure dob is a Date object 


# Compute age at reference date: floor of years between dob and reference date.
# This gives each participant’s last observed age before censoring.
df_plot_1 <- df_plot_1 %>%
  mutate(reference_age = floor(interval(dob, reference_year)/years(1)))

# Reorder columns for readability
df_plot_1 <- df_plot_1 %>%
  relocate('dob', .after = 'eid') %>%
  relocate('stringent_agediag', .after = 'recruit_age') %>%
  relocate('stringent_final', .after = 'stringent_agediag') %>%
  relocate('permissive_agediag', .after = 'stringent_final') %>%
  relocate('permissive_final', .after = 'permissive_agediag') %>%
  relocate('reference_age', .after = 'permissive_final')

# Replace NA values in Age_at_Diagnosis with the last observed age before censoring

# IMPORTANT: Here you are treating NA “age at diagnosis” as “no event observed yet”,
# and setting time-to-event = censoring age. This is a common approach.
df_plot_1$stringent_agediag[is.na(df_plot_1$stringent_agediag)] <- df_plot_1$reference_age[is.na(df_plot_1$stringent_agediag)]
df_plot_1$permissive_agediag[is.na(df_plot_1$permissive_agediag)] <- df_plot_1$reference_age[is.na(df_plot_1$permissive_agediag)]


############################################################
## Survival objects and Kaplan–Meier fits
############################################################


# Surv(time, event) expects:
# - time: numeric follow-up time (here: age at event or age at censoring)
# - event: 1 = event occurred, 0 = censored
surv_object <- Surv(time = df_plot_1$stringent_agediag, event = df_plot_1$stringent_final)
surv_object_3 <- Surv(time = df_plot_1$permissive_agediag, event = df_plot_1$permissive_final)

# Fit Kaplan-Meier model with Genotype as grouping variable

# Fit KM curves stratified by carrier status.
# stringent_carrier must be coded appropriately (often 0/1 or factor).
km_fit <- survfit(surv_object ~ stringent_carrier, data = df_plot_1)
km_fit_3 <- survfit(surv_object_3 ~ stringent_carrier, data = df_plot_1)


survdiff(Surv(time = df_plot_1$stringent_agediag, event = df_plot_1$stringent_final) ~ stringent_carrier, data = df_plot_1)



############################################################
## Extract stats and time-point summaries (Stringent)
############################################################

KM_fit_stats <- survdiff(Surv(time = stringent_agediag, event = stringent_final) ~ stringent_carrier, data = df_plot_1)
chisq <- KM_fit_stats$chisq
df <- length(KM_fit_stats$n) - 1
p_value <- 1 - pchisq(chisq, df)

# Survival estimates at specific ages (18, 40, 60).
# summary(survfit, times=...) gives survival probability S(t) for each stratum at t.
summary_us_18 <- summary(km_fit, times = 18)
summary_us_40 <- summary(km_fit, times = 40)
summary_us_60 <- summary(km_fit, times = 60)

# Convert survival to cumulative incidence-like quantity:
# proportion_with_event = 1 - S(t).
# Also invert CI bounds accordingly.
df_summary_us_18 <- data.frame(
  group = summary_us_18$strata,
  time = summary_us_18$time,
  survival = summary_us_18$surv,
  lower = summary_us_18$lower,
  upper = summary_us_18$upper
) %>%
  mutate(
    proportion_with_event = 1 - survival,
    lower_ci = 1 - upper,
    upper_ci = 1 - lower
  )

df_summary_us_40 <- data.frame(
  group = summary_us_40$strata,
  time = summary_us_40$time,
  survival = summary_us_40$surv,
  lower = summary_us_40$lower,
  upper = summary_us_40$upper
) %>%
  mutate(
    proportion_with_event = 1 - survival,
    lower_ci = 1 - upper,
    upper_ci = 1 - lower
  )

df_summary_us_60 <- data.frame(
  group = summary_us_60$strata,
  time = summary_us_60$time,
  survival = summary_us_60$surv,
  lower = summary_us_60$lower,
  upper = summary_us_60$upper
) %>%
  mutate(
    proportion_with_event = 1 - survival,
    lower_ci = 1 - upper,
    upper_ci = 1 - lower
  )



############################################################
## Extract stats and time-point summaries (Permissive)
############################################################

KM_fit_stats_3 <- survdiff(Surv(time = permissive_agediag, event = permissive_final) ~ stringent_carrier, data = df_plot_1)
chisq_3 <- KM_fit_stats_3$chisq
df_3 <- length(KM_fit_stats_3$n) - 1
p_value_3 <- 1 - pchisq(chisq_3, df_3)

# Survival estimates at specific ages (18, 40, 60).
# summary(survfit, times=...) gives survival probability S(t) for each stratum at t.
summary_p_18 <- summary(km_fit_3, times = 18)
summary_p_40 <- summary(km_fit_3, times = 40)
summary_p_60 <- summary(km_fit_3, times = 60)
summary_p_80 <- summary(km_fit_3, times = 80)


# Convert survival to cumulative incidence-like quantity:
# proportion_with_event = 1 - S(t).
# Also invert CI bounds accordingly.
df_summary_p_18 <- data.frame(
  group = summary_p_18$strata,
  time = summary_p_18$time,
  survival = summary_p_18$surv,
  lower = summary_p_18$lower,
  upper = summary_p_18$upper
) %>%
  mutate(
    proportion_with_event = 1 - survival,
    lower_ci = 1 - upper,
    upper_ci = 1 - lower
  )

df_summary_p_40 <- data.frame(
  group = summary_p_40$strata,
  time = summary_p_40$time,
  survival = summary_p_40$surv,
  lower = summary_p_40$lower,
  upper = summary_p_40$upper
) %>%
  mutate(
    proportion_with_event = 1 - survival,
    lower_ci = 1 - upper,
    upper_ci = 1 - lower
  )

df_summary_p_60 <- data.frame(
  group = summary_p_60$strata,
  time = summary_p_60$time,
  survival = summary_p_60$surv,
  lower = summary_p_60$lower,
  upper = summary_p_60$upper
) %>%
  mutate(
    proportion_with_event = 1 - survival,
    lower_ci = 1 - upper,
    upper_ci = 1 - lower
  )

df_summary_p_80 <- data.frame(
  group = summary_p_80$strata,
  time = summary_p_80$time,
  survival = summary_p_80$surv,
  lower = summary_p_80$lower,
  upper = summary_p_80$upper
) %>%
  mutate(
    proportion_with_event = 1 - survival,
    lower_ci = 1 - upper,
    upper_ci = 1 - lower
  )


