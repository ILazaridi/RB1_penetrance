
library(tidyverse)
library(lubridate)

# ------------------------------------------------------------
# SECTION 1 — Load cohorts: carriers/non-carriers (with EHR)
# ------------------------------------------------------------

# Variant carriers WITH EHR data 
all_variant_carriers_wEHR <- read_csv('person_02917100_000000000000.csv') %>%
  mutate(date_of_birth = ymd_hms(date_of_birth, tz = "UTC")) %>%
  mutate(all_carrier=1)

# Non-carriers WITH EHR data
non_carriers_wEHR <- read_csv('person_81330221_000000000000.csv') %>%
  mutate(date_of_birth = ymd_hms(date_of_birth, tz = "UTC"))


# ------------------------------------------------------------
# SECTION 2 — Carriers after QC filtering
# ------------------------------------------------------------

#Create list of carriers after IGV plot and QC

carriers_after_qc_strict <- read_csv('person_10139945_000000000000.csv') %>% #One individual not included due to lack of EHR data, instead of poor quality read
  mutate(qcd_carrier_strict=1)

# ------------------------------------------------------------
# SECTION 3 — Build phenotype cohort: RB and/or ocular cancer
# ------------------------------------------------------------

#Create phenotype dataframes
RB_phenotype <- read_csv('person_29076054_000000000000.csv') %>%
  mutate(date_of_birth = ymd_hms(date_of_birth, tz = "UTC")) %>%
  mutate(Rb_phenotype = 1)

# Condition records for individuals with Rb and/or ocular cancer
RB_phenotype_condition <- read_csv('condition_29076054_000000000000.csv') %>%
  
  left_join(RB_phenotype %>% select(person_id, date_of_birth), by = "person_id") %>%
  
  # compute age at diagnosis in years
  mutate(condition_start_datetime = ymd_hms(condition_start_datetime, tz = "UTC"),
         condition_end_datetime = ymd_hms(condition_end_datetime, tz = "UTC")) %>%
  
  mutate(age_diag = time_length(interval(date_of_birth, condition_start_datetime), unit="years")) %>%
  
  mutate(age_diag = floor(age_diag)) %>%
  
  relocate('date_of_birth', .after = 'condition_end_datetime') %>%
  relocate('age_diag', .after = 'date_of_birth') %>%
  arrange(age_diag)

# Survey records
RB_phenotype_survey <- read_csv('survey_29076054_000000000000.csv')

# Split a single "answer" field into "question_only" and "answer_only"
RB_phenotype_survey[c('question_only', 'answer_only')] <- str_split_fixed(RB_phenotype_survey$answer, "-", 2)

# ------------------------------------------------------------
# SECTION 4 — Reduce to unique person-level indicators
# ------------------------------------------------------------


# Condition-based indicator: keep one record per person
RB_phenotype_cohort_condition_unique <- RB_phenotype_condition %>%
  distinct(person_id, .keep_all = TRUE) %>%
  select(person_id, age_diag) %>%
  mutate(cohort_condition_Rb = 1)

# Survey-based indicator: filter relevant question regarding to age of diagnosis
RB_phenotype_survey_b <- RB_phenotype_survey %>% 
  filter(question_only == "About how old were you when you were first told you had eye cancer? ")

# ------------------------------------------------------------
# SECTION 5 — Map categorical age bins to numeric ages
# ------------------------------------------------------------


for (i in seq_len(nrow(RB_phenotype_survey_b))) {
  RB_phenotype_survey_b$self_rep_age_max <- ifelse(RB_phenotype_survey_b$answer_only == " Child (0-11)", 
                                                                as.numeric(11), 
                                                                ifelse(RB_phenotype_survey_b$answer_only == " Adolescent (12-17)", 
                                                                       as.numeric(17), 
                                                                       ifelse(RB_phenotype_survey_b$answer_only == " Adult (18-64)", 
                                                                              as.numeric(64), 
                                                                              ifelse(RB_phenotype_survey_b$answer_only == " Older adult (65-74)", 
                                                                                     as.numeric(74), 
                                                                                     ifelse(RB_phenotype_survey_b$answer_only == " Elderly (75+)", 
                                                                                            as.numeric(84), NA_real_)))))}

for (i in seq_len(nrow(RB_phenotype_survey_b))) {
  RB_phenotype_survey_b$self_rep_age_mid <- ifelse(RB_phenotype_survey_b$answer_only == " Child (0-11)", 
                                                                as.numeric(5.5), 
                                                                ifelse(RB_phenotype_survey_b$answer_only == " Adolescent (12-17)", 
                                                                       as.numeric(14.5), 
                                                                       ifelse(RB_phenotype_survey_b$answer_only == " Adult (18-64)", 
                                                                              as.numeric(41), 
                                                                              ifelse(RB_phenotype_survey_b$answer_only == " Older adult (65-74)", 
                                                                                     as.numeric(69.5), 
                                                                                     ifelse(RB_phenotype_survey_b$answer_only == " Elderly (75+)", 
                                                                                            as.numeric(79.5), NA_real_)))))}

for (i in seq_len(nrow(RB_phenotype_survey_b))) {
  RB_phenotype_survey_b$self_rep_age_min <- ifelse(RB_phenotype_survey_b$answer_only == " Child (0-11)", 
                                                                as.numeric(0), 
                                                                ifelse(RB_phenotype_survey_b$answer_only == " Adolescent (12-17)", 
                                                                       as.numeric(12), 
                                                                       ifelse(RB_phenotype_survey_b$answer_only == " Adult (18-64)", 
                                                                              as.numeric(18), 
                                                                              ifelse(RB_phenotype_survey_b$answer_only == " Older adult (65-74)", 
                                                                                     as.numeric(65), 
                                                                                     ifelse(RB_phenotype_survey_b$answer_only == " Elderly (75+)", 
                                                                                            as.numeric(75), NA_real_)))))}

# ------------------------------------------------------------
# SECTION 6 — Merge survey back + person-level survey indicator
# ------------------------------------------------------------


RB_phenotype_survey_c <- full_join(RB_phenotype_survey, RB_phenotype_survey_b, by = 'person_id')

RB_phenotype_cohort_survey_unique <- RB_phenotype_survey_c %>%
  distinct(person_id, .keep_all = TRUE) %>%
  select(person_id, self_rep_age_max, self_rep_age_mid, self_rep_age_min) %>%
  mutate(cohort_survey_Rb = 1)

# ------------------------------------------------------------
# SECTION 7 — Combine phenotype evidence sources
# ------------------------------------------------------------


RB_phenotype_combined <- full_join(RB_phenotype, RB_phenotype_cohort_condition_unique)
RB_phenotype_combined <- full_join(RB_phenotype_combined, RB_phenotype_cohort_survey_unique)

# Replace NA flags with 0 (absence)
RB_phenotype_combined$Rb_phenotype[is.na(RB_phenotype_combined$Rb_phenotype)]=0
RB_phenotype_combined$cohort_condition_Rb[is.na(RB_phenotype_combined$cohort_condition_Rb)]=0
RB_phenotype_combined$cohort_survey_Rb[is.na(RB_phenotype_combined$cohort_survey_Rb)]=0

# ------------------------------------------------------------
# SECTION 8 — Derive diagnosis age variables
# ------------------------------------------------------------

# Compute min across condition-based age_diag and survey-based self-reported age
# min() chooses the earliest age evidence among sources

RB_phenotype_combined$age_diag_max <- 
  apply(RB_phenotype_combined[c('age_diag', 'self_rep_age_max')], 1, function(row) {min(row, na.rm = TRUE)})

RB_phenotype_combined$age_diag_mid <- 
  apply(RB_phenotype_combined[c('age_diag', 'self_rep_age_mid')], 1, function(row) {min(row, na.rm = TRUE)})

RB_phenotype_combined$age_diag_min <- 
  apply(RB_phenotype_combined[c('age_diag', 'self_rep_age_min')], 1, function(row) {min(row, na.rm = TRUE)})

RB_phenotype_combined <- RB_phenotype_combined %>%
  relocate('age_diag_max', .after = 'Rb_phenotype') %>%
  relocate('age_diag_mid', .after = 'age_diag_max') %>%
  relocate('age_diag_min', .after = 'age_diag_mid')

# ------------------------------------------------------------
# SECTION 9 — Create combined cohort: carriers + non-carriers
# ------------------------------------------------------------


#Create dataframe of variant carriers and non_carriers (1 or 0) carriers_after_qc_permissive
combined_cohort_wEHR <- full_join(non_carriers_wEHR, all_variant_carriers_wEHR)
combined_cohort_wEHR <- full_join(combined_cohort_wEHR, carriers_after_qc_strict)

combined_cohort_wEHR$all_carrier[is.na(combined_cohort_wEHR$all_carrier)]=0
combined_cohort_wEHR$qcd_carrier_strict[is.na(combined_cohort_wEHR$qcd_carrier_strict)]=0

# ------------------------------------------------------------
# SECTION 10 — Add phenotype flags to combined cohort
# ------------------------------------------------------------


#Add phenotype columns to dataframe presence or absence (0/1)
combined_cohort_wEHR <- full_join(combined_cohort_wEHR, RB_phenotype_combined)

combined_cohort_wEHR$Rb_phenotype[is.na(combined_cohort_wEHR$Rb_phenotype)]=0
combined_cohort_wEHR$cohort_condition_Rb[is.na(combined_cohort_wEHR$cohort_condition_Rb)]=0
combined_cohort_wEHR$cohort_survey_Rb[is.na(combined_cohort_wEHR$cohort_survey_Rb)]=0

# Drop intermediate columns
combined_cohort_wEHR_final <- combined_cohort_wEHR %>%
  select(-all_of(c('age_diag', 'cohort_condition_Rb', 'self_rep_age_max', 'self_rep_age_mid', 'self_rep_age_min', 'cohort_survey_Rb')))


#write.table(combined_cohort_wEHR_final, file = 'RB1_analysis_total_cohort.tsv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


