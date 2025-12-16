
# ------------------------------------------------------------------------------
# Run phenotype script: creates essential variables required downstream
# ------------------------------------------------------------------------------

# Source code
source("phenotype_script_stringent_permissive.r", echo = FALSE)


# ------------------------------------------------------------------------------
# Load self-reported dataset with age of diagnosis
# Includes participants who self-reported retinoblastoma or eye/adnexal cancer
# and their corresponding age at diagnosis.
#
# NOTE: Requires running "selfrep_eye_cancer_RB_age_diag.r" beforehand.
# ------------------------------------------------------------------------------



# Will require running "selfrep_eye_cancer_RB_age_diag.R" script
RB_eye_cancer_only_agediag <- read_tsv("Selfrep_agediag.tsv")

# Subset to individuals who self reported retinoblastoma specifically
RB_only_agediag <- RB_eye_cancer_only_agediag %>%
  filter(selfrep_RB == 'retinoblastoma')

# ------------------------------------------------------------------------------
# Base dataframe to which age-of-diagnosis variables will be added
# all_data generated from phenotype code: keep core variables + phenotype flags
# ------------------------------------------------------------------------------
all_data_1 <- all_data %>% # Generated in phenotype_script_stringent_permissive.r script
  
  select(all_of(c('eid', 'recruit_age', 'yob', 'dob', 'assess_date', 'sex', 'ethnicity',
                  'stringent_carrier', 'stringent_final', 'permissive_final'))) %>%
  relocate('dob', .after = 'eid') %>%
  relocate('yob', .after = 'dob')


# ==============================================================================
# Stringent phenotype: derive age at diagnosis
# Strategy: merge multiple data sources, extract per-source age-at-diagnosis,
#           then take MIN age across sources per participant.
# ==============================================================================

# ------------------------------------------------------------------------------
# Merge HES records and extract age at diagnosis (Stringent case definition)
# ------------------------------------------------------------------------------

all_data_2 <- full_join(all_data_1, C69_HES_unique, by = 'eid') %>%
  mutate(
    diag_age_HES = ifelse(
      eid %in% C69_HES$eid,
      # Create diagnosis age from the *full* C69_HES table using match on eid.
      #
      # ANNOTATION:
      # - Using `ifelse(..., value, '-')` forces the column to character.
      # - The '-' placeholder later becomes NA when as.numeric() is called.
      C69_HES$diag_age[match(eid, C69_HES$eid)], '-')) %>% 
  select(-all_of('n'))

# ------------------------------------------------------------------------------
# Merge cancer registry records and extract age at diagnosis
# ------------------------------------------------------------------------------

all_data_2 <- full_join(all_data_2, C69_cancer_unique, by = 'eid') %>%
  mutate(
    diag_age_cancer = ifelse(
      eid %in% C69_cancer$eid,
      # Prefer diag_age; if missing, fall back to age
      
      
      coalesce(C69_cancer$diag_age[match(eid, C69_cancer$eid)],
               C69_cancer$age[match(eid, C69_cancer$eid)]),'-')) %>% 
  select(-all_of('n')) 

# ------------------------------------------------------------------------------
# Add self-reported Stringent RB cases and age of diagnosis from self-report
# ------------------------------------------------------------------------------

all_data_2 <- full_join(all_data_2, selfrep_records_stringent, by = 'eid') %>%
  mutate(
    # Only populate self-reported age for those flagged stringent_final == 1
    
    diag_age_selfrep_RB = ifelse(
      stringent_final == 1 & eid %in% RB_eye_cancer_only_agediag$eid,
      RB_eye_cancer_only_agediag$age_diag_3[match(eid, RB_eye_cancer_only_agediag$eid)], NA_real_)) # NA if there is no self-reported age or no ID match in "RB_eye_cancer_only_agediag"

# ------------------------------------------------------------------------------
# Convert age variables from character to numeric
# as.numeric("-") -> NA, which is how the placeholder gets dropped
# ------------------------------------------------------------------------------

all_data_2$diag_age_HES <- as.numeric(all_data_2$diag_age_HES)
all_data_2$diag_age_cancer <- as.numeric(all_data_2$diag_age_cancer)
all_data_2$diag_age_selfrep_RB <- as.numeric(all_data_2$diag_age_selfrep_RB)

# ------------------------------------------------------------------------------
# For each participant: minimum age across sources = Stringent age at dx
# ------------------------------------------------------------------------------

all_data_2$stringent_agediag <- 
  apply(all_data_2[c('diag_age_HES', 'diag_age_cancer',
                     'diag_age_selfrep_RB')], 1, function(row) {
                       min(row, na.rm = TRUE)
                     })

# If all sources were NA, min(..., na.rm=TRUE) returns Inf -> convert back to NA
all_data_2$stringent_agediag[is.infinite(all_data_2$stringent_agediag)] <- NA

# Drop intermediate age columns; keep only final combined age
all_data_2 <- all_data_2 %>%
  select(-all_of(c('diag_age_HES', 'diag_age_cancer', 'diag_age_selfrep_RB')))



# ==============================================================================
# Permissive phenotype: derive age at diagnosis
# Strategy: same as above, but with broader sources (ICD10/ICD9/OPCS/self-report)
# ==============================================================================


# ------------------------------------------------------------------------------
# Merge HES ICD10 and extract age at diagnosis
# ------------------------------------------------------------------------------

all_data_4 <- full_join(all_data_1, ICD10_unique_2, by = 'eid') %>%
  mutate(
    diag_age_HES_2 = ifelse(
      eid %in% ICD10_HES_2$eid,
      ICD10_HES_2$diag_age[match(eid, ICD10_HES_2$eid)], '-')) %>%
  select(-all_of('n'))


# ------------------------------------------------------------------------------
# Merge cancer registry ICD10 and extract age at diagnosis
# ------------------------------------------------------------------------------

all_data_4 <- full_join(all_data_4, cancer_unique_2, by = 'eid') %>%
  mutate(
    diag_age_cancer_2 = ifelse(
      eid %in% ICD10_cancer_2$eid,
      coalesce(ICD10_cancer_2$diag_age[match(eid, ICD10_cancer_2$eid)], 
               ICD10_cancer_2$age[match(eid, ICD10_cancer_2$eid)]),'-')) %>%
  select(-all_of('n'))


# ------------------------------------------------------------------------------
# Merge HES ICD9 and extract age at diagnosis
# ------------------------------------------------------------------------------

all_data_4 <- full_join(all_data_4, ICD9_HES_2_unique, by = 'eid') %>%
  mutate(
    diag_age_HES_ICD9_2 = ifelse(
      eid %in% ICD9_HES_2$eid,
      ICD9_HES_2$diag_age[match(eid, ICD9_HES_2$eid)], '-'))


# ------------------------------------------------------------------------------
# Merge cancer registry ICD9 and extract age at diagnosis
# ------------------------------------------------------------------------------

all_data_4 <- full_join(all_data_4, ICD9_cancer_2_unique, by = 'eid') %>%
  mutate(
    diag_age_cancer_ICD9_2 = ifelse(
      eid %in% ICD9_cancer_2$eid,
      coalesce(ICD9_cancer_2$diag_age[match(eid, ICD9_cancer_2$eid)], 
               ICD9_cancer_2$age[match(eid, ICD9_cancer_2$eid)]),'-'))


# ------------------------------------------------------------------------------
# Merge OPCS procedures (Rb-associated surgeries) and extract age of procedure
# ------------------------------------------------------------------------------

all_data_4 <- full_join(all_data_4, OPCS_proceedure_unique_2, by = 'eid') %>%
  mutate(
    diag_age_OPCS_2 = ifelse(
      eid %in% OPCS_proceedure_2$eid,
      OPCS_proceedure_2$op_age[match(eid, OPCS_proceedure_2$eid)], '-')) %>%
  select(-all_of('n'))


# ------------------------------------------------------------------------------
# Add self-reported Rb / eye / adnexal cancer cases and age of diagnosis
# ------------------------------------------------------------------------------

all_data_4 <- full_join(all_data_4, selfrep_records, by = 'eid') %>%
  mutate(
    diag_age_selfrep_RB_1 = ifelse(
      eid %in% RB_eye_cancer_only_agediag$eid,
      RB_eye_cancer_only_agediag$age_diag_3[match(eid, RB_eye_cancer_only_agediag$eid)], '-')) %>%
  select(-all_of('self_report'))

# ------------------------------------------------------------------------------
# Convert age variables from character to numeric
# ------------------------------------------------------------------------------

all_data_4$diag_age_HES_2 <- as.numeric(all_data_4$diag_age_HES_2)
all_data_4$diag_age_cancer_2 <- as.numeric(all_data_4$diag_age_cancer_2)
all_data_4$diag_age_HES_ICD9_2 <- as.numeric(all_data_4$diag_age_HES_ICD9_2)
all_data_4$diag_age_cancer_ICD9_2 <- as.numeric(all_data_4$diag_age_cancer_ICD9_2)
all_data_4$diag_age_OPCS_2 <- as.numeric(all_data_4$diag_age_OPCS_2)
all_data_4$diag_age_selfrep_RB_1 <- as.numeric(all_data_4$diag_age_selfrep_RB_1)


# ------------------------------------------------------------------------------
# For each participant: minimum age across sources = permissive age at dx
# ------------------------------------------------------------------------------

all_data_4$permissive_agediag <- 
  apply(all_data_4[c('diag_age_HES_2', 'diag_age_cancer_2', 
                     'diag_age_HES_ICD9_2', 'diag_age_cancer_ICD9_2',
                     'diag_age_OPCS_2',
                     'diag_age_selfrep_RB_1')], 1, function(row) {
                       min(row, na.rm = TRUE)
                     })

# Replace Inf (no non-NA age in any source) with NA
all_data_4$permissive_agediag[is.infinite(all_data_4$permissive_agediag)] <- NA

# Drop intermediate columns; keep only permissive_agediag
all_data_4 <- all_data_4 %>%
  select(-all_of(c('diag_age_HES_2', 'diag_age_cancer_2', 
                   'diag_age_HES_ICD9_2', 'diag_age_cancer_ICD9_2',
                   'diag_age_OPCS_2',
                   'diag_age_selfrep_RB_1')))


# ------------------------------------------------------------------------------
# Combine Stringent and permissive age-of-diagnosis
# ------------------------------------------------------------------------------

all_data_with_diag_age <- all_data_2 %>%
  full_join(all_data_4, by = c('eid', 'recruit_age', 'yob', 'dob', 'assess_date', 'sex', 'ethnicity',
                               'stringent_carrier', 'stringent_final', 'permissive_final'))

# Optional file export
# write.table(all_data_with_diag_age,
#             file = 'all_data_with_diag_age_stringent_permissive.tsv',
#             sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

