
library(tidyverse)
library(devtools) 

# Loads a helper function to search through UK Biobank healthcare records
source_url("https://raw.githubusercontent.com/hdg204/UKBB/main/UKBB_Health_Records_New_Project.R")

## ---- Import file with final variant carrier IDs ----
# File contains curated RB1 variant carriers; Assign 1 == carriers.

## ---- Import file with final variant carrier IDs ----
# Reads a curated list of RB1 variant carriers and creates a binary indicator.
# Assumes input file has a column named 'eid' (corresponds to participant ID).
stringent_carriers <- read_tsv("final_RB1_variant_and_carrier_list_27_VAF.tab") %>% 
  select(eid) %>%
  mutate(stringent_carrier = 1) 

## Stringent phenotype criteria > Those with ICD9/10 retinal cancer only OR those with self reported retinoblastoma
#Checked ICD9 records: No participants have ICD9 code for retinal cancer specifically and therefore did not included in the script for stringent criteria

## Stringent phenotype criteria:
## ICD10 retinal cancer only (C69.2) OR self-reported retinoblastoma.
# (You note ICD9 retinal cancer is absent, so excluded.)

## ---- Self-reported data ----
# Loads helper functions for reading self-reported cancer instances.
# NOTE: source() depends on a local file being present in your working directory.
source("selfrep_records_all_instances.R", echo = FALSE)

# Extract self-reported retinoblastoma from baseline + 2 additional instances.
selfrep_records_stringent=read_selfreport_cancer('1075') 
selfrep_records_stringent_b = read_selfreport_cancer_p20001_i1('1075')
selfrep_records_stringent_c = read_selfreport_cancer_p20001_i2('1075')


# Combine sources and remove duplicates
selfrep_records_stringent = selfrep_records_stringent %>%
  bind_rows(selfrep_records_stringent_b) %>%
  bind_rows(selfrep_records_stringent_c) %>%
  unique()



## ---- Hospital episode statistics ----
# Identify HES records with ICD-10 C69.2.

C69_HES =read_ICD10('C69.2')

# Keep unique occurrences; then count number of records per participant.
C69_HES_unique = C69_HES %>% 
  select(eid,diag_icd10,epistart,dob) %>% 
  unique() %>%
  group_by(eid) %>%
  tally()



# ---- Cancer registry data ----
# Identify cancer registry records with ICD-10 C69.2.
C69_cancer=read_cancer('C69.2')

# Count number of registry records per participant.
C69_cancer_unique = C69_cancer %>% 
  select(eid,site,reg_date,dob) %>% 
  unique() %>%
  group_by(eid) %>%
  tally()


# ---- First occurence ----
# Determine earliest diagnosis date and source across available sources.
stringent_phenotype = first_occurence(ICD10 = c('C69.2'), GP = '', OPCS = '', cancer= 'C69.2')

# Get unique participant counts (only count one participant once against all records)
stringent_unique = stringent_phenotype %>% 
  select(eid,source) %>% 
  unique() %>%
  group_by(eid) %>%
  tally()


# Creates a dataset of participant ID with phenotype indicator (1 = case)
stringent = stringent_unique %>%
  full_join(selfrep_records_stringent, by = 'eid') %>%
  mutate(stringent_final=1) %>%
  select(eid, stringent_final)



## ---- Permissive phenotype criteria ----
# Includes broader ICD10/ICD9 cancer codes + procedures + self-report.

# Define all codes in the permissive phenotype criteria
ICD10_codes_2 <- c('H33', 'C69', 'C69.2', 'C75.3', 'Q35',
                   'C34', 
                   'C40', 'C41', 'C43', 'C45', 'C46', 'C47', 'C48', 'C49', 
                   'C50', 'C67', 
                   'C71', 'C72', 
                   'C91', 'C92', 'C93', 'C94', 'C95')
ICD9_codes_2 <- c('749', '190', '1905', '1944', '170', '172', '171', '191', '192', '174', '175', '162', '188', '204', '205', '206', '207', '208')
OPCS4_codes_2 <- c('C01', 'C91', 'F29')
selfrep_codes=c('1075', '1030')



# ---- Extract ICD10 codes from HES/cancer registry ----

# Extract ICD10 codes
ICD10_HES_2=read_ICD10(ICD10_codes_2)
ICD10_unique_2 = ICD10_HES_2 %>% 
  select(eid,diag_icd10) %>% 
  unique() %>%
  group_by(eid) %>%
  tally()

ICD10_cancer_2 = read_cancer(ICD10_codes_2)
cancer_unique_2 = ICD10_cancer_2 %>% 
  select(eid,site) %>% 
  unique() %>%
  group_by(eid) %>%
  tally()


# ---- Extract ICD9 codes ----
# Using read_ICD10() to pull ICD9 codes, then filtering by diag_icd9
ICD9_HES_2 = read_ICD10(c('749', '190', '1905', '1944', '172', '191', '192', '174', '175', '162', '188', '204', '205', '206', '207', '208'))

# Regex filters to keep the exact ICD9 families.
ICD9_HES_2 = ICD9_HES_2 %>% dplyr::filter(grepl('^190|^1905|^1944|^172|^191|^192|^174|^175|^162|^188|^204|^205|^206|^207|^208', diag_icd9))

# This produces a unique list of eids (case flag use-case).
ICD9_HES_2_unique = ICD9_HES_2 %>%
  distinct(eid, .keep_all = TRUE) %>%
  select(eid)

# Cancer registry ICD9 extraction uses read_cancer_2().
ICD9_cancer_2 = read_cancer_2(c('749', '190', '1905', '1944', '170', '172', '171', '191', '192', '174', '175', '162', '188', '204', '205', '206', '207', '208'))
ICD9_cancer_2_unique = ICD9_cancer_2 %>%
  distinct(eid, .keep_all = TRUE) %>%
  select(eid)


## ---- Self-reported data (permissive) ----
# Extract RB or eye/adnexal cancer self-report codes across instances.

selfrep_records=read_selfreport_cancer(selfrep_codes) 
selfrep_records_b = read_selfreport_cancer_p20001_i1(selfrep_codes)
selfrep_records_c = read_selfreport_cancer_p20001_i2(selfrep_codes)

selfrep_records = selfrep_records %>%
  bind_rows(selfrep_records_b) %>%
  bind_rows(selfrep_records_c) %>%
  unique() %>%
  mutate(self_report=1) %>%
  distinct(eid, .keep_all = TRUE)


## ---- OPCS procedures ----
# Extract relevant surgical/procedure codes.
OPCS_proceedure_2 <- read_OPCS(OPCS4_codes_2)
OPCS_proceedure_unique_2 = OPCS_proceedure_2 %>% 
  select(eid,oper4) %>% 
  unique() %>%
  group_by(eid) %>%
  tally()

# Determine earliest diagnosis using first_occurence()
permissive_phenotype = first_occurence(ICD10 = ICD10_codes_2, GP = '', OPCS = OPCS4_codes_2, cancer= ICD10_codes_2)


# Unique participants
permissive_unique = permissive_phenotype %>% 
  select(eid,source) %>% 
  unique() %>%
  group_by(eid) %>%
  tally()

# Creates a dataset of participant ID with phenotype indicator (1 = case)
permissive <- permissive_unique %>%
  full_join(selfrep_records, by = 'eid') %>%
  full_join(ICD9_HES_2_unique, by = 'eid') %>%
  full_join(ICD9_cancer_2_unique, by = 'eid') %>%
  mutate(permissive_final=1) %>%
  select(eid, permissive_final)

# Start with baseline participant table
all_data <- baseline_table

# Remove unused phenotype/biomarker columns
all_data <- all_data %>% 
  select(all_of(c('eid', 'recruit_age', 'yob', 'dob', 'assess_date', 'sex', 'ethnicity')))

# Add RB1 carrier columns to dataframe (0/1) for variant carriers
all_data = full_join(all_data, stringent_carriers, by = 'eid')
all_data$stringent_carrier[is.na(all_data$stringent_carrier)]=0 

# Add phenotype columns to dataframe (0/1) for stringent and permissive categories
all_data = full_join(all_data, stringent, by = 'eid') #stringent
all_data$stringent_final[is.na(all_data$stringent_final)]=0

all_data = full_join(all_data, permissive, by = 'eid') #permissive
all_data$permissive_final[is.na(all_data$permissive_final)]=0


# Optional file export
#write.table(all_data, file = 'RB1_carrier_phenotype_stringent_permissive.tsv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

