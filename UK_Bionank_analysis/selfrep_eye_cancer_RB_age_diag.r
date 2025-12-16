# Load required packages ---------------------------------------------------

install.packages("tidyverse")
library(tidyverse)

# -----------------------------------------------------------------------------
# PURPOSE
# -----------------------------------------------------------------------------
# This script builds a participant-level dataset from UKB Cohort Browser output
# (TSV) generated using the Table Exporter (DNA nexus) containing self-reported cancer codes and the corresponding "year/age first
# occurred" fields across multiple instances/arrays.
#
# It extracts rows with:
#  - "eye and/or adnexal cancer"
#  - "retinoblastoma"
# and derives an "age at diagnosis" variable.
#
# INPUT FILE
#   "RB_eye_cancer_self_rep.tsv" downloaded from cohort browser with fields:
#     - Cancer code, self-reported | Instance {0} | Array {0, 1, 2}
#     - Cancer code, self-reported | Instance {1, 2} | Array {1, 2}
#     - Cancer year/age first occurred | Instance {0} | Array {0, 1, 2}
#     - Cancer year/age first occurred | Instance {1, 2} | Array {1, 2}
#
#
# OUTPUT (optional)
#   "Selfrep_agediag.tsv"
# -----------------------------------------------------------------------------


# Read self-reported dataset ------------------------------------------------
RB_eye_cancer_selfrep <- read_tsv('RB_eye_cancer_self_rep.tsv')

# Standardize key identifiers
RB_eye_cancer_selfrep <- RB_eye_cancer_selfrep %>%
  rename(eid = `Participant ID`) %>%
  rename(YoB = `Year of birth`)


# -----------------------------------------------------------------------------
# Combine each (Cancer code, self-reported) with its corresponding (year/age first occurred)
# into a single string per instance/array such as:
#   "eye and/or adnexal cancer, 35"
# or
#   "retinoblastoma, 1992"
#
# unite(..., remove = FALSE) keeps original columns (useful for QC/debugging).
# na.rm = TRUE drops NAs from the pasted result (so you don't get "NA, 35").
# -----------------------------------------------------------------------------

RB_eye_cancer_diag <- RB_eye_cancer_selfrep %>%
  unite("Instance0_array0", 
        'Cancer code, self-reported | Instance 0 | Array 0',
        'Cancer year/age first occurred | Instance 0 | Array 0',
        sep = ", ", remove = FALSE, na.rm = TRUE) %>%
  unite("Instance0_array1", 
        'Cancer code, self-reported | Instance 0 | Array 1',
        'Cancer year/age first occurred | Instance 0 | Array 1',
        sep = ", ", remove = FALSE, na.rm = TRUE) %>%
  unite("Instance0_array2", 
        'Cancer code, self-reported | Instance 0 | Array 2',
        'Cancer year/age first occurred | Instance 0 | Array 2',
        sep = ", ", remove = FALSE, na.rm = TRUE) %>%
  unite("Instance1_array0", 
        'Cancer code, self-reported | Instance 1 | Array 0',
        'Cancer year/age first occurred | Instance 1 | Array 0',
        sep = ", ", remove = FALSE, na.rm = TRUE) %>%
  unite("Instance1_array1", 
        'Cancer code, self-reported | Instance 1 | Array 1',
        'Cancer year/age first occurred | Instance 1 | Array 1',
        sep = ", ", remove = FALSE, na.rm = TRUE) %>%
  unite("Instance2_array0", 
        'Cancer code, self-reported | Instance 2 | Array 0',
        'Cancer year/age first occurred | Instance 2 | Array 0',
        sep = ", ", remove = FALSE, na.rm = TRUE) %>%
  unite("Instance2_array1", 
        'Cancer code, self-reported | Instance 2 | Array 1',
        'Cancer year/age first occurred | Instance 2 | Array 1',
        sep = ", ", remove = FALSE, na.rm = TRUE)

# Keep only the unified combined columns, plus ID + Year of Birth
RB_eye_cancer_diag <- RB_eye_cancer_diag %>%
  select(all_of(c('eid', 'Instance0_array0', 'Instance0_array1', 'Instance0_array2',
                  'Instance1_array0', 'Instance1_array1', 
                  'Instance2_array0', 'Instance2_array1', 
                  'YoB')))

# -----------------------------------------------------------------------------
# Row-wise extraction:
# For each participant (row), scan across ALL columns and collect those that match
#  - "eye and/or adnexal cancer"
#  - "retinoblastoma"
#
# pmap_chr(across(everything()), ...) builds one string per row.
# values <- c(...) gives the row's values as a vector.
#
# If matches exist: paste them together (comma-separated).
# If no matches: return NA
# -----------------------------------------------------------------------------

RB_eye_cancer_only_agediag <- RB_eye_cancer_diag %>%
  mutate(`selfrep_eyecancer` = pmap_chr(across(everything()), 
                                        ~ {
                                          values <- c(...)  # Get all row values
                                          matches <- values[str_detect(values, "eye and/or adnexal cancer")]
                                          ifelse(length(matches) > 0, paste(matches, collapse = ", "), NA)
                                        })) %>%
  mutate(`selfrep_RB` = pmap_chr(across(everything()), 
                                 ~ {
                                   values <- c(...)  # Get all row values
                                   matches <- values[str_detect(values, "retinoblastoma")]
                                   ifelse(length(matches) > 0, paste(matches, collapse = ", "), NA)
                                 })) %>%
  select(all_of(c('eid', 'selfrep_eyecancer', 'selfrep_RB', 'YoB')))


# -----------------------------------------------------------------------------
# Split the combined strings into:
#  - cancer label column(s)
#  - age/year column(s)
#
# IMPORTANT: This assumes:
#   1) Each match looks like "label, number"
#   2) At most TWO eye cancer occurrences are captured (fields _1 and _2)
#   3) Only ONE retinoblastoma occurrence is captured
#
# sep="," splits on commas. If there are extra commas (or you pasted multiple
# matches together), separation can misalign.
# -----------------------------------------------------------------------------

RB_eye_cancer_only_agediag <- RB_eye_cancer_only_agediag %>%
  separate(selfrep_eyecancer, into = c("selfrep_eyecancer_1", "selfrep_eyecancer_age_1", 
                                       "selfrep_eyecancer_2", 'selfrep_eyecancer_age_2'), sep = ",") %>%
  separate(selfrep_RB, into = c("selfrep_RB", "selfrep_RB_age"), sep = ",")


# Convert age/year strings to numeric for calculations
RB_eye_cancer_only_agediag$selfrep_eyecancer_age_1 <- as.numeric(RB_eye_cancer_only_agediag$selfrep_eyecancer_age_1)
RB_eye_cancer_only_agediag$selfrep_eyecancer_age_2 <- as.numeric(RB_eye_cancer_only_agediag$selfrep_eyecancer_age_2)
RB_eye_cancer_only_agediag$selfrep_RB_age <- as.numeric(RB_eye_cancer_only_agediag$selfrep_RB_age)


# -----------------------------------------------------------------------------
# Derive age at diagnosis from the "year/age first occurred" field.
#
# UKB fields can contain either:
#  - an AGE (e.g., 25)
#  - a YEAR (e.g., 1998)
#
# Heuristic used:
#   if value <= 1000 -> treat as AGE
#   else -> treat as YEAR and convert to AGE via (YEAR - YoB)
#
# This creates a staged variable:
#   age_diag  : from first eye cancer occurrence
#   age_diag_2: update using second eye cancer occurrence if present
#   age_diag_3: update using RB occurrence if present
#
# NOTE: Each step overwrites only if the later field is available.
# -----------------------------------------------------------------------------

RB_eye_cancer_only_agediag <- RB_eye_cancer_only_agediag %>%
  mutate(age_diag = ifelse(selfrep_eyecancer_age_1 <= 1000,
                           selfrep_eyecancer_age_1,
                           selfrep_eyecancer_age_1 - YoB))

RB_eye_cancer_only_agediag <- RB_eye_cancer_only_agediag %>%
  mutate(age_diag_2 = ifelse(is.na(selfrep_eyecancer_age_2), age_diag, ifelse(selfrep_eyecancer_age_2 <= 1000, 
                                                                              selfrep_eyecancer_age_2,
                                                                              selfrep_eyecancer_age_2 - YoB)))

RB_eye_cancer_only_agediag <- RB_eye_cancer_only_agediag %>%
  mutate(age_diag_3 = ifelse(is.na(selfrep_RB_age), age_diag_2, ifelse(selfrep_RB_age <= 1000, 
                                                                       selfrep_RB_age,
                                                                       selfrep_RB_age - YoB)))

# -----------------------------------------------------------------------------
# OUTPUT (optional)
# This file is used downstream to map age of diagnosis to cohort
# (specifically referenced by "age_diag_ustringent_permissive.R")
# -----------------------------------------------------------------------------
# write.table(RB_eye_cancer_only_agediag,
#             file = 'Selfrep_agediag.tsv',
#             sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

