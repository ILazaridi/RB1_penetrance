
library(tidyverse)

# Load data ---------------------------------------------------------------
combined_cohort_wEHR_final <- read_tsv('RB1_analysis_total_cohort.tsv')

#Statistical analysis; QCd carriers (strict) with "relaxed phenotype" 

# Statistical analysis ----------------------------------------------------

# Create a contingency table comparing:
#   - Genotype: Carrier vs non-carrier (qcd carriers)
#   - Phenotype: Rb phenotype (Rb and/or ocular cancer) present vs absent 


contingency_table <- table(combined_cohort_wEHR_final$qcd_carrier_strict, combined_cohort_wEHR_final$relaxed_phenotype)

dimnames(contingency_table) <- list(
  'qcd_carrier (genotype)' = c('non-carrier', 'carrier'),
  'Rb_phenotype (phenotype)' = c('absent', 'present'))

fisher_result <- fisher.test(contingency_table[2:1, 2:1])
