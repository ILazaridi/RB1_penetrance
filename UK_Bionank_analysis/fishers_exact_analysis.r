
library(tidyverse)

# Load data ---------------------------------------------------------------
all_data = read_tsv('RB1_carrier_phenotype_ustringent_permissive_wgs_only.tsv')



# Statistical analysis: stringent phenotype -------------------------------

# Create a contingency table comparing:
#   - Genotype: stringent carrier vs non-carrier
#   - Phenotype: ultra-stringent phenotype present vs absent

contingency_table <- table(all_data$stringent_carrier, all_data$ultra_stringent)

# Add descriptive row and column names to the contingency table
dimnames(contingency_table) <- list(
  'stringent (genotype)' = c('non-carrier', 'carrier'),
  'ultra_stringent (phenotype)' = c('absent', 'present'))


# Perform Fisher's exact test
# Index [2:1, 2:1] reorders rows and columns so that:
#   - carrier and phenotype-present appear first
# This does not change the test result, but standardizes table orientation
fisher_result_1 <- fisher.test(contingency_table[2:1, 2:1])


# Statistical analysis: permissive phenotype ------------------------------


# Create a contingency table comparing:
#   - Genotype: stringent carrier vs non-carrier
#   - Phenotype: permissive phenotype present vs absent
contingency_table_4 <- table(all_data$stringent_carrier, all_data$permissive_final)


# Add descriptive row and column names
dimnames(contingency_table_4) <- list(
  'stringent (genotype)' = c('non-carrier', 'carrier'),
  'permissive (phenotype)' = c('absent', 'present'))


# Perform Fisher's exact test with standardized table orientation
fisher_result_2 <- fisher.test(contingency_table_4[2:1, 2:1])


# Multiple testing correction --------------------------------------------

# Extract p-values from both Fisher's exact tests
p_values <- c(fisher_result_1$p.value, 
              fisher_result_2$p.value)

# Apply Bonferroni correction to control for multiple hypothesis testing
p_corrected <- p.adjust(p_values, method = 'bonferroni')



