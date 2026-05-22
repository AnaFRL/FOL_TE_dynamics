## Author: Ana Rodríguez López
## Date: 07/10/24

# Levene tests based on median to test homocedasticity in conidia-passaged and mycelium-passaged populations

library(car)

# Read the files
YPDA_growth_rate <- read.table("YPDA_growth_rate.txt", sep = "\t", header = T)
YPDA_conidiation_rate <- read.table("YPDA_conidiation_rate.txt", sep = "\t", header = T)
PMM_growth_rate <- read.table("PMM_growth_rate.txt", sep = "\t", header = T)
PMM_conidiation_rate <- read.table("PMM_conidiation_rate.txt", sep = "\t", header = T)

# Change the structure of the data: 
  # create one column with the numerical values, other column with the population name
  data_long_YPDA_growth <- stack(YPDA_growth_rate)
  data_long_PMM_growth <- stack(PMM_growth_rate)
  data_long_YPDA_conidiation <- stack(YPDA_conidiation_rate)
  data_long_PMM_conidiation <- stack(PMM_conidiation_rate)
  # remove empty rows
  data_long_YPDA_growth <- data_long_YPDA_growth[-c(82:102),]
  data_long_PMM_growth <- data_long_PMM_growth[-c(82:102),]
  data_long_YPDA_conidiation <- data_long_YPDA_conidiation[-c(82:102),]
  data_long_PMM_conidiation <- data_long_PMM_conidiation[-c(82:102),]
  # change decimal punctuation and set values column as numeric
  data_long_YPDA_growth$values <- as.numeric(gsub(",", ".", data_long_YPDA_growth$values))
  data_long_PMM_growth$values <- as.numeric(gsub(",", ".", data_long_PMM_growth$values))
  data_long_YPDA_conidiation$values <- as.numeric(gsub(",", ".", data_long_YPDA_conidiation$values))
  data_long_PMM_conidiation$values <- as.numeric(gsub(",", ".", data_long_PMM_conidiation$values))
  
# Run Levene tests
  
# YPDA growth rate
levene_YPDA_growth <- leveneTest(values ~ ind, data = data_long_YPDA_growth, center = "median")

# PMM growth rate
levene_PMM_growth <- leveneTest(values ~ ind, data = data_long_PMM_growth, center = "median")

# YPDA conidiation rate
levene_YPDA_conidiation <- leveneTest(values ~ ind, data = data_long_YPDA_conidiation, center = "median")

# PMM conidiation rate
levene_PMM_conidiation <- leveneTest(values ~ ind, data = data_long_PMM_conidiation, center = "median")

