## data-raw/DATASET.R

library(readr)
library(dplyr)
library(tidyr)

message("Processing Human Protein Atlas Data...")

# 1. Load Data
hpa_raw <- read_tsv("data-raw/proteinatlas.tsv", col_types = cols())

# 2. Filter Categories
target_categories <- c("Tissue enriched", "Group enriched", "Tissue enhanced")

# 3. Process
# We ONLY use the 'Gene' column since you confirmed it contains the Symbols.
hpa_dataset <- hpa_raw %>%
  filter(`RNA tissue specificity` %in% target_categories) %>%
  select(
    Gene = Gene,  # This contains the Symbols
    Specificity = `RNA tissue specificity`,
    Tissue_String = `RNA tissue specific nTPM`
  ) %>%
  # Split the "tissue: score; tissue: score" string
  separate_rows(Tissue_String, sep = ";") %>%
  separate(Tissue_String, into = c("Tissue", "nTPM"), sep = ":") %>%
  mutate(
    Tissue = trimws(Tissue),
    nTPM = as.numeric(nTPM)
  ) %>%
  select(Gene, Specificity, Tissue) %>%
  distinct()

# 4. Save to package
usethis::use_data(hpa_dataset, overwrite = TRUE)