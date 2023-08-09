# setwd, source libraries, functions
source("code/00_libraries.R")
source("code/00_functions.R")

################################################################################

# import data 
# **isolate missense, single res indels**
data <- read_csv("data/variant/raw/bihler_2023_mod.csv") %>% 
  filter(variant_type %in% c("Missense", "In-frame in/del"))
  # multiple residue indels and "complex" (ie multiple missense) variant types
  # are ignored here; consider separately annotating these by labeling their "resno" 
  # with first residue involved and appending these to the existing df

# cleanup treatment labels - make compact
data <- data %>% mutate(treatment = case_when(
  treatment == "DMSO" ~ "DMSO",
  treatment == "ELX + TEZ" ~ "ELXTEZ"))

# concatenate variants to single row for different treatments
data <- data %>% 
  pivot_wider(names_from = treatment, 
              values_from = c(fsk_pct_wt:total_protein_pct_wt_sem))

# extract resno, resid, resmut
data <- data %>% 
  extract(col = hgvs_protein, 
          into=c('p','resid','resno','resmut'),
          regex = "(p\\.)([A-Z][a-z][a-z])([0-9]+)([A-Za-z][a-z][a-z])", 
          convert=TRUE, remove = FALSE) %>%
  select(-p)

# import lookup table, map common numbering to resno
lookup <- read_csv("data/sequence/lookup.csv") %>% select(-resid)
lookup <- lookup %>% unite(col = "dom_sse", c(dom, sse), sep = ".", remove = FALSE)
colnames(lookup)[1] <- c("resno")
data <- left_join(data, lookup)
rm(lookup)

# label with deficit type (""folding", "conductance", "dual", none")
data <- data %>%
  mutate(deficit_class = case_when(
    band_c_pct_wt_DMSO >= 10 & fsk_pct_wt_DMSO >= 10 ~ "none",
    band_c_pct_wt_DMSO < 10 & fsk_pct_wt_DMSO < 10 ~ "dual",
    band_c_pct_wt_DMSO < 10  & fsk_pct_wt_DMSO >= 10 ~ "fold",
    band_c_pct_wt_DMSO >= 10 & fsk_pct_wt_DMSO < 10 ~ "conductance"
  ))

# label with responsive/refractory to trikafta ("trikafta_response")
data <- data %>%
  mutate(trikafta_response = case_when(
    fsk_plus_iva_pct_wt_ELXTEZ >= 10 ~ "yes",
    fsk_plus_iva_pct_wt_ELXTEZ < 10 ~ "no",
  ))

# add numeric term capturing primary deficit type ("cond_vs_fold")
data <- data %>%
  mutate(cond_vs_fold = log2(fsk_pct_wt_DMSO/band_c_pct_wt_DMSO))

# add numeric term capturing primary modulator subtype response ("iva_vs_elxtez")
data <- data %>%
  mutate(iva_vs_elxtez = log2(fsk_plus_iva_pct_wt_DMSO/fsk_pct_wt_ELXTEZ) )

# write csv
# write_csv(data, "data/variant/processed/bihler_2023_clean.csv")
# write_csv(data, "docs/bihler_2023_clean.csv")
# last written 20230807