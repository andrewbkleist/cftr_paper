# setwd, source libraries, functions
source("code/00_libraries.R")
source("code/00_functions.R")

################################################################################

# import data
data <- read_csv("data/variant/processed/bihler_2023_clean.csv")

# isolate res/ref
res <- data %>%
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>% #pathogenic only
  filter(trikafta_response == "yes") %>% # RESPONSIVE
  select(resno, hgvs_protein, resid, resmut, CAN, dom, sse, dom_sse) %>%
  select(resno) %>% unique() # 232 variants, 185 unique residues 
res <- paste(c(res$resno), collapse='+')
paste0("select responsive, resi ",res)
# writeLines(paste0("select responsive, resi ",res), con ="output/01_variant_analysis/for_pymol/responsive.txt", useBytes = FALSE)
rm(res)
# written 20230807

ref <- data %>%
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>% #pathogenic only
  filter(trikafta_response == "no") %>% # RESPONSIVE
  select(resno, hgvs_protein, resid, resmut, CAN, dom, sse, dom_sse) %>%
  select(resno) %>% unique() # 83 variants, 65 unique residues 
ref <- paste(c(ref$resno), collapse='+')
paste0("select refractory, resi ",ref)
# writeLines(paste0("select refractory, resi ",ref), con ="output/01_variant_analysis/for_pymol/refractory.txt", useBytes = FALSE)
rm(ref)
# written 20230807

# pymol session 20230806_responsive_vs_refractory.pse made from these atom
# selections
# updated 20230806