# setwd, source libraries, functions
source("code/00_libraries.R")
source("code/00_functions.R")

################################################################################

# (1) IMPORT, CLEAN, ANNOTATE DATA ---------------------------------------------
# import, isolate missense, single res indels
data <- read_csv("data/variant/raw/bihler_2023_mod.csv") %>% 
  filter(variant_type %in% c("Missense", "In-frame in/del"))

# cleanup treatment labels
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

# import lookup, map to resno
lookup <- read_csv("data/sequence/lookup.csv") %>% select(-resid)
lookup <- lookup %>% unite(col = "dom_sse", c(dom, sse), sep = ".", remove = FALSE)
colnames(lookup)[1] <- c("resno")
data <- left_join(data, lookup)
rm(lookup)

# consider separately annotating multi-residue indel and complex with first
# residue involved and appending these to existing list

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
# write_csv(data, "data/variant/raw/bihler_2023_clean.csv")
# write_csv(data, "docs/bihler_2023_clean.csv")

# (2) PLOTS --------------------------------------------------------------------
# PLOT 1: FOLDING VS CONDUCTANCE SCATTER PLOT
# Folding ("band_c_pct_wt_DMSO") vs. conductance ("fsk_pct_wt_DMSO") among
# missense and indels (in-frame) disease-causing variants
p <- data %>%
  # filter(variant %in% c("WT-ctrl", "G551D-ctrl", "G551D", "F508del","F508del-ctrl", "I507del")) %>%
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>%
  ggplot(aes(band_c_pct_wt_DMSO, fsk_pct_wt_DMSO, fill = trikafta_response, text  = variant)) +
  geom_point(shape = 21, colour = "black", size = 5, stroke = 0.5) +
  scale_fill_manual(values=c("red4", "grey70")) +
  # xlim(0,20) + 
  # ylim(0,20) + 
  geom_vline(xintercept = c(10)) +
  geom_hline(yintercept = c(10)) +
  theme_minimal()
p
ggplotly(p)
rm(p)

# PLOT 2: FUNCTIONAL BIAS VS. MODULATOR BIAS SCATTER PLOT
# Variant deficit bias ("cond_vs_fold") versus variant modulator response 
# ("iva_vs_elxtez") bias 
p <- data %>%
  # filter(variant %in% c("WT-ctrl", "G551D-ctrl", "G551D", "F508del","F508del-ctrl", "I507del")) %>%
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>%
  filter(trikafta_response == "yes") %>%
  ggplot(aes(cond_vs_fold, iva_vs_elxtez, fill = trikafta_response, text  = variant)) +
  geom_point(shape = 21, colour = "black", size = 5, stroke = 0.5) +
  scale_fill_manual(values=c("grey70")) +
  geom_vline(xintercept = c(0)) +
  geom_hline(yintercept = c(0)) +
  theme_minimal()
p
ggplotly(p)
rm(p)

# PLOT 3-4: FUNCTIONAL BIAS AMONG RESPONSIVE VS. REFRACTORY DONUT PLOTS
# TRIKAFTA RESPONSIVE
pie <- data %>% 
  filter(variant_type %in% c("Missense", "In-frame in/del")) %>%
  # filter missense/in-frame indel
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>%
  filter(trikafta_response == "yes") %>%
  select(variant, deficit_class) %>%
  unique() %>%
  count(deficit_class)
temp <- as.numeric(pie %>% select(n) %>% colSums())
pie <- pie %>% 
  mutate(pct = n / temp)
rm(temp)
  
hsize <- 1
pie %>%
  ggplot(aes(x = 1, y = pct, fill = deficit_class)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  # facet_grid( ~ modulator_approved) +
  theme_minimal()
rm(pie, hsize)

# TRIKAFTA REFRACTORY
pie <- data %>% 
  filter(variant_type %in% c("Missense", "In-frame in/del")) %>%
  # filter missense/in-frame indel
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>%
  filter(trikafta_response == "no") %>%
  select(variant, deficit_class) %>%
  unique() %>%
  count(deficit_class)
temp <- as.numeric(pie %>% select(n) %>% colSums())
pie <- pie %>% 
  mutate(pct = n / temp)
rm(temp)

hsize <- 1
pie %>%
  ggplot(aes(x = 1, y = pct, fill = deficit_class)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  # facet_grid( ~ modulator_approved) +
  theme_minimal()
rm(pie, hsize)

# PLOT 5: FOLDING AND CONDUCTANCE FOR EACH RESNO BARPLOT
# devise ordering variables
data <- data %>% mutate(ave_fold_cond =(band_c_pct_wt_DMSO*fsk_pct_wt_DMSO/2))
data <- data %>% mutate(rownum = row_number()) 
data <- data %>% unite(ord, c(ave_fold_cond,rownum), sep = "", remove = FALSE)
data$ord <- as.numeric(data$ord)
data <- data %>% unite(variant2, c(variant,rownum), sep = ".", remove = FALSE)
data$variant2 <- factor(data$variant2, levels = data$variant2[order(data$ord)])

# order by resno
p <- data %>%
  ggplot(aes(fill = trikafta_response, text  = variant)) +
  geom_bar(aes(resno, fsk_pct_wt_DMSO), stat = "identity") +
  geom_bar(aes(resno, (band_c_pct_wt_DMSO)*(-1)), stat = "identity") +
  scale_fill_manual(values=c("red4", "grey70")) +
  theme_minimal()
p
ggplotly(p)
rm(p)

# order by mean fold/conductance
p <- data %>%
  ggplot(aes(fill = trikafta_response, text  = variant)) +
  geom_bar(aes(variant2, fsk_pct_wt_DMSO), stat = "identity") +
  geom_bar(aes(variant2, (band_c_pct_wt_DMSO)*(-1)), stat = "identity") +
  scale_fill_manual(values=c("red4", "grey70")) +
  theme_minimal()
p
ggplotly(p)
rm(p)

# now try lolipop plot
p <- data %>%
  ggplot(aes(fill = trikafta_response, color = trikafta_response, text  = variant)) +
  geom_segment(aes(x=variant2, xend=variant2, y=fsk_pct_wt_DMSO, yend=band_c_pct_wt_DMSO)) +
  geom_point( aes(x=variant2, y=fsk_pct_wt_DMSO), size=1 ) +
  geom_point( aes(x=variant2, y=band_c_pct_wt_DMSO), size=1) +
  
  scale_color_manual(values=c("red4", "grey70")) +
  scale_fill_manual(values=c("red4", "grey70")) +
  
  theme_minimal()
p
ggplotly(p)

# plot composite number
# p <- data %>%
#   ggplot(aes(fill = trikafta_response, text  = variant)) +
#   geom_bar(aes(resno, cond_vs_fold), stat = "identity") +
#   scale_fill_manual(values=c("red4", "grey70")) +
#   theme_minimal()
# p
# ggplotly(p)
# rm(p)

# PLOT 6: POTENTIATOR AND CORRECTOR FOR EACH RESNO BARPLOT
# devise ordering variables
data <- data %>% mutate(ave_iva_elxtez =(fsk_plus_iva_pct_wt_DMSO*fsk_pct_wt_ELXTEZ/2))
data <- data %>% mutate(rownum = row_number()) 
data <- data %>% unite(ord, c(ave_iva_elxtez,rownum), sep = "", remove = FALSE)
data$ord <- as.numeric(data$ord)
data <- data %>% unite(variant2, c(variant,rownum), sep = ".", remove = FALSE)
data$variant2 <- factor(data$variant2, levels = data$variant2[order(data$ord)])

# order by resno
p <- data %>%
  ggplot(aes(fill = trikafta_response, text  = variant)) +
  geom_bar(aes(resno, fsk_plus_iva_pct_wt_DMSO), stat = "identity") +
  geom_bar(aes(resno, (fsk_pct_wt_ELXTEZ)*(-1)), stat = "identity") +
  scale_fill_manual(values=c("red4", "grey70")) +
  theme_minimal()
p
ggplotly(p)
rm(p)

# order by mean response to corrector/potentiators
p <- data %>%
  ggplot(aes(fill = trikafta_response, text  = variant)) +
  geom_bar(aes(variant2, fsk_plus_iva_pct_wt_DMSO), stat = "identity") +
  geom_bar(aes(variant2, (fsk_pct_wt_ELXTEZ)*(-1)), stat = "identity") +
  scale_fill_manual(values=c("red4", "grey70")) +
  theme_minimal()
p
ggplotly(p)
rm(p)

# plot composite number
# p <- data %>%
#   ggplot(aes(fill = trikafta_response, text  = variant)) +
#   geom_bar(aes(resno, iva_vs_elxtez), stat = "identity") +
#   scale_fill_manual(values=c("red4", "grey70")) +
#   theme_minimal()
# p
# ggplotly(p)
# rm(p)

# now try lolipop plot (trying for elx/tez relative to baseline)
p <- data %>%
  ggplot(aes(fill = trikafta_response, color = trikafta_response, text  = variant)) +
  geom_segment(aes(x=variant2, xend=variant2, y=fsk_pct_wt_DMSO, yend=fsk_pct_wt_ELXTEZ)) +
  geom_point( aes(x=variant2, y=fsk_pct_wt_DMSO), size=1, shape=19) +
  geom_point( aes(x=variant2, y=fsk_pct_wt_ELXTEZ), size=1, shape=17) +
  
  scale_color_manual(values=c("red4", "grey70")) +
  scale_fill_manual(values=c("red4", "grey70")) +
  geom_hline(yintercept = c(10), linetype="dashed", alpha = 0.5) +
  theme_minimal()
p
ggplotly(p)

# (3) IDENTIFY RESPONSIVE & NON-RESPONSIVE -------------------------------------
# characterize position by domain
data <- data %>% mutate(dom_mod = case_when(
  dom == "N1" ~ "NBD1",
  dom == "N2" ~ "NBD2",
  dom == "L" ~ "lasso",
  dom == "R" ~ "R",
  dom == "T1" & !(dom_sse %in% c("T1.4", "T1.5", "T1.45")) ~ "TMD1",
  dom == "T2" & !(dom_sse %in% c("T2.4", "T2.5", "T2.45")) ~ "TMD2",
  dom_sse %in% c("T2.4", "T2.5", "T2.45") ~ "TMD1",
  dom_sse %in% c("T1.4", "T1.5", "T1.45") ~ "TMD2"
))

# get percentage by region
temp <- data %>% filter(!(is.na(trikafta_response))) %>% 
  group_by(trikafta_response) %>% 
  count(dom_mod) %>% ungroup() %>%
  group_by(trikafta_response) %>% 
  mutate(total = sum(n)) %>% ungroup() %>%
  mutate(pct = n/total)

order <- c("lasso", "TMD1", "NBD1", "R", "TMD2", "NBD2")
temp$dom_mod <- factor(temp$dom_mod, levels = order)
temp$trikafta_response <- factor(temp$trikafta_response, levels = c("yes", "no"))

temp %>%
  ggplot(aes(dom_mod, pct, fill = trikafta_response)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=c("grey50", "red4")) +
  # coord_flip() +
  theme_minimal()

# isolate res/ref
res <- data %>%
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>% #pathogenic only
  filter(trikafta_response == "yes") %>% # RESPONSIVE
  select(resno, hgvs_protein, resid, resmut, CAN, dom, sse, dom_sse) %>%
  select(resno) %>% unique() # 232 variants, 185 unique residues 
res <- paste(c(res$resno), collapse='+')
paste0("select responsive, resi ",res)

ref <- data %>%
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>% #pathogenic only
  filter(trikafta_response == "no") %>% # RESPONSIVE
  select(resno, hgvs_protein, resid, resmut, CAN, dom, sse, dom_sse) %>%
  select(resno) %>% unique() # 83 variants, 65 unique residues 
ref <- paste(c(ref$resno), collapse='+')
paste0("select refractory, resi ",ref)

# pymol session 20230806_responsive_vs_refractory.pse made from these atom
# selections
# updated 20230806



#########
# install.packages("shiny")
# install.packages("remotes")
library(shiny)
library(remotes)
# remotes::install_github("Appsilon/shiny.molstar")
library(shiny.molstar)
# install.packages("glue")
library(glue)

pdbId <- "6msm"

shinyApp(
  ui = basicPage(
    tags$main(
      tags$div(
        class = "box",
        ##
        Molstar(
          pdbId = pdbId,
          dimensions = c(300, 300),
          showAxes = TRUE
        ),
        ##
        tags$hr(),
        tags$span(
          "Molecular visualization of pdbID:",
          tags$a(
            href = glue::glue("https://www.ebi.ac.uk/pdbe/entry/pdb/{pdbId}"),
            pdbId
          )
        )
      )
    )
  ),
  server = function(input, output) {
  }
)


