# setwd, source libraries, functions
source("code/00_libraries.R")
source("code/00_functions.R")

################################################################################

# (1) SCATTERPLOTS -------------------------------------------------------------
# import data
data <- read_csv("data/variant/processed/bihler_2023_clean.csv")

# PLOT: FOLDING VS CONDUCTANCE SCATTER PLOT
# Folding ("band_c_pct_wt_DMSO") vs. conductance ("fsk_pct_wt_DMSO") among
# missense and indels (in-frame) disease-causing variants
# **constrained x- and y-lim for better interpretability, missing values beyond axes **
p <- data %>%
  # filter(variant %in% c("WT-ctrl", "G551D-ctrl", "G551D", "F508del","F508del-ctrl", "I507del", "R117H")) %>%
  filter(variant %in% c("WT-ctrl", "G551D-ctrl", "G551D", "F508del","F508del-ctrl", "I507del", "R117H", "Y517C","Y563D", "Y563H", "Y563N")) %>%
  # filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>%
  filter(!(is.na(trikafta_response))) %>%
  ggplot(aes(band_c_pct_wt_DMSO, fsk_pct_wt_DMSO, fill = trikafta_response, text  = variant)) +
  geom_point(shape = 21, colour = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values=c("red4", "grey70")) +
  # xlim(0,50) +
  # ylim(0,20) +
  geom_vline(xintercept = c(10)) +
  geom_hline(yintercept = c(10)) +
  theme_minimal()
p
ggplotly(p)
# ggsave(filename = "folding_vs_conductance.pdf",
#        plot = p,
#        path = "output/01_variant_analysis/plots/",
#        width = 5, height = 3)
# written 20230807
rm(p)

# PLOT: FUNCTIONAL BIAS VS. MODULATOR BIAS SCATTER PLOT
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


# (2) DONUT/PIE CHART PLOTS ----------------------------------------------------

# PLOT: FUNCTIONAL BIAS AMONG RESPONSIVE VS. REFRACTORY DONUT PLOTS

data <- read_csv("data/variant/processed/bihler_2023_clean.csv") %>%
  filter( (band_c_pct_wt_DMSO < 10) | (fsk_pct_wt_DMSO < 10)) %>%
  filter(!(is.na(trikafta_response))) %>%
  group_by(trikafta_response) %>%
  select(variant, deficit_class, trikafta_response) %>%
  unique() %>%
  count(deficit_class) %>%
  mutate(total = sum(n)) %>%
  mutate(pct = n/total) %>%
  ungroup()

hsize <- 1
p <- data %>% 
  ggplot(aes(x = 1, y = pct, fill = deficit_class)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  facet_grid(trikafta_response ~ .) +
  theme_minimal()
p
# ggsave(filename = "deficit_class_donut.pdf",
#        plot = p,
#        path = "output/01_variant_analysis/plots/",
#        width = 5, height = 4)
# written 20230807
rm(data, hsize,p)



# (3) BAR & LOLLIPOP PLOTS -----------------------------------------------------
# import data
data <- read_csv("data/variant/processed/bihler_2023_clean.csv")

# PLOT: FOLDING AND CONDUCTANCE FOR EACH RESNO BARPLOT
# (devise ordering variables)
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

# PLOT: FOLDING AND CONDUCTANCE FOR EACH RESNO BARPLOT
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

# PLOT: FOLDING AND CONDUCTANCE FOR EACH RESNO
# now try lollipop plot
# order by mean fold/conductance
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

# PLOT: FOLDING AND CONDUCTANCE FOR EACH RESNO BARPLOT
# plot composite numberby resno
# p <- data %>%
#   ggplot(aes(fill = trikafta_response, text  = variant)) +
#   geom_bar(aes(resno, cond_vs_fold), stat = "identity") +
#   scale_fill_manual(values=c("red4", "grey70")) +
#   theme_minimal()
# p
# ggplotly(p)
# rm(p)

# PLOT: POTENTIATOR AND CORRECTOR FOR EACH RESNO BARPLOT
# (devise ordering variables)
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

# PLOT: POTENTIATOR AND CORRECTOR FOR EACH RESNO BARPLOT
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

# PLOT: POTENTIATOR AND CORRECTOR FOR EACH RESNO BARPLOT
# plot composite number by resno
# p <- data %>%
#   ggplot(aes(fill = trikafta_response, text  = variant)) +
#   geom_bar(aes(resno, iva_vs_elxtez), stat = "identity") +
#   scale_fill_manual(values=c("red4", "grey70")) +
#   theme_minimal()
# p
# ggplotly(p)
# rm(p)

# PLOT: POTENTIATOR AND CORRECTOR FOR EACH RESNO LOLLIPOP
# now try lolipop plot (elx/tez relative to baseline)
p <- data %>%
  # filter(variant %in% c("WT-ctrl", "G551D-ctrl", "G551D", "F508del","F508del-ctrl", "I507del", "R117H", "Y517C","Y563D", "Y563H", "Y563N")) %>%
  
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