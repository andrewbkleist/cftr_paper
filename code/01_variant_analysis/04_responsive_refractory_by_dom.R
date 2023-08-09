# setwd, source libraries, functions
source("code/00_libraries.R")
source("code/00_functions.R")

################################################################################

# import data
data <- read_csv("data/variant/processed/bihler_2023_clean.csv")

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

# order for plot
order <- c("lasso", "TMD1", "NBD1", "R", "TMD2", "NBD2")
temp$dom_mod <- factor(temp$dom_mod, levels = rev(order))
temp$trikafta_response <- factor(temp$trikafta_response, levels = rev(c("yes", "no")))

# plot
p <- temp %>%
  ggplot(aes(dom_mod, pct, fill = trikafta_response)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=c("red4", "grey50")) +
  coord_flip() +
  theme_minimal()

# ggsave(filename = "responsiveness_by_domain.pdf", 
#        plot = p,
#        path = "output/01_variant_analysis/plots/",
#        width = 4, height = 3)
# written 20230807

rm(temp,p)
