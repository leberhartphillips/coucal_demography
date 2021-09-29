# load libraries
source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# load functions
function.sources = list.files(path = "R/functions", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

# load capture histories
data.sources = list.files(path = "data/cooked", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

# # load output
# BC_hazard_rate_boot <- 
#   readRDS("output/bootstraps/hazard/cooked/BC_hazard_ASR_bootstrap_result_one.rds")
# 
# # load output
# WBC_hazard_rate_boot <- 
#   readRDS("output/bootstraps/hazard/cooked/WBC_hazard_ASR_bootstrap_result_one.rds")

# # clean up the output from the bootstrap procedure and save as rds
# BC_hazard_rate_boot_tidy <- 
#   hazard_boot_out_wrangle(species = "BC", niter = 1000, 
#                           output_dir = "output/bootstraps/hazard/cooked/",
#                           rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc_trans_no_imm")
# 
# # clean up the output from the bootstrap procedure and save as rds
# WBC_hazard_rate_boot_tidy <- 
#   hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
#                           output_dir = "output/bootstraps/hazard/cooked/",
#                           rds_file = "_hazard_ASR_bootstrap_result_w_WBC_ad_surv_stoc_no_imm")

# load tidy output
BC_hazard_rate_boot_tidy <-
  readRDS("output/bootstraps/hazard/cooked/BC_haz_sur_ASR_boot_tidy_stoc_trans_no_imm.rds")

# load tidy output
WBC_hazard_rate_boot_tidy <-
  readRDS("output/bootstraps/hazard/cooked/WBC_haz_sur_ASR_boot_tidy_stoc_no_imm.rds")

ASR_boot <- 
  bind_rows(BC_hazard_rate_boot_tidy$ASR_ests_boot,
            WBC_hazard_rate_boot_tidy$ASR_ests_boot) %>%
  mutate(species = factor(species, levels = c("BC", "WBC")))

CI <- 0.95

ASR_boot_summary <- 
  ASR_boot %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(lcl_ASR = stats::quantile(M_Adult, (1 - CI)/2, na.rm = TRUE),
                   ucl_ASR = stats::quantile(M_Adult, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_ASR = mean(M_Adult),
                   med_ASR = median(M_Adult),
                   max_ASR = max(M_Adult),
                   min_ASR = min(M_Adult))

# define the factor levels of the population variable so that the populations
# are in an order that reflects the ASR (male biased to female biased)
# ASR_boot <- 
#   ASR_boot %>% 
#   mutate(species = factor(species, levels = c("BC", "WBC")))

Figure_2b <- 
  ggplot() +
  annotate("rect", xmin=0, xmax=0.5, ymin=0, ymax=600, alpha=0.7,
           fill= pull(ggthemes_data$wsj$palettes$colors6[3,2])) +
  annotate("rect", xmin=0.5, xmax=1, ymin=0, ymax=600, alpha=0.7,
           fill= pull(ggthemes_data$wsj$palettes$colors6[2,2])) +
  annotate("text", x = c(0.025), y = c(300),
           label = c("\u2640"), size = 4, colour = "grey10",
           family="Menlo", vjust = c(0), hjust = c(0.5)) +
  annotate("text", x = c(0.975), y = c(300),
           label = c("\u2642"), size = 4, colour = "grey10",
           family="Menlo", vjust = c(1), hjust = c(0.5)) +
  #   annotate("text", x = c(0.2), y = c(80),
  #          label = c("female"), size = 2,
  #          vjust = c(0), hjust = c(0.5)) +
  # annotate("text", x = c(0.8), y = c(80),
  #          label = c("male"), size = 2,
  #          vjust = c(1), hjust = c(0.5)) +
  geom_histogram(binwidth = 0.01, data = ASR_boot, aes(x = M_Adult), fill = "grey30") +
  geom_errorbarh(data = ASR_boot_summary, aes(y = 575, x = lcl_ASR, xmin = lcl_ASR, xmax = ucl_ASR), 
                 color = "black", size = 0.3, linetype = "solid") +
  coord_flip() +
  facet_grid(. ~ species) +
  theme_bw() +
  theme(text = element_text(family="Franklin Gothic Book"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title.x = element_text(size=7, vjust=-0.1),
        axis.text.x  = element_text(size=6, angle = 45, hjust = 1, colour = "black"),
        axis.title.y = element_text(size=7, hjust=0.5, vjust = -2.75, margin = margin(0, 11, 0, 0)),
        axis.text.y  = element_text(size=6, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.2, colour = "grey40"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks.x = element_line(size = 0.2, colour = "grey40"),
        panel.border = element_blank(),
        plot.margin = unit(c(0.2, 1.35, 0.405, 0.1), "cm"),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_blank(),
        # strip.text = element_text(size=11)) +
        strip.text = element_blank()) +
  ylab("Bootstrap frequency") +
  xlab("Adult sex ratio\n(proportion \u2642, 95% CI)") +
  # xlab("Adult sex ratio\n(proportion male)") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 600), expand = c(0, 0), breaks = seq(from = 0, to = 500, 100))

# check colorblindness issues with colors
cvdPlot(Figure_2b)

ggsave(Figure_2b,
       filename = "products/figures/jpg/ASR_plot.jpeg",
       width = 4.75,
       height = 1.59, units = "in",
       dpi = 600,
       # compression="lzw",
       scale = 1)
