# load packages
source("scripts/01_libraries.R")

# load wrangled bootstrap output
load("output/wrangled/Black_Coucal_bootstrap_result_clean.rds")
load("output/wrangled/White-browed_Coucal_bootstrap_result_clean.rds")

ASR_boot <- 
  bind_rows(BC_boot_out$ASR_boot,
            WBC_boot_out$ASR_boot) %>%
  mutate(species = factor(species, levels = c("BC", "WBC")))

CI <- 0.95

ASR_boot_summary <- 
  ASR_boot %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(ucl_ASR = stats::quantile(ASR_estimate, (1 - CI)/2, na.rm = TRUE),
                   lcl_ASR = stats::quantile(ASR_estimate, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_ASR = mean(ASR_estimate),
                   med_ASR = median(ASR_estimate),
                   max_ASR = max(ASR_estimate),
                   min_ASR = min(ASR_estimate))

# define the factor levels of the population variable so that the populations
# are in an order that reflects the ASR (male biased to female biased)
# ASR_boot <- 
#   ASR_boot %>% 
#   mutate(species = factor(species, levels = c("BC", "WBC")))

Figure_2b <- 
  ggplot() +
  annotate("rect", xmin=0, xmax=0.5, ymin=0, ymax=160, alpha=0.6,
           fill= brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=0.5, xmax=1, ymin=0, ymax=160, alpha=0.6,
           fill= brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(0.025), y = c(80),
           label = c("\u2640"), size = 4, colour = "grey10",
           family="Menlo", vjust = c(0), hjust = c(0.5)) +
  annotate("text", x = c(0.975), y = c(80),
           label = c("\u2642"), size = 4, colour = "grey10",
           family="Menlo", vjust = c(1), hjust = c(0.5)) +
  #   annotate("text", x = c(0.2), y = c(80),
  #          label = c("female"), size = 2,
  #          vjust = c(0), hjust = c(0.5)) +
  # annotate("text", x = c(0.8), y = c(80),
  #          label = c("male"), size = 2,
  #          vjust = c(1), hjust = c(0.5)) +
  geom_histogram(binwidth = 0.025, data = ASR_boot, aes(x = ASR_estimate), fill = "grey30") +
  geom_errorbarh(data = ASR_boot_summary, aes(y = 150, x = lcl_ASR, xmin = lcl_ASR, xmax = ucl_ASR), 
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
        plot.margin = unit(c(0.2, 1.35, 0.405, 0.08), "cm"),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_blank(),
        # strip.text = element_text(size=11)) +
        strip.text = element_blank()) +
  ylab("Bootstrap frequency") +
  xlab("Adult sex ratio\n(proportion \u2642, 95% CI)") +
  # xlab("Adult sex ratio\n(proportion male)") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 160), expand = c(0, 0), breaks=c(0, 80, 150))
Figure_2b

ggsave(Figure_2b,
       filename = "figs_and_tabs/ASR_plot_v2.tiff",
       width = 4.75,
       height = 1.59, units = "in",
       dpi = 1200,
       compression="lzw",
       scale = 1)
