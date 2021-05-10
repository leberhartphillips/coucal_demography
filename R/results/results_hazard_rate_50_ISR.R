# load libraries
source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# load functions
function.sources = list.files(path = "R/functions", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source, .GlobalEnv)

# load capture histories
data.sources = list.files(path = "data/cooked", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

# load output
BC_hazard_rate_boot_50_ISR <- 
  readRDS("output/bootstraps/hazard/cooked/BC_hazard_ASR_bootstrap_result_50_ISR.rds")

# load output
WBC_hazard_rate_boot <- 
  readRDS("output/bootstraps/hazard/cooked/WBC_hazard_ASR_bootstrap_result_one.rds")

# clean up the output from the bootstrap procedure and save as rds
BC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "BC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_50_ISR")

# clean up the output from the bootstrap procedure and save as rds
WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

# plot each iteration's hazard function
# surv_plot <-
  ggplot() +
  luke_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic")) +
  geom_vline(xintercept = 36, linetype = "dashed", alpha = 0.5, color = "grey20") +
  geom_vline(xintercept = 15, linetype = "dashed", alpha = 0.5, color = "grey20") +
  geom_line(data = bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot,
                             WBC_hazard_rate_boot_tidy$hazard_rates_boot),
            aes(x = age, y = estimate, 
                group = interaction(iter, sex), 
                color = sex),
            alpha = 0.05) +
  scale_colour_manual(values = rev(plot_palette_sex)) +
  ylab("Estimated daily survival rate") +
  xlab("Age (Days since hatching)") + 
  scale_y_continuous(limits = c(0.9, 1),
                     breaks = seq(from = 0.9, to = 1, by = 0.025)) +
  scale_x_continuous(limits = c(0, 70), 
                     breaks = seq(0, 70, by = 5), 
                     expand = c(0.01, 0.01),
                     labels = as.character(seq(0, 70, by = 5))) +
  facet_grid(species ~ ., labeller = as_labeller(species_names)) +
  annotate(geom = "text", y = 0.9, x = 15/2,
           label = "nestling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = 0.9, x = (36 - 15)/2 + 15,
           label = "groundling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = 0.9, x = (70 - 36)/2 + 36,
           label = "fledgling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5)
surv_plot

CI <- 0.95

sex_diff_background <-
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot) %>% 
  select(species, age, sex, iter, fit) %>% 
  pivot_wider(names_from = c(sex), values_from = c(fit)) %>% 
  mutate(sex_diff = Male - Female,
         age_f = as.factor(age),
         species = factor(species, levels = c("BC", "WBC"))) %>% 
  group_by(species, age_f) %>% 
  ggplot(.) +
  geom_line(aes(x = age, y = sex_diff,
                group = iter), color = "black",
            alpha = 0.05) +
  annotate("rect", ymin = 0, ymax = 0.125, xmin = 0, xmax = 70, alpha = 0.6,
           fill= brewer.pal(8, "Dark2")[c(1)]) +
  annotate("rect", ymin = -0.125, ymax = 0, xmin = 0, xmax = 70, alpha = 0.6,
           fill= brewer.pal(8, "Dark2")[c(2)]) +
  annotate("text", x = c(1), y = c(0.05),
           label = c("\u2642"), size = 4, colour = "grey10",
           family="Menlo",
           vjust = c(0), hjust = c(0)) +
  annotate("text", x = c(1), y = c(-0.05),
           label = c("\u2640"), size = 4, colour = "grey10",
           family="Menlo",
           vjust = c(1), hjust = c(0)) +
  luke_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        axis.title.x = element_text(colour = "white"),
        axis.text.x  = element_text(colour = "white"),
        axis.ticks.x = element_line(colour = "white"),
        axis.title.y = element_text(colour = "white"),
        axis.text.y  = element_text(colour = "white"),
        axis.ticks.y = element_line(colour = "white")) +
  ylab("Estimated sex bias in daily survival rate") +
  xlab("Age (Days since hatching)") + 
  scale_y_continuous(limits = c(-0.125, 0.125),
                     breaks = seq(from = -0.1, to = 0.1, by = 0.05),
                     expand = c(0, 0)) +
  facet_grid(species ~ ., labeller = as_labeller(species_names)) +
  scale_x_continuous(breaks = seq(0, 70, by = 5),
                     expand = c(0, 0),
                     labels = as.character(seq(0, 70, by = 5))) +
  annotate(geom = "text", y = -0.115, x = 15/2,
           label = "nestling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = -0.115, x = (36 - 15)/2 + 15,
           label = "groundling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = -0.115, x = (70 - 36)/2 + 36,
           label = "fledgling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5)

sex_diff_foreground <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot) %>% 
  select(species, age, sex, iter, fit) %>% 
  pivot_wider(names_from = c(sex), values_from = c(fit)) %>% 
  mutate(sex_diff = Male - Female,
         age_f = as.factor(age),
         species = factor(species, levels = c("BC", "WBC"))) %>% 
  group_by(species, age_f) %>% 
  dplyr::summarise(lcl_diff = stats::quantile(sex_diff, (1 - CI)/2, na.rm = TRUE),
                   ucl_diff = stats::quantile(sex_diff, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_diff = mean(sex_diff),
                   med_diff = median(sex_diff),
                   max_diff = max(sex_diff),
                   min_diff = min(sex_diff)) %>% 
  ggplot(.) +
  luke_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  geom_vline(xintercept = "36", linetype = "dashed", alpha = 0.5, color = "black") +
  geom_vline(xintercept = "15", linetype = "dashed", alpha = 0.5, color = "black") +
  geom_errorbar(aes(ymin = lcl_diff, ymax = ucl_diff, y = avg_diff,
                    x = age_f), 
                color = "white", size = 0.3, linetype = "solid") +
  geom_point(aes(y = med_diff, x = age_f), size = 0.3, color = "white") +
  # geom_hline(yintercept = 0, alpha = 0.5, color = "black") +
  facet_grid(species ~ ., labeller = as_labeller(species_names)) +
  scale_colour_manual(values = rev(plot_palette_sex)) +
  ylab("Estimated sex bias in daily survival rate (± 95% CI)") +
  xlab("Age (Days since hatching)") + 
  scale_y_continuous(limits = c(-0.125, 0.125),
                     breaks = seq(from = -0.1, to = 0.1, by = 0.05),
                     expand = c(0, 0)) +
  scale_x_discrete(breaks = as.character(seq(0, 70, by = 5)),
                   expand = c(0.01, 0.01),
                   labels = as.character(seq(0, 70, by = 5)))

# jpeg(filename = "products/figures/sex_differences_plot.jpeg",
#      width = 4.75,
#      height = 6,
#      units = "in",
#      res = 600)

grid.newpage()
grid::pushViewport( grid::viewport(
  layout = grid::grid.layout(1, 1, widths = unit(1, "npc")))) 
print(sex_diff_background, newpage = FALSE)
print(sex_diff_foreground, newpage = FALSE)
# print(HSR_plot, newpage = FALSE)
grid::popViewport()
# dev.off()

ASR_boot <- 
  bind_rows(BC_hazard_rate_boot_tidy$ASR_ests_boot,
            WBC_hazard_rate_boot_tidy$ASR_ests_boot) %>%
  mutate(species = factor(species, levels = c("BC", "WBC")))

CI <- 0.95

ASR_boot_summary <- 
  ASR_boot %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(ucl_ASR = stats::quantile(M_Adult, (1 - CI)/2, na.rm = TRUE),
                   lcl_ASR = stats::quantile(M_Adult, 1 - (1 - CI)/2, na.rm = TRUE),
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
  annotate("rect", xmin=0.25, xmax=0.5, ymin=0, ymax = 200, alpha=0.6,
           fill= brewer.pal(8, "Dark2")[c(2)]) +
  annotate("rect", xmin=0.5, xmax=0.75, ymin=0, ymax = 200, alpha=0.6,
           fill= brewer.pal(8, "Dark2")[c(1)]) +
  annotate("text", x = c(0.3), y = c(100),
           label = c("\u2640"), size = 4, colour = "grey10",
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  annotate("text", x = c(0.7), y = c(100),
           label = c("\u2642"), size = 4, colour = "grey10",
           family="Menlo", vjust = c(0.5), hjust = c(0.5)) +
  geom_histogram(binwidth = 0.01, data = ASR_boot, aes(x = M_Adult), fill = "grey30") +
  geom_errorbarh(data = ASR_boot_summary, aes(y = 175, x = lcl_ASR, xmin = lcl_ASR, xmax = ucl_ASR), 
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
  scale_x_continuous(limits = c(0.25, 0.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 200), expand = c(0, 0), breaks=c(0, 50, 100, 150))
Figure_2b