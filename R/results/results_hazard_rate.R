# load libraries
source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# load parameter distributions
source("R/analysis/analysis_fledge_age.R")
source("R/analysis/analysis_age_first_flight.R")

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
# 
# # clean up the output from the bootstrap procedure and save as rds
# BC_hazard_rate_boot_tidy <- 
#   hazard_boot_out_wrangle(species = "BC", niter = 1000, 
#                           output_dir = "output/bootstraps/hazard/cooked/",
#                           rds_file = "_hazard_ASR_bootstrap_result_one")
# 
# # clean up the output from the bootstrap procedure and save as rds
# WBC_hazard_rate_boot_tidy <- 
#   hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
#                           output_dir = "output/bootstraps/hazard/cooked/",
#                           rds_file = "_hazard_ASR_bootstrap_result_one")

# load tidy output
BC_hazard_rate_boot_tidy <-
  readRDS("output/bootstraps/hazard/cooked/BC_haz_sur_ASR_boot_tidy_stoc_trans_no_imm.rds")

# load tidy output
WBC_hazard_rate_boot_tidy <-
  readRDS("output/bootstraps/hazard/cooked/WBC_haz_sur_ASR_boot_tidy_stoc_no_imm.rds")

flight_dat <-
  bind_rows(coucal_fledge_age, coucal_flight_age) %>% 
  dplyr::select(trait, species, sex, mean) %>% 
  pivot_wider(names_from = trait, values+)
  group_by(species, sex) %>% 
  mutate(mean = ifelse(trait == "flight_age", )) %>% 
  mutate(mean = round(mean))
  
flight_dat <- 
  data.frame(species = c("BC","WBC"),
             end_nestling = c(13, 14),
             end_nestling_lower = c(12, 13),
             end_nestling_upper = c(13, 15),
             end_groundling = c(36, 32),
             end_groundling_lower = c(34, 29),
             end_groundling_upper = c(38, 35))

grp_means <- 
  bind_rows(BC_hazard_rate_boot_tidy$vital_rate_ests_boot,
          WBC_hazard_rate_boot_tidy$vital_rate_ests_boot) %>% 
  filter(rate == "development") %>% 
  pivot_wider(names_from = stage, values_from = value) %>% 
  mutate(flight_age = fledge_age + flight_age) %>% 
  pivot_longer(c(fledge_age, flight_age), names_to = "stage", values_to = "value") %>% 
  mutate(stage_sex = paste(sex, stage, sep = "_")) %>% 
  group_by(species, stage_sex, sex, stage) %>% 
  dplyr::summarise(grp.mean = mean(value))

BC_density_plot <- 
  bind_rows(BC_hazard_rate_boot_tidy$vital_rate_ests_boot,
          WBC_hazard_rate_boot_tidy$vital_rate_ests_boot) %>% 
  filter(rate == "development" & species == "BC") %>% 
  pivot_wider(names_from = stage, values_from = value) %>% 
  mutate(flight_age = fledge_age + flight_age) %>% 
  pivot_longer(c(fledge_age, flight_age), names_to = "stage", values_to = "value") %>% 
  mutate(stage_sex = paste(sex, stage, sep = "_")) %>% 
  ggplot() +
  geom_density(aes(value, y=..scaled.., fill = fct_rev(stage_sex)),
               alpha = 0.7, color = "grey40") + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(values = rev(sex_pal3)) +
  scale_x_continuous(limits = c(0, 70),
                     expand = c(0.01, 0.01))

WBC_density_plot <-
  bind_rows(BC_hazard_rate_boot_tidy$vital_rate_ests_boot,
            WBC_hazard_rate_boot_tidy$vital_rate_ests_boot) %>% 
  filter(rate == "development" & species == "WBC") %>% 
  pivot_wider(names_from = stage, values_from = value) %>% 
  mutate(flight_age = fledge_age + flight_age) %>% 
  pivot_longer(c(fledge_age, flight_age), names_to = "stage", values_to = "value") %>% 
  mutate(stage_sex = paste(sex, stage, sep = "_")) %>% 
  ggplot() +
  geom_density(aes(value, y=..scaled.., fill = fct_rev(stage_sex)),
               alpha = 0.7, color = "grey40") + 
  theme_void() +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(values = rev(sex_pal3)) +
  scale_x_continuous(limits = c(0, 70),
                     expand = c(0.01, 0.01))

BC_surv_plot <-
  ggplot() +
  luke_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0.5,0), "cm")) +
  geom_vline(data = filter(grp_means, species == "BC" & stage == "fledge_age"),
             aes(xintercept = grp.mean, color = sex), linetype = "dashed", alpha = 1) +
  geom_vline(data = filter(grp_means, species == "BC" & stage == "flight_age"),
             aes(xintercept = grp.mean, color = sex), linetype = "dashed", alpha = 1) +
  geom_line(data = filter(BC_hazard_rate_boot_tidy$hazard_rates_boot, 
                          iter %in% c(as.character(seq(from = 1, to = 500, by = 1)))),
            aes(x = age, y = estimate, 
                group = interaction(iter, sex), 
                color = sex),
            alpha = 0.05) +
  facet_grid(species ~ ., labeller = as_labeller(species_names)) +
  scale_colour_manual(values = sex_pal2) +
  ylab("Estimated daily survival rate") +
  xlab("Age (Days since hatching)") + 
  scale_y_continuous(limits = c(0.9, 1),
                     breaks = seq(from = 0.9, to = 1, by = 0.025)) +
  scale_x_continuous(limits = c(0, 70), 
                     breaks = seq(0, 70, by = 5), 
                     expand = c(0.01, 0.01),
                     labels = as.character(seq(0, 70, by = 5))) +
  annotate(geom = "text", y = 0.91, 
           x = mean(c(pull(filter(grp_means, species == "BC" & stage == "fledge_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "BC" & stage == "fledge_age")[2,5], grp.mean)))/2,
           label = "nestling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = 0.91, 
           x = mean(c(pull(filter(grp_means, species == "BC" & stage == "fledge_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "BC" & stage == "fledge_age")[2,5], grp.mean))) +
             ((mean(c(pull(filter(grp_means, species == "BC" & stage == "flight_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "BC" & stage == "flight_age")[2,5], grp.mean))) -
                 mean(c(pull(filter(grp_means, species == "BC" & stage == "fledge_age")[1,5], grp.mean),
                        pull(filter(grp_means, species == "BC" & stage == "fledge_age")[2,5], grp.mean))))/2),
           label = "groundling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = 0.91, 
           x = mean(c(pull(filter(grp_means, species == "BC" & stage == "flight_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "BC" & stage == "flight_age")[2,5], grp.mean))) +
             ((70 - mean(c(pull(filter(grp_means, species == "BC" & stage == "flight_age")[1,5], grp.mean), 
                           pull(filter(grp_means, species == "BC" & stage == "flight_age")[2,5], grp.mean))))/2),
           label = "fledgling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5)

WBC_surv_plot <-
  ggplot() +
  luke_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        plot.margin = unit(c(0,0,0.5,0), "cm")) +
  geom_vline(data = filter(grp_means, species == "WBC" & stage == "fledge_age"),
             aes(xintercept = grp.mean, color = sex), linetype = "dashed", alpha = 1) +
  geom_vline(data = filter(grp_means, species == "WBC" & stage == "flight_age"),
             aes(xintercept = grp.mean, color = sex), linetype = "dashed", alpha = 1) +
  geom_line(data = filter(WBC_hazard_rate_boot_tidy$hazard_rates_boot, 
                          iter %in% c(as.character(seq(from = 1, to = 500, by = 1)))),
            aes(x = age, y = estimate, 
                group = interaction(iter, sex), 
                color = sex),
            alpha = 0.05) +
  facet_grid(species ~ ., labeller = as_labeller(species_names)) +
  scale_colour_manual(values = sex_pal2) +
  ylab("                                                         Estimated daily survival rate") +
  xlab("Age (days since hatching)") + 
  scale_y_continuous(limits = c(0.9, 1),
                     breaks = seq(from = 0.9, to = 1, by = 0.025)) +
  scale_x_continuous(limits = c(0, 70), 
                     breaks = seq(0, 70, by = 5), 
                     expand = c(0.01, 0.01),
                     labels = as.character(seq(0, 70, by = 5))) +
  annotate(geom = "text", y = 0.91, 
           x = mean(c(pull(filter(grp_means, species == "WBC" & stage == "fledge_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "WBC" & stage == "fledge_age")[2,5], grp.mean)))/2,
           label = "nestling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = 0.91, 
           x = mean(c(pull(filter(grp_means, species == "WBC" & stage == "fledge_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "WBC" & stage == "fledge_age")[2,5], grp.mean))) +
             ((mean(c(pull(filter(grp_means, species == "WBC" & stage == "flight_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "WBC" & stage == "flight_age")[2,5], grp.mean))) -
                 mean(c(pull(filter(grp_means, species == "WBC" & stage == "fledge_age")[1,5], grp.mean),
                        pull(filter(grp_means, species == "WBC" & stage == "fledge_age")[2,5], grp.mean))))/2),
           label = "groundling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = 0.91, 
           x = mean(c(pull(filter(grp_means, species == "WBC" & stage == "flight_age")[1,5], grp.mean),
                      pull(filter(grp_means, species == "WBC" & stage == "flight_age")[2,5], grp.mean))) +
             ((70 - mean(c(pull(filter(grp_means, species == "WBC" & stage == "flight_age")[1,5], grp.mean), 
                           pull(filter(grp_means, species == "WBC" & stage == "flight_age")[2,5], grp.mean))))/2),
           label = "fledgling",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5)
  
coucal_juvenile_plot <-
  (BC_density_plot / BC_surv_plot / WBC_density_plot / WBC_surv_plot) + 
  plot_layout(heights = unit(c(0.5, 4, 0.5, 4), c('cm', 'cm', 'cm', 'cm')),
              widths = unit(c(7, 7, 7, 7), c('cm', 'cm', 'cm', 'cm')))

ggsave(plot = coucal_juvenile_plot,
       filename = "products/figures/coucal_juvenile_plot.jpg",
       width = 7 * 1.4,
       height = sum(c(0.5, 4, 0.5, 4)) * 1.3, units = "cm", dpi = 600)
  
# filter(BC_hazard_rate_boot_tidy$hazard_rates_boot, age == 0) %>% 
#   ggplot() +
#   geom_histogram(aes(fit), binwidth = 0.001) +
#   facet_grid(sex~ .)
# 
# filter(WBC_hazard_rate_boot_tidy$hazard_rates_boot, age == 12) %>% 
#   ggplot() +
#   geom_histogram(aes(fit), binwidth = 0.001) +
#   facet_grid(sex ~ .)
# 
# # plot each iteration's hazard function
# surv_plot <-
#   ggplot() +
#   luke_theme +
#   theme(legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(size = 12, face = "italic")) +
#   geom_vline(data = filter(flight_dat, trait == "fledge_age"),
#              aes(xintercept = mean), linetype = "dashed", alpha = 1) +
#   geom_vline(data = filter(flight_dat, trait == "flight_age"),
#              aes(xintercept = mean), linetype = "dashed", alpha = 1) +
#   # facet_grid(species ~ ., labeller = as_labeller(species_names)) +
#   # geom_vline(data = flight_dat,
#   #            aes(xintercept = end_groundling),
#   #            linetype = "dashed", alpha = 0.5, color = "grey20") +
#   # geom_rect(data = filter(flight_dat, trait == "fledge_age" & sex == "F"),
#   #          geom = "rect", 
#   #          xmin = round(pull(filter(flight_dat, trait == "fledge_age" & sex == "F"), CI_low)), xmax = round(CI_high), 
#   #          ymin = -Inf, ymax = Inf,
#   #          alpha = 0.2, 
#   #          fill = sex_pal2[1]) +
#   # geom_rect(data = filter(flight_dat, trait == "fledge_age" & sex == "M"),
#   #           aes(xmin = round(CI_low), xmax = round(CI_high),
#   #               ymin = -Inf, ymax = Inf, group = species), 
#   #               alpha = 1, fill = sex_pal2[2]) +
#   # annotate(data = filter(flight_dat, trait == "flight_age" & sex == "F"),
#   #          aes(xmin = round(CI_low), xmax = round(CI_high),
#   #          ymin = -Inf, ymax = Inf), alpha = 0.2, fill = sex_pal2[1], geom = "rect") +
#   # annotate(data = filter(flight_dat, trait == "flight_age" & sex == "M"),
#   #          aes(xmin = round(CI_low), xmax = round(CI_high),
#   #          ymin = -Inf, ymax = Inf), alpha = 0.2, fill = sex_pal2[2], geom = "rect") +
#   geom_line(data = bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot,
#                              WBC_hazard_rate_boot_tidy$hazard_rates_boot)[c(1:2000), ],
#             aes(x = age, y = estimate, 
#                 group = interaction(iter, sex), 
#                 color = sex),
#             alpha = 0.01) +
#   scale_colour_manual(values = rev(sex_pal2)) +
#   ylab("Estimated daily survival rate") +
#   xlab("Age (Days since hatching)") + 
#   scale_y_continuous(limits = c(0.9, 1),
#                      breaks = seq(from = 0.9, to = 1, by = 0.025)) +
#   scale_x_continuous(limits = c(0, 70), 
#                      breaks = seq(0, 70, by = 5), 
#                      expand = c(0.01, 0.01),
#                      labels = as.character(seq(0, 70, by = 5))) +
#   facet_grid(species ~ ., labeller = as_labeller(species_names)) +
#   annotate(geom = "text", y = 0.9, x = 15/2,
#            label = "nestling",
#            color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
#   annotate(geom = "text", y = 0.9, x = (36 - 15)/2 + 15,
#            label = "groundling",
#            color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
#   annotate(geom = "text", y = 0.9, x = (70 - 36)/2 + 36,
#            label = "fledgling",
#            color = "black", size = 3, fontface = 'italic', hjust = 0.5)
# surv_plot
# 
# ggsave(surv_plot,
#        filename = "products/figures/offsring_hazard_functions.jpeg",
#        # compression = "none",
#        width = 6,
#        height = 6,
#        units = "in",
#        dpi = 600)

CI <- 0.95

sex_diff_background <-
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot) %>% 
  dplyr::select(species, age, sex, iter, estimate) %>% 
  pivot_wider(names_from = c(sex), values_from = c(estimate)) %>% 
  mutate(sex_diff = Male - Female,
         age_f = as.factor(age),
         species = factor(species, levels = c("BC", "WBC"))) %>% 
  group_by(species, age_f) %>% 
  ggplot(.) +
  geom_line(aes(x = age, y = sex_diff,
                group = iter), color = "black",
            alpha = 0.05) +
  annotate("rect", ymin = 0, ymax = 0.125, xmin = 0, xmax = 70, alpha = 0.6,
           fill= sex_pal2[c(2)]) +
  annotate("rect", ymin = -0.125, ymax = 0, xmin = 0, xmax = 70, alpha = 0.6,
           fill= sex_pal2[c(1)]) +
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
        strip.text = element_text(size = 12, face = "italic", color = "white"),
        axis.title.x = element_text(colour = "white"),
        axis.text.x  = element_text(colour = "white"),
        axis.ticks.x = element_line(colour = "white"),
        axis.title.y = element_text(colour = "white"),
        axis.text.y  = element_text(colour = "white"),
        axis.ticks.y = element_line(colour = "white")) +
  ylab("Estimated sex bias in daily mortality rate") +
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

CI = 0.95

sex_diff_foreground <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  bind_rows(BC_hazard_rate_boot_tidy$hazard_rates_boot) %>% 
  dplyr::select(species, age, sex, iter, estimate) %>% 
  pivot_wider(names_from = c(sex), values_from = c(estimate)) %>% 
  mutate(sex_diff = Male - Female,
         age_f = as.factor(age),
         species = factor(species, levels = c("BC", "WBC"))) %>% 
  group_by(species, age_f) %>% 
  dplyr::summarise(lcl_diff = stats::quantile(sex_diff, (1 - CI)/2, na.rm = TRUE),
                   ucl_diff = stats::quantile(sex_diff, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_diff = mean(sex_diff),
                   med_diff = median(sex_diff),
                   # max_diff = max(sex_diff),
                   # min_diff = min(sex_diff),
                   lcl_Male = stats::quantile(Male, (1 - CI)/2, na.rm = TRUE),
                   ucl_Male = stats::quantile(Male, 1 - (1 - CI)/2, na.rm = TRUE),
                   # avg_Male = mean(Male),
                   med_Male = median(Male),
                   # max_Male = max(Male),
                   # min_Male = min(Male),
                   lcl_Female = stats::quantile(Female, (1 - CI)/2, na.rm = TRUE),
                   ucl_Female = stats::quantile(Female, 1 - (1 - CI)/2, na.rm = TRUE),
                   # avg_Female = mean(Female),
                   med_Female = median(Female),
                   # max_Female = max(Female),
                   # min_Female = min(Female)
                   ) %>% 
  ggplot(.) +
  luke_theme +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  # geom_vline(data = flight_dat,
  #            aes(xintercept = end_nestling), linetype = "dashed", alpha = 0.5, color = "grey20") +
  # geom_vline(data = flight_dat,
  #            aes(xintercept = end_groundling),
  #            linetype = "dashed", alpha = 0.5, color = "grey20") +
  geom_errorbar(aes(ymin = lcl_diff, ymax = ucl_diff, y = avg_diff,
                    x = age_f), 
                 color = "white", size = 0.3, linetype = "solid", alpha = 0.5) +
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
grid::popViewport()
# dev.off()
