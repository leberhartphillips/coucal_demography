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
BC_hazard_rate_boot <- 
  readRDS("output/bootstraps/hazard/cooked/BC_hazard_ASR_bootstrap_result_one.rds")

# load output
WBC_hazard_rate_boot <- 
  readRDS("output/bootstraps/hazard/cooked/WBC_hazard_ASR_bootstrap_result_one.rds")

# clean up the output from the bootstrap procedure and save as rds
BC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "BC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

# clean up the output from the bootstrap procedure and save as rds
WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

# calculate the sex differences in stage specific rates
BC_sex_diff_hazard_output <- 
  sex_diff_hazard(BC_hazard_rate_boot_tidy, niter = 1000) %>% 
  mutate(species = "BC")

WBC_sex_diff_hazard_output <- 
  sex_diff_hazard(WBC_hazard_rate_boot_tidy, niter = 1000) %>% 
  mutate(species = "WBC")

# consolidate results
sex_diff_survival_output <- 
  bind_rows(BC_sex_diff_hazard_output,
            WBC_sex_diff_hazard_output) %>% 
  mutate(stage = factor(stage, 
                        levels = c("Nestling", "Groundling", 
                                   "Fledgling", "Adult")),
         species = factor(species, 
                          levels = c("BC", "WBC")))

# calculate some summary statistics
sex_diff_survival_summary <- 
  sex_diff_survival_output %>%
  dplyr::group_by(stage, species) %>%
  dplyr::summarise(avg = mean(difference, na.rm = TRUE),
                   median = median(difference, na.rm = TRUE),
                   var = var(difference, na.rm = TRUE),
                   max = max(difference, na.rm = TRUE),
                   min = min(difference, na.rm = TRUE))

# specify custom color palette to distingush first-year stages 
# (i.e. chicks and juveniles) from adults
cbPalette <- c("#D9D9D9", "#D9D9D9", "#D9D9D9", "#A6A6A6", "#A6A6A6")

species_names <- 
  c('BC' = "Black Coucal",
    'WBC' = "White-browed Coucal")

junk0 <- 
  data.frame(species = c("BC", "WBC"),
             stage = "HSR",
             difference = c(NA,NA))

sex_diff_survival_output2 <- 
  bind_rows(sex_diff_survival_output, junk0) %>% 
  mutate(difference = as.numeric(difference),
         stage = factor(stage, 
                        levels = c("HSR", "Nestling", "Groundling", 
                                   "Fledgling", "Adult")),
         species = factor(species, 
                          levels = c("BC", "WBC")))

# from Safari's email on February 7th, 2020, HSR values are:
# Black coucal = 0.4955[0.4577 - 0.5334]
# White-browed coucal = 0.5198 [0.4751 - 0.5643]
HSR_df <- 
  data.frame(species = c("BC", "WBC"),
             mean = c(0.4955, 0.5198),
             upper_CI = c(0.5334, 0.5643),
             lower_CI = c(0.4577, 0.4751),
             variable = "HSR")

junk1 <- 
  data.frame(species = as.character(HSR_df$species), 
             mean = c(NA,NA), 
             lower_CI = c(NA,NA), 
             upper_CI = c(NA,NA),
             variable = "Nestling")
junk2 <- 
  data.frame(species = as.character(HSR_df$species), 
             mean = c(NA,NA), 
             lower_CI = c(NA,NA), 
             upper_CI = c(NA,NA),
             variable = "Groundling")
junk3 <- 
  data.frame(species = as.character(HSR_df$species), 
             mean = c(NA,NA), 
             lower_CI = c(NA,NA), 
             upper_CI = c(NA,NA),
             variable = "Fledgling")
junk4 <- 
  data.frame(species = as.character(HSR_df$species), 
             mean = c(NA,NA), 
             lower_CI = c(NA,NA), 
             upper_CI = c(NA,NA),
             variable = "Adult")
HSR_df2 <- 
  bind_rows(HSR_df, junk1, junk2, junk3, junk4) %>% 
  mutate(variable = factor(variable, 
                           levels = c("HSR", "Nestling", "Groundling", 
                                      "Fledgling", "Adult")),
         mean = as.numeric(mean),
         lower_CI = as.numeric(lower_CI),
         upper_CI = as.numeric(upper_CI),
         species = factor(species, 
                          levels = c("BC", "WBC")))

# Figure 2a: plot the sex-biases in survival across the three stages
vital_rate_theme <- 
  theme_bw() +
  theme(
    text = element_text(family = "Franklin Gothic Book"),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.1, "cm"),
    panel.border = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_blank()
  )

theme_set(vital_rate_theme)

surv_diff_plot <-
  ggplot(aes(y = difference, x = stage, fill = stage), 
         data = sex_diff_survival_output2) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975), color = "grey10",
              scale = "width", trim = TRUE, adjust = 1, size = 0.25) +
  # geom_jitter()+#draw_quantiles = c(0.025, 0.5, 0.975), 
  #             #scale = "width", trim = TRUE, adjust = 1, size = 0.25) +
  # coord_flip() +
  facet_grid(. ~ species, labeller = as_labeller(species_names)) +
  theme(axis.title.x = element_text(size = 7, colour = "white"),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 1, 
                                    colour = "white"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 7, hjust = 0.5, vjust = -3),
        axis.text.y  = element_text(size = 6, colour = "black"),
        axis.ticks.y = element_line(size = 0.2),
        # panel.border = element_rect(colour = "red"),
        plot.margin = unit(c(0.2, 0.39, 0.405, 1.55), "cm"),
        strip.text = element_text(size = 6, colour = "white")
  ) +
  scale_fill_manual(values = cbPalette) +
  # scale_color_manual(values = cbPalette) +
  scale_y_continuous(limits = c(-0.65, 0.65), 
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5), 
                     expand = c(0, 0), position = "right",
                     labels = c("-0.50","-0.25", 
                                expression(phantom("-")*"0.00"), 
                                expression(phantom("-")*"0.25"),
                                expression(phantom("-")*"0.50"))) +
  xlab("Life stage") +
  ylab("Sex bias in survival") +
  scale_x_discrete(labels = c("HSR" = expression(italic("\u03C1")),
                              "Nestling" = expression(S["n"]),
                              "Groundling" = expression(S["g"]),
                              "Fledgling" = expression(S["f"]),
                              "Adult" = expression(phi["ad"])))

background <-
  ggplot(aes(y = difference, x = stage, fill = stage),
         data = sex_diff_survival_output2) +
  # coord_flip() +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = -Inf, ymax = 0, alpha = 0.6,
           fill = brewer.pal(8, "Dark2")[c(1)]) +
  annotate("rect", xmin = -Inf, xmax = Inf, 
           ymin = 0, ymax = Inf, alpha = 0.6,
           fill = brewer.pal(8, "Dark2")[c(2)]) +
  annotate("text", x = c(5), y = c(-0.5),
           label = c("\u2640"), size = 4,
           family = "Menlo",
           vjust = c(0), hjust = c(0.5), colour = "grey10") +
  annotate("text", x = c(1), y = c(0.5),
           label = c("\u2642"), size = 4,
           family = "Menlo",
           vjust = c(1), hjust = c(0.5), colour = "grey10") +
  # annotate("text", x = c(2), y = c(-0.45),
  #          label = c("female"), size = 2,
  #          vjust = c(0), hjust = c(0.5)) +
  # annotate("text", x = c(2), y = c(0.45),
  #          label = c("male"), size = 2,
  #          vjust = c(1), hjust = c(0.5)) +
  facet_grid(. ~ species, labeller = as_labeller(species_names)) +
  theme(axis.title.x = element_text(size = 7, colour = "white"),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 0.9, 
                                    vjust = 1, colour = "white"),
        axis.ticks.x = element_line(size = 0.2, colour = "white"),
        axis.title.y = element_text(size = 7, hjust = 0.5, vjust = -1, 
                                    colour = "white"),
        axis.text.y  = element_text(size = 6, colour = "white"),
        axis.ticks.y = element_line(size = 0.2, colour = "white"),
        # panel.border = element_rect(colour = "black"),
        plot.margin = unit(c(0.2, 1.35, 1.05, 0.65), "cm"),
        strip.text = element_text(size = 6, colour = "white")
  ) +
  scale_x_discrete(labels = c("HSR" = expression(italic("\u03C1")),
                              "Nestling" = expression(S["n"]),
                              "Groundling" = expression(S["g"]),
                              "Fledgling" = expression(S["f"]),
                              "Adult" = expression(phi["ad"]))) +
  scale_y_continuous(limits = c(-0.65, 0.65), expand = c(0, 0)) +
  xlab("Life stage") +
  ylab(expression(paste("Hatching sex ratio\n(",italic("\u03C1"),", 95% CI)"), 
                  sep = ""))

HSR_plot <- 
  ggplot(aes(fill = variable), data = HSR_df2) +
  geom_pointrange(data = HSR_df2, aes(y = mean, x = variable, ymin = lower_CI,
                                      ymax = upper_CI), size = 0.2, 
                  fatten = 0.1) +
  geom_blank(data = HSR_df2, aes(y = mean, x = variable)) +
  # coord_flip() +
  facet_grid(. ~ species, labeller = as_labeller(species_names)) +
  theme(axis.title.x = element_text(size = 7, colour = "black"),
        axis.text.x  = element_text(size = 8, angle = 45, hjust = 0.9, 
                                    vjust = 1, colour = "black"),
        axis.ticks.x = element_line(size = 0.2, colour = "grey40"),
        axis.title.y = element_text(size = 10, hjust = 0.5, vjust = 1.5, 
                                    colour = "black"),
        axis.text.y  = element_text(size = 6, colour = "black"),
        axis.ticks.y = element_line(size = 0.2, colour = "grey40"),
        # panel.border = element_rect(colour = "blue"),
        plot.margin = unit(c(0.2, 1.35, 0.405, -0.07), "cm"),
        strip.text = element_text(size = 6, colour = "black")
  ) +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  xlab("Life stage") +
  ylab(expression(atop(x = "", 
                       y = atop(x = "Hatching sex ratio", 
                                y = paste("(", italic("\u03C1"), ", 95% CI)", 
                                          sep = ""))))) +
  # ylab(expression(paste("Hatching sex ratio\n(",italic("\u03C1"),", 95% CI)"), 
  #                 sep = "")) +
  scale_x_discrete(labels = c("HSR" = expression(italic("\u03C1")),
                              "Nestling" = expression(S["n"]),
                              "Groundling" = expression(S["g"]),
                              "Fledgling" = expression(S["f"]),
                              "Adult" = expression(phi["ad"])))

# tiff(filename = "figs_and_tabs/demographic_differences_plot.tiff",
#      compression = "none",
#      width = 4.75,
#      height = 2,
#      units = "in",
#      res = 1200)

grid.newpage()
grid::pushViewport( grid::viewport(
  layout = grid::grid.layout(1, 1, widths = unit(1, "npc")))) 
print(background, newpage = FALSE)
print(surv_diff_plot, newpage = FALSE)
print(HSR_plot, newpage = FALSE)
grid::popViewport()
# dev.off()