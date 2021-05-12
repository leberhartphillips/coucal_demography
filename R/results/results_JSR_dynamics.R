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
BC_JSR_out <- 
  readRDS("output/bootstraps/hazard/cooked/BC_hazard_JSR_bootstrap_result.rds")

# load output
WBC_JSR_out <- 
  readRDS("output/bootstraps/hazard/cooked/WBC_hazard_JSR_bootstrap_result.rds")

JSR_boot <- 
  bind_rows(BC_JSR_out,
            WBC_JSR_out) %>%
  mutate(species = factor(species, levels = c("BC", "WBC")))

CI <- 0.95

end_JSR_boot <- 
  JSR_boot %>%
  filter(!is.na(JSR)) %>% 
  mutate(JSR = ifelse(JSR < 0, 0, JSR)) %>% 
  arrange(species, iteration, desc(age)) %>% 
  dplyr::group_by(species, iteration) %>%
  slice(1)

JSR_boot_summary <- 
  end_JSR_boot %>% 
  dplyr::group_by(species) %>%
  dplyr::summarise(ucl_JSR = stats::quantile(JSR, (1 - CI)/2, na.rm = TRUE),
                   lcl_JSR = stats::quantile(JSR, 1 - (1 - CI)/2, na.rm = TRUE),
                   avg_JSR = mean(JSR, na.rm = TRUE),
                   med_JSR = median(JSR, na.rm = TRUE),
                   max_JSR = max(JSR, na.rm = TRUE),
                   min_JSR = min(JSR, na.rm = TRUE))

#### plot JSR dynamics and distribution ----
ggplot(data = bind_rows(BC_JSR_out,
                        BC_JSR_out)) +
  geom_line(aes(y = JSR, x = as.numeric(age), group = iteration),
            alpha = 0.05) +
  geom_hline(yintercept = 0.5, color = "white") +
  facet_grid(species ~ .)

# Figure_2b <- 
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
  geom_histogram(binwidth = 0.01, data = end_JSR_boot, aes(x = JSR), fill = "grey30") +
  geom_errorbarh(data = JSR_boot_summary, aes(y = 175, x = lcl_JSR, xmin = lcl_JSR, xmax = ucl_JSR), 
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
  xlab("Juvenile sex ratio\n(proportion \u2642, 95% CI)") +
  # xlab("Adult sex ratio\n(proportion male)") +
  scale_x_continuous(limits = c(0.25, 0.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 200), expand = c(0, 0), breaks=c(0, 50, 100, 150))
# Figure_2b