# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)

# custom color palette for the plotting of Juvenile and Adult stats
cbPalette <- c("#D9D9D9", "#D9D9D9", "#D9D9D9", 
               "#D9D9D9", "#A6A6A6", "#A6A6A6",
               "#A6A6A6")

species_names <- 
  c('BC' = "Black Coucal",
    'WBC' = "White-browed Coucal")

analysis_names <- c(
  'male' = "Male Mo scenario",
  'female' = "Female Mo scenario"
)


# plot the comparative LTRE results
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
    panel.spacing.x = unit(0.3, "lines"),
    panel.spacing.y = unit(0.7, "lines"),
    strip.background = element_blank()
  )

theme_set(vital_rate_theme)



LTRE_ASR_plot_background <-
  ggplot2::ggplot(data = filter(LTRE_coucal_ASR, parameter != "Immigrant sex ratio"),
                  aes(x = parameter, y = contribution, fill = parameter)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.15, ymax=0, alpha=0.6,
           fill=brewer.pal(8, "Set1")[c(1)]) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=0.15, alpha=0.6,
           fill=brewer.pal(8, "Set1")[c(2)]) +
  annotate("text", x = c(2), y = c(-0.14),
           label = c("\u2640"), size = 4, family = "Menlo",
           vjust = c(0), hjust = c(0.5), colour = "grey10") +
  annotate("text", x = c(2), y = c(0.14),
           label = c("\u2642"), size = 4, family = "Menlo",
           vjust = c(1), hjust = c(0.5), colour = "grey10") +
  facet_grid(sex ~ species, 
             labeller = labeller(.cols = species_names, .rows = analysis_names)) +
  theme(axis.title.x = element_text(size = 7, vjust = -0.1, colour = "white"),
        axis.text.x  = element_text(size = 7, angle = 45, hjust = 1, colour = "white"),
        axis.ticks.x = element_line(size = 0.2, colour = "white"),
        strip.text.x = element_text(size = 6, colour = "white"),
        
        axis.title.y = element_text(size = 7, hjust=0.5, vjust = 3.5, colour = "white"),
        axis.text.y  = element_text(size = 6, colour = "white"),
        axis.ticks.y = element_blank(),
        strip.text.y = element_text(size = 6, colour = "white"),
        
        # panel.border = element_rect(colour = "red"),
        plot.margin = unit(c(0.2, 0.85, #0.78, 
                             1.15, 0.6), "cm")) +
  
  # scale_x_discrete(labels = c("Hatching sex ratio" = expression(italic("\u03C1")),
  #                             "Nestling survival" = expression(italic(S["n"])),
  #                             "Groundling survival" = expression(italic(S["g"])),
  #                             "Fledgling survival" = expression(italic(S["f"])),
  #                             "Adult survival" = expression(italic(phi["ad"])),
  #                             "Mating system" = expression(italic("h")),
  #                             "Immigrant sex ratio" = "ISR")) +
  scale_y_continuous(limits = c(-0.15 ,0.15), expand = c(0, 0),
                     breaks = seq(from = -0.15, to = 0.15, by = 0.03)) +
  ylab("Absolute contribution to adult sex ratio bias") +
  xlab("Sex bias in parameter")

LTRE_ASR_plot <-
  ggplot2::ggplot() +
  geom_bar(data = filter(LTRE_coucal_ASR, parameter != "Immigrant sex ratio"),
           aes(x = parameter, y = contribution, fill = parameter), 
           color = "black", stat = "identity", size = 0.2) + 
  facet_grid(sex ~ species, 
             labeller = labeller(.cols = species_names, .rows = analysis_names)) +
  theme(axis.title.x = element_text(size=7, vjust=-0.1),
        axis.text.x  = element_text(size=7, angle = 45, hjust = 1),
        axis.ticks.x = element_line(size = 0.2, colour = "grey40"),
        strip.text.x = element_text(size = 6),
        
        axis.title.y = element_text(size=7, hjust=0.5, vjust = 3.5),
        axis.text.y  = element_text(size=6),
        axis.ticks.y = element_line(size = 0.2, colour = "grey40"),
        strip.text.y = element_text(size = 6),
        # panel.background = element_rect(fill = "transparent"),
        
        # panel.border = element_rect(colour = "blue"),
        plot.margin = unit(c(0.2, 0.85, 0.2, 0.6), "cm")) +
  
  scale_fill_manual(values = cbPalette) +
  
  scale_y_continuous(limits = c(-0.15 ,0.15), expand = c(0, 0),
                     breaks = seq(from = -0.150, to = 0.150, by = 0.03)) +
  
  # scale_x_discrete(labels = c("Hatching sex ratio" = expression(italic("\u03C1")),
  #                             "Nestling survival" = expression(italic(S["n"])),
  #                             "Groundling survival" = expression(italic(S["g"])),
  #                             "Fledgling survival" = expression(italic(S["f"])),
  #                             "Adult survival" = expression(italic(phi["ad"])),
  #                             "Mating system" = expression(italic("h")),
  #                             "Immigrant sex ratio" = "ISR")) +
  ylab("Absolute contribution to adult sex ratio bias") +
  xlab("Sex bias in parameter")

# jpeg(filename = "products/figures/LTRE_ASR_plot.jpeg",
#      # compression = "none",
#      width = 4.75,
#      height = 3,
#      units = "in",
#      res = 600)

grid.newpage()
pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) )
print(LTRE_ASR_plot_background, newpage = FALSE)
print(LTRE_ASR_plot, newpage = FALSE)
grid::popViewport()
# dev.off()

# Determine how much larger the contribution of each vital rates is compared to juvenile survival
# groundling vs nestling:
LTRE_coucal_ASR <- 
  filter(LTRE_coucal_ASR) %>% 
  dplyr::select(-model) %>% 
  mutate(parameter2 = as.character(parameter)) %>% 
  dplyr::rename(parameter1 = parameter)

prop_contribution_table <- 
  LTRE_coucal_ASR %>% 
  tidyr::expand(., parameter1 = parameter1, parameter2 = parameter1) %>%
  left_join(., dplyr::select(LTRE_coucal_ASR, parameter1, contribution, sex, species), 
            by = c("parameter1")) %>% 
  dplyr::rename(contribution1 = contribution) %>% 
  left_join(., dplyr::select(LTRE_coucal_ASR, parameter2, contribution, sex, species), 
            by = c("parameter2", "species", "sex")) %>% 
  dplyr::rename(contribution2 = contribution) %>% 
  mutate(prop_diff = abs(contribution1) / abs(contribution2),
         parameter1 = factor(parameter1, 
                             levels = c("Immigrant sex ratio",
                                        "Mating system",
                                        "Adult survival",
                                        "Fledgling survival",
                                        "Groundling survival",
                                        "Nestling survival",
                                        "Hatching sex ratio")),
         parameter2 = factor(parameter2, 
                             levels = c("Immigrant sex ratio",
                                        "Mating system",
                                        "Adult survival",
                                        "Fledgling survival",
                                        "Groundling survival",
                                        "Nestling survival",
                                        "Hatching sex ratio")))

LTRE_heatmap <- 
  ggplot(filter(prop_contribution_table, parameter1 %in% c(#"Immigrant sex ratio",
                                                           "Mating system",
                                                           "Adult survival",
                                                           "Fledgling survival",
                                                           "Groundling survival",
                                                           "Nestling survival",
                                                           "Hatching sex ratio") & 
                  parameter2 %in% c(#"Immigrant sex ratio",
                                    "Mating system",
                                    "Adult survival",
                                    "Fledgling survival",
                                    "Groundling survival",
                                    "Nestling survival",
                                    "Hatching sex ratio")),
         aes(parameter1, parameter2, fill = prop_diff)) + 
  geom_tile() +
  facet_grid(sex ~ species,
             labeller = labeller(.cols = species_names, .rows = analysis_names)) +
  theme(
    axis.text.x  = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 7),
    axis.title.x  = element_text(size = 7),
    axis.title.y  = element_text(size = 7),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "grey40", fill = NA),
    
  ) +
  scale_x_discrete(labels = c("Hatching sex ratio" = expression(italic("\u03C1")),
                              "Nestling survival" = expression(italic(S["n"])),
                              "Groundling survival" = expression(italic(S["g"])),
                              "Fledgling survival" = expression(italic(S["f"])),
                              "Adult survival" = expression(italic(phi["ad"])),
                              "Mating system" = expression(italic("h")),
                              "Immigrant sex ratio" = "ISR")) +
  scale_y_discrete(labels = c("Hatching sex ratio" = expression(italic("\u03C1")),
                              "Nestling survival" = expression(italic(S["n"])),
                              "Groundling survival" = expression(italic(S["g"])),
                              "Fledgling survival" = expression(italic(S["f"])),
                              "Adult survival" = expression(italic(phi["ad"])),
                              "Mating system" = expression(italic("h")),
                              "Immigrant sex ratio" = "ISR")) +
  scale_fill_gradient(low = "white", 
                      high = brewer.pal(8, "Set1")[c(2)],
                      name = "Proportional contribution of A compared to B") +#,
  # guide = guide_legend(
  #   # direction = "horizontal",
  #   title.position = "top")) +
  #   # label.position = "bottom")) +
  xlab("Parameter A") +
  ylab("Parameter B")

ggsave(LTRE_heatmap,
       filename = "products/figures/LTRE_heatmap.jpeg",
       # compression = "none",
       width = 4,
       height = 4,
       units = "in",
       dpi = 600)