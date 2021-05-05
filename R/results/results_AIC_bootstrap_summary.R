source("scripts/01_libraries.R")

# load wrangled bootstrap output
load("output/wrangled/Black_Coucal_bootstrap_result_clean.rds")
load("output/wrangled/White-browed_Coucal_bootstrap_result_clean.rds")

# secify the levels of the species variable
BC_boot_out$nestling_AIC_tables_boot$species <- 
  factor(BC_boot_out$nestling_AIC_tables_boot$species,
         levels = c("BC",
                    "WBC"))

roles <- function(x) sub("[^_]*_","",x ) 

species_names_AIC <- c(
  'BC'="Black Coucal",
  'WBC'="White-browed Coucal")

parameter_names_AIC <- c(
  "nestling_S" = "Nestling\nsurvival",
  "fledgling_S" = "Fledgling\nsurvial",
  "fledgling_p" = "Fledgling\nlive-encounter",
  "fledgling_r" = "Fledgling\ndead-recovery",
  "fledgling_F" = "Fledgling\nsite-fidelity",
  "adult_Phi" = "Adult\napparent survival",
  "adult_p" = "Adult\nlive-encounter\nprobability"
)

# bind all the AIC results together
AIC_tables_boot <- 
  bind_rows(BC_boot_out$nestling_AIC_tables_boot,
            WBC_boot_out$nestling_AIC_tables_boot,
            BC_boot_out$fledgling_AIC_tables_boot,
            WBC_boot_out$fledgling_AIC_tables_boot,
            BC_boot_out$adult_AIC_tables_boot,
            WBC_boot_out$adult_AIC_tables_boot) %>% 
  gather(., parameter, structure, S, r, Phi, F, p) %>%
  mutate(par_stage = paste(stage, parameter, sep = "_")) %>% 
  group_by(species, par_stage) %>% 
  mutate(structure = factor(structure, levels = unique(structure[order(DeltaAICc)])))

# draw the figure
Figure_S3_Delta_AIC <-
  AIC_tables_boot %>% 
  as.data.frame() %>% 
  mutate(structure = str_replace(structure, "Time", "age")) %>% 
  mutate(structure = str_replace(structure, "time", "year")) %>%
  filter(!is.na(structure) & parameter != "Phi") %>% 
  # filter(as.data.frame(AIC_tables_boot), parameter != "Phi")
  ggplot(.,
         # ggplot(cbind(as.data.frame(AIC_tables_boot), 
         #              V4 = paste(AIC_tables_boot$species, 
         #                         AIC_tables_boot$structure, sep = "_")), 
         aes(x = structure, y = DeltaAICc)) + #reorder(V4, DeltaAICc), y = DeltaAICc) ) + 
  coord_flip() +
  theme_bw() +
  geom_boxplot(width = 0.6, fill = "grey70", outlier.size = 0.000, size = 0.2) +
  facet_grid(par_stage ~ species, labeller = labeller(species = species_names_AIC, 
                                                      par_stage = parameter_names_AIC), 
             scales = "free_y") +
  # theme(#text = element_text(family="Arial"),
  #   legend.position = "none",
  #   axis.title.x = element_blank(),
  #   axis.text.x  = element_blank(), 
  #   axis.title.y = element_text(size=10, margin = margin(0, 10, 0, 0), angle = 0),
  #   axis.text.y  = element_text(size=9, angle = 0),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   axis.ticks.y = element_line(size = 0.2, colour = "grey40"),
  #   axis.ticks.length = unit(0.1, "cm"),
  #   axis.ticks.x = element_line(size = 0.2, colour = "grey40"),
#   plot.margin = unit(c(0.2,0.3,0.2,0.2), "cm"),
#   panel.spacing = unit(0.3, "lines"),
#   strip.background = element_blank(), 
#   strip.text = element_text(size=10)) +
# scale_x_discrete(labels=roles) +
scale_y_continuous(limits=c(0,5)) + #, position = "right") +
  xlab("Model structure") + 
  ylab("Delta AIC") +
  scale_x_discrete(labels = c("~year" = "~year",
                              "~sex + year" = "~sex + year",
                              "~sex" = "~sex",
                              "~1" = "~1",
                              "~sex * Quadratic" = "~sex * age^2",
                              "~sex * Cubic" = "~sex * age^3",
                              "~sex * age" = "~sex * age"))

# save the figure
ggsave(Figure_S3_Delta_AIC,
       filename = "figs_and_tabs/AIC_breakdown.tiff",
       compression = "none",
       width = 5,
       height = 8,
       units = "in",
       dpi = 300)