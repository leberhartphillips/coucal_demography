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
BC_known_fate_single_mod_boot_nest <- 
  readRDS("output/bootstraps/single_models/cooked/BC_known_fate_single_mod_boot_nest.rds")

BC_known_fate_single_mod_boot_fled <- 
  readRDS("output/bootstraps/single_models/cooked/BC_known_fate_single_mod_boot_fled.rds")

# clean up the output from the bootstrap procedure and save as rds
BC_nest_single_mod_boot_tidy <- 
  single_model_boot_out_wrangle(species = "BC", stage = "nest", niter = 1000, 
                                output_dir = "output/bootstraps/single_models/cooked/",
                                rds_file = "_known_fate_single_mod_boot_")

# clean up the output from the bootstrap procedure and save as rds
BC_feld_single_mod_boot_tidy <- 
  single_model_boot_out_wrangle(species = "BC", stage = "fled", niter = 1000, 
                                output_dir = "output/bootstraps/single_models/cooked/",
                                rds_file = "_known_fate_single_mod_boot_")

BC_nest_surv_plot <- 
  ggplot() +
  luke_theme +
  theme(legend.position = "none") +
  # geom_vline(xintercept = 36.85 - 15, linetype = "dashed", alpha = 0.5, color = "orange") +
  # annotate("rect", xmin = 34.71 - 15, xmax = 38.99 - 15, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "orange") +
  # annotate("rect", xmin = 32.73 - 15, xmax = 36.94 - 15, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "green") +
  # geom_vline(xintercept = 34.85 - 15, linetype = "dashed", alpha = 0.5, color = "green") +
  geom_line(data = filter(BC_nest_single_mod_boot_tidy$reals_known_fate_boot, parameter == "S"), 
            aes(x = age, y = estimate, 
                group = interaction(iter, sex), 
                color = sex),
            alpha = 0.05) +
  scale_colour_brewer(palette = "Dark2", direction = -1) + 
  ylab("Estimated daily survival rate") +
  xlab("Age (Days since hatching)") + 
  scale_y_continuous(limits = c(0.9, 1),
                     breaks = seq(from = 0.9, to = 1, by = 0.025)) +
  scale_x_continuous(limits = c(0, 55), 
                     breaks = seq(0, 55, by = 5), 
                     expand = c(0.01, 0.01),
                     labels = as.character(seq(15, 70, by = 5)))

BC_fled_surv_plot <- 
  ggplot() +
  luke_theme +
  theme(legend.position = "none") +
  geom_vline(xintercept = 36.85, linetype = "dashed", alpha = 0.5, color = "orange") +
  annotate("rect", xmin = 34.71, xmax = 38.99, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "orange") +
  annotate("rect", xmin = 32.73, xmax = 36.94, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "green") +
  geom_vline(xintercept = 34.85, linetype = "dashed", alpha = 0.5, color = "green") +
  geom_line(data = filter(BC_feld_single_mod_boot_tidy$reals_known_fate_boot, parameter == "S"), 
            aes(x = age + 15, y = estimate, 
                group = interaction(iter, sex), 
                color = sex),
            alpha = 0.05) +
  geom_line(data = filter(BC_nest_single_mod_boot_tidy$reals_known_fate_boot, parameter == "S"), 
            aes(x = age, y = estimate, 
                group = interaction(iter, sex), 
                color = sex),
            alpha = 0.05) +
  scale_colour_brewer(palette = "Dark2", direction = -1) + 
  ylab("Estimated daily survival rate") +
  xlab("Age (Days since hatching)") + 
  scale_y_continuous(limits = c(0.9, 1),
                     breaks = seq(from = 0.9, to = 1, by = 0.025)) +
  scale_x_continuous(limits = c(0, 70), 
                     breaks = seq(0, 70, by = 5), 
                     expand = c(0.01, 0.01),
                     labels = as.character(seq(0, 70, by = 5)))