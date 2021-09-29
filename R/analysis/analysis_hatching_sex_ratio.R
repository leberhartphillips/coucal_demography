# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

# import raw csv data into R
egg_data <-
  read_xls("data/raw/Egg_surv_data_2001_2019_20210524.xls", na = "NA", col_types = "text") %>% 
  dplyr::select(species, nest_ID, year, `Total clutch size`, all_chicks, 
                nst_suc_sho, nest_success, `male_chi+eggs`, `female_chi+eggs`,
                `unkn_sex_eggs+chiks`) %>% 
  filter(species != "CTC") %>% 
  filter(!is.na(nest_success)) %>% 
  dplyr::rename(n_eggs = `Total clutch size`,
                n_chicks = all_chicks,
                male_eggs = `male_chi+eggs`, 
                female_eggs = `female_chi+eggs`,
                unknown_eggs = `unkn_sex_eggs+chiks`) %>% 
  dplyr::mutate(n_eggs = na_if(n_eggs, "?"),
                n_chicks = na_if(n_chicks, "?")) %>% 
  dplyr::mutate(n_eggs = as.numeric(n_eggs),
                n_chicks = as.numeric(n_chicks),
                male_eggs = as.numeric(male_eggs),
                female_eggs = as.numeric(female_eggs),
                unknown_eggs = as.numeric(unknown_eggs)) %>% 
  filter(n_eggs >= n_chicks)

#### hatching sex ratio summary ----
egg_HSR_data <- 
  egg_data %>% 
  dplyr::filter(n_eggs == (male_eggs + female_eggs + unknown_eggs)) %>% 
  dplyr::filter(unknown_eggs == 0) %>% 
  mutate(HSR = male_eggs/(male_eggs + female_eggs + unknown_eggs)) %>%
  mutate(species_plot = ifelse(species == "WBC", 2.2, 0.8),
         species_plot2 = ifelse(species == "WBC", 2, 1))

egg_HSR_data %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(n_distinct(nest_ID))

mod_HSR_BC <- 
  lme4::glmer(cbind(male_eggs, female_eggs) ~ 
                1 + 
                (1| nest_ID), 
              data = filter(egg_HSR_data, species == "BC"), 
              family = binomial)

mod_HSR_WBC <- 
  lme4::glmer(cbind(male_eggs, female_eggs) ~ 
                1 + 
                (1| nest_ID), 
              data = filter(egg_HSR_data, species == "WBC"), 
              family = binomial)

# extract model coefficients
mod_HSR_BC_coefs <- 
  model_parameters(mod_HSR_BC) %>%
  as.data.frame(.)

mod_HSR_WBC_coefs <- 
  model_parameters(mod_HSR_WBC) %>%
  as.data.frame(.)

# # run tidy bootstrap to obtain model diagnostics
# tidy_mod_HSR_BC <-
#   tidy(mod_HSR_BC, conf.int = TRUE, conf.method = "boot", nsim = 1000)
# 
# tidy_mod_HSR_WBC <-
#   tidy(mod_HSR_WBC, conf.int = TRUE, conf.method = "boot", nsim = 1000)

coucal_HSR <- 
  data.frame(trait = "hatching_sex_ratio",
             species = c("BC", "WBC"),
             mean = c(invlogit(model_parameters(BC_HSR)$Coefficient),
                          invlogit(model_parameters(WBC_HSR)$Coefficient)),
             CI_low = c(invlogit(model_parameters(BC_HSR)$CI_low),
                        invlogit(model_parameters(WBC_HSR)$CI_low)),
             CI_high = c(invlogit(model_parameters(BC_HSR)$CI_high),
                         invlogit(model_parameters(WBC_HSR)$CI_high)),
             n_nests = c(filter(egg_HSR_data, species == "BC") %>% 
                           summarise(n_nests = n_distinct(nest_ID)) %>% 
                           pull(n_nests),
                         filter(egg_HSR_data, species == "WBC") %>% 
                           summarise(n_nests = n_distinct(nest_ID)) %>% 
                           pull(n_nests)),
             n_years = c(filter(egg_HSR_data, species == "BC") %>% 
                           summarise(n_nests = n_distinct(year)) %>% 
                           pull(n_nests),
                         filter(egg_HSR_data, species == "WBC") %>% 
                           summarise(n_nests = n_distinct(year)) %>% 
                           pull(n_nests))) %>% 
  mutate(sd = ifelse(!is.na(CI_low), 
                         approx_sd(x1 = CI_low, x2 = CI_high),
                         CI_low),
         species_plot = ifelse(species == "WBC", 1.8, 1.2))

# HSR_plot <-
#   ggplot2::ggplot() + 
#   geom_boxplot(data = egg_HSR_data,
#                aes(x = species_plot, y = HSR,
#                    group = species, fill = species),
#                color = "grey50",
#                width = 0.05, alpha = 0.5,
#                position = position_dodge(width = 0)) +
#   geom_errorbar(data = coucal_HSR, 
#                 aes(x = species_plot, ymax = CI_high, ymin = CI_low),
#                 alpha = 1, color = "black", width = 0.05, lwd = 0.5) + 
#   geom_point(data = coucal_HSR, 
#              aes(x = species_plot, y = mean, fill = species),
#              lwd = 3, shape = 21, color= "black") +
#   geom_jitter(data = egg_HSR_data, 
#               aes(x = species_plot2, y = HSR, 
#                   group = species, 
#                   fill = species, color = species), 
#               width = 0.02, alpha = 0.2, shape = 19) +
#   # coord_flip() +
#   luke_theme +
#   theme(legend.position = "none",
#         panel.border = element_blank(),
#         strip.background = element_blank(),
#         # strip.text = element_text(size = 12, face = "italic"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(12),
#         # axis.title.y = element_blank(),
#         # axis.text.y = element_blank(),
#         # axis.ticks.x = element_blank(),
#         # axis.ticks.y = element_blank(),
#         legend.background = element_blank(),
#         # panel.grid = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey70", size = 0.1),
#         plot.title = element_text(hjust = 0.5, face = "italic")) +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
#   scale_x_discrete(labels = c("1" = "Black Coucal",
#                               "2" = "White-browed Coucal")) +
#   ylab(expression(paste("Hatching sex ratio (prop. male)" %+-%  "95% CI", sep = ""))) +
#   xlab("Species") +
#   scale_color_manual(values = c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Set1")[1])) +
#   scale_fill_manual(values = c(brewer.pal(6, "Dark2")[1], brewer.pal(6, "Set1")[1])) +
#   ggtitle("Black Coucal           White-browed\n                                Coucal")
# 
# HSR_plot
# 
# ggsave(plot = HSR_plot,
#        filename = "products/figures/jpg/HSR_plot.jpg",
#        width = 4,
#        height = 4, units = "in")

rm(egg_data, egg_HSR_data, BC_HSR, WBC_HSR)