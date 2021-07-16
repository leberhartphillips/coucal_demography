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
  dplyr::filter(unknown_eggs == 0)

BC_HSR <- 
  lme4::glmer(cbind(male_eggs, female_eggs) ~ 
                1 + 
                (1| nest_ID), 
              data = filter(egg_HSR_data, species == "BC"), 
              family = binomial)

WBC_HSR <- 
  lme4::glmer(cbind(male_eggs, female_eggs) ~ 
                1 + 
                (1| nest_ID), 
              data = filter(egg_HSR_data, species == "WBC"), 
              family = binomial)

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
                         CI_low))

rm(egg_data, egg_HSR_data, BC_HSR, WBC_HSR)