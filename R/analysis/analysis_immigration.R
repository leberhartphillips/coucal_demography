# load libraries
source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# load functions
function.sources = list.files(path = "R/functions", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source, .GlobalEnv)

mating_dat <- 
  # read raw data
  read.csv("data/raw/Coucal_Nr_mates_2001_2019_20200123.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>% 
  # str()
  
  # make all entries lower case for consistency
  mutate(site = tolower(site),
         age_status = tolower(age_status),
         sex = tolower(sex)) %>% 
  
  # select variables of interest
  # select(species, ring_ID, lab_no, sex, year, site, nest_ID, pref_age, 
  #        Fledged_status, postf_age, postf_status, ageC, lay_date, hatch_order) %>% 
  
  # remove all white space from data
  mutate(across(everything(), ~str_trim(~.x))) %>% 
  # mutate(across(where(is.character), str_remove_all, pattern = fixed(" ")))
  mutate(across(.cols = everything(), 
                .fns = str_replace_all(
                  string = ..1, 
                  pattern = " ", 
                  replacement = ""))) %>% 
  
  # specify empty data as NA
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  
  # exclude all individuals that died in the nest
  # filter(Fledged_status == "yes") %>% 
  
  # classify columns
  mutate(sex = as.factor(sex),
         species = as.factor(species),
         Nr_partners = as.numeric(Nr_partners),
         age_status = as.factor(age_status),
         site = as.factor(site)) %>% 
  
  # remove rows with missing sex and age, and post-fledged status info
  filter(!is.na(sex) & !is.na(Nr_partners) & !is.na(species) & species != "CTC") %>% 
  
  # make a unique id for each individual
  mutate(ind_ID = str_replace_all(ring_ID, fixed(" "), "")) %>% 
  mutate(ID_partner.s. = str_replace_all(ID_partner.s., fixed(" "), "")) %>% 
  
  # make a unique id for each individual
  mutate(ind_ID = str_replace_all(ind_ID, fixed("_"), "")) %>% 
  mutate(ID_partner.s. = str_replace_all(ID_partner.s., fixed("_"), "")) %>% 
  
  
  # consolidate to variables of interest
  dplyr::select(year, species, ind_ID, site, sex, age_status, Nr_partners, ID_partner.s.)

CI <- 0.95

mating_dat %>% 
  group_by(species, year, site, sex) %>% 
  dplyr::summarise(n_ind = n_distinct(ind_ID)) %>% 
  pivot_wider(names_from = sex, values_from = n_ind) %>% 
  mutate(ISR = male / (female + male)) %>% 
  filter(!is.na(ISR) & site == "kapunga") %>%
  group_by(species, site) %>% 
  dplyr::summarise(mean_annual_ISR = mean(ISR, na.rm = TRUE),
                   var_annual_ISR = var(ISR, na.rm = TRUE),
                   median_annual_ISR = median(ISR, na.rm = TRUE),
                   sd_annual_ISR = sd(ISR, na.rm = TRUE),
                   lcl_annual_ISR = stats::quantile(ISR, (1 - CI)/2, na.rm = TRUE),
                   ucl_annual_ISR = stats::quantile(ISR, 1 - (1 - CI)/2, na.rm = TRUE),
                   max_annual_ISR = max(ISR),
                   min_annual_ISR = min(ISR),
                   n_years = n_distinct(year),
                   mean_f = mean(female, na.rm = TRUE),
                   mean_M = mean(male, na.rm = TRUE),
                   min_f = min(female, na.rm = TRUE),
                   max_f = max(female, na.rm = TRUE))

BC_n_years = 16
BC_ISR_mean <- 0.738
BC_ISR_CI <- 0.037
BC_sd <- sqrt(BC_n_years) * ((BC_ISR_mean + BC_ISR_CI) - (BC_ISR_mean - BC_ISR_CI))/3.92

rnorm(1, ISR_distribution[1], ISR_distribution[2])

WBC_n_years = 12
WBC_ISR_mean <- 0.524
WBC_ISR_CI <- 0.027
WBC_sd <- sqrt(WBC_n_years) * ((WBC_ISR_mean + WBC_ISR_CI) - (WBC_ISR_mean - WBC_ISR_CI))/3.92

                   
                   
                   
                   
