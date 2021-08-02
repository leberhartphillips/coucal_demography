# load libraries
source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# load functions
function.sources = list.files(path = "R/functions", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

# read raw data
mating_dat <- 
  read.csv("data/raw/Coucal_Nr_mates_2001_2019_20200123.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>% 
  
  # make all entries lower case for consistency
  mutate(site = tolower(site),
         age_status = tolower(age_status),
         sex = tolower(sex)) %>% 
  
  # remove all white space from data
  mutate(across(.cols = everything(), 
                str_remove_all, pattern = fixed(" "))) %>% 
  
  # specify empty data as NA
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  
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
  
  # consolidate to variables of interest
  dplyr::select(year, species, ind_ID, site, sex, age_status, Nr_partners, ID_partner.s.)

# CI <- 0.95

coucal_ISR <- 
  mating_dat %>% 
  filter((sex == "female" | (sex == "male" & `ID_partner.s.` == "unringedf"))  & site == "kapunga") %>% 
  group_by(species, year, sex) %>% 
  dplyr::summarise(n_partners = sum(Nr_partners),
                   n_focals = n_distinct(ind_ID),
                   .groups = "drop") %>% 
  mutate(n_focals = ifelse(sex == "male", n_partners, n_focals)) %>% 
  group_by(species, year) %>% 
  dplyr::summarise(n_males = sum(n_partners),
                   n_females = sum(n_focals),
                   .groups = "drop") %>% 
  mutate(ISR = n_males / (n_females + n_males),
         total = n_females + n_males) %>%
  filter(!is.na(ISR) & total >= 20) %>%
  arrange(total) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  dplyr::summarise(trait = "immigration_sex_ratio", 
                   mean = mean(ISR, na.rm = TRUE),
                   # var_annual_ISR = var(ISR, na.rm = TRUE),
                   # median_annual_ISR = median(ISR, na.rm = TRUE),
                   sd = sd(ISR, na.rm = TRUE),
                   # lcl_annual_ISR = stats::quantile(ISR, (1 - CI)/2, na.rm = TRUE),
                   # ucl_annual_ISR = stats::quantile(ISR, 1 - (1 - CI)/2, na.rm = TRUE),
                   max_annual_ISR = max(ISR),
                   min_annual_ISR = min(ISR),
                   n_years = n_distinct(year),
                   mean_f = mean(n_females, na.rm = TRUE),
                   mean_M = mean(n_males, na.rm = TRUE),
                   min_f = min(n_females, na.rm = TRUE),
                   max_f = max(n_females, na.rm = TRUE),
                   min_m = min(n_males, na.rm = TRUE),
                   max_m = max(n_males, na.rm = TRUE),
                   .groups = "drop")

rm(mating_dat)