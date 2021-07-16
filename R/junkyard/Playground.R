# playground

#### wrangle ----
status_dat_all <- 
  # read raw data
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>% 
  
  # rename ring_ID column
  dplyr::rename(ring_ID = Ring_ID) %>% 
  
  # make all entries lower case for consistency
  mutate(Fledged_status = tolower(Fledged.),
         site = tolower(site)) %>% 
  
  # select variables of interest
  select(species, ring_ID, lab_no, sex, year, site, nest_ID, pref_age, 
         Fledged_status, postf_age, postf_status, ageC, statusC, 
         lay_date, hatch_order) %>% 
  
  # remove all white space from data
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  
  # specify empty data as NA
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  
  # exclude all individuals that died in the nest
  # filter(Fledged_status == "yes") %>% 
  
  # classify columns
  mutate(sex = as.factor(sex),
         ageC = as.numeric(ageC),
         postf_age = as.numeric(postf_age),
         postf_status = as.numeric(postf_status),
         hatch_order = as.numeric(hatch_order),
         pref_age = as.numeric(pref_age)) %>% 
  
  # remove rows with missing sex, age, and status info
  filter(!is.na(sex) & !is.na(ageC) & !is.na(statusC)) %>% 
  
  # make a unique id for each individual
  mutate(ind_ID = paste(nest_ID, lab_no, ring_ID, sep = "_"),
         
         # create the age of entry into the data (all at age 15)
         entry = 0,
         
         # specify the age of death or censoring
         exit = ageC,
         
         # make the event numeric and specify if 
         # the individual died (1) or was censored (0)
         event = as.numeric(statusC)) %>% 
  
  # consolidate to variables of interest
  dplyr::select(species, ind_ID, nest_ID, year, sex, entry, exit, event)

BC_all_dat <- 
  status_dat_all %>% 
  filter(species == "BC")

WBC_all_dat <- 
  status_dat_all %>% 
  filter(species == "WBC")

#### specs ----
offspring = BC_all_dat
k = 4
HSR = 0.4955
h = 1/2.9
egg_survival = 0.32
# ISR = 0.738
ISR = 0.5
immigrant_pop_size = 100
fledge_age = 15
flight_age = 36
bootstrap_name = "BC_boot_one"
species = "BC"
iter_add = 1
prefix_number = "boot_one_"
niter = 2
adult_survival_rate = 0.3
alpha_value = 1.4
time_vector = seq(0, 70, 1)

num_boot = 1
offspring = WBC_all_dat
k = 4
HSR = 0.5198
h = 1/1.1
egg_survival = 0.18
ISR = 0.524
adult_survival_rate = 0.3
immigrant_pop_size = 100
fledge_age = 15
flight_age = 32
bootstrap_name = "WBC_boot_one"
species = "WBC"
iter_add = 1
prefix_number = "boot_one_"
max_time = 70

rm(hzd_curve_try_M, hzd_curve_try_F, F_haz_ss_curve, M_haz_ss_curve, M_haz_ss, F_haz_ss, bootstrap_data_list, time_vector)

set.seed(2)
#### model runs ----
# run the sampling function and specify the datasets
bootstrap_data_list <- 
  bootstrap_hazard_data(offspring = offspring, 
                        num_boot = num_boot,
                        species = species, 
                        iter_add = iter_add,
                        alpha_value = alpha_value,
                        max_time = max_time)

########
coucal_boot_list = bootstrap_data_list
coucal_boot_list
# specify the bootstrapped data samples (from the previous function)
offspring_data <- coucal_boot_list[["offspring_boot"]]

# clean up capture histories
offspring_data <-    
  offspring_data %>% 
  ungroup() %>% 
  as.data.frame()

# fit smoothed spline of hazard function for either sex
M_haz_ss <- sshzd(Surv(exit, event, entry) ~ exit, 
                  data = filter(offspring_data, sex == "M"), 
                  alpha = alpha_value)

F_haz_ss <- sshzd(Surv(exit, event, entry) ~ exit, 
                  data = filter(offspring_data, sex == "F"), 
                  alpha = alpha_value)

haz_ss_function <- list(Male_haz_ss = M_haz_ss,
                        Female_haz_ss = F_haz_ss)

# extract fitted estimates from the spline function
M_haz_ss_curve <- 
  hzdcurve.sshzd(object = M_haz_ss, time = coucal_boot_list[["time_vector"]], se = TRUE)

F_haz_ss_curve <- 
  hzdcurve.sshzd(object = F_haz_ss, time = coucal_boot_list[["time_vector"]], se = TRUE)

haz_ss_curve <- 
  expand.grid(species = species, 
              age = coucal_boot_list[["time_vector"]],
              sex = c("Male", "Female")) %>% 
  mutate(fit = c(M_haz_ss_curve$fit, F_haz_ss_curve$fit),
         se = c(M_haz_ss_curve$se, F_haz_ss_curve$se)) %>% 
  mutate(estimate = 1 - fit,
         upper = 1 - fit * exp(1.96 * se),
         lower = 1 - fit / exp(1.96 * se),
         iter = num_boot)

# transform the daily nestling survival (DCS) to apparent fledgling success
# by calculating the product of all DCS estimates:
coucal_nestling_survival <-
  haz_ss_curve %>% 
  filter(age <= fledge_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>%
  mutate(stage = "nestling",
         rate = "survival")

coucal_groundling_survival <-
  haz_ss_curve %>% 
  filter(age < flight_age & age > fledge_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "groundling",
         rate = "survival")

coucal_fledgling_survival <-
  haz_ss_curve %>% 
  filter(age >= flight_age) %>% 
  group_by(sex) %>% 
  dplyr::summarise(value = prod(estimate), .groups = 'drop') %>% 
  mutate(stage = "fledgling",
         rate = "survival")

coucal_adult_survival <-
  expand.grid(sex = c("Male", "Female"),
              value = adult_survival_rate,
              stage = c("adult"),
              rate = c("survival"))

adult_F_immigrants <- immigrant_pop_size * (1 - ISR)
adult_M_immigrants <- immigrant_pop_size * ISR

coucal_adult_immigration <- 
  data.frame(sex = c("Female", "Male"),
             value = c(adult_F_immigrants, adult_M_immigrants),
             stage = c("adult"),
             rate = c("immigration"))

coucal_egg_survival <- 
  data.frame(sex = NA,
             value = egg_survival,
             stage = c("egg"),
             rate = c("survival"))

# Bind the juvenile and adult dataframe with the nestlings
coucal_vital_rates <- 
  bind_rows(coucal_egg_survival,
            coucal_nestling_survival,
            coucal_groundling_survival,
            coucal_fledgling_survival,
            coucal_adult_survival,
            coucal_adult_immigration) %>% 
  as.data.frame() %>% 
  mutate(iter = num_boot, #+ ((iter_add - 1) * niter),
         species = species) %>% 
  arrange(sex, stage, rate)

# Create a list of demographic rates from the survival analyses above
demographic_rates <- list(Egg_survival = coucal_vital_rates[11, 2],
                          F_Nestling_survival = coucal_vital_rates[5, 2],
                          F_Groundling_survival = coucal_vital_rates[4, 2],
                          F_Fledgling_survival = coucal_vital_rates[3, 2],
                          F_Adult_survival = coucal_vital_rates[2, 2],
                          F_Adult_immigration = coucal_vital_rates[1, 2],
                          M_Nestling_survival = coucal_vital_rates[10, 2],
                          M_Groundling_survival = coucal_vital_rates[9, 2],
                          M_Fledgling_survival = coucal_vital_rates[8, 2],
                          M_Adult_survival = coucal_vital_rates[7, 2],
                          M_Adult_immigration = coucal_vital_rates[6, 2],
                          
                          # Define hatching sex ratio
                          HSR = HSR,
                          
                          # Define the mating system (h), and clutch size (k)
                          h = h,
                          k = k)

# Build matrix based on rates specified in the list above
demographic_matrix <- coucal_matrix(demographic_rates)

# populate sex-specific adult immigration matrix
immigration_vector <- matrix(data = c(0, demographic_rates$F_Adult_immigration, 
                                      0, demographic_rates$M_Adult_immigration), 
                             nrow = 4, ncol = 1)

# Determine the ASR at the stable stage distribution
ASR_SSD <- matrix_ASR(M = demographic_matrix,
                      h = demographic_rates$h,
                      HSR = demographic_rates$HSR, iterations = 100,
                      num_boot = num_boot,
                      species = species,
                      immigrant_pop_size = immigrant_pop_size,
                      ISR = ISR, iter_add = 1)

# Extract ASR
ASR_estimate <- ASR_SSD$ASR

# make a list of all the results from this iteration
bootstrap_results_list <- 
  list(#offspring_hazard_function = haz_ss_function, 
    offspring_hazard_rates = haz_ss_curve,
    coucal_vital_rates = coucal_vital_rates,
    ASR_SSD = ASR_SSD,
    bootstrapped_data = coucal_boot_list)

# run the survival analysis and ASR deduction on the sampled data
result <- bootstrap_hazard_ASR(coucal_boot_list = bootstrap_data_list, 
                               num_boot = num_boot, 
                               egg_survival = egg_survival, 
                               ISR = ISR, 
                               immigrant_pop_size = immigrant_pop_size, 
                               HSR = HSR, h = h, k = k, 
                               flight_age = flight_age, 
                               bootstrap_name = bootstrap_name,
                               species = species,
                               iter_add = iter_add,
                               prefix_number = prefix_number)
result

three <- bootstrap_results_list
two <- bootstrap_results_list
one <- bootstrap_results_list
load("output/bootstraps/hazard/raw/BC_boot_one_2.Rds")
bootstrap_results_list$offspring_hazard_function$Male_haz_ss$


lengths(two)
lengths(three)

test <- two$offspring_hazard_function$M_haz_ss["terms"] <- NULL

project.sshzd(bootstrap_results_list$offspring_hazard_function$Male_haz_ss)

length(bootstrap_results_list$offspring_hazard_function$M_haz_ss)
length(two$offspring_hazard_function$M_haz_ss$quad$pt)
length(one$offspring_hazard_function$M_haz_ss$se.aux)

three$offspring_hazard_function$M_haz_ss$
bootstrap_results_list$offspring_hazard_rates
bootstrap_results_list$coucal_vital_rates
bootstrap_results_list$ASR_SSD$stage.vectors

#########
num_boot = 1
offspring = WBC_all_dat
k = 4
HSR = 0.5198
h = 1/1.1
egg_survival = 0.18
ISR = 0.524
adult_survival_rate = 0.3
immigrant_pop_size = 100
fledge_age = 15
flight_age = 32
bootstrap_name = "WBC_boot_one"
species = "WBC"
iter_add = 1
prefix_number = "boot_one_"
time_vector = seq(0, 70, 1)

rm(hzd_curve_try_M, hzd_curve_try_F)

# set.seed(2)
# set attempt to 0 at start of each loop
attempt <- 0

max_time <- 70

# sample a new offspring dataset containing only one nest 
# member per draw
offspring_boot <- 
  offspring %>% 
  group_by(nest_ID) %>% 
  sample_n(1)

# store simulated estimates only if peak >= 1 and <= 10 and it's less than
# 100 attempts
while( (exists("hzd_curve_try_M") == FALSE | exists("hzd_curve_try_F") == FALSE) && attempt <= max_time) {
  
  time_vector <- seq(0, max_time - attempt, 1)
  
  # next attempt
  attempt <- attempt + 1
  
  try(
    hzd_ss_try_M <- sshzd(Surv(exit, event, entry) ~ exit, 
                          data = filter(offspring_boot, sex == "M"), 
                          alpha = alpha_value)
  )
  
  try(
    hzd_ss_try_F <- sshzd(Surv(exit, event, entry) ~ exit, 
                          data = filter(offspring_boot, sex == "F"), 
                          alpha = alpha_value)
  )
  
  # simulate an estimate
  try(
    hzd_curve_try_M <- 
      hzdcurve.sshzd(object = hzd_ss_try_M, time = time_vector, se = TRUE)
  )
  try(
    hzd_curve_try_F <- 
      hzdcurve.sshzd(object = hzd_ss_try_F, time = time_vector, se = TRUE)
  )
}

time_vector

# make a list of these two datasets, which will be used in the next function
out <- list(offspring_boot = offspring_boot, 
            iter = num_boot + ((iter_add - 1) * niter),
            species = species,
            time_vector = time_vector)
###########
detect_dat <- 
  read_xlsx("data/raw/All_coucal_waypoints_2001_2019_20200202.xlsx", na = "NA", col_types = "text") %>% 
  dplyr::select(species, ring_ID, sex, year, site, age_status, date_dec) %>% 
  mutate(month = str_sub(date_dec, start = 5, end = 6),
         day = str_sub(date_dec, start = 7, end = 8)) %>% 
  mutate(date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d")) %>% 
  dplyr::select(-month, -day, -date_dec) %>% 
  mutate(across(c("sex", "site", "age_status"), tolower)) %>%
  mutate(sex = ifelse(sex == "Female", "F", ifelse(sex == "Male", "M", "XXX")),
         age_status = ifelse(age_status == "adult", "A", ifelse(age_status == "juvenile", "J", "XXX")),
         ring_ID = str_replace_all(string = ring_ID, fixed(" "), "")) %>% 
  mutate(across(c("species", "ring_ID", "sex", "site", "age_status"), as.factor)) %>% 
  mutate(sex = ifelse(str_detect(ring_ID, pattern = "[Ff]emale"), "F",
                      ifelse(str_detect(ring_ID, pattern = "[Mm]ale"), "M", sex)))

detect_dat %>% 
  dplyr::select(-date) %>% 
  filter(site == "Kapunga" & species != "CTC") %>%
  distinct() %>% 
  group_by(species, year, site) %>% 
  dplyr::summarise(n_ind = n_distinct(ring_ID)) %>% 
  filter(n_ind > 20)


detect_dat %>% 
  dplyr::select(-date) %>% 
  distinct() %>% 
  group_by(species, year, site, sex) %>% 
  dplyr::summarise(n_ind = n_distinct(ring_ID)) %>% 
  pivot_wider(names_from = sex, values_from = n_ind) %>% 
  mutate(ISR = M / (F + M),
         total = F + M) %>% 
  filter(!is.na(ISR) & site == "Kapunga" & species != "CTC" & total > 20) %>%
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
                   mean_f = mean(F, na.rm = TRUE),
                   mean_M = mean(M, na.rm = TRUE),
                   min_f = min(F, na.rm = TRUE),
                   max_f = max(F, na.rm = TRUE))

################
ISR_mod <- 
  glm(cbind(n_males, n_females) ~ species, 
      family = binomial,
      data = ISR_data)

invlogit(model_parameters(ISR_mod)$Coefficient[1])
invlogit(model_parameters(ISR_mod)$CI_low[1])
invlogit(model_parameters(ISR_mod)$CI_high[1])

coucal_HSR <- 
  data.frame(species = c("BC", "WBC"),
             mean_HSR = c(invlogit(model_parameters(ISR_mod)$Coefficient),
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
  mutate(sd_est = ifelse(!is.na(CI_low), 
                         approx_sd(x1 = CI_low, x2 = CI_high),
                         CI_low))
ISR_data %>% 
  group_by(species) %>% 
  dplyr::summarise(mean_annual_ISR = mean(ISR, na.rm = TRUE),
                   var_annual_ISR = var(ISR, na.rm = TRUE),
                   median_annual_ISR = median(ISR, na.rm = TRUE),
                   sd_annual_ISR = sd(ISR, na.rm = TRUE),
                   lcl_annual_ISR = stats::quantile(ISR, (1 - CI)/2, na.rm = TRUE),
                   ucl_annual_ISR = stats::quantile(ISR, 1 - (1 - CI)/2, na.rm = TRUE),
                   max_annual_ISR = max(ISR),
                   min_annual_ISR = min(ISR),
                   n_years = n_distinct(year),
                   mean_f = mean(n_females, na.rm = TRUE),
                   mean_M = mean(n_males, na.rm = TRUE),
                   min_f = min(n_females, na.rm = TRUE),
                   max_f = max(n_females, na.rm = TRUE))

mating_dat %>% 
  filter((sex == "female" | (sex == "male" & `ID_partner.s.` == "unringedf"))  & site == "kapunga") %>% 
  group_by(species, year, sex) %>% 
  dplyr::summarise(n_partners = sum(Nr_partners),
                   n_focals = n_distinct(ind_ID)) %>% 
  mutate(ISR = n_males / (n_females + n_males),
         total = n_females + n_males) %>% 
  filter(!is.na(ISR)) %>%
  arrange(desc(total)) %>% 
  group_by(species) %>% 
  dplyr::summarise(mean_annual_ISR = mean(ISR, na.rm = TRUE),
                   var_annual_ISR = var(ISR, na.rm = TRUE),
                   median_annual_ISR = median(ISR, na.rm = TRUE),
                   sd_annual_ISR = sd(ISR, na.rm = TRUE),
                   lcl_annual_ISR = stats::quantile(ISR, (1 - CI)/2, na.rm = TRUE),
                   ucl_annual_ISR = stats::quantile(ISR, 1 - (1 - CI)/2, na.rm = TRUE),
                   max_annual_ISR = max(ISR),
                   min_annual_ISR = min(ISR),
                   n_years = n_distinct(year),
                   mean_f = mean(n_females, na.rm = TRUE),
                   mean_M = mean(n_males, na.rm = TRUE),
                   min_f = min(n_females, na.rm = TRUE),
                   max_f = max(n_females, na.rm = TRUE))

mating_dat %>% 
  filter(year == 2001) %>% 
  filter(ind_ID %in% c("GN41707", "GN41713", "GN41703"))

mating_dat %>% 
  filter((sex == "female" & `ID_partner.s.` == "unringedf") | (sex == "male" & `ID_partner.s.` == "unringedm"))

mating_dat %>% 
  filter((sex == "male" & `ID_partner.s.` == "unringedf"))

mating_dat %>% 
  group_by(species, year, site, sex) %>% 
  dplyr::summarise(n_partners = sum(Nr_partners)) %>% 
  pivot_wider(names_from = sex, values_from = n_partners) %>% 
  dplyr::rename(male_partners = female,
                female_partners = male) %>% 
  mutate(ISR = male_partners / (female_partners + male_partners),
         total = female_partners + male_partners) %>% 
  filter(!is.na(ISR) & site == "kapunga") %>%
  arrange(desc(total))
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
                   mean_f = mean(male_partners, na.rm = TRUE),
                   mean_M = mean(female_partners, na.rm = TRUE),
                   min_f = min(male_partners, na.rm = TRUE),
                   max_f = max(male_partners, na.rm = TRUE))

mating_dat %>% 
  filter(species == "BC", year == 2014, site == "kapunga")

mating_dat %>% 
  group_by(species, year, site, sex) %>% 
  dplyr::summarise(n_partners = sum(Nr_partners)) %>% 
  pivot_wider(names_from = sex, values_from = n_partners) %>% 
  dplyr::rename(male_partners = female,
                female_partners = male) %>% 
  mutate(ISR = male_partners / (female_partners + male_partners),
         total = female_partners + male_partners) %>% 
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
                   mean_f = mean(male_partners, na.rm = TRUE),
                   mean_M = mean(female_partners, na.rm = TRUE),
                   min_f = min(male_partners, na.rm = TRUE),
                   max_f = max(male_partners, na.rm = TRUE))

mating_dat %>% 
  group_by(species, year, site, sex) %>% 
  dplyr::summarise(n_ind = n_distinct(ind_ID)) %>% 
  pivot_wider(names_from = sex, values_from = n_ind) %>% 
  mutate(ISR = male / (female + male),
         total = female + male) %>% 
  filter(!is.na(ISR) & site == "kapunga" & total > 15) %>%
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


######
#### one step matrix model ----

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
WBC_hazard_rate_boot <- 
  readRDS("output/bootstraps/hazard/cooked/WBC_hazard_ASR_bootstrap_result_one.rds")

# clean up the output from the bootstrap procedure and save as rds
WBC_hazard_rate_boot_tidy <- 
  hazard_boot_out_wrangle(species = "WBC", niter = 1000, 
                          output_dir = "output/bootstraps/hazard/cooked/",
                          rds_file = "_hazard_ASR_bootstrap_result_one")

average_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  group_by(sex, age) %>% 
  dplyr::summarise(med_surv = median(estimate))

F_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Female")

M_rates <- 
  WBC_hazard_rate_boot_tidy$hazard_rates_boot %>% 
  filter(sex == "Male")

average_rates
ages <- seq(0, max(average_rates$age), 1)

juv_pop_storage <- 
  matrix(numeric(2 * length(ages)),
         ncol = 2, nrow = 70, 
         dimnames = list(ages, c("nF", "nM")))

F_juv_pop_storage_boot <- 
  matrix(, nrow = 70, ncol = niter)

M_juv_pop_storage_boot <- 
  matrix(, nrow = 70, ncol = niter)

n_eggs = 100
k = 4
HSR = 0.5198
h = 1/1.1
WBC_egg_survival = 0.18
adult_survival_rate = 0.3
ISR = 0.524
immigrant_pop_size = 1000
fledge_age = 15
flight_age = 32

A <- matrix(sapply(matrix_str, eval, vr, NULL), 
            nrow = sqrt(length(matrix_str)), byrow=TRUE, 
            dimnames = list(stages, stages))

N_F_immigrants <- 
  immigrant_pop_size * (1 - ISR)

N_M_immigrants <- 
  immigrant_pop_size * ISR

# # reset the starting stage distribution for simulation (all with 10 individuals)
# m <- rep(10, no_stages) 
# 
# # number of male adults at time t
# M2 <- stage[4, i] * A["M_Adult", "M_Adult"] + I[4, 1]
# 
# # number of female adults at time t
# F2 <- stage[2, i] * A["F_Adult", "F_Adult"] + I[2, 1]

# Female freq-dep fecundity of Female chicks
F_F_eggs <- 
  ((k * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * (1 - HSR))

# Female freq-dep fecundity of Male chicks
F_M_eggs <- 
  ((k * N_M_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * HSR)

# Male freq-dep fecundity of Female chicks
M_F_eggs <- 
  ((k * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * (1 - HSR))

# Male freq-dep fecundity of Male chicks
M_M_eggs <- 
  ((k * N_F_immigrants) / (N_M_immigrants + (N_F_immigrants / h)) * HSR)

N_F_eggs <- F_F_eggs * N_F_immigrants + M_F_eggs * N_M_immigrants
N_M_eggs <- F_M_eggs * N_F_immigrants + M_M_eggs * N_M_immigrants
juv_pop_storage[1, "nF"] <- N_F_eggs * egg_survival
juv_pop_storage[1, "nM"] <- N_M_eggs * egg_survival

for(i in 2:(max(average_rates$age) + 1)){
  juv_pop_storage[i, "nF"] <- 
    juv_pop_storage[i - 1, "nF"] * average_rates[which(average_rates$sex == "Female"), ][["med_surv"]][i - 1]
  juv_pop_storage[i, "nM"] <- 
    juv_pop_storage[i - 1, "nM"] * average_rates[which(average_rates$sex == "Male"), ][["med_surv"]][i - 1]
}

for(j in 1:niter){
  F_juv_pop_storage_boot[1, j] <- N_F_eggs * egg_survival
  M_juv_pop_storage_boot[1, j] <- N_M_eggs * egg_survival
  for(i in 2:70){
    F_juv_pop_storage_boot[i, j] <- 
      F_juv_pop_storage_boot[i - 1, j] * F_rates[which(F_rates$iter == j), ][["estimate"]][i - 1]
    M_juv_pop_storage_boot[i, j] <- 
      M_juv_pop_storage_boot[i - 1, j] * M_rates[which(M_rates$iter == j), ][["estimate"]][i - 1]
  }
}

juv_pop_storage_clean <- 
  juv_pop_storage %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "age") %>% 
  pivot_longer(!age, names_to = "sex", values_to = "n")

F_juv_pop_storage_clean <- 
  data.frame(F_juv_pop_storage_boot) %>% 
  mutate(age = c(0:69)) %>% 
  # as.data.frame() %>% 
  # rownames_to_column(var = "age") %>% 
  pivot_longer(!age, names_to = "iteration", values_to = "n") %>% 
  mutate(sex = "F") %>% 
  arrange(iteration)

M_juv_pop_storage_clean <- 
  data.frame(M_juv_pop_storage_boot) %>% 
  mutate(age = c(0:69)) %>% 
  # as.data.frame() %>% 
  # rownames_to_column(var = "age") %>% 
  pivot_longer(!age, names_to = "iteration", values_to = "n") %>% 
  mutate(sex = "M") %>% 
  arrange(iteration)

ggplot(data = bind_rows(F_juv_pop_storage_clean,
                        M_juv_pop_storage_clean)) +
  geom_line(aes(y = n, x = as.numeric(age), group = iteration)) +
  facet_grid(sex ~ .)

JSR <- 
  juv_pop_storage %>% 
  as.data.frame() %>% 
  mutate(JSR = nM / (nM + nF)) %>% 
  rownames_to_column(var = "age")

JSR_boot_WBC <- 
  left_join(F_juv_pop_storage_clean, M_juv_pop_storage_clean, 
            by = c("age", "iteration"))  %>% 
  mutate(JSR = n.y / (n.y + n.x))

filter(JSR_boot_WBC, age == 69) %>% 
  summary(mean(JSR))

filter(JSR_boot, age == 69) %>% 
  summary(mean(JSR))

ggplot(data = JSR_boot) +
  geom_line(aes(y = JSR, x = as.numeric(age), group = iteration),
            alpha = 0.05) +
  geom_hline(yintercept = 0.5, color = "white")


#############
#### estimated from the chick data ----
chick_data <-
  read.csv("data/raw/Chicks_BC_WBC_sex_alloc_data_all_20200925.csv") %>% 
  dplyr::filter(Clutch_size == Brood_size) %>% 
  dplyr::select(Spp, Year, Nest_ID, Clutch_size, Brood_size, Status, Sex, 
                Prob_sex_male, Prob_sex_female, Absolute_hatchingorder, EPY) %>% 
  dplyr::group_by(Nest_ID) %>% 
  dplyr::mutate(n_ind_sexed = n(),
                n_males = sum(Prob_sex_male),
                n_females = sum(Prob_sex_female)) %>% 
  dplyr::filter(n_ind_sexed == Clutch_size)

BC_HSR_2 <- 
  lme4::glmer(cbind(Prob_sex_male, Prob_sex_female) ~ 
                1 + 
                (1| Nest_ID), 
              data = filter(chick_data, Spp == "BC"), 
              family = binomial)

WBC_HSR_2 <- 
  lme4::glmer(cbind(Prob_sex_male, Prob_sex_female) ~ 
                1 + 
                (1| Nest_ID), 
              data = filter(chick_data, Spp == "WBC"), 
              family = binomial)

coucal_HSR_2 <- 
  data.frame(species = c("BC", "WBC"),
             mean_HSR = c(invlogit(model_parameters(BC_HSR_2)$Coefficient),
                          invlogit(model_parameters(WBC_HSR_2)$Coefficient)),
             CI_low = c(invlogit(model_parameters(BC_HSR_2)$CI_low),
                        invlogit(model_parameters(WBC_HSR_2)$CI_low)),
             CI_high = c(invlogit(model_parameters(BC_HSR_2)$CI_high),
                         invlogit(model_parameters(WBC_HSR_2)$CI_high)),
             n_nests = c(filter(chick_data, Spp == "BC") %>% 
                           summarise(n_nests = n_distinct(Nest_ID)) %>% 
                           pull(n_nests),
                         filter(chick_data, Spp == "WBC") %>% 
                           summarise(n_nests = n_distinct(Nest_ID)) %>% 
                           pull(n_nests)),
             n_years = c(filter(chick_data, Spp == "BC") %>% 
                           summarise(n_nests = n_distinct(Year)) %>% 
                           pull(n_nests),
                         filter(chick_data, Spp == "WBC") %>% 
                           summarise(n_nests = n_distinct(Year)) %>% 
                           pull(n_nests))) %>% 
  mutate(sd_est = ifelse(!is.na(CI_low), 
                         approx_sd(x1 = CI_low, x2 = CI_high),
                         CI_low))