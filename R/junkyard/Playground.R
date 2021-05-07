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