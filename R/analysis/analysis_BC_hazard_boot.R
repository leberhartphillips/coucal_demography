# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
sapply(function.sources, source)

# load capture histories
data.sources = list.files(path = "cooked_data", 
                          pattern="*ch.rds$", full.names = TRUE, 
                          ignore.case = TRUE)
sapply(data.sources, load, .GlobalEnv)

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
  dplyr::select(species, ind_ID, nest_ID, year, sex, entry, exit, event, hatch_order)

BC_all_dat <- 
  status_dat_all %>% 
  filter(species == "BC")

#### Bootstrap run ----
niter = 1000
set.seed(14)

# run bootstrap procedure on Black Coucals
BC_hazard_ASR_bootstrap_result_one <-
  pbsapply(1:niter, run_bootstrap_hazard_ASR,
           offspring = BC_all_dat,
           k = 4,
           HSR = 0.4955,
           h = 1/2.9,
           egg_survival = 0.32,
           adult_survival_rate = 0.3,
           ISR = 0.738,
           immigrant_pop_size = 100,
           fledge_age = 15,
           flight_age = 36,
           bootstrap_name = "BC_boot_one",
           species = "BC",
           iter_add = 1,
           prefix_number = "boot_one_",
           time_vector = seq(0, 70, 1))

# save model output
saveRDS(object = BC_hazard_ASR_bootstrap_result_one, 
        file = "output/bootstraps/hazard/cooked/BC_hazard_ASR_bootstrap_result_one.rds")