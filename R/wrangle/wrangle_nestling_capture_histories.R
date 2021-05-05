source("scripts/01_libraries.R")
source("scripts/02_func_known_fate_df_convert().R")

# import raw csv data
dat <-  
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", 
           header = TRUE, na.strings = c("", " ", "NA")) %>% 
  mutate(sex = as.factor(sex),
         ageC = as.numeric(ageC)) %>% 
  filter(!is.na(sex) & !is.na(ageC)) %>% 
  mutate(Final_status = tolower(Final_status),
         Fledged_status = tolower(Fledged.),
         ring_ID = str_replace_all(string = Ring_ID, fixed(" "), ""))

# create species-specific dataframes
BC_dat <- 
  dat %>% 
  filter(species == "BC")

WBC_dat <- 
  dat %>% 
  filter(species == "WBC")

# convert the raw data to the known-fate format for nestlings
Black_Coucal_nestling_Known_ch <- 
  known_fate_df_convert(df = filter(BC_dat, !is.na(Fledged_status)), 
                        interval_col = 5, status_col = 44, 
                        n_intervals = 15, stage = "nestling") %>% 
  select(ch, sex, year, nest_ID, lab_no, ring_ID, 
         lay_date, clutch_size, brood_size, hatch_order, age_diff) %>% 
  rename(clutch_s = clutch_size,
         brood_s = brood_size,
         hatch_ord = hatch_order)

White_browed_Coucal_nestling_Known_ch <- 
  known_fate_df_convert(df = filter(WBC_dat, !is.na(Fledged_status)), 
                        interval_col = 5, status_col = 44, 
                        n_intervals = 15, stage = "nestling") %>% 
  select(ch, sex, year, nest_ID, lab_no, ring_ID,
         lay_date, clutch_size, brood_size, hatch_order, age_diff) %>% 
  rename(clutch_s = clutch_size,
         brood_s = brood_size,
         hatch_ord = hatch_order)

# save cleaned capture history data
save(Black_Coucal_nestling_Known_ch, 
     file = "cooked_data/Black_Coucal_nestling_Known_ch.rds")
save(White_browed_Coucal_nestling_Known_ch, 
     file = "cooked_data/White-browed_Coucal_nestling_Known_ch.rds")