source("scripts/01_libraries.R")

# import raw csv data into R
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
  distinct()
  

# assess the first and last year of study for each species
detect_dat %>% 
  group_by(species) %>% 
  summarise(first_year = min(year),
            last_year = max(year))

# make annual capture histories for adults
BC_detect_dat_A <-
  detect_dat %>% 
  filter(age_status == "A" & species == "BC")

# use the BaSTA function "CensusToCaptHist()" function to convert long format
# encounter histories of each individual, to wide format with 1's and 0's for 
# each year of encounter
BC_detect_dat_A_ch <- 
  CensusToCaptHist(ID = BC_detect_dat_A$ring_ID,
                   d = BC_detect_dat_A$year) %>% 
  mutate(ring_ID = rownames(.),
         ID = as.character(ID)) %>% 
  left_join(., select(BC_detect_dat_A, ring_ID, sex, year), by = "ring_ID") %>% 
  distinct()

Black_Coucal_adult_CJS_ch <- 
  data.frame(ch = apply(BC_detect_dat_A_ch[, 2:19 ] , 1, paste, collapse = "")) %>% 
  bind_cols(., select(BC_detect_dat_A_ch, sex, ring_ID)) %>% 
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  distinct()

WBC_detect_dat_A <-
  detect_dat %>% 
  filter(age_status == "A" & species == "WBC")

WBC_detect_dat_A_ch <- 
  CensusToCaptHist(ID = WBC_detect_dat_A$ring_ID,
                   d = WBC_detect_dat_A$year) %>% 
  mutate(ring_ID = rownames(.),
         ID = as.character(ID)) %>% 
  left_join(., select(WBC_detect_dat_A, ring_ID, sex, year), by = "ring_ID") %>% 
  distinct()

White_browed_Coucal_adult_CJS_ch <- 
  data.frame(ch = apply(WBC_detect_dat_A_ch[, 2:16 ] , 1, paste, collapse = "")) %>% 
  bind_cols(., select(WBC_detect_dat_A_ch, sex, ring_ID)) %>% 
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  distinct()

# save cleaned capture history data
save(Black_Coucal_adult_CJS_ch, 
     file = "cooked_data/Black_Coucal_adult_CJS_ch.rds")
save(White_browed_Coucal_adult_CJS_ch, 
     file = "cooked_data/White-browed_Coucal_adult_CJS_ch.rds")
