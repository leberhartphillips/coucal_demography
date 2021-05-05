source("R/project/project_libraries.R")
# source("R/scripts/02_func_rearrange_ch_for_age.R")

# import and consolidate data
detect_dat <- 
  read_xlsx("data/raw/All_coucal_waypoints_2001_2019_20200202.xlsx", 
            na = "NA", col_types = "text") %>% 
  select(species, ring_ID, sex, year, site, age_status, date_dec, age_days) %>% 
  mutate(month = str_sub(date_dec, start = 5, end = 6),
         day = str_sub(date_dec, start = 7, end = 8)) %>% 
  mutate(date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d")) %>% 
  select(-month, -day, -date_dec) %>% 
  mutate(across(c("sex", "site", "age_status"), tolower)) %>%
  mutate(sex = ifelse(sex == "female", "F", ifelse(sex == "male", "M", "XXX")),
         age_status = ifelse(age_status == "adult", "A", 
                             ifelse(age_status == "juvenile", "J", "XXX")),
         ring_ID = str_replace_all(string = ring_ID, fixed(" "), "")) %>% 
  mutate(across(c("species", "ring_ID", "sex", "site", "age_status"), as.factor),
         age_days = as.numeric(age_days))

# make daily capture histories for Juveniles for each year radio-tracking was done
detect_dat_daily_2014 <-
  detect_dat %>% 
  filter(age_status == "J" & year == 2014) %>% 
  na.omit() %>% 
  mutate(date2 = as.character(date))

detect_dat_daily_2014_ch <- 
  CensusToCaptHist(ID = detect_dat_daily_2014$ring_ID,
                   d = detect_dat_daily_2014$date, 
                   timeInt = "D",
                   dformat = "yyyy-mm-dd") %>% 
  mutate(ring_ID = rownames(.),
         ID = as.character(ID)) %>% 
  left_join(., select(detect_dat, species, ring_ID, sex, site, age_status, year), 
            by = "ring_ID") %>% 
  distinct()

detect_dat_daily_2015 <-
  detect_dat %>% 
  filter(age_status == "J" & year == 2015) %>% 
  na.omit() %>% 
  mutate(date2 = as.character(date))

detect_dat_daily_2015_ch <- 
  CensusToCaptHist(ID = detect_dat_daily_2015$ring_ID,
                   d = detect_dat_daily_2015$date, 
                   timeInt = "D",
                   dformat = "yyyy-mm-dd") %>% 
  mutate(ring_ID = rownames(.),
         ID = as.character(ID)) %>% 
  left_join(., select(detect_dat, species, ring_ID, sex, site, age_status, year), 
            by = "ring_ID") %>% 
  distinct()

detect_dat_daily_2016 <-
  detect_dat %>% 
  filter(age_status == "J" & year == 2016) %>% 
  na.omit() %>% 
  mutate(date2 = as.character(date))

detect_dat_daily_2016_ch <- 
  CensusToCaptHist(ID = detect_dat_daily_2016$ring_ID,
                   d = detect_dat_daily_2016$date, 
                   timeInt = "D",
                   dformat = "yyyy-mm-dd") %>% 
  mutate(ring_ID = rownames(.),
         ID = as.character(ID)) %>% 
  left_join(., select(detect_dat, species, ring_ID, sex, site, age_status, year), 
            by = "ring_ID") %>% 
  distinct()

detect_dat_daily_2017 <-
  detect_dat %>% 
  filter(age_status == "J" & year == 2017) %>% 
  na.omit() %>% 
  mutate(date2 = as.character(date))

detect_dat_daily_2017_ch <- 
  CensusToCaptHist(ID = detect_dat_daily_2017$ring_ID,
                   d = detect_dat_daily_2017$date, 
                   timeInt = "D",
                   dformat = "yyyy-mm-dd") %>% 
  mutate(ring_ID = rownames(.),
         ID = as.character(ID)) %>% 
  left_join(., select(detect_dat, species, ring_ID, sex, site, age_status, year), 
            by = "ring_ID") %>% 
  distinct()

detect_dat_daily_2018 <-
  detect_dat %>% 
  filter(age_status == "J" & year == 2018) %>% 
  na.omit() %>% 
  mutate(date2 = as.character(date))

detect_dat_daily_2018_ch <- 
  CensusToCaptHist(ID = detect_dat_daily_2018$ring_ID,
                   d = detect_dat_daily_2018$date, 
                   timeInt = "D",
                   dformat = "yyyy-mm-dd") %>% 
  mutate(ring_ID = rownames(.),
         ID = as.character(ID)) %>% 
  left_join(., select(detect_dat, species, ring_ID, sex, site, age_status, year), 
            by = "ring_ID") %>% 
  distinct()

# bind all annual capture histories together
detect_dat_daily_ch <- 
  bind_rows(detect_dat_daily_2015_ch,
            detect_dat_daily_2014_ch,
            detect_dat_daily_2016_ch,
            detect_dat_daily_2017_ch,
            detect_dat_daily_2018_ch)

# determine the age at first and last observation for each individual
age_first_obs <- 
  detect_dat %>% 
  group_by(ring_ID) %>% 
  summarise(age_first_obs = min(age_days),
            age_last_obs = max(age_days))

# create a two-character string for each encounter and clean the output
detect_dat_daily_ch <- 
  sapply(detect_dat_daily_ch[2:139], function(x) paste(x, "0", sep = "")) %>% 
  cbind(detect_dat_daily_ch[c(1,140:length(detect_dat_daily_ch))]) %>% 
  mutate(across(everything(), ~str_replace(string = .x, pattern = "NA0", "00"))) %>% 
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x)))

# Import status data
status_dat_chicks <- 
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::rename(ring_ID = Ring_ID) %>% 
  mutate(Fledged_status = tolower(Fledged.),
         site = tolower(site)) %>% 
  select(species, ring_ID, sex, year, site, nest_ID, pref_age, Fledged_status, postf_age, postf_status, lay_date) %>% 
  mutate(across(everything(), ~str_trim(.x))) %>% 
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>% 
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  mutate(age_status = "J")

# join detection history with status data
detect_status_join <-
  left_join(detect_dat_daily_ch, status_dat_chicks, by = "ring_ID") %>% 
  filter(Fledged_status == "yes") %>% 
  filter(year.x == year.y) %>%
  filter(sex.x == sex.y) %>%
  filter(species.x == species.y) %>% 
  filter(!is.na(postf_status)) %>% 
  select(-species.y, -sex.y, -site.y, -year.y, -age_status.y) %>% 
  rename(species = species.x,
         sex = sex.x,
         site = site.x,
         age_status = age_status.x,
         year = year.x)

# determine the last and first detection for each individual 
max_age_index <- 
  apply(detect_status_join[, c(1:138)], 1, function(x) which(x == "10")) %>% 
  lapply(., function(x) x[which.max(x)]) %>% 
  unlist(.)

min_age_index <- 
  apply(detect_status_join[, c(1:138)], 1, function(x) which(x == "10")) %>% 
  lapply(., function(x) x[which.min(x)]) %>% 
  unlist(.)

# put "11" at the last detection if the status was a "1" (i.e., dead)
for(i in 1:nrow(detect_status_join)){
  if(detect_status_join$postf_status[i] == "1"){
    detect_status_join[i, max_age_index[i]] <- "11"
  }
}

# bind the min and max age info to the capture history data
detect_status <-
  bind_cols(detect_status_join, max_age_index) %>% 
  bind_cols(., min_age_index) %>% 
  rename(date_last_obs = '...152',
         date_first_obs = '...153') %>% 
  left_join(., age_first_obs, by = "ring_ID") %>% 
  arrange((age_first_obs)) %>% 
  mutate(hatch_date = as.numeric(lay_date) + as.numeric(age_first_obs)) %>% 
  drop_na()

# consolidate capture histories for each sex
Black_Coucal_fledgling_Burnham_ch <-
  apply(detect_status[, c(1:138)], 1, paste, collapse = "") %>% 
  as.data.frame() %>% 
  rename(ch = '.') %>% 
  mutate(ch = as.character(ch)) %>% 
  bind_cols(., detect_status[, c(139:length(detect_status))]) %>% 
  filter(species == "BC") %>%
  rearrange_ch_for_age(df = ., ch_col = 1, 
                       age_last_obs_col = 18,
                       age_first_obs_col = 17, 
                       max_age = 70,
                       min_age = 15) %>% 
  filter(str_detect(ch, "11", negate = TRUE) | str_count(ch, "1") > 2) %>% 
  select(everything(), -c(ID, age_status, pref_age, Fledged_status,
                          postf_age, postf_status, lay_date, date_last_obs,
                          date_first_obs, age_first_obs, age_last_obs,
                          hatch_date, species)) %>%
  drop_na() %>%
  mutate_at(vars(ring_ID, sex, site, year, nest_ID), as.factor)

White_browed_Coucal_fledgling_Burnham_ch <-
  apply(detect_status[, c(1:138)], 1, paste, collapse = "") %>% 
  as.data.frame() %>% 
  rename(ch = '.') %>% 
  mutate(ch = as.character(ch)) %>% 
  bind_cols(., detect_status[, c(139:length(detect_status))]) %>% 
  filter(species == "WBC") %>%
  rearrange_ch_for_age(df = ., ch_col = 1, 
                       age_last_obs_col = 18,
                       age_first_obs_col = 17, 
                       max_age = 70,
                       min_age = 15) %>% 
  filter(str_detect(ch, "11", negate = TRUE) | str_count(ch, "1") > 2) %>% 
  select(everything(), -c(ID, age_status, pref_age, Fledged_status,
                          postf_age, postf_status, lay_date, date_last_obs,
                          date_first_obs, age_first_obs, age_last_obs,
                          hatch_date, species)) %>%
  drop_na() %>%
  mutate_at(vars(ring_ID, sex, site, year, nest_ID), as.factor)

# save cleaned capture history data
save(Black_Coucal_fledgling_Burnham_ch, 
     file = "cooked_data/Black_Coucal_fledgling_Burnham_ch.rds")
save(White_browed_Coucal_fledgling_Burnham_ch, 
     file = "cooked_data/White-browed_Coucal_fledgling_Burnham_ch.rds")