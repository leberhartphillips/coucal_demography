# load libraries
source("R/project/project_libraries.R")
source("R/project/project_plotting.R")

# load functions
function.sources = list.files(path = "R/functions", 
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

# age_plumage_dist_data <-
  read_xlsx("data/raw/age_classes_2001-2020_2.xlsx", na = "NA", col_types = "text") %>%
  dplyr::rename_all(~str_replace_all(., "\\s+", "")) %>%
  dplyr::mutate(species = tolower(species)) %>%
  dplyr::filter(str_detect(string = species, pattern = "black")) %>% 
  dplyr::mutate(month = str_sub(`date(yymmdd)`, start = 5, end = 6),
         day = str_sub(`date(yymmdd)`, start = 7, end = 8)) %>% 
  dplyr::mutate(date = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d")) %>% 
  # remove all white space from data
  dplyr::mutate(dplyr::across(.cols = dplyr::everything(), 
                str_remove_all, pattern = fixed(" "))) %>% 
  dplyr::mutate(species = "BC",
                sex = ifelse(sex == "female", "F", ifelse(sex == "male", "M", "XXX"))) %>% 
  dplyr::mutate(dplyr::across(.cols = dplyr::everything(), 
                              str_replace_all, pattern = fixed("onebarred"), replacement = fixed("barred"))) %>%
    
  dplyr::mutate(Secondarycoverts = ifelse(SN == "462", "rufous", Secondarycoverts)) %>% 
  dplyr::mutate(Secondarycoverts = ifelse(SN == "409", "barred", Secondarycoverts)) %>% 
    
  dplyr::mutate(age_LEP = ifelse(str_detect(primarycoverts, "barred") & 
                                   str_detect(primaries, "barred") & 
                                   str_detect(Secondarycoverts, "barred") & 
                                   str_detect(secondaries, "barred"), 1,
                                 
                                 ifelse(str_detect(primaries, "barred") | 
                                        str_detect(secondaries, "barred"), 2,
                                        
                                        ifelse((str_detect(primarycoverts, "barred") |
                                               str_detect(Secondarycoverts, "barred")) &
                                               (str_detect(primaries, "rufous") |
                                               str_detect(secondaries, "rufous")), 2,
                                               
                                               ifelse(str_detect(primarycoverts, "rufous") & 
                                                         str_detect(primaries, "rufous") & 
                                                         str_detect(Secondarycoverts, "rufous") & 
                                                         str_detect(secondaries, "rufous"), 3, 0))))) %>%
  dplyr::mutate(age = ifelse(age == "juvenile", 1, age)) %>%
  dplyr::mutate(CHECK = ifelse(age != age_LEP, "CHECK", "")) %>% 
  # dplyr::filter(age != age_LEP) %>%
  # dplyr::select(SN, year, species, Alu, sex, date, primarycoverts, primaries,
  #               Secondarycoverts, secondaries, tail, legcoverts, face, age, age_LEP) #%>%
  write.csv(file = "data/raw/age_classes_2001-2020_LEP3.csv")
  # View()
  dplyr::mutate(age_LEP = age) %>% 
  dplyr::filter(!is.na(age_LEP)) %>% 
  dplyr::filter(year != "2020") %>% 
  dplyr::group_by(year, sex, age_LEP) %>%
  dplyr::summarise(n_ages = n_distinct(Alu)) %>% 
  dplyr::group_by(year, sex) %>% 
  dplyr::mutate(prop = n_ages/sum(n_ages)) %>%
  dplyr::mutate(age_LEP = as.character(age_LEP)) %>% 
  dplyr::mutate(age_LEP = ifelse(age_LEP == "3", "3+", age_LEP)) %>% 
  dplyr::mutate(age_LEP = fct_rev(as.factor(age_LEP)))

age_plumage_dist_data2 <- 
  age_plumage_dist_data %>% 
  dplyr::group_by(sex, age_LEP) %>% 
  dplyr::summarise(mean_n_ages = mean(n_ages)) %>% 
  dplyr::group_by(sex) %>% 
  dplyr::mutate(prop = mean_n_ages/sum(mean_n_ages),
                year = "Grand\naverage") %>% 
  dplyr::bind_rows(age_plumage_dist_data, .) %>% 
  dplyr::mutate(sex = ifelse(sex == "F", "Female", "Male"))

ggplot(age_plumage_dist_data2, 
       aes(fill = sex, alpha = age_LEP, 
           y = prop, x = sex)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(. ~ year) +
  theme_bw() +
  theme(text = element_text(family = "Franklin Gothic Book"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.1, "cm"),
    panel.border = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    strip.background = element_blank(),
    legend.position = "top") +
  scale_fill_manual(values = sex_pal2, name = "Sex", guide = FALSE) +
  scale_alpha_manual(values = c(0.3, 0.65, 1), name = "Age (years)",
                     guide = guide_legend(title.position = "top", title.hjust = 0.5, label.position = "bottom")) +
  ylab("Proportion of breeding adults") +
  xlab("Sex")

age_plumage_dist_data2 %>% 
  dplyr::filter(str_detect(year, "Grand")) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(sex, age_LEP, mean_n_ages, prop) %>% 
  dplyr::arrange(sex, desc(age_LEP)) %>% 
  gt(rowname_col = "row",
     groupname_col = "sex") %>%
  cols_label(age_LEP = html("<i>Age class</i>"), 
             mean_n_ages = "Annual\nmean no. ind.",
             prop = "Proportion") %>% 
  fmt_number(columns = vars(mean_n_ages, prop),
             decimals = 2,
             use_seps = FALSE) %>% 
  cols_align(align = "right",
             columns = vars(mean_n_ages)) %>%
  cols_align(align = "left",
             columns = vars(prop)) %>%
  tab_options(row_group.font.weight = "bold",
              row_group.background.color = brewer.pal(9,"Greys")[3],
              table.font.size = 12,
              data_row.padding = 3,
              row_group.padding = 4,
              summary_row.padding = 2,
              column_labels.font.size = 14,
              row_group.font.size = 12,
              table.width = pct(60))
  

  # dplyr::mutate(age = ifelse(age == "adult", 1, age)) %>%
  # 
  # dplyr::filter(age_LEP != age & year == 2018) %>%
  # arrange(age)# %>%
  # # pull(SN)
  # # dplyr::select(Alu, sex   primarycoverts primaries Secondarycoverts secondaries tail)
  # # dplyr::filter(primaries == "onebarred")
  # # dplyr::filter(secondaries == "onebarred")
  # # dplyr::filter(age %in% c("adult", "juvenile"))
  # # dplyr::filter(is.na(age_LEP))
  # 
  # dplyr::group_by(year, sex, age_LEP, age) %>% 
  # dplyr::summarise(n_ages = n()) %>% 
  # dplyr::filter(age_LEP != age) %>% 
  # dplyr::select(SN, year, species, Alu, sex, date, primarycoverts, primaries, 
  #               Secondarycoverts, secondaries, tail, legcoverts, face, age, age_LEP) %>% 
  # write.csv(file = "data/raw/age_classes_2001-2020_LEP.csv")
  # View()
  
  dplyr::pull(age_LEP) %>% 
  as.factor() %>% 
  levels
  dplyr::select(species, ring_ID, sex, year, site, age_status, date_dec) %>% 
  # dplyr::select(-month, -day, -date_dec) %>% 
  mutate(across(c("sex", "site", "age_status"), tolower)) %>%
  mutate(sex = ifelse(sex == "Female", "F", ifelse(sex == "Male", "M", "XXX")),
         age_status = ifelse(age_status == "adult", "A", ifelse(age_status == "juvenile", "J", "XXX")),
         ring_ID = str_replace_all(string = ring_ID, fixed(" "), "")) %>% 
  mutate(across(c("species", "ring_ID", "sex", "site", "age_status"), as.factor)) %>% 
  mutate(sex = ifelse(str_detect(ring_ID, pattern = "[Ff]emale"), "F",
                      ifelse(str_detect(ring_ID, pattern = "[Mm]ale"), "M", sex)),
         month_year = format(date, "%Y-%m"),
         month_numeric = as.numeric(month),
         year_numeric = as.numeric(year),
         week_numeric = as.numeric(strftime(date, format = "%V"))) %>% 
  mutate(month_std = round(scale_by(month_numeric ~ year_numeric, ., scale = 0), digits = 0)) %>% 
  mutate(month_std = month_std + abs(min(month_std, na.rm = TRUE))) %>% 
  mutate(week_std = round(scale_by(week_numeric ~ year_numeric, ., scale = 0), digits = 0)) %>% 
  mutate(week_std = week_std + abs(min(week_std, na.rm = TRUE)))
  