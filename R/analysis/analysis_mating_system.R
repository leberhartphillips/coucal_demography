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

# # ***
# #   ## Quantifying mating system
# #   
# #   To put our estimate of ASR in the context of breeding behavior, we quantified sex bias in mating system based on behavioral obersvations from the field. Females of *Charadrius* species are more likely to desert broods and seek serial mates than males.  Thus, we expected that females would have more mates per year than males. 
# # 
# # **Step one: wrangle the data**
# #   remove any cases in which one mate was not identified (i.e., "NA")
# # ```{r}
# mating_df <- 
#   breeding_data[which(!is.na(breeding_data$female) & !is.na(breeding_data$male)),]
# # ```
# # determine the number of families used in the mating system analysis (i.e. the sample size)
# # ```{r}
# length(unique(mating_df$family_ID))
# # ```
# # bind the two mates together to make a unique pair
# # ```{r}
# mating_df$pair <- as.factor(paste(mating_df$female, mating_df$male, sep = "-"))
# # ```
# # determine how many mating attempts each individual had each year
# # ```{r}
# females <- reshape2::dcast(mating_df, female + population ~ year)
# males <- reshape2::dcast(mating_df, male + population ~ year)
# # ```
# # determine how many different mates each individual had over their lifetime in the population
# # ```{r}
# number_males_p_female <- 
#   stats::aggregate(male ~ female, mating_df, function(x) length(unique(x)))
# number_females_p_male <- 
#   stats::aggregate(female ~ male, mating_df, function(x) length(unique(x)))
# # ```
# # join these two dataframes together and define as numeric
# # ```{r}
# females <- dplyr::inner_join(females, number_males_p_female)
# females[,c(3:16)] <- 
#   lapply(females[,c(3:16)], as.numeric)
# males <- dplyr::inner_join(males, number_females_p_male)
# males[,c(3:16)] <- 
#   lapply(males[,c(3:16)], as.numeric)
# # ```
# # calculate the total number of mating attempts over each individual's lifetime
# # ```{r}
# females$attempts <- rowSums(females[, c(3:16)])
# males$attempts <- rowSums(males[, c(3:16)])
# # ```
# # calculate the number of years breeding
# # ```{r}
# females$years <- rowSums(females[, c(3:16)] > 0)
# males$years <- rowSums(males[, c(3:16)] > 0)
# # ```
# # filter out all individuals that only had one mating attempt
# # ```{r}
# females_no_1 <- dplyr::filter(females, male  != 1 | years != 1 | attempts != 1)
# males_no_1 <- dplyr::filter(males, female  != 1 | years != 1 | attempts != 1)
# # ```
# # tidy up dataframes then bind them together
# # ```{r}
# females_no_1$sex <- "Female"
# females_no_1$sex <- as.factor(females_no_1$sex)
# colnames(females_no_1)[c(1,17)] <- c("focal", "mate")
# males_no_1$sex <- "Male"
# males_no_1$sex <- as.factor(males_no_1$sex)
# colnames(males_no_1)[c(1,17)] <- c("focal", "mate")
# mating <- rbind(females_no_1, males_no_1)
# # ```
# # calculate the number of mates per year
# # ```{r}
# mating$no_mates_per_year <- mating$mate/mating$years
# # corrected for long-term monogamy (i.e., mu in Eq. 4), whereby 1 is given to 
# # individuals that have <1 mate per year
# 
# mating_dat %>% 
#   group_by(species, sex)
# 
# mating$no_mates_per_year_mono <- 
#   ifelse(mating$no_mates_per_year < 1, 1, mating$no_mates_per_year)
# # ```
# # summarise the matings by sex and determine "h", the average annual number of mates per female
# # ```{r}
sex_specific_mating_system <- 
  mating_dat %>% 
  group_by(species, sex) %>% 
  dplyr::summarise(mean_annual_no_mates = mean(Nr_partners, na.rm = TRUE),
                   var_annual_no_mates = var(Nr_partners, na.rm = TRUE),
                   median_annual_no_mates = median(Nr_partners, na.rm = TRUE),
                   sd_annual_no_mates = sd(Nr_partners, na.rm = TRUE),
                   n = n_distinct(ind_ID)) %>% 
  mutate(sample_size = paste("n = ", n, sep = ""),
         species_plot = ifelse(species == "WBC", 1.9, 1.1),
         species_lab = ifelse(species == "WBC", 1, 2))

# To obtain a female-based h index for each population the inverse of the mean 
# mu is calculated (i.e., Eq. 5)
BC_h <- 
  1/as.numeric(sex_specific_mating_system[which(
    sex_specific_mating_system$sex == "female" &
    sex_specific_mating_system$species == "BC"), "mean_annual_no_mates"])
WBC_h <- 
  1/as.numeric(sex_specific_mating_system[which(
    sex_specific_mating_system$sex == "female" &
    sex_specific_mating_system$species == "WBC"), "mean_annual_no_mates"])

# display the h values for each population (these are used in the mating function 
# of the matrix model)
BC_h
WBC_h

# ```
# 
# Figure S3: plot the sex-specific distributions of mating system
# ```{r, fig.width=6, fig.height=3, fig.align="center"}
# define the factor levels of the population variable so that the populations are in an 
# order that reflects the ASR (male biased to female biased)
mating_dat$species <- 
  factor(mating_dat$species ,
         levels = c("BC",
                    "WBC"))

mating_dat_plotting <-
  mating_dat %>% 
  mutate(species_dummy = ifelse(species == "WBC", 1, 0)) %>%
  mutate(species_plot = ifelse(species == "WBC", 2.1, 0.9))

mating_system_plot <-
  ggplot2::ggplot() +
  theme_bw() +
  geom_jitter(aes(y = Nr_partners, x = species_plot, 
                  color = sex),
              data = mating_dat_plotting, 
              width = 0.1, height = 0.1,
              alpha = 0.75, 
              fill = "white", shape = 16) +
  geom_pointrange(data = sex_specific_mating_system, 
    aes(y = mean_annual_no_mates, x = species_plot, 
        ymin = (mean_annual_no_mates-sd_annual_no_mates), 
        ymax = (mean_annual_no_mates+sd_annual_no_mates),
        fill = sex),
    size = 0.8, shape = 21, color = "black") +
  geom_text(aes(y = 0.2, x = species_lab, label = sample_size), 
            data = sex_specific_mating_system,
            size = 3) +
  geom_hline(yintercept = 1.5, linetype = "dashed",
             alpha = 0.5, color = "grey20") +
  facet_grid(. ~ sex, labeller = as_labeller(sex_names)) +
  luke_theme +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic")) +
  scale_x_continuous(limits = c(0.5, 2.5), 
                   breaks = c(1, 2), labels = species_names) +
  scale_y_continuous(limits = c(0, 5.5), expand = c(0, 0), 
                     breaks = c(1, 2, 3, 4, 5),
                     labels = c(1, 2, 3, 4, 5)) +
  ylab("Number of unique mates\nper year (± 1 SD)") +
  scale_color_manual(values = plot_palette_sex) +
  scale_fill_manual(values = plot_palette_sex) +
  annotate(geom = "text", y = 1.25, x = 1.5,
           label = "monogamy",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5) +
  annotate(geom = "text", y = 1.75, x = 1.5,
           label = "polygamy",
           color = "black", size = 3, fontface = 'italic', hjust = 0.5)

ggsave(mating_system_plot,
       filename = "products/figures/mating_system.jpeg",
       width = 6,
       height = 4,
       units = "in",
       dpi = 600)