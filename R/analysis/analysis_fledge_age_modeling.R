#### Script Preparation ----
# load packages
source("R/project/project_libraries.R")

# load functions
function.sources = list.files(path = "R/functions",
                              pattern = "*\\().R$", full.names = TRUE, 
                              ignore.case = TRUE)
try (sapply(function.sources, source), silent = TRUE)

#### Data Wrangle ----
# load raw data
fledge_dat <- 
  read.csv("data/raw/Coucal_chick_survival_2001-2019_20200129.csv", 
           header = TRUE, stringsAsFactors = FALSE, na.strings = c("", " ", "NA")) %>%
  dplyr::select(species, sex, Fledge_age, nest_ID, year, lab_no, Ring_ID, site) %>% 
  filter(site == "Kapunga" & species != "CTC") %>% 
  mutate(Nst_No = str_replace_all(string = nest_ID, fixed(" "), ""),
         Ring_ID = str_replace_all(string = Ring_ID, fixed(" "), ""),
         lab_no = str_replace_all(string = lab_no, fixed(" "), "")) %>% 
  filter(!is.na(Fledge_age)) %>% 
  mutate(Ind_ID = paste(Nst_No, Ring_ID, lab_no, sep = "_")) %>%
  mutate(sex_plot = ifelse(sex == "M", 2.2, 0.8))

# # assess normality of "days_since_fledging" variable
# ggplot(data = fledge_dat) +
#   geom_histogram(aes(Fledge_age), binwidth = 1) +
#   facet_grid(sex ~ species)

# check for repeated measures within nest
# fledge_dat %>% 
#   group_by(nest_ID) %>% 
#   dplyr::summarise(n_ = n())

#### Modeling (BC) ----
# "Fledge_age" as dependent variable, interaction with sex and species
mod_fledge_age_BC <- 
  lmer(Fledge_age ~ sex + 
         (1 | nest_ID) + (1 | year), 
       data = filter(fledge_dat, species == "BC"))

# extract model coefficients
mod_fledge_age_BC_coefs <- 
  model_parameters(mod_fledge_age_BC) %>%
  as.data.frame(.)
# run tidy bootstrap to obtain model diagnostics
tidy_mod_fledge_age_BC <-
  tidy(mod_fledge_age_BC, conf.int = TRUE, conf.method = "boot", nsim = 1000)

# run rptR to obtain repeatabilities of random effects
rpt_mod_fledge_age_BC <-
  rpt(Fledge_age ~ sex + 
        (1 | nest_ID) + (1 | year),
      grname = c("nest_ID", "year", "Fixed"),
      data = filter(fledge_dat, species == "BC"),
      datatype = "Gaussian",
      nboot = 1000, npermut = 1000, ratio = TRUE,
      adjusted = TRUE, ncores = 4, parallel = TRUE)

# run partR2 to obtain marginal R2, parameter estimates, and beta weights
R2m_mod_fledge_age_BC <-
  partR2(mod_fledge_age_BC,
         partvars = c("sex"),
         R2_type = "marginal",
         nboot = 1000, CI = 0.95, max_level = 1)

R2c_mod_fledge_age_BC <-
  partR2(mod_fledge_age_BC,
         partvars = c("sex"),
         R2_type = "conditional",
         nboot = 1000, CI = 0.95, max_level = 1)

# save model, tidy, rptR, and partR2 output as a list
stats_mod_fledge_age_BC <-
  list(mod = mod_fledge_age_BC,
       tidy = tidy_mod_fledge_age_BC,
       rptR = rpt_mod_fledge_age_BC,
       partR2m = R2m_mod_fledge_age_BC,
       partR2c = R2c_mod_fledge_age_BC)

save(stats_mod_fledge_age_BC,
     file = "output/mod_stats/Stats_mod_fledge_age_BC.rds")

load(file = "output/mod_stats/stats_mod_fledge_age_BC.rds")

#### Table of effect sizes (BC) ----
# Retrieve sample sizes
sample_sizes <-
  fledge_dat %>% 
  filter(species == "BC") %>% 
  summarise(Year = n_distinct(year),
            Individual = n_distinct(Ring_ID),
            Nests = n_distinct(nest_ID))

sample_sizes <- 
  as.data.frame(t(as.data.frame(sample_sizes))) %>%
  rownames_to_column("term") %>% 
  dplyr::rename(estimate = V1) %>% 
  mutate(stat = "n")

# clean model component names
mod_comp_names <- 
  data.frame(comp_name = c("Sex (Male)",
                           "Total Marginal \U1D479\U00B2",
                           "Sex",
                           "Total Conditional \U1D479\U00B2",
                           "Nest",
                           "Year",
                           "Residual",
                           "Nest",
                           "Year",
                           "Residual",
                           "Years",
                           "Individuals",
                           "Observations (i.e., Nests)"))

# Fixed effect sizes (non-standardized)
fixefTable <- 
  stats_mod_fledge_age_BC$tidy %>% 
  dplyr::filter(effect == "fixed") %>% 
  dplyr::select(term, estimate, conf.low, conf.high) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed")

# Fixed effect sizes (standardized)
fixef_bw_Table <- 
  stats_mod_fledge_age_BC$partR2m$BW %>% 
  # dplyr::select(term, estimate, CI_lower, CI_upper) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed_bw") %>% 
  dplyr::rename(conf.low = CI_lower,
                conf.high = CI_upper)

# Semi-partial R2 estimates
fixef_bw_Table <-
  model_parameters(stats_mod_fledge_age_BC$mod, standardize = "refit") %>%
  dplyr::select(Parameter, Coefficient, CI_low, CI_high) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed_bw") %>% 
  dplyr::rename(conf.low = CI_low,
                conf.high = CI_high,
                term = Parameter,
                estimate = Coefficient) %>% 
  filter(term != "(Intercept)")

R2Table <- 
  bind_rows(stats_mod_fledge_age_BC$partR2m$R2,
            stats_mod_fledge_age_BC$partR2c$R2[1,]) %>% 
  dplyr::select(term, estimate, CI_lower, CI_upper) %>% 
  as.data.frame() %>% 
  mutate(stat = "partR2") %>% 
  dplyr::rename(conf.low = CI_lower,
                conf.high = CI_upper)

# Random effects variances
ranefTable <- 
  stats_mod_fledge_age_BC$tidy %>% 
  dplyr::filter(effect == "ran_pars") %>% 
  dplyr::select(group, estimate, conf.low, conf.high) %>% 
  as.data.frame() %>% 
  mutate(stat = "rand") %>% 
  dplyr::rename(term = group) %>% 
  mutate(estimate = estimate^2,
         conf.high = conf.high^2,
         conf.low = conf.low^2)

# Adjusted repeatabilities
coefRptTable <- 
  stats_mod_fledge_age_BC$rptR$R_boot %>% 
  dplyr::select(-Fixed) %>% 
  mutate(residual = 1 - rowSums(.)) %>% 
  apply(., 2, 
        function(x) c(mean (x), quantile (x, prob = c(0.025, 0.975)))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("term") %>% 
  dplyr::rename(estimate = V1,
                conf.low = `2.5%`,
                conf.high = `97.5%`) %>% 
  mutate(stat = "RptR")

# Store all parameters into a single table and clean it up
allCoefs_mod <- 
  bind_rows(fixef_bw_Table,
            R2Table,
            ranefTable, 
            coefRptTable, 
            sample_sizes) %>% 
  bind_cols(.,
            mod_comp_names) %>%
  mutate(coefString = ifelse(!is.na(conf.low),
                             paste0("[", 
                                    round(conf.low, 2), ", ", 
                                    round(conf.high, 2), "]"),
                             NA),
         effect = c(rep("Fixed effects \U1D6FD (standardized)", nrow(fixef_bw_Table)),
                    rep("Partitioned \U1D479\U00B2", nrow(R2Table)),
                    rep("Random effects \U1D70E\U00B2", nrow(ranefTable)),
                    rep("Adjusted repeatability \U1D45F", nrow(coefRptTable)),
                    rep("Sample sizes \U1D45B", nrow(sample_sizes)))) %>%
  dplyr::select(effect, everything())

# re-organize model components for table
# allCoefs_mod <-
#   allCoefs_mod[c(5, 6, 1:4, 7, 13, 11, 12, 8:10, 14:22), ]

# draw gt table
table_mod_fledge_age_BC <- 
  allCoefs_mod %>% 
  dplyr::select(effect, comp_name, estimate, coefString) %>% 
  gt(rowname_col = "row",
     groupname_col = "effect") %>% 
  cols_width(vars(comp_name) ~ pct(50),
             vars(estimate) ~ pct(20),
             vars(coefString) ~ pct(30)) %>% 
  cols_label(comp_name = html("<i>Black Coucals <br> age at leaving nest</i>"),
             estimate = "Mean estimate",
             coefString = "95% confidence interval") %>% 
  fmt_number(columns = vars(estimate),
             rows = 1:10,
             decimals = 2,
             use_seps = FALSE) %>% 
  fmt_number(columns = vars(estimate),
             rows = 11:13,
             decimals = 0,
             use_seps = FALSE) %>% 
  fmt_missing(columns = 1:4,
              missing_text = "") %>% 
  cols_align(align = "left",
             columns = vars(comp_name)) %>% 
  tab_options(row_group.font.weight = "bold",
              row_group.background.color = brewer.pal(9,"Greys")[3],
              table.font.size = 12,
              data_row.padding = 3,
              row_group.padding = 4,
              summary_row.padding = 2,
              column_labels.font.size = 14,
              row_group.font.size = 12,
              table.width = pct(60))

table_mod_fledge_age_BC

#### Modeling (WBC) ----
# "Fledge_age" as dependent variable, interaction with sex and species
mod_fledge_age_WBC <- 
  lmer(Fledge_age ~ sex + 
         (1 | nest_ID) + (1 | year), 
       data = filter(fledge_dat, species == "WBC"))

# extract model coefficients
mod_fledge_age_WBC_coefs <- 
  model_parameters(mod_fledge_age_WBC) %>%
  as.data.frame(.)

# run tidy bootstrap to obtain model diagnostics
tidy_mod_fledge_age_WBC <-
  tidy(mod_fledge_age_WBC, conf.int = TRUE, conf.method = "boot", nsim = 1000)

# run rptR to obtain repeatabilities of random effects
rpt_mod_fledge_age_WBC <-
  rpt(Fledge_age ~ sex + 
        (1 | nest_ID) + (1 | year),
      grname = c("nest_ID", "year", "Fixed"),
      data = filter(fledge_dat, species == "WBC"),
      datatype = "Gaussian",
      nboot = 1000, npermut = 1000, ratio = TRUE,
      adjusted = TRUE, ncores = 4, parallel = TRUE)

# run partR2 to obtain marginal R2, parameter estimates, and beta weights
R2m_mod_fledge_age_WBC <-
  partR2(mod_fledge_age_WBC,
         partvars = c("sex"),
         R2_type = "marginal",
         nboot = 1000, CI = 0.95, max_level = 1)

R2c_mod_fledge_age_WBC <-
  partR2(mod_fledge_age_WBC,
         partvars = c("sex"),
         R2_type = "conditional",
         nboot = 1000, CI = 0.95, max_level = 1)

# save model, tidy, rptR, and partR2 output as a list
stats_mod_fledge_age_WBC <-
  list(mod = mod_fledge_age_WBC,
       tidy = tidy_mod_fledge_age_WBC,
       rptR = rpt_mod_fledge_age_WBC,
       partR2m = R2m_mod_fledge_age_WBC,
       partR2c = R2c_mod_fledge_age_WBC)

save(stats_mod_fledge_age_WBC,
     file = "output/mod_stats/Stats_mod_fledge_age_WBC.rds")

load(file = "output/mod_stats/stats_mod_fledge_age_WBC.rds")

#### Table of effect sizes (WBC) ----
# Retrieve sample sizes
sample_sizes <-
  fledge_dat %>% 
  filter(species == "WBC") %>% 
  summarise(Year = n_distinct(year),
            Individual = n_distinct(Ring_ID),
            Nests = n_distinct(nest_ID))

sample_sizes <- 
  as.data.frame(t(as.data.frame(sample_sizes))) %>%
  rownames_to_column("term") %>% 
  dplyr::rename(estimate = V1) %>% 
  mutate(stat = "n")

# clean model component names
mod_comp_names <- 
  data.frame(comp_name = c("Sex (Male)",
                           "Total Marginal \U1D479\U00B2",
                           "Sex",
                           "Total Conditional \U1D479\U00B2",
                           "Nest",
                           "Year",
                           "Residual",
                           "Nest",
                           "Year",
                           "Residual",
                           "Years",
                           "Individuals",
                           "Observations (i.e., Nests)"))

# Fixed effect sizes (non-standardized)
fixefTable <- 
  stats_mod_fledge_age_WBC$tidy %>% 
  dplyr::filter(effect == "fixed") %>% 
  dplyr::select(term, estimate, conf.low, conf.high) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed")

# Fixed effect sizes (standardized)
fixef_bw_Table <- 
  stats_mod_fledge_age_WBC$partR2m$BW %>% 
  # dplyr::select(term, estimate, CI_lower, CI_upper) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed_bw") %>% 
  dplyr::rename(conf.low = CI_lower,
                conf.high = CI_upper)

# Semi-partial R2 estimates
fixef_bw_Table <-
  model_parameters(stats_mod_fledge_age_WBC$mod, standardize = "refit") %>%
  dplyr::select(Parameter, Coefficient, CI_low, CI_high) %>% 
  as.data.frame() %>% 
  mutate(stat = "fixed_bw") %>% 
  dplyr::rename(conf.low = CI_low,
                conf.high = CI_high,
                term = Parameter,
                estimate = Coefficient) %>% 
  filter(term != "(Intercept)")

R2Table <- 
  bind_rows(stats_mod_fledge_age_WBC$partR2m$R2,
            stats_mod_fledge_age_WBC$partR2c$R2[1,]) %>% 
  dplyr::select(term, estimate, CI_lower, CI_upper) %>% 
  as.data.frame() %>% 
  mutate(stat = "partR2") %>% 
  dplyr::rename(conf.low = CI_lower,
                conf.high = CI_upper)

# Random effects variances
ranefTable <- 
  stats_mod_fledge_age_WBC$tidy %>% 
  dplyr::filter(effect == "ran_pars") %>% 
  dplyr::select(group, estimate, conf.low, conf.high) %>% 
  as.data.frame() %>% 
  mutate(stat = "rand") %>% 
  dplyr::rename(term = group) %>% 
  mutate(estimate = estimate^2,
         conf.high = conf.high^2,
         conf.low = conf.low^2)

# Adjusted repeatabilities
coefRptTable <- 
  stats_mod_fledge_age_WBC$rptR$R_boot %>% 
  dplyr::select(-Fixed) %>% 
  mutate(residual = 1 - rowSums(.)) %>% 
  apply(., 2, 
        function(x) c(mean (x), quantile (x, prob = c(0.025, 0.975)))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("term") %>% 
  dplyr::rename(estimate = V1,
                conf.low = `2.5%`,
                conf.high = `97.5%`) %>% 
  mutate(stat = "RptR")

# Store all parameters into a single table and clean it up
allCoefs_mod <- 
  bind_rows(fixef_bw_Table,
            R2Table,
            ranefTable, 
            coefRptTable, 
            sample_sizes) %>% 
  bind_cols(.,
            mod_comp_names) %>%
  mutate(coefString = ifelse(!is.na(conf.low),
                             paste0("[", 
                                    round(conf.low, 2), ", ", 
                                    round(conf.high, 2), "]"),
                             NA),
         effect = c(rep("Fixed effects \U1D6FD (standardized)", nrow(fixef_bw_Table)),
                    rep("Partitioned \U1D479\U00B2", nrow(R2Table)),
                    rep("Random effects \U1D70E\U00B2", nrow(ranefTable)),
                    rep("Adjusted repeatability \U1D45F", nrow(coefRptTable)),
                    rep("Sample sizes \U1D45B", nrow(sample_sizes)))) %>%
  dplyr::select(effect, everything())

# re-organize model components for table
# allCoefs_mod <-
#   allCoefs_mod[c(5, 6, 1:4, 7, 13, 11, 12, 8:10, 14:22), ]

# draw gt table
table_mod_fledge_age_WBC <- 
  allCoefs_mod %>% 
  dplyr::select(effect, comp_name, estimate, coefString) %>% 
  gt(rowname_col = "row",
     groupname_col = "effect") %>% 
  cols_width(vars(comp_name) ~ pct(50),
             vars(estimate) ~ pct(20),
             vars(coefString) ~ pct(30)) %>%
  cols_label(comp_name = html("<i>White-browed Coucals <br> age at leaving nest</i>"),
             estimate = "Mean estimate",
             coefString = "95% confidence interval") %>% 
  fmt_number(columns = vars(estimate),
             rows = 1:10,
             decimals = 2,
             use_seps = FALSE) %>% 
  fmt_number(columns = vars(estimate),
             rows = 11:13,
             decimals = 0,
             use_seps = FALSE) %>% 
  fmt_missing(columns = 1:4,
              missing_text = "") %>% 
  cols_align(align = "left",
             columns = vars(comp_name)) %>% 
  tab_options(row_group.font.weight = "bold",
              row_group.background.color = brewer.pal(9,"Greys")[3],
              table.font.size = 12,
              data_row.padding = 3,
              row_group.padding = 4,
              summary_row.padding = 2,
              column_labels.font.size = 14,
              row_group.font.size = 12,
              table.width = pct(40)) 

table_mod_fledge_age_WBC

# export tables to disk
table_mod_fledge_age_BC %>% 
  gtsave("table_mod_fledge_age_BC.rtf", path = "products/tables/rtf/")

table_mod_fledge_age_BC %>% 
  gtsave("table_mod_fledge_age_BC.png", path = "products/tables/png/")

table_mod_fledge_age_WBC %>% 
  gtsave("table_mod_fledge_age_WBC.rtf", path = "products/tables/rtf/")

table_mod_fledge_age_WBC %>% 
  gtsave("table_mod_fledge_age_WBC.png", path = "products/tables/png/")

#### Forest plot of model predictions ----

# extract model coefficients
mod_fledge_age_BC_coefs <- 
  model_parameters(mod_fledge_age_BC) %>%
  as.data.frame(.)

mod_fledge_age_WBC_coefs <- 
  model_parameters(mod_fledge_age_BC) %>%
  as.data.frame(.)

mod_fledge_age_BC_fits <- 
  as.data.frame(effect(term = "sex", mod = stats_mod_fledge_age_BC$mod, 
                       xlevels = list(sex = c("M", "F")))) %>%
  mutate(sex_plot = ifelse(sex == "M", 1.8, 1.2),
         species = "BC")

mod_fledge_age_WBC_fits <- 
  as.data.frame(effect(term = "sex", mod = stats_mod_fledge_age_WBC$mod, 
                       xlevels = list(sex = c("M", "F")))) %>%
  mutate(sex_plot = ifelse(sex == "M", 1.8, 1.2),
         species = "WBC")

mod_fledge_age_fits <- 
  bind_rows(mod_fledge_age_BC_fits, 
            mod_fledge_age_WBC_fits) %>% 
  `rownames<-`( NULL )

# plot_palette_recruit <- brewer.pal(6, "Dark2")[c(2,3)]

mod_fledge_age_fits[1, 2] - mod_fledge_age_fits[2, 2] 
mod_fledge_age_fits[1, 4] - mod_fledge_age_fits[2, 4]
mod_fledge_age_fits[1, 5] - mod_fledge_age_fits[2, 5] 

fledge_age_plot <-
  ggplot2::ggplot() + 
  geom_boxplot(data = fledge_dat,
               aes(x = sex_plot, y = Fledge_age,
                   group = sex, fill = sex),
               color = "grey50",
               width = 0.05, alpha = 0.5,
               position = position_dodge(width = 0)) +
  geom_errorbar(data = mod_fledge_age_fits, 
                aes(x = sex_plot, ymax = upper, ymin = lower),
                alpha = 1, color = "black", width = 0.05, lwd = 0.5) + 
  geom_point(data = mod_fledge_age_fits, 
             aes(x = sex_plot, y = fit, fill = sex),
             lwd = 3, shape = 21, color= "black") +
  geom_jitter(data = fledge_dat, 
              aes(x = sex, y = Fledge_age, 
                  group = sex, 
                  fill = sex, color = sex), 
              width = 0.02, alpha = 0.2, shape = 19) +
  coord_flip() +
  luke_theme +
  theme(legend.position = "none",
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "italic"),
        axis.title.x = element_text(size = 12),
        # axis.title.y = element_blank(),
        # axis.text.y = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        # panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "grey70", size = 0.1),
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(7.5, 20.5), breaks = seq(8, 20, 2)) +#c(-60, -30, 0, 30, 60)) +
  scale_x_discrete(labels = c("M" = "Male",
                              "F" = "Female")) +
  ylab(expression(paste("Age at leaving nest (days)" %+-%  "95% CI", sep = ""))) +
  # xlab("Origin") +
  scale_color_manual(values = sex_pal2) +
  scale_fill_manual(values = sex_pal2) +
  facet_grid(species ~ ., labeller = as_labeller(species_names))# +
# annotate(geom = "text", x = 0.5, y = 58,
#          label = "First nests of the season",
#          color = "black", size = 3, fontface = 'italic', hjust = 0)

fledge_age_plot



# ggsave(plot = fledge_age_plot,
#        filename = "products/figures/svg/fledge_age_plot.svg",
#        width = 6,
#        height = 4, units = "in")

ggsave(plot = fledge_age_plot,
       filename = "products/figures/jpg/fledge_age_plot.jpg",
       width = 6,
       height = 4, units = "in")

#### Forest plot of effect sizes ----
# Standardized fixed effects
laydate_mod_forest_plot_fixef <-
  allCoefs_mod %>%
  filter(str_detect(effect, "Fixed") & 
           term != "(Intercept)") %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 "Local recruit",
                                 "Between ind. last breeding age", 
                                 "Between ind. first breeding age", 
                                 "Within ind. quadratic age", "Within ind. linear age",
                                 "Mother tarsus length")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  # theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.x = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(face = 'italic', hjust = 0.5)) +
  ylab("Fixed effects") +
  xlab(expression(italic(paste("Standardized effect size (", beta,")" %+-% "95% CI", sep = "")))) +
  scale_y_discrete(position = "right") +
  ggtitle('Lay date model')

# Semi-partial R2 estimates
laydate_mod_forest_plot_partR2 <-
  allCoefs_mod %>%
  filter(str_detect(effect, "Partitioned") & str_detect(comp_name, "Conditional", negate = TRUE)) %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 # "Seasonality",
                                 "Local recruit",
                                 "Selective disappearance",
                                 "Selective appearance",
                                 "Senescence",
                                 "Mother tarsus length",
                                 "Total Marginal \U1D479\U00B2")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  # theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.x = element_blank(),
        plot.background = element_blank()) +
  scale_y_discrete(labels = c(#"Seasonality" = expression("Seasonality"),
    "Selective disappearance" = expression("Selective disappearance"),
    "Selective appearance" = expression("Selective appearance"),
    "Senescence" = expression("Senescence"),
    "Mother tarsus length" = expression("Mother tarsus length"),
    "Recruit status" = expression("Recruit status"),
    "Total Marginal \U1D479\U00B2" = expression(paste("Total marginal ", italic("R"), ''^{2}, sep = ""))),
    position = "right") +
  ylab(expression(paste("Semi-partial ", italic("R"),''^{2}, sep = ""))) +
  xlab(expression(italic(paste("Variance explained (R", ''^{2}, ")" %+-% "95% CI", sep = ""))))

# Random effect variances
laydate_mod_forest_plot_randef <-
  allCoefs_mod %>%
  filter(str_detect(effect, "Random")) %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 "Residual",
                                 "Year",
                                 "Individual")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  # theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.x = element_blank(),
        plot.background = element_blank()) +
  ylab("Random\neffects") +
  xlab(expression(italic(paste("Variance (", sigma, ''^{2}, ")" %+-% "95% CI", sep = "")))) +
  scale_y_discrete(position = "right")

# Adjusted repeatabilities
laydate_mod_forest_plot_rptR <-
  allCoefs_mod %>%
  filter(str_detect(effect, "repeat")) %>%
  mutate(comp_name = fct_relevel(comp_name,
                                 "Residual",
                                 "Year",
                                 "Individual")) %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     y = comp_name),
                 alpha = 1, color = col_all, 
                 size = 0.5,
                 height = 0) +
  geom_point(aes(y = comp_name, x = estimate),
             size = 3, shape = 21, 
             fill = "#ECEFF4", col = col_all, 
             alpha = 1, stroke = 0.5) +
  luke_theme +
  # theme(axis.title.x = element_text(size = 10)) +
  theme(axis.title.x = element_blank(),
        plot.background = element_blank()) +
  
  ylab("Intra-class\ncorrelation") +
  xlab(expression(italic(paste("Adjusted repeatability (r)" %+-% "95% CI", sep = "")))) +
  scale_y_discrete(position = "right")

# Patchwork plot
laydate_mod_forest_plot_combo <-
  (laydate_mod_forest_plot_fixef / laydate_mod_forest_plot_partR2 / 
     # laydate_mod_forest_plot_randef / 
     laydate_mod_forest_plot_rptR) + 
  plot_annotation(tag_levels = 'A', title = 'Lay date model', 
                  theme = theme(plot.title = element_text(face = 'italic', hjust = 0.2))) +
  plot_layout(heights = unit(c(4, 4, 
                               # 2.5, 
                               2.5), c('cm', 'cm', 
                                       # 'cm',
                                       'cm')))

laydate_mod_forest_plot_combo

forest_plot_combo <-
  ((eggv_mod_forest_plot_fixef + laydate_mod_forest_plot_fixef) / 
     (eggv_mod_forest_plot_partR2 + laydate_mod_forest_plot_partR2) / 
     # eggv_mod_forest_plot_randef / 
     (eggv_mod_forest_plot_rptR + laydate_mod_forest_plot_rptR)) + 
  # plot_annotation(tag_levels = 'A', title = 'Egg volume model', 
  #                 theme = theme(plot.title = element_text(face = 'italic', hjust = 0.2))) +
  plot_layout(heights = unit(c(4.5, 4, 
                               # 2.5, 
                               2.5), c('cm', 'cm', 
                                       # 'cm', 
                                       'cm')))
# forest_plot_combo

ggsave(plot = forest_plot_combo,
       filename = "products/figures/jpg/Figure_3.jpg",
       width = 8,
       height = 6.5, units = "in")
