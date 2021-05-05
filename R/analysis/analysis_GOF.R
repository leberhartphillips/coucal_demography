# # load packages
# source("scripts/01_libraries.R")
# 
# # load functions
# function.sources = list.files(path = "scripts", 
#                               pattern = "*\\().R$", full.names = TRUE, 
#                               ignore.case = TRUE)
# sapply(function.sources, source, .GlobalEnv)
# 
# # load capture histories
# data.sources = list.files(path = "cooked_data", 
#                           pattern="*ch.rds$", full.names = TRUE, 
#                           ignore.case = TRUE)
# sapply(data.sources, load, .GlobalEnv)

## ----loadlibraries---------------------
# install.packages("R2ucare") # first time only
library(R2ucare) # For Goodness of fit tests
library(dplyr) # for tidy data
library(magrittr) # for pipes
library(marked)

## ----reshapedata---------------------------------------------------------
# Load dipper data (in "marked", but also in R2ucare package in different format)
# I'm loading `marked` version to be consistent with other tutorials/workshops
# and for fitting models below in "marked".
data(dipper, package = "marked")

# Full data
dipper.ch.gof <- dipper$ch %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper))

BC.ch.gof <- 
  Black_Coucal_fledgling_CJS_ch$ch %>% 
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(Black_Coucal_fledgling_CJS_ch))

# Females only
dipper.fem.ch.gof <- dipper$ch[dipper$sex == "Female"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper[dipper$sex == "Female",]))

BC.fem.ch.gof <- 
  Black_Coucal_fledgling_CJS_ch$ch[Black_Coucal_fledgling_CJS_ch$sex == "F"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "F",]))

# Males only
dipper.mal.ch.gof <- dipper$ch[dipper$sex == "Male"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(dipper[dipper$sex == "Male",]))

BC.mal.ch.gof <- 
  Black_Coucal_fledgling_CJS_ch$ch[Black_Coucal_fledgling_CJS_ch$sex == "M"] %>%
  strsplit('') %>%
  sapply(`[`) %>%
  t() %>%
  unlist() %>%
  as.numeric %>%
  matrix(nrow = nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "M",]))

## ----overallcjsgoffem----------------------------------------------------
# first argument = capture history matrix, second argument = capture history frequency (vector of 1's for our example)
R2ucare::overall_CJS(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",]))) # Females only
R2ucare::overall_CJS(BC.fem.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "F",]))) # Females only
R2ucare::overall_CJS(BC.mal.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "M",]))) # Females only
R2ucare::overall_CJS(BC.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch))) # Females only

## ----test2ct-------------------------------------------------------------
# first argument = capture history matrix, second argument is frequency of each capture history (1 for example)
test2ct_fem_dipper <- test2ct(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",]))) 
test2ct_fem_dipper

test2ct_fem <- test2ct(BC.fem.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "F",]))) 
test2ct_fem

test2ct_mal <- test2ct(BC.mal.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "M",]))) 
test2ct_mal

test2ct_all <- test2ct(BC.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch))) 
test2ct_all

## ----test2cl-------------------------------------------------------------
test2cl_fem_dipper <- test2cl(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
test2cl_fem_dipper

test2cl_fem <- test2cl(BC.fem.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "F",])))
test2cl_fem

test2cl_mal <- test2cl(BC.mal.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "M",])))
test2cl_mal

test2cl_all <- test2cl(BC.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch)))
test2cl_all

## ----test3sr-------------------------------------------------------------
test3sr_fem_dipper <- test3sr(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
test3sr_fem_dipper

test3sr_fem <- test3sr(BC.fem.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "F",])))
test3sr_fem

test3sr_mal <- test3sr(BC.mal.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "M",])))
test3sr_mal

test3sr_all <- test3sr(BC.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch)))
test3sr_all

## ----test3sm-------------------------------------------------------------
test3sm_fem_dipper <- test3sm(dipper.fem.ch.gof, rep(1,nrow(dipper[dipper$sex == "Female",])))
test3sm_fem_dipper

test3sm_fem <- test3sm(BC.fem.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "F",])))
test3sm_fem

test3sm_mal <- test3sm(BC.mal.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch[Black_Coucal_fledgling_CJS_ch$sex == "M",])))
test3sm_mal

test3sm_all <- test3sm(BC.ch.gof, rep(1,nrow(Black_Coucal_fledgling_CJS_ch)))
test3sm_all

## ----test2test3addition--------------------------------------------------
test2ct_fem_dipper$test2ct[1] + test2cl_fem_dipper$test2cl[1] + test3sr_fem_dipper$test3sr[1] + test3sm_fem_dipper$test3sm[1]

test2ct_fem$test2ct[1] + test2cl_fem$test2cl[1] + test3sr_fem$test3sr[1] + test3sm_fem$test3sm[1]

test2ct_mal$test2ct[1] + test2cl_mal$test2cl[1] + test3sr_mal$test3sr[1] + test3sm_mal$test3sm[1]

test2ct_all$test2ct[1] + test2cl_all$test2cl[1] + test3sr_all$test3sr[1] + test3sm_all$test3sm[1]