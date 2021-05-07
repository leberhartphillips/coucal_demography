source("R/project/project_libraries.R")

# bootstrap_data() randomly samples with replacement from the nestling and 
# fledgling datasets, while making sure that if an individual existing in both
# datasets was sampled from the nestling data it was also sampled in the 
# fledgling data. Each bootstrapped sample has the same length as the original
# data.
bootstrap_data <- function(adult, fledgling, nestling, 
                           num_boot, species, iter_add) {
  
  # sample a new nestling mark-recapture dataset containing only one nest 
  # member per draw
  nestling_boot <- 
    nestling %>% 
    group_by(nest_ID) %>% 
    sample_n(1)
  
  # sample a new fledgling mark-recapture dataset containing only one nest 
  # member per draw
  fledgling_boot <- 
    fledgling %>% 
    group_by(nest_ID) %>%
    sample_n(1)
  
  # sample a new adult mark-recapture dataset the same size as the original, 
  # with replacement
  adult_boot <-
    adult %>%
    sample_n(nrow(.), replace = TRUE)
  
  # make a list of these two datasets, which will be used in the next function
  out <- list(nestling_boot = nestling_boot, 
              fledgling_boot = fledgling_boot,
              adult_boot = adult_boot,
              iter = num_boot + ((iter_add - 1) * niter),
              species = species)
}
