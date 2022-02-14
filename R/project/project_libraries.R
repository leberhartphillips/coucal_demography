## load (and install if necessary) packages need for project

# a vector of all the packages needed in the project's scipts
packages_required_in_project <- c("RMark",
                                  "tidyverse",
                                  "readxl",
                                  "BaSTA",
                                  "pbapply",
                                  "RColorBrewer",
                                  "grid",
                                  "Rmisc",
                                  "gss",
                                  "arm",
                                  "partR2",
                                  "parameters",
                                  "standardize",
                                  "colorBlindness",
                                  "ggthemes",
                                  "patchwork",
                                  "gt",
                                  "rptR",
                                  "tidybayes",
                                  "broom.mixed",
                                  "effects",
                                  "patchwork",
                                  "devtools",
                                  "unmarked",
                                  "R2ucare",
                                  "marked")

# of the required packages, check if some need to be installed
new.packages <- 
  packages_required_in_project[!(packages_required_in_project %in% 
                                   installed.packages()[,"Package"])]

# install all packages that are not locally available
if(length(new.packages)) install.packages(new.packages)

# load all the packages into the current R session
lapply(packages_required_in_project, require, character.only = TRUE)