source("R/project/project_libraries.R")

# rearrange_ch_for_age() takes a capture history in known fate format that has
# intervals specifying a date of the breeding season and rearranges it so that
# the first interval is the first day of age for a given individual
rearrange_ch_for_age <- 
  function(df, ch_col, age_first_obs_col, age_last_obs_col, max_age, min_age){
    
    df[, ch_col] <- gsub("(^|[^0-9])0+", "\\1", df[, ch_col], perl = TRUE)
    
    for(i in 1:nrow(df)){
      df[i, ch_col] <- paste(paste(rep("0", 
                                       times = ((df[i, age_first_obs_col] - 1) * 2)), 
                                   collapse = ""), 
                             df[i, ch_col], 
                             sep = "")
    }
    
    df[, ch_col] <- str_sub(df[, ch_col], start = 1, max_age * 2)
    
    for(j in 1:nrow(df)){
      
      if(nchar(df[j, ch_col] < max_age * 2)){
        df[j, ch_col] <- paste(df[j, ch_col], 
                               paste(rep("0", 
                                         times = ((max_age * 2) - nchar(df[j, ch_col]))), 
                                     collapse = ""), 
                               sep = "")
      }
      else{
        df[j, ch_col] <- df[j, ch_col]
      }
    }
    df[, ch_col] <- str_sub(df[, ch_col], start = (min_age * 2) - 1, end = max_age * 2)
    df
  }