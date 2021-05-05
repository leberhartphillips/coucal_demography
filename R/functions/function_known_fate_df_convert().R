source("R/project/project_libraries.R")

# known_fate_df_convert() converts the census data to a capture history that can
# be used in the nestling and fledgling survival analysese
# df : raw data
# interval_col
# status_col : column in the data that specifies the status of the individual
# n_intervals : number of intervals to produce in the capture history
# stage : life-history stage (nestling or both)
# nestling_duration : time in days between hatching and leaving the nest
# id_col : column in the data that specifies the identity of the individual
# sex_col : column in the data that specifies the sex of the individual

known_fate_df_convert <- 
  function(df, interval_col, status_col, n_intervals, stage, nestling_duration,
           id_col, sex_col){
    
    if(missing(n_intervals)) {
      # determine the number of intervals to include
      n_intervals <- max(df[, interval_col])
    }
    
    df$age_interval_adj <- ifelse(df[, interval_col] > n_intervals, 
                                  n_intervals, df[, interval_col])
    
    # define as a dataframe
    df <- as.data.frame(df)
    
    if(stage == "fledgling"){
      # for-loop to specify the last encounter as a "11" if an individual was 
      # detected dead or "00" if an individual was censored then fill the rest
      # of the interval with "00"
      df$age_interval_adj <- df$age_interval_adj - nestling_duration
      df <- filter(df, age_interval_adj > 0)
      n_intervals <- n_intervals - nestling_duration
      for(i in 1:nrow(df)){
        ch <- c(rep("10", n_intervals))
        
        ch[df[, "age_interval_adj"][i]] <- ifelse(df[, status_col][i] == 1, "11", "00")
        
        if(df[, "age_interval_adj"][i] == n_intervals & df[, status_col][i] == 0){
          ch[df[, "age_interval_adj"][i]] <- "00"
        }
        else if(df[, "age_interval_adj"][i] == n_intervals & df[, status_col][i] == 1){
          ch[df[, "age_interval_adj"][i]] <- "11"
        }
        else{
          ch[c((df[, "age_interval_adj"][i] + 1) : n_intervals)] <- "00"
        }
        
        df$ch[i] <- paste(ch, collapse = "")
      }
      
      return(df)
    }
    if(stage == "nestling"){
      # for-loop to specify the last encounter as a "11" if an individual was 
      # detected dead or "00" if an individual was censored then fill the rest
      # of the interval with "00"
      for(i in 1:nrow(df)){
        ch <- c(rep("10", n_intervals))
        
        ch[df[, "age_interval_adj"][i]] <- ifelse(df[, status_col][i] == "no", "11", "10")
        
        if(df[, "age_interval_adj"][i] < n_intervals & df[, status_col][i] == "no"){
          ch[c((df[, "age_interval_adj"][i] + 1) : n_intervals)] <- "00"
        }
        else if(df[, "age_interval_adj"][i] == n_intervals & df[, status_col][i] == "no"){
          ch[df[, "age_interval_adj"][i]] <- "11"
        }
        df$ch[i] <- paste(ch, collapse = "")
      }
      
      return(df)
    }
    if(stage == "both"){
      # for-loop to specify the last encounter as a "0" if an individual was
      # detected dead then fill the rest of the interval with "00"
      
      life_table <- matrix(data = numeric(nrow(df) * (max(df[, interval_col]) + 4)), 
                           nrow = nrow(df), ncol = (max(df[, interval_col]) + 4))
      
      for(i in 1:nrow(df)){
        
        life_table[i,] <- c(rep("1", max(df[, interval_col])), "1", NA, NA, NA)
        
        life_table[i, c(df[, interval_col][i] : max(df[, interval_col]))] <- "0"
        
        life_table[i, (max(df[, interval_col]) + 2)] <- 
          ifelse(df[, status_col][i] == 0, "0", 
                 ifelse(df[, status_col][i] == 1, as.character(df[, interval_col][i])))
        
        life_table[i, (max(df[, interval_col]) + 3)] <- as.character(df[, id_col][i])
        
        life_table[i, (max(df[, interval_col]) + 4)] <- as.character(df[, sex_col][i])
        
        
      }
      life_table <- as.data.frame(life_table)
      colnames(life_table) <- as.character(c(1:max(df[, interval_col]), "birth", "death", "ID", "sex"))
      return(life_table)
    }
  }