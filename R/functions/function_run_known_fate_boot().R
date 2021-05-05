run_known_fate_boot <- 
  function(num_boot, ch_input,
           bootstrap_name,
           species, iter_add, prefix_number, stage_name){
    # run the survival analysis and ASR deduction on the sampled data
    result <- known_fate_boot(ch_input = ch_input, 
                              num_boot = num_boot, 
                              bootstrap_name = bootstrap_name,
                              species = species,
                              iter_add = iter_add,
                              prefix_number = prefix_number,
                              stage_name = stage_name)
    result
  }
