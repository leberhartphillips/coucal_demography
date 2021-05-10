BC_ISR_results <- 
  BC_treat_sensitivity_analysis$ISR_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value")

BC_HSR_results <- 
  BC_treat_sensitivity_analysis$HSR_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value")

BC_h_results <- 
  BC_treat_sensitivity_analysis$h_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value")

BC_k_results <- 
  BC_treat_sensitivity_analysis$k_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value")

BC_treat_sensitivity_analysis$vr_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value") %>% 
  left_join(., BC_ISR_results, by = "value") %>% 
  left_join(., BC_HSR_results, by = "value") %>% 
  pivot_longer(!value, names_to = "vr", values_to = "ASR") %>% 
  mutate(value = as.numeric(value),
         ASR = as.numeric(ASR)) %>% 
  ggplot() +
  geom_line(aes(x = value, y = ASR, color = vr)) +
  # geom_vline(xintercept = 0.5) +
  # geom_vline(xintercept = BC_HSR) +
  geom_hline(yintercept = 0.5) +
  theme(legend.position = "right")


WBC_ISR_results <- 
  WBC_treat_sensitivity_analysis$ISR_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value")

WBC_HSR_results <- 
  WBC_treat_sensitivity_analysis$HSR_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value")

WBC_treat_sensitivity_analysis$vr_pert_ASR %>%
  as.data.frame() %>% 
  rownames_to_column(var = "value") %>% 
  left_join(., WBC_ISR_results, by = "value") %>% 
  left_join(., WBC_HSR_results, by = "value") %>% 
  pivot_longer(!value, names_to = "vr", values_to = "ASR") %>% 
  mutate(value = as.numeric(value),
         ASR = as.numeric(ASR)) %>% 
  ggplot() +
  geom_line(aes(x = value, y = ASR, color = vr)) +
  # geom_vline(xintercept = 0.5) +
  # geom_vline(xintercept = WBC_HSR) +
  geom_hline(yintercept = 0.5) +
  theme(legend.position = "right")
  
            
            BC_treat_sensitivity_analysis$ISR_pert_ASR)

row.names(BC_treat_sensitivity_analysis$ISR_pert_ASR)
