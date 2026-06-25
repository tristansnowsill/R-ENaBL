#########################################################################################################################################################
#########################################################################################################################################################
#
# 4: Natural history
#
#########################################################################################################################################################
#########################################################################################################################################################

# Natural history parameters are samples
# One row for each sample, one column for each variable
# Save as csv or equivalent



sample_initial_state <- function(States_list, PLCO, All_parameters) {
  
  Cen_logit_PLCO <- log(PLCO / (1 - PLCO)) + 3.66649
  
  NSCLC_PLCO_coef <- All_parameters %>% filter(Parameter_label == "nh_entry_nsclc_plco") %>% pull(Modelled_value)
  SCLC_PLCO_coef <- All_parameters %>% filter(Parameter_label == "nh_entry_sclc_plco") %>% pull(Modelled_value)
  
  Int_data <- All_parameters %>%
    filter(startsWith(Parameter_label, "nh_entry"),
           endsWith(Parameter_label, "int")) %>%
    select(Parameter_label, Modelled_value) %>%
    mutate(
      State = case_when(
        Parameter_label == "nh_entry_nsclc_IA1_int" ~ "Preclinical_NSCLC_IA1",
        Parameter_label == "nh_entry_nsclc_IA2_int" ~ "Preclinical_NSCLC_IA2",
        Parameter_label == "nh_entry_nsclc_IA3_int" ~ "Preclinical_NSCLC_IA3",
        Parameter_label == "nh_entry_nsclc_IB_int"  ~ "Preclinical_NSCLC_IB",
        Parameter_label == "nh_entry_nsclc_II_int"  ~ "Preclinical_NSCLC_II",
        Parameter_label == "nh_entry_nsclc_III_int" ~ "Preclinical_NSCLC_III",
        Parameter_label == "nh_entry_nsclc_IV_int"  ~ "Preclinical_NSCLC_IV",
        Parameter_label == "nh_entry_sclc_lim_int"  ~ "Preclinical_SCLC_Limited",
        Parameter_label == "nh_entry_sclc_ext_int"  ~ "Preclinical_SCLC_Extensive"
      ),
      PLCO_coef = if_else(
        grepl("nsclc", Parameter_label),
        NSCLC_PLCO_coef,
        SCLC_PLCO_coef
      ),
      LP = Modelled_value + PLCO_coef * Cen_logit_PLCO,
      ExpLP = exp(LP)
    )
  
  Prob_data <- bind_rows(
    tibble(State = "No_cancer", LP = 0, ExpLP = 1), 
    Int_data %>% select(State, LP, ExpLP)) %>%
    mutate(Prob = ExpLP / sum(ExpLP),
           CumProb = cumsum(Prob))
  
  U <- runif(1)
  Start_state <- Prob_data$State[which(U <= Prob_data$CumProb)[1]]
  Start_state_code <- States_list[[Start_state]]
  
  return(as.numeric(Start_state_code))
}









NH_sample_next_state <- function(Current_state_patient, States_list, Event_table, 
                                 PLCO, Smoking_status, All_parameters) {
  
  
  if (Current_state == "No_cancer") {
    
    if (Smoking_status == "Former") {
      Beta <- 0.0470296
    } else if (Smoking_status == "Current") {
      Beta <- 0.1096189
    } else {
      Beta <- NA
    }
    
    Cen_logit_PLCO <- log(PLCO / (1 - PLCO)) + 3.66649
    
    # Time to NSCLC incidence
    M_NSCLC <- All_parameters %>% filter(Parameter_label == "nh_post_inc_nsclc_m") %>% pull(Modelled_value)
    C_NSCLC <- All_parameters %>% filter(Parameter_label == "nh_post_inc_nsclc_c") %>% pull(Modelled_value)
    Shape_NSCLC <- Beta * M_NSCLC
    Rate_NSCLC <- exp(C_NSCLC + (M_NSCLC * Cen_logit_PLCO))
    U <- runif(1)
    
    Time_to_NSCLC <- log(-(Shape_NSCLC / Rate_NSCLC) * log(U) + 1) / Shape_NSCLC
    
    # Time to SCLC incidence
    M_SCLC <- All_parameters %>% filter(Parameter_label == "nh_post_inc_sclc_m") %>% pull(Modelled_value)
    C_SCLC <- All_parameters %>% filter(Parameter_label == "nh_post_inc_sclc_c") %>% pull(Modelled_value)
    Shape_SCLC <- Beta * M_SCLC
    Rate_SCLC <- exp(C_SCLC + (M_SCLC * Cen_logit_PLCO))
    U <- runif(1)
    
    Time_to_SCLC <- log(-(Shape_SCLC / Rate_SCLC) * log(U) + 1) / Shape_SCLC
    
    
    # Next state and time to state
    Time_to_next_state <- min(Time_to_NSCLC, Time_to_SCLC)
    Next_state <- c("Preclinical_NSCLC_IA1", "Preclinical_SCLC_limited")[which.min(c(Time_to_NSCLC, Time_to_SCLC))]
    
  } else if (Current_state %in% c("Preclinical_NSCLC_IA1", "Preclinical_NSCLC_IA2",
                                  "Preclinical_NSCLC_IA3", "Preclinical_NSCLC_IB",
                                  "Preclinical_NSCLC_II", "Preclinical_NSCLC_III",
                                  "Preclinical_NSCLC_IV", "Preclinical_SCLC_limited",
                                  "Preclinical_SCLC_extensive")) {
    
    Sigma <- All_parameters %>% filter(Parameter_label == "nh_prog_re_sigma") %>% pull(Modelled_value)
    NH_acc_factor <- rlnorm(1, meanlog = -0.5 * Sigma^2, sdlog = Sigma)
    
    Constant_rate <- All_parameters %>% filter(Parameter_label == "r_IF_NLST_mid_screening") %>% pull(Modelled_value)
    
    Current_events <- Event_table %>%
      filter(Current_state == Current_state_patient, !is.na(Parameter_label))
    
    Current_events <- Current_events %>%
      left_join(All_parameters %>% select(Parameter_label, Modelled_value),
                by = "Parameter_label") %>%
      mutate(Rate = case_when(
        startsWith(Next_state, "Preclinical") ~ Modelled_value * NH_acc_factor,
        startsWith(Next_state, "Clinical") ~ Modelled_value + Constant_rate,
        TRUE ~ Modelled_value)) 
    
    Current_events <- Current_events %>%
      #mutate(Time_to = rexp(n(), rate = Rate))
      mutate(U = runif(n()), Time_to = -log(U) / Rate)
    
    Next_state <- Current_events$Next_state[which.min(Current_events$Time_to)]
    Time_to_next_state <- min(Current_events$Time_to)
    
  } else {
    
    # Next state and time to state
    Time_to_next_state <- NA
    Next_state <- NA
    
  }
  
  
  return(list(Next_state = Next_state, Time = Time_to_next_state))
}








































###############################################################################################
# Old code


Generate_nh_preclinical_progression <- function(All_parameters, Event_table) {
  
  nh_prog_re_sigma <- All_parameters %>% filter(Parameter_label == "nh_prog_re_sigma") %>% pull(Modelled_value)
  NH_acc_factor <- rlnorm(1, meanlog = -0.5 * nh_prog_re_sigma^2, sdlog = nh_prog_re_sigma)
  
  Event_table_prog <- Event_table %>%
    filter(startsWith(Current_state, "Preclinical")) %>%
    filter(!(Next_state %in% c("Screening", "Other_cause_death"))) %>%
    mutate(Parameter_label = c(
      "nh_prog_lambda_1",
      "nh_prog_xi_1",
      "nh_prog_lambda_2",
      "nh_prog_xi_2",
      "nh_prog_lambda_3",
      "nh_prog_xi_3",
      "nh_prog_lambda_4",
      "nh_prog_xi_4",
      "nh_prog_lambda_5",
      "nh_prog_xi_5",
      "nh_prog_lambda_6",
      "nh_prog_xi_6",
      "nh_prog_xi_7",
      "nh_prog_lambda_7",
      "nh_prog_mu_1",
      "nh_prog_phi_1",
      "nh_prog_phi_2",
      "nh_prog_mu_2"
    )
    ) %>%
    left_join(
      All_parameters %>%
        select(Parameter_label, Modelled_value),
      by = "Parameter_label"
    )
  
}






Generate_nh_preclinical_progression <- function(All_parameters, Event_table, Current_state_patient) {
  
  nh_prog_re_sigma <- All_parameters %>% filter(Parameter_label == "nh_prog_re_sigma") %>% pull(Modelled_value)
  NH_acc_factor <- rlnorm(1, meanlog = -0.5 * nh_prog_re_sigma^2, sdlog = nh_prog_re_sigma)
  
  Current_events <- Event_table %>%
    filter(Current_state == Current_state_patient)
  
}




