#########################################################################################################################################################
#########################################################################################################################################################
#
# 3: Population - characteristics and mortality
#
#########################################################################################################################################################
#########################################################################################################################################################



# Sample one individual with 
# - Age at entry
# - Sex
# - Smoking status
# - PLCO risk
# - Health state at entry
# - Age of other-cause death



sample_plco_current <- function(Age, All_parameters, All_user_inputs, All_model_inputs, i=NULL) {
  
  PLCO_data <- All_parameters %>% filter(startsWith(Parameter_label, "plco"))
  
  Age_minus_60 <- PLCO_data %>% filter(Parameter_label == "plco_cs_age_minus_60") %>% pull(Modelled_value)
  Cons <- PLCO_data %>% filter(Parameter_label == "plco_cs_cons") %>% pull(Modelled_value)
  Alpha <- PLCO_data %>% filter(Parameter_label == "plco_cs_alpha") %>% pull(Modelled_value)
  Omega <- PLCO_data %>% filter(Parameter_label == "plco_cs_omega") %>% pull(Modelled_value)
  DF <- PLCO_data %>% filter(Parameter_label == "plco_cs_df") %>% pull(Modelled_value)
  
  if (All_user_inputs$Model_method == "Fixed") {
    UPLCO <- All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
      dplyr::pull(U_PLCO)
  } else {
    UPLCO <- runif(1)
  }
  
  U1 <- qt(UPLCO, df = DF + 1)
  X0 <- qt(UPLCO, df = DF)
  
  Threshold <- (Alpha * X0 * sqrt((DF + 1) / (DF + X0^2)))
  
  Z <- ifelse(U1 < Threshold, X0, -X0)
  
  LP <- Omega * Z + Cons + Age_minus_60 * (Age - 60)
  
  return(plogis(LP))
}





sample_plco_former <- function(Age, All_parameters, All_user_inputs, All_model_inputs, i=NULL) {
  
  PLCO_data <- All_parameters %>% filter(startsWith(Parameter_label, "plco"))
  
  Age_minus_60 <- PLCO_data %>% filter(Parameter_label == "plco_fs_age_minus_60") %>% pull(Modelled_value)
  Cons <- PLCO_data %>% filter(Parameter_label == "plco_fs_cons") %>% pull(Modelled_value)
  Alpha <- PLCO_data %>% filter(Parameter_label == "plco_fs_alpha") %>% pull(Modelled_value)
  Omega <- PLCO_data %>% filter(Parameter_label == "plco_fs_omega") %>% pull(Modelled_value)
  
  if (All_user_inputs$Model_method == "Fixed") {
    UPLCO <- All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
      dplyr::pull(U_PLCO)
  } else {
    UPLCO <- runif(1)
  }
  
  U1 <- qnorm(UPLCO)
  X0 <- qnorm(UPLCO)
  
  Threshold <- Alpha * X0
  
  Z <- ifelse(U1 < Threshold, X0, -X0)
  
  LP <- Omega * Z + Cons + Age_minus_60 * (Age - 60)
  
  return(plogis(LP))
}






sample_patient_characteristics <- function(All_parameters, All_user_inputs, All_model_inputs, Age_in=FALSE, i=NULL) {
  
  Pmale <- All_parameters %>% filter(Parameter_label == "p_male") %>% pull(Modelled_value)
  Age_mean <- All_parameters %>% filter(Parameter_label == "pop_age_mean") %>% pull(Modelled_value)
  Age_sd <- All_parameters %>% filter(Parameter_label == "pop_age_sd") %>% pull(Modelled_value)
  Age_LL <- All_parameters %>% filter(Parameter_label == "pop_age_LL") %>% pull(Modelled_value)
  Age_UL <- All_parameters %>% filter(Parameter_label == "pop_age_UL") %>% pull(Modelled_value)
  
  
  # Age
  if (Age_in == TRUE) {
    Patient_age <- 55.1
  } else {
    if (All_user_inputs$Model_method == "Fixed") {
      UAge <- Age_LL + (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
        dplyr::pull(U_Age)) * (Age_UL - Age_LL)
    } else {
      UAge <- runif(n = 1, min = Age_LL, max = Age_UL)
    }
    Patient_age <- round(qnorm(p = UAge, mean = Age_mean, sd = Age_sd), 2)
  }
  
  # Gender
  if (All_user_inputs$Model_method == "Fixed") {
    Ugender <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
       dplyr::pull(U_Gender))
  } else {
    Ugender <- runif(1)
  }
  Patient_gender <- if (Ugender < Pmale) {1} else {2}
  
  # Age gender score
  if (Patient_gender == 1 & Patient_age < 65) {
    Age_gen_score <- 1
  } else if (Patient_gender == 1 & Patient_age >= 65) {
    Age_gen_score  <- 2
  } else if (Patient_gender == 2 & Patient_age < 65) {
    Age_gen_score <- 3
  } else {
    Age_gen_score <- 4
  }
  
  # Smoking status -> 1:Former, 2:Current
  Smoking_status <- smoking_status(All_parameters = All_parameters,
                                   Age_gen_score = Age_gen_score,
                                   All_user_inputs = All_user_inputs,
                                   All_model_inputs = All_model_inputs,
                                   i = i) 
  
  # PLCO
  if (Smoking_status == 1) {
    PLCO <- sample_plco_former(Age = Patient_age, 
                               All_parameters = All_parameters, 
                               All_user_inputs = All_user_inputs,
                               All_model_inputs = All_model_inputs,
                               i = i)
  } else {
    PLCO <- sample_plco_current(Age = Patient_age, 
                                All_parameters = All_parameters,
                                All_user_inputs = All_user_inputs,
                                All_model_inputs = All_model_inputs,
                                i = i)
  }
  
  
  
  Comb_method <- All_user_inputs$Mortality_comb_method
  Age_step <- All_model_inputs$Age_step
  General_mortality_fits <- All_model_inputs$General_mortality_fits
  States_list <- All_model_inputs$States_list
  Event_table <- All_model_inputs$Event_table
  
  # Age/time to death from other cause
  OCM <- other_cause_mortality(All_parameters = All_parameters,
                               Comb_method = Comb_method, 
                               Age_patient = Patient_age,
                               Gender = Patient_gender, 
                               PLCO = PLCO,
                               Smoking_status = Smoking_status,
                               Age_step = Age_step, 
                               General_mortality_fits = General_mortality_fits,
                               All_user_inputs = All_user_inputs,
                               All_model_inputs = All_model_inputs,
                               i = i)
  
  UP_output <- Uptake_persist(All_parameters = All_parameters,
                              All_model_inputs = All_model_inputs,
                              Age_entry = Patient_age,
                              Smoking_status = Smoking_status-1,
                              PLCO = PLCO,
                              Gender = Patient_gender-1,
                              i = i)
  
  # Sample initial state
  Samp_state <- sample_initial_state(States_list = States_list, 
                                     PLCO = PLCO, 
                                     All_parameters = All_parameters,
                                     All_user_inputs = All_user_inputs,
                                     All_model_inputs = All_model_inputs,
                                     i = i)
  Start_state <- names(States_list[Samp_state])
  Stage <- Samp_state
  
  if (All_user_inputs$Respond_or_join) {
    Respond_or_join <- Respond_join(All_parameters = All_parameters,
                                 All_user_inputs = All_user_inputs,
                                 Age_patient = Patient_age,
                                 PLCO = PLCO,
                                 i=i)
  } else {
    Respond_or_join <- c(1,1,1)
  }
  
  # If initial state is no cancer, when would they get cancer and what type?
  Next_event <- NH_sample_next_state(Current_state_patient = Start_state, 
                                     States_list = States_list,
                                     Event_table = Event_table, 
                                     PLCO = PLCO,
                                     Smoking_status = Smoking_status,
                                     All_parameters = All_parameters,
                                     All_user_inputs = All_user_inputs,
                                     All_model_inputs = All_model_inputs,
                                     i = i,
                                     j = 1)
  
  Next_state <- Next_event$Next_state
  Next_state <- as.numeric(All_model_inputs$States_list[Next_state])
  Next_time_to_state <- Next_event$Time
  
  
  Patient_data <- c(Patient_age, Patient_gender,
                    Age_gen_score, Smoking_status,
                    PLCO, UP_output, OCM$Age_at_death, OCM$Time_to_death,
                    Samp_state, Stage, 1, Respond_or_join,
                    Next_state, Next_time_to_state, 0)
  
  
  return(Patient_data)
}


patient_characteristics <- function(All_parameters, All_model_inputs, All_user_inputs, i_patient1) {

  
  if (All_user_inputs$Model_route == "Risk_in") {
    
    if (All_user_inputs$Model_method == "Fixed") {
      UAgein <- (All_model_inputs$Random_patients %>% filter(i_patient == i_patient1) %>% 
                   dplyr::pull(U_Age_in))
    } else {
      UAgein <- runif(1)
    }
    
    if (UAgein < All_model_inputs$Prop_age_in) {
      Patient_data <- sample_patient_characteristics(All_parameters = All_parameters,
                                                     All_user_inputs = All_user_inputs,
                                                     All_model_inputs = All_model_inputs,
                                                     Age_in = TRUE,
                                                     i = i_patient1)
    } else {
      Patient_data <- sample_patient_risk_in(All_parameters = All_parameters, 
                                             All_user_inputs = All_user_inputs,
                                             All_model_inputs = All_model_inputs,
                                             i = i_patient1)
    }
    
    
  } else {
    
    # Sample age, gender, smoking status, PLCO, etc.
    Patient_data <- sample_patient_characteristics(All_parameters = All_parameters,
                                                   All_user_inputs = All_user_inputs,
                                                   All_model_inputs = All_model_inputs,
                                                   i = i_patient1)
  }
  
  
  return(Patient_data)
}


smoking_status <- function(All_parameters, Age_gen_score, All_user_inputs, All_model_inputs, i=NULL) {
  
  Smoking_data <- All_parameters %>% filter(startsWith(Parameter_label, "p_smk"))
  
  P_1_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_m_55-64") %>% pull(Modelled_value)
  P_1_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_m_55-64") %>% pull(Modelled_value)
  P_2_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_m_65+") %>% pull(Modelled_value)
  P_2_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_m_65+") %>% pull(Modelled_value)
  P_3_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_f_55-64") %>% pull(Modelled_value)
  P_3_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_f_55-64") %>% pull(Modelled_value)
  P_4_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_f_65+") %>% pull(Modelled_value)
  P_4_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_f_65+") %>% pull(Modelled_value)
  
  Smoking_status_df <- data.frame(
    Age_gender_score = c(1,1,2,2,3,3,4,4),
    Smoker = c("Former", "Current", "Former", "Current", 
               "Former", "Current", "Former", "Current"),
    Prob = c(P_1_F/(P_1_C+P_1_F), P_1_C/(P_1_C+P_1_F), 
             P_2_F/(P_2_C+P_2_F), P_2_C/(P_2_C+P_2_F), 
             P_3_F/(P_3_C+P_3_F), P_3_C/(P_3_C+P_3_F), 
             P_4_F/(P_4_C+P_4_F), P_4_C/(P_4_C+P_4_F)))
  

  Smoking_status_Prob <- Smoking_status_df %>% filter(Age_gender_score == Age_gen_score, Smoker == 'Former') %>% pull(Prob)
  
  if (All_user_inputs$Model_method == "Fixed") {
    Usmoking <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                  dplyr::pull(U_Smoking))
  } else {
    Usmoking <- runif(1)
  }
  
  Smoking_status <- if (Usmoking < Smoking_status_Prob) {1} else {2}
  
  return(Smoking_status)
}





other_cause_mortality <- function(All_parameters, Comb_method, Age_patient, Gender, PLCO,
                                  Smoking_status, Age_step, General_mortality_fits,
                                  All_user_inputs, All_model_inputs, i=NULL) {
  
  Cen_logit_PLCO <- log(PLCO / (1 - PLCO)) + 3.66649
  
  # NLST component
  NLST_component <- All_parameters %>% filter(startsWith(Parameter_label, "ocm_coef")) %>%
                dplyr::select(Parameter_label, Modelled_value)
  #NLST_component <- data.frame(
   # Label = c("OCM_coef_intercept", "OCM_coef_age", "OCM_coef_plco",
    #          "OCM_coef_curr_smok", "OCM_coef_plco_curr_smok"),

    #Coefficients = c(-5.802614, 0.045320, 0.520297, 0.312519, -0.166465),
    #Value = c(1, Age_patient-61.41793, Cen_logit_PLCO, 1, Cen_logit_PLCO)
  #)
  Smoking <- Smoking_status - 1

  NLST_component <- NLST_component %>% mutate(
    #Value = c(1, (Age_patient-61.41793), Cen_logit_PLCO, Smoking, (Cen_logit_PLCO*Smoking)),
    Value = c(1, 10.582, -30.882, 1, -30.882),
    LP = Modelled_value * Value
  )
  
  Lambda <- exp(sum(NLST_component$LP))
  Gamma <- 1.4166
  
  
  # General mortality component
  # Either Fit_male or Fit_female
  if (Gender == 1) {
    Fit_gm <- General_mortality_fits$Fit_male
  } else if (Gender == 2) {
    Fit_gm <- General_mortality_fits$Fit_female
  } else {
    Fit_gm <- NA
  }
  Fit_gm_coef <- coef(Fit_gm)
  
  # Other cause mortality
  OCM <- data.frame(
    Age = seq(Age_patient, 110, Age_step)
  )
  
  OCM <- OCM %>% mutate(NLST = Lambda * Gamma * (Age - Age_patient)^(Gamma - 1),
                        General_pop = exp(Fit_gm_coef[1] + 
                                            Fit_gm_coef[2] * (Age - 0.5))
  ) %>%
    mutate(Combination = if (Comb_method == "MAX") {
      pmax(NLST, General_pop)
    } else {
      NLST + General_pop
    }) %>%
    mutate(
      Hazard_increment = 0.5 * (Combination + lag(Combination)) * Age_step,
      Hazard_increment = if_else(is.na(Hazard_increment), 0, Hazard_increment),
      Cum_hazard = cumsum(Hazard_increment)
    )
  
  if (All_user_inputs$Model_method == "Fixed") {
    UOCM <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                   dplyr::pull(U_OCM))
  } else {
    UOCM <- runif(1)
  }
  
  Transformed <- -log(UOCM)
  
  Index <- max(which(OCM$Cum_hazard <= Transformed))
  L1 <- OCM$Cum_hazard[Index]
  L2 <- if (Index < nrow(OCM)) {
    OCM$Cum_hazard[Index + 1]
  } else {
    L1
  }
  
  Age_at_death <- Age_patient + Age_step * (Index - 1) +
    ifelse(L2 > L1, Age_step * ((Transformed - L1) / (L2 - L1)), 0)
  
  Time_to_death <- Age_at_death - Age_patient
  
  return(list(Age_at_death=Age_at_death, Time_to_death=Time_to_death))
}







lung_cancer_mortality <- function(Age, Stage, Route, All_parameters, All_user_inputs, 
                                  All_model_inputs, i=NULL) {
  

  if (startsWith(Stage, "Clinical_NSCLC")) {
    # Extract model parameters
    Intercept <- All_parameters %>% filter(Parameter_label == "s_nsclc_cons") %>% pull(Modelled_value) 
    Age_coefficient <- All_parameters %>% filter(Parameter_label == "s_nsclc_age_minus_64") %>% pull(Modelled_value)
    Log_sigma <- All_parameters %>% filter(Parameter_label == "s_nsclc_ln_sigma") %>% pull(Modelled_value)
    Kappa <- All_parameters %>% filter(Parameter_label == "s_nsclc_kappa") %>% pull(Modelled_value)
    
    # Select the route coefficient
    Route_labels <- c(Screening = "s_nsclc_screening",
                      Interval = "s_nsclc_inter_canc",
                      Outside = "s_nsclc_post_scrn_canc")
    
    Route_coefficient <- All_parameters %>% filter(Parameter_label == Route_labels[Route]) %>% pull(Modelled_value)
    
    # Select the stage coefficient
    Stage_labels <- c(Clinical_NSCLC_IA1 = "s_nsclc_stage IA1",
                      Clinical_NSCLC_IA2 = "s_nsclc_stage IA2",
                      Clinical_NSCLC_IA3 = "s_nsclc_stage IA3",
                      Clinical_NSCLC_IB = "s_nsclc_stage IB",
                      Clinical_NSCLC_II = "s_nsclc_stage II",
                      Clinical_NSCLC_III = "s_nsclc_stage III",
                      Clinical_NSCLC_IV = "s_nsclc_stage IV")
    
    Stage_coefficient <- All_parameters %>% filter(Parameter_label == Stage_labels[Stage]) %>% pull(Modelled_value)
    
  } else {
    
    # Extract model parameters
    Intercept <- All_parameters %>% filter(Parameter_label == "s_sclc_cons") %>% pull(Modelled_value) 
    Age_coefficient <- All_parameters %>% filter(Parameter_label == "s_sclc_age_minus_64") %>% pull(Modelled_value)
    Log_sigma <- All_parameters %>% filter(Parameter_label == "s_sclc_ln_sigma") %>% pull(Modelled_value)
    Kappa <- All_parameters %>% filter(Parameter_label == "s_sclc_kappa") %>% pull(Modelled_value)
    
    # Select the route coefficient
    Route_labels <- c(Screening = "s_sclc_screening",
                      Interval = "s_sclc_inter_canc",
                      Outside = "s_sclc_post_scrn_canc")
    
    Route_coefficient <- All_parameters %>% filter(Parameter_label == Route_labels[Route]) %>% pull(Modelled_value)
    
    # Select the stage coefficient
    Stage_labels <- c(Clinical_SCLC_limited = "s_sclc_limited",
                      Clinical_SCLC_extensive = "s_sclc_extensive")
    
    Stage_coefficient <- All_parameters %>% filter(Parameter_label == Stage_labels[Stage]) %>% pull(Modelled_value)
    
  }
  
  # Linear predictor
  Mu <- Intercept + Route_coefficient + Stage_coefficient + Age_coefficient * (Age - 64)
  
  Sigma <- exp(Log_sigma)
  
  # Generalised gamma calculations
  Gamma_shape <- abs(Kappa)^(-2)
  if (All_user_inputs$Model_method == "Fixed") {
    ULCM <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                       dplyr::pull(U_LCM))
  } else {
    ULCM <- runif(1)
  }
  Log_normal_k <- exp(Mu + Sigma * qnorm(ULCM))
  Gamma_probability <- ULCM #if (Kappa > 0) {1 - ULCM} else {ULCM}
  Gamma_untransformed <- qgamma(p = Gamma_probability, shape = Gamma_shape, scale = 1)
  #Gamma_transformed <- exp(Mu + (Sigma / Kappa) * log(Kappa^2 * Gamma_untransformed))
  Gamma_transformed <- exp(Mu + (Sigma / Kappa) * log(Gamma_untransformed / Gamma_shape))
  Survival_time <- if (Kappa == 0) {Log_normal_k} else {Gamma_transformed}
  
  
  return(Survival_time)
}








sample_patient_risk_in <- function(All_parameters, All_user_inputs, All_model_inputs,
                                   i = i_patient) {

  
  Risk_in_results <- build_risk_in_tables(Starting_counts = All_model_inputs$Starting_counts,
                                          Transition_probabilities = All_model_inputs$Transition_probabilities)
  
  
  Age_probabilities <- Risk_in_results$Risk_in_population %>%
    filter(Age > 55) %>%
    group_by(Age) %>%
    summarise(
      Total = sum(Count),
      .groups = "drop"
    ) %>%
    arrange(Age) %>%
    mutate(
      Probability = Total / sum(Total),
      Cumulative_probability = lag(
        cumsum(Probability),
        default = 0
      )
    )
  
  if (All_user_inputs$Model_method == "Fixed") {
    UAge <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                       dplyr::pull(U_Age))
  } else {
    UAge <- runif(1)
  }
  
  Sampled_age <- Age_probabilities %>%
    filter(Cumulative_probability >= UAge) %>%
    slice(1) %>%
    pull(Age)
  
  if (All_user_inputs$Model_method == "Fixed") {
    UComponent <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
               dplyr::pull(U_Component))
  } else {
    UComponent <- runif(1)
  }
  
  Selected_component <- Risk_in_results$Previous_state_cumulative %>%
    filter(Target_age == Sampled_age, Cumulative_probability_mass <= UComponent) %>%
    slice_tail(n = 1)
  
  Smoking_status <- as.character(Selected_component$Smoking_status)
  Smoking_status2 <- as.integer(Smoking_status == "Current") + 1
  
  Previous_risk_band <- as.character(Selected_component$Previous_risk_band)
  
  if (All_user_inputs$Model_method == "Fixed") {
    URisk_sample <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                     dplyr::pull(U_Risk_sample))
  } else {
    URisk_sample <- runif(1)
  }
  
  PLCO <- 1 / (1 + exp(qnorm(p = URisk_sample, mean = -Selected_component$E1, sd = sqrt(Selected_component$V1))))
  
  Gender <- sample_gender(All_parameters = All_parameters, 
                          Sampled_age = Sampled_age, 
                          Smoking_status = Smoking_status,
                          PLCO = PLCO) 
  
  # Age gender score
  if (Gender == 1 & Sampled_age < 65) {
    Age_gen_score <- 1
  } else if (Gender == 1 & Sampled_age >= 65) {
    Age_gen_score  <- 2
  } else if (Gender == 2 & Sampled_age < 65) {
    Age_gen_score <- 3
  } else {
    Age_gen_score <- 4
  }
  
  Comb_method <- All_user_inputs$Mortality_comb_method
  Age_step <- All_model_inputs$Age_step
  General_mortality_fits <- All_model_inputs$General_mortality_fits
  States_list <- All_model_inputs$States_list
  Event_table <- All_model_inputs$Event_table
  
  # Age/time to death from other cause
  OCM <- other_cause_mortality(All_parameters = All_parameters,
                               Comb_method = Comb_method, 
                               Age_patient = Sampled_age,
                               Gender = Gender, 
                               PLCO = PLCO,
                               Smoking_status = Smoking_status2,
                               Age_step = Age_step, 
                               General_mortality_fits = General_mortality_fits,
                               All_user_inputs = All_user_inputs,
                               All_model_inputs = All_model_inputs,
                               i = i)
  
  UP_output <- Uptake_persist(All_parameters = All_parameters,
                              All_model_inputs = All_model_inputs,
                              Age_entry = Sampled_age,
                              Smoking_status = as.integer(Smoking_status == "Current"),
                              PLCO = PLCO,
                              Gender = Gender-1,
                              i = i)
  
  # Sample initial state
  Samp_state <- sample_initial_state(States_list = States_list, 
                                     PLCO = PLCO, 
                                     All_parameters = All_parameters,
                                     All_user_inputs = All_user_inputs,
                                     All_model_inputs = All_model_inputs,
                                     i = i)
  Start_state <- names(States_list[Samp_state])
  Stage <- Samp_state
  
  if (All_user_inputs$Respond_or_join) {
    Respond_or_join <- Respond_join(All_parameters = All_parameters,
                                    All_user_inputs = All_user_inputs,
                                    Age_patient = Patient_age,
                                    PLCO = PLCO,
                                    i=i)
  } else {
    Respond_or_join <- c(1,1,1)
  }
  
  # If initial state is no cancer, when would they get cancer and what type?
  Next_event <- NH_sample_next_state(Current_state_patient = Start_state, 
                                     States_list = States_list,
                                     Event_table = Event_table, 
                                     PLCO = PLCO,
                                     Smoking_status = Smoking_status2,
                                     All_parameters = All_parameters,
                                     All_user_inputs = All_user_inputs,
                                     All_model_inputs = All_model_inputs,
                                     i = i,
                                     j = 1)
  
  Next_state <- Next_event$Next_state
  Next_state <- as.numeric(All_model_inputs$States_list[Next_state])
  Next_time_to_state <- Next_event$Time
  
  
  Patient_data <- c(Sampled_age, Gender, Age_gen_score, Smoking_status2,
                    PLCO, UP_output, OCM$Age_at_death, OCM$Time_to_death,
                    Samp_state, Stage, 1, Respond_or_join, Next_state, Next_time_to_state, 0)
  
  
  return(Patient_data)
}








calculate_transition_row <- function(Mean_logit, Variance_logit, Risk_bands) {

  SD_logit <- sqrt(Variance_logit)
  
  Result <- Risk_bands %>%
    transmute(
      Future_risk_band = Risk_band,
      Probability = pnorm(Upper_logit, mean = Mean_logit, sd = SD_logit) -
                    pnorm(Lower_logit, mean = Mean_logit, sd = SD_logit))
  
  return(Result)
}




build_transition_probabilities <- function(Starting_risk_moments, Beta_current,
                                           Beta_former, Sigma_eta, Risk_bands) {
  
  Starting_risk_moments %>% mutate(Smoking_effect = case_when(
        Smoking_status == "Current" ~ Beta_current,
        Smoking_status == "Former"  ~ Beta_former),
      
        E1 = E0 + Smoking_effect * Wait_time,
        V1 = V0 + Sigma_eta * sqrt(Wait_time)) %>%
    rowwise() %>% 
    mutate(Transition_row = list(calculate_transition_row(Mean_logit = E1,
                                                          Variance_logit = V1,
                                                          Risk_bands = Risk_bands))) %>%
    ungroup() %>%
    unnest(Transition_row) %>%
    dplyr::select(Wait_time,Smoking_status, Starting_risk_band, Start_lower, 
                  Start_upper, E0, V0, Smoking_effect, E1, V1, Future_risk_band, Probability)
}




build_risk_in_tables <- function(Starting_counts, Transition_probabilities, Minimum_age = 55,
                                 Maximum_age = 74, Risk_band_levels = c("<0.5%", "0.5%-1.0%",
                                    "1.0%-1.2%", "1.2%-1.4%", "1.4%-1.5%", ">1.5%")) {
  
  Smoking_levels <- c("Current", "Former")
  
  Risk_in_population <- crossing(Age = Minimum_age:Maximum_age,
                                 Smoking_status = factor(Smoking_levels,levels = Smoking_levels),
                                 Risk_band = factor(Risk_band_levels,levels = Risk_band_levels)) %>%
    left_join(Starting_counts %>% mutate(
          Smoking_status = factor(Smoking_status, levels = Smoking_levels),
          Risk_band = factor(Risk_band, levels = Risk_band_levels)),
          by = c("Age", "Smoking_status", "Risk_band")) %>%
    mutate(Count = coalesce(Count, 0)) %>%
    arrange(Age, Smoking_status, Risk_band)
  
  # Store the previous-state masses for each target age
  Previous_state_mass_list <- list()
  
  for (Target_age in seq((Minimum_age + 1), to = Maximum_age)) {
    
    Age_contributions <- Transition_probabilities %>%
      mutate(Target_age = Target_age,
             Source_age = Target_age - Wait_time) %>%
      filter(Source_age >= Minimum_age) %>%
      left_join(Risk_in_population %>%
          transmute(Source_age = Age, Smoking_status, Starting_risk_band = Risk_band,
            Source_count = Count), by = c("Source_age", "Smoking_status", "Starting_risk_band")) %>%
      mutate(Source_count = coalesce(Source_count, 0),
             Contribution = Source_count * Probability)
    
    # Sum the contributions to get the population at the target age
    New_age_population <- Age_contributions %>%
      group_by(Age = Target_age, Smoking_status, Risk_band = Future_risk_band) %>%
      summarise(Count = sum(Contribution), .groups = "drop") %>%
      mutate(Smoking_status = factor(Smoking_status, levels = Smoking_levels),
             Risk_band = factor(Risk_band,levels = Risk_band_levels))
    
    Risk_in_population <- Risk_in_population %>% rows_update(New_age_population,
        by = c("Age", "Smoking_status", "Risk_band"))
    
    # Sum over future risk bands to recover the source-state mass
    Previous_state_mass_age <- Age_contributions %>%
      group_by(Target_age, Wait_time, Source_age, Smoking_status,
        Previous_risk_band = Starting_risk_band, E1, V1) %>%
      summarise(Source_mass = sum(Contribution), .groups = "drop")
    
    Previous_state_mass_list[[as.character(Target_age)]] <- Previous_state_mass_age
  }
  
  Previous_state_mass_raw <- bind_rows(Previous_state_mass_list)
  
  Previous_state_grid <- crossing(Target_age = seq((Minimum_age + 1), Maximum_age),
    Smoking_status = factor(Smoking_levels, levels = Smoking_levels),
    Previous_risk_band = factor(Risk_band_levels, levels = Risk_band_levels))

  Previous_state_mass <- Previous_state_grid %>% 
    left_join(Previous_state_mass_raw %>%
        mutate(Smoking_status = factor(Smoking_status, levels = Smoking_levels),
          Previous_risk_band = factor(Previous_risk_band, levels = Risk_band_levels)),
        by = c("Target_age", "Smoking_status", "Previous_risk_band")) %>%
    mutate(Source_mass = coalesce(Source_mass, 0)) %>%
    group_by(Target_age) %>% mutate(Total_mass = sum(Source_mass),
             Probability_mass = if_else(Total_mass > 0, Source_mass / Total_mass,0)) %>%
    arrange(Smoking_status, Previous_risk_band, .by_group = TRUE) %>%
    mutate(Sampled_col = row_number()) %>%
    ungroup()
  
  Previous_state_cumulative <- Previous_state_mass %>%
    group_by(Target_age) %>%
    arrange(Smoking_status, Previous_risk_band, .by_group = TRUE) %>%
    mutate(Cumulative_probability_mass = lag(cumsum(Probability_mass), default = 0)) %>%
    ungroup()

  return(list(Risk_in_population = Risk_in_population,
              Previous_state_mass = Previous_state_mass,
              Previous_state_cumulative = Previous_state_cumulative))
}




sample_gender <- function(All_parameters, Sampled_age, Smoking_status, PLCO) {
  
  Gender_knots <- c(All_parameters %>% filter(Parameter_label == "Knot_k1") %>% pull(Modelled_value),
                    All_parameters %>% filter(Parameter_label == "Knot_k2") %>% pull(Modelled_value),
                    All_parameters %>% filter(Parameter_label == "Knot_k3") %>% pull(Modelled_value),
                    All_parameters %>% filter(Parameter_label == "Knot_k4") %>% pull(Modelled_value))
  
  Logit_plco <- qlogis(PLCO)
  
  Spline_basis <- Hmisc::rcspline.eval(x = Logit_plco,
                                       knots = Gender_knots,
                                       inclx = TRUE, norm = 2)
  
  Plco_spl1 <- Spline_basis[1]
  Plco_spl2 <- Spline_basis[2]
  Plco_spl3 <- Spline_basis[3]
  
  Age_minus_60 <- Sampled_age - 60
  
  Current_smoker <- as.integer(Smoking_status == "Current")
  
  Eta <- All_parameters %>% filter(Parameter_label == "_cons") %>% pull(Modelled_value) +
    All_parameters %>% filter(Parameter_label == "age_minus_60") %>% pull(Modelled_value) * Age_minus_60 +
    All_parameters %>% filter(Parameter_label == "1.currentsmok") %>% pull(Modelled_value) * Current_smoker +
    All_parameters %>% filter(Parameter_label == "currentsmokage_minus_60") %>% pull(Modelled_value) * Current_smoker * Age_minus_60 +
    All_parameters %>% filter(Parameter_label == "plco_spl1") %>% pull(Modelled_value) * Plco_spl1 +
    All_parameters %>% filter(Parameter_label == "plco_spl2") %>% pull(Modelled_value) * Plco_spl2 +
    All_parameters %>% filter(Parameter_label == "plco_spl3") %>% pull(Modelled_value) * Plco_spl3
  
  P_female <- plogis(Eta)

  if (All_user_inputs$Model_method == "Fixed") {
    UGender <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                       dplyr::pull(U_Gender))
  } else {
    UGender <- runif(1)
  }
  
  Gender <- ifelse(UGender < P_female,2,1)
  
  return(Gender)
}





Uptake_persist <- function(All_parameters, All_model_inputs, Age_entry, Smoking_status, PLCO, Gender, i=NULL) {


  Gender_knots <- c(All_parameters %>% filter(Parameter_label == "Knot_k1") %>% pull(Modelled_value),
                    All_parameters %>% filter(Parameter_label == "Knot_k2") %>% pull(Modelled_value),
                    All_parameters %>% filter(Parameter_label == "Knot_k3") %>% pull(Modelled_value),
                    All_parameters %>% filter(Parameter_label == "Knot_k4") %>% pull(Modelled_value))

  Logit_plco <- qlogis(PLCO)
  
  Spline_basis <- Hmisc::rcspline.eval(x = Logit_plco,
                                     knots = Gender_knots,
                                     inclx = TRUE, norm = 2)

  Plco_spl1 <- Spline_basis[1]
  Plco_spl2 <- Spline_basis[2]
  Plco_spl3 <- Spline_basis[3]

  Riskin_predictors <- c(
    age_60_64 = as.integer(Age_entry >= 60 & Age_entry < 65),
    age_65_69 = as.integer(Age_entry >= 65 & Age_entry < 70),
    age_70_74 = as.integer(Age_entry >= 70 & Age_entry < 75),
    age_75_79 = as.integer(Age_entry >= 75 & Age_entry < 80),
    
    Current_smoker = Smoking_status,
    
    Q1_IMD = 0,
    Q2_IMD = 0,
    Q3_IMD = 0,
    Q4_IMD = 0,
    Q5_IMD = 0,
    
    Plco_spl1 = Plco_spl1,
    Plco_spl2 = Plco_spl2,
    Plco_spl3 = Plco_spl3,
    
    Gender = Gender,
    ethnicity_other = 0,
    
    intercept = 1
  )
  
  
  Coef_matrix1 <- as.matrix(All_model_inputs$Coef_uptake[, (2:6)])
  Contributions1 <- Riskin_predictors * Coef_matrix1
  Log_odds1 <- exp(apply(Contributions1, 2, sum))
  Probs1 <- Log_odds1/sum(Log_odds1) 
  IMD_cumulative <- c(0, cumsum(Probs1)[-length(Probs1)])
  
  if (All_user_inputs$Model_method == "Fixed") {
    UIMD <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
               dplyr::pull(U_IMD))
  } else {
    UIMD <- runif(1)
  }
  
  IMD <- max(which(IMD_cumulative <= UIMD))
  Riskin_predictors[paste0("Q", IMD, "_IMD")] <- 1
  
  Coef_matrix2 <- as.matrix(All_model_inputs$Coef_uptake[, (7:10)])
  Contributions2 <- Riskin_predictors * Coef_matrix2
  Log_odds2 <- exp(-apply(Contributions2, 2, sum))
  Probs2 <- 1/(1 + Log_odds2)
  
  if (All_user_inputs$Model_method == "Fixed") {
    UWhite <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
               dplyr::pull(U_White))
    ULHC <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
               dplyr::pull(U_LHC))
    ULDCT <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
               dplyr::pull(U_LDCT))
    UPersist <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
               dplyr::pull(U_Persist))
  } else {
    UWhite <- runif(1)
    ULHC <- runif(1)
    ULDCT <- runif(1)
    UPersist <- runif(1)
  }
  OtherEth <- as.integer(UWhite > as.numeric(Probs2["White"]))
  Attend_lhc <- as.integer(ULHC < as.numeric(Probs2["Attend_lhc"]))
  Attend_ldct <- as.integer(ULDCT < as.numeric(Probs2["Attend_ldct"]))
  Persist <- as.integer(UPersist < as.numeric(Probs2["Persist"]))
  

  return(c(IMD, OtherEth, Attend_lhc, Attend_ldct, Persist))
}





Respond_join <- function(All_parameters, All_user_inputs, Age_patient, PLCO, i=NULL) {
  
  
  if (All_user_inputs$Model_method == "Fixed") {
    URespond <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                   dplyr::pull(U_Respond))
    UJoin <- (All_model_inputs$Random_patients %>% filter(i_patient == i) %>% 
                dplyr::pull(U_Join))
  } else {
    URespond <- runif(1)
    UJoin <- runif(1)
  }
  PRespond <- All_parameters %>% filter(Parameter_label == "p_respond") %>% pull(Modelled_value)
  PJoin <- All_parameters %>% filter(Parameter_label == "p_join") %>% pull(Modelled_value)
  
  Respond <- as.integer(URespond < PRespond)
  
  if (Respond == 1) {
    
    # Is the patient eligible for screening strategies?
    if (All_user_inputs$Screen_design == "Design_1"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 55 & Age_patient <= 80 & PLCO >= 0.015) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_2"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 60 & Age_patient <= 80 & PLCO >= 0.015) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_3"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 55 & Age_patient <= 75 & PLCO >= 0.015) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_4"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 60 & Age_patient <= 75 & PLCO >= 0.015) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_5"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 55 & Age_patient <= 80 & PLCO >= 0.025) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_6"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 60 & Age_patient <= 80 & PLCO >= 0.025) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_7"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 55 & Age_patient <= 75 & PLCO >= 0.025) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_8"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 60 & Age_patient <= 75 & PLCO >= 0.025) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_9"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 55 & Age_patient <= 80 & PLCO >= 0.05) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_10"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 60 & Age_patient <= 80 & PLCO >= 0.05) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_11"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 55 & Age_patient <= 75 & PLCO >= 0.05) {1} else {0}
    } else if (All_user_inputs$Screen_design == "Design_12"& All_user_inputs$Model_route != "Risk_in") {
      Eligibility <- if (Age_patient >= 60 & Age_patient <= 75 & PLCO >= 0.05) {1} else {0}
    } else {
      Eligibility <- 1
    }
    
    if (Eligibility == 1) {
      Join_screening <- as.integer(UJoin < PJoin)
    } else {
      Join_screening <- 0
    }
    
  } else {
    Eligibility <- 0
    Join_screening <- 0
  }
  
  
  return(c(Respond, Eligibility, Join_screening))
}











