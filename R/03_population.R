#########################################################################################################################################################
#########################################################################################################################################################
#
# 3: Population - characteristics
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



sample_plco_current <- function(Age, Params) {
  
  Age_minus_60 <- Params %>% filter(Parameter_label == "plco_cs_age_minus_60") %>% pull(Modelled_value)
  Cons <- Params %>% filter(Parameter_label == "plco_cs_cons") %>% pull(Modelled_value)
  Alpha <- Params %>% filter(Parameter_label == "plco_cs_alpha") %>% pull(Modelled_value)
  Omega <- Params %>% filter(Parameter_label == "plco_cs_omega") %>% pull(Modelled_value)
  DF <- Params %>% filter(Parameter_label == "plco_cs_df") %>% pull(Modelled_value)
  
  U1 <- qt(runif(1), df = DF + 1)
  X0 <- qt(runif(1), df = DF)
  
  Threshold <- (Alpha * X0 * sqrt((DF + 1) / (DF + X0^2)))
  
  Z <- ifelse(U1 < Threshold, X0, -X0)
  
  LP <- Omega * Z + Cons + Age_minus_60 * (Age - 60)
  
  return(plogis(LP))
}





sample_plco_former_never <- function(Age, Params, Smoking_status) {
  
  if (Smoking_status == 0) {
    Age_minus_60 <- Params %>% filter(Parameter_label == "plco_ns_age_minus_60") %>% pull(Modelled_value)
    Cons <- Params %>% filter(Parameter_label == "plco_ns_cons") %>% pull(Modelled_value)
    Alpha <- Params %>% filter(Parameter_label == "plco_ns_alpha") %>% pull(Modelled_value)
    Omega <- Params %>% filter(Parameter_label == "plco_ns_omega") %>% pull(Modelled_value)
  } else {
    Age_minus_60 <- Params %>% filter(Parameter_label == "plco_fs_age_minus_60") %>% pull(Modelled_value)
    Cons <- Params %>% filter(Parameter_label == "plco_fs_cons") %>% pull(Modelled_value)
    Alpha <- Params %>% filter(Parameter_label == "plco_fs_alpha") %>% pull(Modelled_value)
    Omega <- Params %>% filter(Parameter_label == "plco_fs_omega") %>% pull(Modelled_value)
  }
  
  U1 <- qnorm(runif(1))
  X0 <- qnorm(runif(1))
  
  Threshold <- Alpha * X0
  
  Z <- ifelse(U1 < Threshold, X0, -X0)
  
  LP <- Omega * Z + Cons + Age_minus_60 * (Age - 60)
  
  return(plogis(LP))
}






sample_patient_characteristics <- function(n, Pmale, Age_mean, Age_sd, Age_LL, Age_UL, Never_switch, All_parameters) {
  
  # Age
  u <- runif(n = n, min = Age_LL, max = Age_UL)
  Patient_age <- qnorm(p = u, mean = Age_mean, sd = Age_sd)
  
  # Gender
  Patient_gender <- sample(x = c(1,2), size = n, prob=c(Pmale, 1-Pmale))
  
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
  
  # Smoking status -> 1:Current, 2:Former, 0:Never
  Smoking_data <- All_parameters %>% filter(startsWith(Parameter_label, "p_smk"))
  Smoking_status_df <- smoking_status(Never_switch = Never_switch, Smoking_data = Smoking_data) 
  Smoking_status_Prob <- Smoking_status_df %>% filter(Age_gender_score == Age_gen_score) %>% pull(Prob)
  Smoking_status <- sample(x = c(1,2,0), size = 1, prob = Smoking_status_Prob)
  
  # PLCO
  PLCO_data <- All_parameters %>% filter(startsWith(Parameter_label, "plco"))
  if (Smoking_status == 1) {
    PLCO <- sample_plco_current(Age = Patient_age, Params = PLCO_data)
  } else {
    PLCO <- sample_plco_former_never(Age = Patient_age, Params = PLCO_data, Smoking_status = Smoking_status)
  }
  
  return(c(Patient_age, Patient_gender, Age_gen_score, Smoking_status, PLCO))
}







smoking_status <- function(Never_switch, Smoking_data) {
  
  P_1_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_m_55-64") %>% pull(Modelled_value)
  P_1_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_m_55-64") %>% pull(Modelled_value)
  P_2_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_m_65+") %>% pull(Modelled_value)
  P_2_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_m_65+") %>% pull(Modelled_value)
  P_3_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_f_55-64") %>% pull(Modelled_value)
  P_3_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_f_55-64") %>% pull(Modelled_value)
  P_4_C <- Smoking_data %>% filter(Parameter_label == "p_smk_c_f_65+") %>% pull(Modelled_value)
  P_4_F <- Smoking_data %>% filter(Parameter_label == "p_smk_f_f_65+") %>% pull(Modelled_value)
  
  Smoking_status <- data.frame(
    Age_gender_score = c(1,1,1,2,2,2,3,3,3,4,4,4),
    Smoker = c("Current", "Former", "Never",
               "Current", "Former", "Never",
               "Current", "Former", "Never",
               "Current", "Former", "Never"),
    Prob = c(P_1_C, P_1_F, (1 - P_1_C - P_1_F),
             P_2_C, P_2_F, (1 - P_2_C - P_2_F),
             P_3_C, P_3_F, (1 - P_3_C - P_3_F),
             P_4_C, P_4_F, (1 - P_4_C - P_4_F))
  )
  
  if (Never_switch == 0) {
    Smoking_status <- Smoking_status %>%
      group_by(Age_gender_score) %>%
      mutate(Prob = case_when(
        Smoker == "Current" ~ Prob[Smoker == "Current"] /
          sum(Prob[Smoker %in% c("Current", "Former")]),
        
        Smoker == "Former" ~ 1 -
          Prob[Smoker == "Current"] /
          sum(Prob[Smoker %in% c("Current", "Former")]),
        
        Smoker == "Never" ~ 0
      )
      ) %>%
      ungroup()
  }
  return(Smoking_status)
}



General_mortality_df <- data.frame(
  Age = seq(55, 100),
  Males = c(0.00477,0.00540,0.00587,0.00642,0.00695,0.00762,0.00835,0.00927,
            0.01023,0.01101,0.01210,0.01342,0.01454,0.01586,0.01744,0.01846,
            0.02048,0.02258,0.02583,0.02853,0.03190,0.03574,0.03961,0.04449,
            0.04928,0.05548,0.06193,0.06880,0.07835,0.08851,0.09937,0.11288,
            0.12606,0.14289,0.16199,0.17308,0.19612,0.21760,0.24099,0.27045,
            0.30179,0.33712,0.36166,0.39346,0.45109,0.47637),
  Females = c(0.00320,0.00355,0.00382,0.00424,0.00463,0.00507,0.00551,0.00628,
              0.00674,0.00729,0.00802,0.00861,0.00942,0.01038,0.01136,0.01252,
              0.01343,0.01536,0.01750,0.01935,0.02167,0.02451,0.02769,0.03145,
              0.03509,0.03922,0.04461,0.05019,0.05790,0.06605,0.07519,0.08670,
              0.09804,0.11242,0.12646,0.14304,0.16244,0.18204,0.20329,0.22811,
              0.25801,0.28753,0.31021,0.34518,0.37451,0.42623)
)


# General mortality models

Fit_male <- lm(log(Males) ~ Age, data = General_mortality_df)
#coef(Fit_male)
Fit_female <- lm(log(Females) ~ Age, data = General_mortality_df)
#coef(Fit_female)

General_mortality_fits <- list(Fit_male = Fit_male, 
                               Fit_female = Fit_female)




other_cause_mortality <- function(Comb_method, Age_patient, Gender, PLCO,
                                  Age_step, General_mortality_fits) {
  
  # CHECK OCM_COEF_CURR_SMOK VALUE
  Cen_logit_PLCO <- log(PLCO / (1 - PLCO)) + 3.66649
  
  # NLST component
  NLST_component <- data.frame(
    Label = c("OCM_coef_intercept", "OCM_coef_age", "OCM_coef_plco",
              "OCM_coef_curr_smok", "OCM_coef_plco_curr_smok"),
    # WILL CHANGE TO PARAMETER VALUES
    Coefficients = c(-5.802614, 0.045320, 0.520297, 0.312519, -0.166465),
    Value = c(1, Age_patient-61.41793, Cen_logit_PLCO, 1, Cen_logit_PLCO)
  )
  
  NLST_component <- NLST_component %>% mutate(
    LP = Coefficients * Value
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
                                            Fit_gm_coef[2] * (Age_patient - 0.5))
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
  
  U <- runif(1)
  Transformed <- -log(U)
  
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
  
  return(c(Age_at_death, Time_to_death))
}


