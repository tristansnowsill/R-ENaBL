#########################################################################################################################################################
#########################################################################################################################################################
#
# 2: Parameters
#
#########################################################################################################################################################
#########################################################################################################################################################


setwd("~/HEPD/LCS/R-ENaBL")
All_parameters <- read.csv("Parameters.csv", header=TRUE)
setwd("~/HEPD/LCS/R-ENaBL/R")

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


States_list <- data.frame(
  No_cancer = 1,
  Preclinical_NSCLC_IA1 = 2,
  Preclinical_NSCLC_IA2 = 3,
  Preclinical_NSCLC_IA3 = 4,
  Preclinical_NSCLC_IB = 5,
  Preclinical_NSCLC_II = 6,
  Preclinical_NSCLC_III = 7,
  Preclinical_NSCLC_IV = 8,
  Preclinical_SCLC_limited = 9,
  Preclinical_SCLC_extensive = 10,
  Clinical_NSCLC_IA1 = 11,
  Clinical_NSCLC_IA2 = 12,
  Clinical_NSCLC_IA3 = 13,
  Clinical_NSCLC_IB = 14,
  Clinical_NSCLC_II = 15,
  Clinical_NSCLC_III = 16,
  Clinical_NSCLC_IV = 17,
  Clinical_SCLC_limited = 18,
  Clinical_SCLC_extensive = 19,
  Screening = 20,
  Other_cause_death = 21,
  Lung_cancer_death = 22)

State_lookup <- States_list %>%
  pivot_longer(
    cols = everything(),
    names_to = "State_name",
    values_to = "State_code")

# Event table
# ASSUMING NO 'RECURRENCE' AFTER DIAGNOSIS ATM

Event_table <- data.frame(
  Current_state = c(
    "No_cancer",
    "No_cancer",
    "No_cancer",
    "No_cancer",
    
    "Preclinical_NSCLC_IA1",
    "Preclinical_NSCLC_IA1",
    "Preclinical_NSCLC_IA1",
    "Preclinical_NSCLC_IA1",
    
    "Preclinical_NSCLC_IA2",
    "Preclinical_NSCLC_IA2",
    "Preclinical_NSCLC_IA2",
    "Preclinical_NSCLC_IA2",
    
    "Preclinical_NSCLC_IA3",
    "Preclinical_NSCLC_IA3",
    "Preclinical_NSCLC_IA3",
    "Preclinical_NSCLC_IA3",
    
    "Preclinical_NSCLC_IB",
    "Preclinical_NSCLC_IB",
    "Preclinical_NSCLC_IB",
    "Preclinical_NSCLC_IB",
    
    "Preclinical_NSCLC_II",
    "Preclinical_NSCLC_II",
    "Preclinical_NSCLC_II",
    "Preclinical_NSCLC_II",
    
    "Preclinical_NSCLC_III",
    "Preclinical_NSCLC_III",
    "Preclinical_NSCLC_III",
    "Preclinical_NSCLC_III",
    
    "Preclinical_NSCLC_IV",
    "Preclinical_NSCLC_IV",
    "Preclinical_NSCLC_IV",
    "Preclinical_NSCLC_IV",
    
    "Preclinical_SCLC_limited",
    "Preclinical_SCLC_limited",
    "Preclinical_SCLC_limited",
    "Preclinical_SCLC_limited",
    
    "Preclinical_SCLC_extensive",
    "Preclinical_SCLC_extensive",
    "Preclinical_SCLC_extensive",
    "Preclinical_SCLC_extensive",
    
    "Clinical_NSCLC_IA1",
    "Clinical_NSCLC_IA1",
    
    "Clinical_NSCLC_IA2",
    "Clinical_NSCLC_IA2",
    
    "Clinical_NSCLC_IA3",
    "Clinical_NSCLC_IA3",
    
    "Clinical_NSCLC_IB",
    "Clinical_NSCLC_IB",
    
    "Clinical_NSCLC_II",
    "Clinical_NSCLC_II",
    
    "Clinical_NSCLC_III",
    "Clinical_NSCLC_III",
    
    "Clinical_NSCLC_IV",
    "Clinical_NSCLC_IV",
    
    "Clinical_SCLC_limited",
    "Clinical_SCLC_limited",
    
    "Clinical_SCLC_extensive",
    "Clinical_SCLC_extensive",
    
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening",
    "Screening"
  ),
  
  Next_state = c(
    "Preclinical_NSCLC_IA1",
    "Preclinical_SCLC_limited",
    "Screening",
    "Other_cause_death",
    
    "Preclinical_NSCLC_IA2",
    "Clinical_NSCLC_IA1",
    "Screening",
    "Other_cause_death",
    
    "Preclinical_NSCLC_IA3",
    "Clinical_NSCLC_IA2",
    "Screening",
    "Other_cause_death",
    
    "Preclinical_NSCLC_IB",
    "Clinical_NSCLC_IA3",
    "Screening",
    "Other_cause_death",
    
    "Preclinical_NSCLC_II",
    "Clinical_NSCLC_IB",
    "Screening",
    "Other_cause_death",
    
    "Preclinical_NSCLC_III",
    "Clinical_NSCLC_II",
    "Screening",
    "Other_cause_death",
    
    "Preclinical_NSCLC_IV",
    "Clinical_NSCLC_III",
    "Screening",
    "Other_cause_death",
    
    "Clinical_NSCLC_IV",
    "Screening",
    "Other_cause_death",
    "Undiagnosed_cancer_death",
    
    "Preclinical_SCLC_extensive",
    "Clinical_SCLC_limited",
    "Screening",
    "Other_cause_death",
    
    "Clinical_SCLC_extensive",
    "Screening",
    "Other_cause_death",
    "Undiagnosed_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "Other_cause_death",
    "Lung_cancer_death",
    
    "No_cancer",
    "Preclinical_NSCLC_IA1",
    "Preclinical_NSCLC_IA2",
    "Preclinical_NSCLC_IA3",
    "Preclinical_NSCLC_IB",
    "Preclinical_NSCLC_II",
    "Preclinical_NSCLC_III",
    "Preclinical_NSCLC_IV",
    "Preclinical_SCLC_limited",
    "Preclinical_SCLC_extensive",
    "Clinical_NSCLC_IA1",
    "Clinical_NSCLC_IA2",
    "Clinical_NSCLC_IA3",
    "Clinical_NSCLC_IB",
    "Clinical_NSCLC_II",
    "Clinical_NSCLC_III",
    "Clinical_NSCLC_IV",
    "Clinical_SCLC_limited",
    "Clinical_SCLC_extensive"
  ),
  
  Parameter_label = c(
    NA,NA,NA,NA,
    
    "nh_prog_lambda_1", "nh_prog_xi_1", NA, NA,
    "nh_prog_lambda_2", "nh_prog_xi_2", NA, NA,
    "nh_prog_lambda_3", "nh_prog_xi_3", NA, NA,
    "nh_prog_lambda_4", "nh_prog_xi_4", NA, NA,
    "nh_prog_lambda_5", "nh_prog_xi_5", NA, NA,
    "nh_prog_lambda_6", "nh_prog_xi_6", NA, NA,
    "nh_prog_xi_7", NA, NA, "nh_prog_lambda_7", 
    "nh_prog_mu_1", "nh_prog_phi_1", NA, NA,
    "nh_prog_phi_2", NA, NA, "nh_prog_mu_2", 

    rep(NA, 37)
  )
)


Screening_list <- data.frame(
  Negative = 1,
  Positive = 2)

Follow_up_list <- data.frame(
  None = 1,
  Three = 2,
  Twelve = 3
)

Smoking_list <- c(
  Former = 1,
  Current = 2
)


Random_patients <- read.csv("Random_patients.csv", header=TRUE)
Random_events <- read.csv("Random_events.csv", header=TRUE)
Random_screening <- read.csv("Random_screening.csv", header=TRUE)


Risk_bands <- data.frame(
  Risk_band = c("<0.5%","0.5%-1.0%","1.0%-1.2%","1.2%-1.4%","1.4%-1.5%",">1.5%"),
  Lower = c(0, 0.005, 0.010, 0.012, 0.014, 0.015),
  Upper = c(0.005, 0.010, 0.012, 0.014, 0.015, 1)) %>%
  mutate(Lower_logit = qlogis(Lower),
         Upper_logit = qlogis(Upper)) %>%
  mutate(E0 = case_when(
    is.finite(Lower_logit) & is.finite(Upper_logit) ~
      (Lower_logit + Upper_logit) / 2,
    TRUE ~ NA_real_),
    V0 = case_when(
      is.finite(Lower_logit) & is.finite(Upper_logit) ~
        (Upper_logit - Lower_logit)^2 / 12,
      TRUE ~ NA_real_))


Starting_risk_moments <- data.frame(
  
  Wait_time = c(1,1,1,3,3,7,7,10,10,10),
  Smoking_status = c("Current","Current","Former","Current","Former","Current",
                     "Former","Current","Former","Former"),
  Starting_risk_band = c("1.2%-1.4%","1.4%-1.5%", "1.4%-1.5%",
                         "1.0%-1.2%","1.2%-1.4%",
                         "0.5%-1.0%", "1.0%-1.2%",
                         "<0.5%", "<0.5%", "0.5%-1.0%"),
  Start_lower = c(0.012,0.014,0.014,0.01,0.012,0.005,0.01,0.0025,0.0025,0.005),
  Start_upper = c(0.014,0.015,0.015,0.012,0.014,0.01,0.012,0.005,0.005,0.01)) %>%
  mutate(
    Lower_logit = qlogis(Start_lower),
    Upper_logit = qlogis(Start_upper),
    
    E0 = (Lower_logit + Upper_logit) / 2,
    
    V0 = (Upper_logit - Lower_logit)^2 / 12)





Transition_probabilities <- build_transition_probabilities(
  Starting_risk_moments = All_model_inputs$Starting_risk_moments,
  Beta_current = All_parameters %>% filter(Parameter_label == "Beta_current") %>% pull(Modelled_value),
  Beta_former = All_parameters %>% filter(Parameter_label == "Beta_former") %>% pull(Modelled_value),
  Sigma_eta = All_parameters %>% filter(Parameter_label == "Sigma_eta") %>% pull(Modelled_value),
  Risk_bands = All_model_inputs$Risk_bands
)

Starting_counts <- data.frame(
  
  Age = rep(55, length(Risk_band_levels)*2),
  Smoking_status = c(rep("Current", length(Risk_band_levels)), rep("Former", length(Risk_band_levels))),
  Risk_band = rep(Risk_band_levels, 2),
  Count = c(1000 * All_parameters %>% filter(Parameter_label == "p_plco_00_05_c") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_05_10_c") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_10_12_c") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_12_14_c") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_14_15_c") %>% pull(Modelled_value),
            1000 * 0.0441,
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_00_05_f") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_05_10_f") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_10_12_f") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_12_14_f") %>% pull(Modelled_value),
            1000 * All_parameters %>% filter(Parameter_label == "p_plco_14_15_f") %>% pull(Modelled_value),
            1000 * 0.0287)
)


Risk_in_results <- build_risk_in_tables(Starting_counts = All_model_inputs$Starting_counts,
                                        Transition_probabilities = All_model_inputs$Transition_probabilities)


Age_probabilities <- Risk_in_results$Risk_in_population %>%
  filter(Age > 55) %>%
  group_by(Age) %>%
  summarise(
    Total = sum(Count),
    .groups = "drop"
  ) 

Prop_age_in <- sum(Starting_counts$Count)/(sum(Starting_counts$Count) + sum(Age_probabilities$Total))


Coef_uptake <- data.frame(
  
  predictor = c(
    "age_60_64",
    "age_65_69",
    "age_70_74",
    "age_75_79",
    "current_smoker",
    "Q1_IMD",
    "Q2_IMD",
    "Q3_IMD",
    "Q4_IMD",
    "Q5_IMD",
    "plco_spl1",
    "plco_spl2",
    "plco_spl3",
    "female",
    "ethnicity_other",
    "intercept"
  ),
  
  IMD_Q1 = c(All_parameters %>% filter(Parameter_label == "Q1_60_64") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q1_65_69") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q1_70_74") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q1_75_79") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q1_Current_smoker") %>% pull(Modelled_value),
             0,0,0,0,0,
             All_parameters %>% filter(Parameter_label == "Q1_PLCO_spline_basis_1") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q1_PLCO_spline_basis_2") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q1_PLCO_spline_basis_3") %>% pull(Modelled_value),
             0,0,
             All_parameters %>% filter(Parameter_label == "Q1_cons") %>% pull(Modelled_value)),
  IMD_Q2 = c(All_parameters %>% filter(Parameter_label == "Q2_60_64") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q2_65_69") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q2_70_74") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q2_75_79") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q2_Current_smoker") %>% pull(Modelled_value),
             0,0,0,0,0,
             All_parameters %>% filter(Parameter_label == "Q2_PLCO_spline_basis_1") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q2_PLCO_spline_basis_2") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q2_PLCO_spline_basis_3") %>% pull(Modelled_value),
             0,0,
             All_parameters %>% filter(Parameter_label == "Q2_cons") %>% pull(Modelled_value)),
  IMD_Q3 = c(All_parameters %>% filter(Parameter_label == "Q3_60_64") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q3_65_69") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q3_70_74") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q3_75_79") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q3_Current_smoker") %>% pull(Modelled_value),
             0,0,0,0,0,
             All_parameters %>% filter(Parameter_label == "Q3_PLCO_spline_basis_1") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q3_PLCO_spline_basis_2") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q3_PLCO_spline_basis_3") %>% pull(Modelled_value),
             0,0,
             All_parameters %>% filter(Parameter_label == "Q3_cons") %>% pull(Modelled_value)),
  IMD_Q4 = c(All_parameters %>% filter(Parameter_label == "Q4_60_64") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q4_65_69") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q4_70_74") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q4_75_79") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q4_Current_smoker") %>% pull(Modelled_value),
             0,0,0,0,0,
             All_parameters %>% filter(Parameter_label == "Q4_PLCO_spline_basis_1") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q4_PLCO_spline_basis_2") %>% pull(Modelled_value),
             All_parameters %>% filter(Parameter_label == "Q4_PLCO_spline_basis_3") %>% pull(Modelled_value),
             0,0,
             All_parameters %>% filter(Parameter_label == "Q4_cons") %>% pull(Modelled_value)),
  IMD_Q5 = rep(0,16),
  
  White = c(All_parameters %>% filter(Parameter_label == "White_60_64") %>% pull(Modelled_value),
            All_parameters %>% filter(Parameter_label == "White_65_69") %>% pull(Modelled_value),
            All_parameters %>% filter(Parameter_label == "White_70_74") %>% pull(Modelled_value),
            All_parameters %>% filter(Parameter_label == "White_75_79") %>% pull(Modelled_value),
            All_parameters %>% filter(Parameter_label == "White_current_smoker") %>% pull(Modelled_value),
            0,
            All_parameters %>% filter(Parameter_label == "Q2_IMD") %>% pull(Modelled_value),
            All_parameters %>% filter(Parameter_label == "Q3_IMD") %>% pull(Modelled_value),
            All_parameters %>% filter(Parameter_label == "Q4_IMD") %>% pull(Modelled_value),
            All_parameters %>% filter(Parameter_label == "Q5_IMD") %>% pull(Modelled_value),
            0,0,0,0,0,
            All_parameters %>% filter(Parameter_label == "White_cons") %>% pull(Modelled_value)),
  Attend_lhc = c(0,
                 All_parameters %>% filter(Parameter_label == "LHC_age65_74") %>% pull(Modelled_value),
                 All_parameters %>% filter(Parameter_label == "LHC_age65_74") %>% pull(Modelled_value),
                 All_parameters %>% filter(Parameter_label == "LHC_age75") %>% pull(Modelled_value),
                 0,
                 All_parameters %>% filter(Parameter_label == "LHC_imdQ5") %>% pull(Modelled_value),
                 All_parameters %>% filter(Parameter_label == "LHC_imdQ4") %>% pull(Modelled_value),
                 All_parameters %>% filter(Parameter_label == "LHC_imdQ3") %>% pull(Modelled_value),
                 All_parameters %>% filter(Parameter_label == "LHC_imdQ2") %>% pull(Modelled_value),
                 0,0,0,0,
                 All_parameters %>% filter(Parameter_label == "LHC_genderFemale") %>% pull(Modelled_value),
                 All_parameters %>% filter(Parameter_label == "LHC_ethnicityOtherethnicities") %>% pull(Modelled_value),
                 All_parameters %>% filter(Parameter_label == "LHC_Intercept") %>% pull(Modelled_value)),
  Attend_ldct = c(0,
                  All_parameters %>% filter(Parameter_label == "LDCT_age65_74") %>% pull(Modelled_value),
                  All_parameters %>% filter(Parameter_label == "LDCT_age65_74") %>% pull(Modelled_value),
                  All_parameters %>% filter(Parameter_label == "LDCT_age75") %>% pull(Modelled_value),
                  0,
                  All_parameters %>% filter(Parameter_label == "LDCT_imdQ5") %>% pull(Modelled_value),
                  All_parameters %>% filter(Parameter_label == "LDCT_imdQ4") %>% pull(Modelled_value),
                  All_parameters %>% filter(Parameter_label == "LDCT_imdQ3") %>% pull(Modelled_value),
                  All_parameters %>% filter(Parameter_label == "LDCT_imdQ2") %>% pull(Modelled_value),
                  0,0,0,0,
                  All_parameters %>% filter(Parameter_label == "LDCT_genderFemale") %>% pull(Modelled_value),
                  All_parameters %>% filter(Parameter_label == "LDCT_ethnicityOtherethnicities") %>% pull(Modelled_value),
                  All_parameters %>% filter(Parameter_label == "LDCT_Intercept") %>% pull(Modelled_value)),
  Persist = c(0,0,0,0,
              All_parameters %>% filter(Parameter_label == "2ndLDCT_currentsmoker") %>% pull(Modelled_value),
              0,
              All_parameters %>% filter(Parameter_label == "2ndLDCT_imdQ4") %>% pull(Modelled_value),
              All_parameters %>% filter(Parameter_label == "2ndLDCT_imdQ3") %>% pull(Modelled_value),
              All_parameters %>% filter(Parameter_label == "2ndLDCT_imdQ2") %>% pull(Modelled_value),
              All_parameters %>% filter(Parameter_label == "2ndLDCT_imdQ1") %>% pull(Modelled_value),
              0,0,0,
              All_parameters %>% filter(Parameter_label == "2ndLDCT_genderFemale") %>% pull(Modelled_value),
              All_parameters %>% filter(Parameter_label == "2ndLDCT_ethnicityOtherethnicities") %>% pull(Modelled_value),
              All_parameters %>% filter(Parameter_label == "2ndLDCT_Intercept") %>% pull(Modelled_value))
)




All_model_inputs <- list(General_mortality_fits = General_mortality_fits,
                         States_list = States_list,
                         Event_table = Event_table,
                         Screening_list = Screening_list,
                         Smoking_list = Smoking_list,
                         Age_step = 0.5,
                         Random_patients = Random_patients,
                         Random_events = Random_events,
                         Random_screening = Random_screening,
                         Risk_bands = Risk_bands,
                         Starting_risk_moments = Starting_risk_moments,
                         Transition_probabilities = Transition_probabilities,
                         Starting_counts = Starting_counts,
                         Prop_age_in = Prop_age_in,
                         Coef_uptake = Coef_uptake)




run_model <- function(All_parametersM = All_parameters,
                      All_model_inputsM = All_model_inputs,
                      All_user_inputsM = All_user_inputs,
                      Mortality_comb_methodM = "SUM",
                      Patient_charsM = NULL,
                      Screening_strategyM = "No_screening",
                      Initial_time_to_screenM = 0.08,
                      Patient_numberM = 10,
                      Respond_or_joinM = TRUE,
                      Model_methodM = "Random",
                      Model_routeM = "Risk_in",
                      Screen_designM = "None") {

  
  All_user_inputs <- list(Mortality_comb_method=Mortality_comb_methodM, 
                          Patient_chars=Patient_charsM,
                          Screening_strategy=Screening_strategyM,
                          Initial_time_to_screen=Initial_time_to_screenM,
                          Patient_number=Patient_numberM,
                          Respond_or_join=Respond_or_joinM,
                          Model_method=Model_methodM,
                          Model_route=Model_routeM,
                          Screen_design=Screen_designM)
  
  
  
  # Run model
  Simmodel <- simmer("Lung_cancer_screening_DES")  
  
  Traj_main <- make_natural_history_traj(Env = Simmodel, 
                                         All_parameters = All_parameters, #All_costs, All_utilities, 
                                         All_model_inputs = All_model_inputs,
                                         All_user_inputs = All_user_inputs)
  
  
  Simmodel %>%
    add_generator(
      name_prefix = "patient",
      trajectory = Traj_main,
      distribution = at(rep(x = 0, times = All_user_inputs$Patient_number)),
      mon = 2)
  
  
  Simmodel %>% reset() %>% run()
  
  
  # Collect attributes
  Atts <- get_mon_attributes(Simmodel)
  
  # Useful keys
  Atts_results <- Atts %>% filter(key %in% c("Age", "Gender", "Smoking_status", "PLCO", "IMD",
                                             "OtherEth", "Attend_LHC", "Attend_LDCT", "Persist",
                                             "Age_of_death", "Time_of_death", "State", "Join_screening",
                                             # "Next_cancer_state", "Next_cancer_state_time",
                                             "Screening_no", "Screen_detected", "diagnosed"))#, "Diagnosed_stage"))
  
  # Just keep maximum cost for each patient
  Atts_results <- Atts_results %>%
    group_by(name, key) %>%
    arrange(name, time) %>%
    ungroup()
  
  Atts_labelled <- Atts_results %>%
    left_join(
      State_lookup,
      by = c("value" = "State_code")
    ) %>%
    mutate(
      Value_label = ifelse(
        key == "State",
        State_name,
        as.character(value)
      )
    ) %>% dplyr::select(-State_name, -replication, -value) %>%
    mutate(
      Value_label = case_when(
        key == "Gender" & Value_label == 1 ~ "Male",
        key == "Gender" & Value_label == 2 ~ "Female",
        key == "Smoking_status" & Value_label == 1 ~ "Former",
        key == "Smoking_status" & Value_label == 2 ~ "Current",
        key == "Screen_detected" & Value_label == 1 ~ "Negative",
        key == "Screen_detected" & Value_label == 2 ~ "Positive",
        key == "Join_screening" & Value_label == 0 ~ "No",
        key == "Join_screening" & Value_label == 1 ~ "Yes",
        key == "OtherEth" & Value_label == 0 ~ "White",
        key == "OtherEth" & Value_label == 1 ~ "Other",
        key == "Attend_LHC" & Value_label == 0 ~ "No",
        key == "Attend_LHC" & Value_label == 1 ~ "Yes",
        key == "Attend_LDCT" & Value_label == 0 ~ "No",
        key == "Attend_LDCT" & Value_label == 1 ~ "Yes",
        key == "Persist" & Value_label == 0 ~ "No",
        key == "Persist" & Value_label == 1 ~ "Yes",
        TRUE ~ as.character(Value_label)
      )
    )
  
  Patient_list <- list()
  
  for (i in unique(Atts_labelled$name)) {
    Patient_list[[i]] <- Atts_labelled %>%
      filter(name == i)
  }
  
  return(Patient_list)
}


