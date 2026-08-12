#########################################################################################################################################################
#########################################################################################################################################################
#
# 5: Screening
#
#########################################################################################################################################################
#########################################################################################################################################################


# Screening scenarios
#
# 1. Age at programme qualification
#     - 55 to 80
#     - 60 to 80
#     - 55 to 75
#     - 60 to 75
#
# 2. 6-year lung cancer risk
#     - 1.5%
#     - 2.5%
#     - 5%
#
# 3. Programme designs
#     - Single one-off LDCT screen at programme entry
#     - Three consecutive annual LDCT screens from programme entry
#     - Annual LDCT screens from entry to 80 years
#     - Biennial LDCT screens from entry to 80 years

# First screening will occur after N time (0.08 in excel)
# Single screening happens after this time then that is it
# Triple screening happens after this time and then at exactly 1 and 2 years - 
#     there may also be an option for triple every 2 years
# Annual screening happens after this time and then exactly every year to age 80
# Biennial screening happens after this time and then exactly every two years to age 80








# make_screening_traj

# Creates a simmer trajectory for a screening event. The trajectory applies screening costs 
# and temporary screening disutility, simulates the screening result, and branches to 
# either continued natural history or screen-detected diagnosis
#
# @param Env A simmer simulation environment
# @param All_parameters List containing model parameters
# @param All_costs List containing cost data frames
# @param All_utilities List containing utility data frames
#
# @return A simmer trajectory object

make_screening_traj <- function(Env, All_parameters, All_model_inputs,
                                All_user_inputs) { #All_costs, All_utilities) {
  
  traj_screening <- trajectory("screening") %>%
    
    # Set current state to screening
    set_attribute(key = "State", value = 20) %>% 
    
    # Total number of screens
    set_attribute(key = 'Screening_no', value = 1, mod = "+") %>%
    
    # Schedule next screening time
    set_attribute(keys = c('Next_screening_time', "Last_screening_time"), value = function() {
      
      Time_of_death <- get_attribute(.env = Env, key = "Time_of_death")
      Current_time <- now(.env = Env)
      Screening_no <- get_attribute(.env = Env, key = "Screening_no")
      Age_current <- get_attribute(.env = Env, key = "Age") #+ Current_time
      
      if (All_user_inputs$Screening_strategy == "Single_screen") {
        Next_screening_time <- Time_of_death + 1
        
      } else if (All_user_inputs$Screening_strategy == "Triple_screen") {
        if (Screening_no < 3) {
          Next_screening_time <- floor(Current_time) + 1
        } else {
          Next_screening_time <- Time_of_death + 1
        }
      } else if (All_user_inputs$Screening_strategy == "Annual_screen") {
        if (Age_current <= 79) {
          Next_screening_time <- floor(Current_time) + 1
        } else {
          Next_screening_time <- Time_of_death + 1
        }
      } else if (All_user_inputs$Screening_strategy == "Biennial_screen") {
        if (Age_current <= 78) {
          Next_screening_time <- floor(Current_time) + 2
        } else {
          Next_screening_time <- Time_of_death + 1
        }
      }
    
      Last_screening_time <- Current_time
      as.numeric(c(Next_screening_time, Last_screening_time))}) %>%
    
    # Add discounted screening test costs
    #set_attribute(key = "Total_cost", value = function() {
      
      #Current_time <- simmer::now(.env = Env)
      #Cost <- All_costs$Diagnosis_costs %>% filter(State == "Screening", Step == 1) %>% pull(Cost)
      #Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      #Cost <- disc_cost(Cost = Cost, Time_in_days = Current_time, Discount_rate = Disc_rate)
      
      #as.numeric(Cost)
    #}, mod = "+") %>%
    
    # Decide whether screening detects existing cancer
    set_attribute("Screen_detected", function() {
      
      Stage <- names(All_model_inputs$States_list[get_attribute(.env = Env, key = "Stage")])
      S_result <- screening_result(Stage = Stage, 
                                   Sensitivity = All_parameters %>% filter(Input_param_category == "LDCT sensitivity") %>% dplyr::select(Parameter_label, Modelled_value), 
                                   Specificity = All_parameters %>% filter(Parameter_label == "spec_LDCT") %>% dplyr::pull(Modelled_value),
                                   All_user_inputs = All_user_inputs,
                                   All_model_inputs = All_model_inputs,
                                   i = as.integer(sub("patient", "", get_name(.env = Env))) + 1,
                                   j = get_attribute(.env = Env, key = 'Screening_no'))
    
      as.numeric(All_model_inputs$Screening_list[S_result])
    }) %>%
    
    # QALYs during screening-result waiting period
    #set_attribute(keys = "Total_QALYs", values = function() {
      
      #Q_utility <- get_attribute(.env = Env, keys = "Current_utility") 
      #Q_time <- 7
      #Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      #QALYs <- disc_QALYs(Utility = Q_utility, Time_total = simmer::now(.env = Env), Current_time = Q_time, Discount_rate = Disc_rate)
      
      #as.numeric(QALYs)
    #}, mod = "+") %>%
    
    # BRANCH
    
    branch(option = function() get_attribute(.env = Env, "Screen_detected"), continue = TRUE,
           
           # 1. Screen negative: return to natural history (includes false negative)
           trajectory("screen_negative") %>%
             
             # Restore state to underlying cancer state
             set_attribute(key = "State", value = function() {
               
               Stage <- get_attribute(.env = Env, key = "Stage")
               State <- Stage 
               
               as.numeric(State)
             }) %>%
             
             simmer::rollback(target = "No_or_Undiag_state", times = Inf),
           
           # 2. Screen positive: proceed to diagnosis or rollback
           trajectory("screen_positive") %>%
             
             # Determine whether this is a true positive or a false positive
             set_attribute(key = "Diagnosed", value = function() {
               
               Stage <- names(All_model_inputs$States_list[get_attribute(.env = Env, key = "Stage")])
               if (Stage == 'No_cancer') {
                 Diagnosed <- 1
               } else {
                 Diagnosed <- 2
               }
               
               as.numeric(Diagnosed)
             }) %>%
             
             #set_attribute(keys = c("Follow_up", "Time_to_FU"), value = function() {
               
              # if (All_user_inputs$Model_method == "Fixed") {
               #  UFollow_up <- (All_model_inputs$Random_screening %>% 
                #                  filter(i_patient == as.integer(sub("patient", "", get_name(.env = Env))) + 1,
                 #                        j_no == get_attribute(.env = Env, key = 'Screening_no')) %>% 
                  #                dplyr::pull(U_Follow_up))
               #} else {
                # UFollow_up <- runif(1)
               #}
               
               #Prob_3m <- All_parameters %>% filter(Parameter_label == "pr_indeter_3m") %>% dplyr::pull(Modelled_value)
               #Prob_12m <- All_parameters %>% filter(Parameter_label == "pr_indeter_12m") %>% dplyr::pull(Modelled_value)
                 
               #if (UFollow_up < Prob_3m) {
                # Follow_up <- "Three"
                 #FU_time <- 3
               #} else if (UFollow_up < (Prob_3m + Prob_12m)) {
                # Follow_up <- "Twelve"
                 #FU_time <- 12
               #} else {
                # Follow_up <- "None"
                 #FU_time <- get_attribute(.env = Env, key = "Time_of_death") + 1
               #}
               
               #as.numeric(c(Follow_up_list[Follow_up], FU_time))
             #}) %>%

             
             # BRANCH
             
             branch(option = function() get_attribute(.env = Env, "Diagnosed"), continue = TRUE,
                    
                    # 1. False positive: return to no cancer state
                    trajectory("False_positive") %>%
                      set_attribute(key = "State", value = 1) %>%
                      simmer::rollback(target = "No_or_Undiag_state", times = Inf),
                    
                    # 2. True positive: follow screen-detected diagnosis pathway
                    make_diagnosis_traj(
                      Env = Env,
                      All_parameters = All_parameters,
                      All_model_inputs = All_model_inputs,
                      All_user_inputs = All_user_inputs,
                      Traj_name = "Traj_screen_diagnosed",
                      Pathway = as.numeric(All_model_inputs$States_list[sub("^Preclinical", "Clinical", 
                                             names(All_model_inputs$States_list[get_attribute(.env = Env, key = "Stage")]))]),
                      Route_detection = "Screening")
                    
             )
    )
  
  return(traj_screening)
}





screening_result <- function(Stage, Sensitivity, Specificity, All_user_inputs, 
                             All_model_inputs, i=NULL, j=NULL) {
  
  Valid_stages <- c("No_cancer", "Preclinical_NSCLC_IA1", "Preclinical_NSCLC_IA2",
                    "Preclinical_NSCLC_IA3", "Preclinical_NSCLC_IB", "Preclinical_NSCLC_II",
                    "Preclinical_NSCLC_III", "Preclinical_NSCLC_IV", "Preclinical_SCLC_limited",
                    "Preclinical_SCLC_extensive")
  if (!Stage %in% Valid_stages) {
    stop(paste("Stage must be one of:", paste(Valid_stages, collapse = ", ")), call. = FALSE)
  }
  
  Probabilities <- c(Sensitivity$Modelled_value, Specificity)
  if (any(!is.numeric(Probabilities)) || any(is.na(Probabilities)) ||
      any(Probabilities < 0) || any(Probabilities > 1)) {
    stop("Sensitivities and specificity must be probabilities between 0 and 1", call. = FALSE)
  }
  
  # Simulate screening result
  if (Stage == "No_cancer") {
    # False positive - extra cost
    #Test_positive <- rbinom(1, 1, (1 - Specificity))
    if (All_user_inputs$Model_method == "Fixed") {
      UFalse_positive <- (All_model_inputs$Random_screening %>% filter(i_patient == i,
                                                                       j_no == j) %>% 
                           dplyr::pull(U_False_positive))
    } else {
      UFalse_positive <- runif(1)
    }

    Test_positive <- as.integer(UFalse_positive < (1 - Specificity))
    
  } else if (Stage %in% c("Preclinical_NSCLC_IV", "Preclinical_SCLC_extensive")) {
    # Late stage always gives positive result
    Test_positive <- 1
    
  } else {
    # True positive with sensitivity
    Parameter <- sub("^Preclinical_", "sens_", Stage)
    Sensitivity_stage <- Sensitivity %>% filter(Parameter_label == Parameter) %>% pull(Modelled_value)
    #Test_positive <- rbinom(1, 1, Sensitivity_stage)
    
    if (All_user_inputs$Model_method == "Fixed") {
      UPositive <- (All_model_inputs$Random_screening %>% filter(i_patient == i,
                                                                 j_no == j) %>% 
                            dplyr::pull(U_Positive))
    } else {
      UPositive <- runif(1)
    }
    
    Test_positive <- as.integer(UPositive < Sensitivity_stage)
  }

  Result <- ifelse(Test_positive == 1, "Positive", "Negative")
  
  return(Result)
}





