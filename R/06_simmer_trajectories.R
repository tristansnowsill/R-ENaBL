#########################################################################################################################################################
#########################################################################################################################################################
#
# 6: Simmer trajectories
#
#########################################################################################################################################################
#########################################################################################################################################################





# make_natural_history_traj

# Creates the main simmer trajectory for the natural history component of the model
# Patients are sampled into an initial state, accrue QALYs while in no-cancer or undiagnosed states,
# and transition to cancer progression, diagnosis, screening, or death
#
# @param Env A simmer simulation environment
# @param Patient_data Data frame containing patient IDs
# @param All_parameters List containing model parameters
# @param All_costs List containing model costs
# @param All_utilities List containing model utilities
#
# @return A simmer trajectory object

make_natural_history_traj <- function(Env, All_parameters, All_costs, All_utilities, All_model_inputs, All_user_inputs) {
  
  NoCancer_Undiagnosed_traj <- trajectory("NoCancer_Undiagnosed") %>%
    
    
    # Sample initial patient characteristics
    set_attribute(keys = c("Age", "Gender", "Age_gender_score", "Smoking_status", 
                           "PLCO", "IMD", "OtherEth", "Attend_LHC", "Attend_LDCT", "Persist",
                           "Age_of_death", "Time_of_death", "State","Stage", "No_cancer_states",
                           "Respond", "Eligibility", "Join_screening", "Next_cancer_state",
                           "Next_cancer_state_time", "Last_screening_time"), values = function() {
      
      if (is.null(All_user_inputs$Patient_chars)) {
        
        
        if (All_user_inputs$Model_method == "Fixed") {
          Sim_name <- get_name(.env = Env)
          i_patient1 <- as.integer(sub("patient", "", Sim_name)) + 1
        } else {
          i_patient1 <- NULL
        }
        
        Patient_data <- patient_characteristics(All_parameters = All_parameters, 
                                                All_model_inputs = All_model_inputs, 
                                                All_user_inputs = All_user_inputs,
                                                i_patient1 = i_patient1)
        
      } else {
        
        Sim_name <- get_name(.env = Env)
        Sim_index <- as.integer(sub("patient", "", Sim_name)) + 1
        Patient_data <- All_user_inputs$Patient_chars[Sim_index,]
      }
                         
      as.numeric(Patient_data)  
    }) %>%
  
    
    # Initialise cost, utility and QALY attributes
    #set_attribute(keys = c("Total_cost", "Current_utility", "Total_QALYs"), value = function() {
      
    #  Utility_start <- All_utilities$Baseline_utility %>% pull(Utility)
    #  c(0, Utility_start, 0)
    #}) %>%
    
    # Initialise first screening time
    set_attribute(keys = c("Next_screening_time", "Screening_no"), value = function() {
      
      Join_screening <- get_attribute(.env = Env, key = "Join_screening")
      Time_of_death <- get_attribute(.env = Env, key = "Time_of_death")
      
      if (All_user_inputs$Screening_strategy == "No_screening") {
        Initial_time <- Time_of_death + 1
      } else {
        if (Join_screening == 1) {
          Initial_time <- All_user_inputs$Initial_time_to_screen
        } else if (Join_screening == 0) {
          Initial_time <- Time_of_death + 1
        } # else - eligibility error
      }
      
      as.numeric(c(Initial_time, 0))    
    }) %>%
    
    # Mark return point for natural history loop
    log_(message = "", tag = "No_or_Undiag_state", level = 2) %>%
    
    # Check if next state is progression, death or screening
    set_attribute(keys = c("Next_state", "Time_in_state"), value = function() {
      
      Next_state <- names(States_list[get_attribute(.env = Env, key = "Next_cancer_state")])
      Next_state_time <- get_attribute(.env = Env, key = "Next_cancer_state_time")
      
      #Age_current <- get_attribute(.env = Env, key = "Age")
      #Age_of_death <- get_attribute(.env = Env, key = "Age_of_death")
      
      Time_of_death <- get_attribute(.env = Env, key = "Time_of_death")
      
      Possible_states <- c(Next_state, "Other_cause_death")
      #Age_of_event <- c((Age_current + Next_time_to_state), Age_of_death)
      Time_of_events <- c(Next_state_time, Time_of_death)
      
      Next_state <- Possible_states[which.min(Time_of_events)]
      Next_state <- as.numeric(All_model_inputs$States_list[Next_state])
      Current_time <- now(.env = Env)
      Time_in_state <- min(Time_of_events) - Current_time
      
      
      
      # If screening occurs before the next natural-history transition,
      # move to the screening state instead
      if (Time_in_state + Current_time > get_attribute(.env = Env, key = 'Next_screening_time')) {
        Next_state <- 20
        Time_in_state <- get_attribute(.env = Env, key = 'Next_screening_time') - Current_time
      } else {
        Next_state <- Next_state
        Time_in_state <- Time_in_state
      }
      
      c(Next_state, Time_in_state)
    }) %>%
    
    timeout_from_attribute(key = "Time_in_state") %>%
    
    set_attribute(key = "Age", value = function() {
      Time_age <- get_attribute(.env = Env, key = 'Time_in_state')
      as.numeric(Time_age)
      }, mod = "+") %>%
    
    # QALYs during time spent in current state
    #set_attribute(keys = "Total_QALYs", values = function() {
      
     # Q_utility <- get_attribute(.env = Env, keys = "Current_utility") 
      #Q_time <- get_attribute(.env = Env, keys = "Time_in_state") 
      #Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      #QALYs <- disc_QALYs(Utility=Q_utility, Time_total=simmer::now(.env = Env), Current_time=Q_time, Discount_rate=Disc_rate)
      
      #as.numeric(QALYs)
    #}, mod = "+") %>%
    
    # BRANCH
    
    branch(option = function() get_attribute(.env = Env, keys = "Next_state"), continue = TRUE,
           
           # 1. No cancer - 0% probability
           trajectory(name = "No_cancer") %>%
             log_(message = "Problem with no cancer state!", level = 2),
           
           # 2. Preclinical NSCLC IA1
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_NSCLC_IA1", 
                                 Pathway_no = 2),
           
           # 3. Preclinical NSCLC IA2
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_NSCLC_IA2", 
                                 Pathway_no = 3),
           
           # 4. Preclinical NSCLC IA3
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_NSCLC_IA3", 
                                 Pathway_no = 4),
           
           # 5. Preclinical NSCLC IB
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_NSCLC_IB", 
                                 Pathway_no = 5),
           
           # 6. Preclinical NSCLC II
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_NSCLC_II", 
                                 Pathway_no = 6),
           
           # 7. Preclinical NSCLC III
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_NSCLC_III", 
                                 Pathway_no = 7),
           
           # 8. Preclinical NSCLC IV
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_NSCLC_IV", 
                                 Pathway_no = 8),
           
           # 9. Preclinical SCLC limited
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_SCLC_limited", 
                                 Pathway_no = 9),
           
           # 10. Preclinical SCLC extensive
           make_progression_traj(Env = Env,
                                 All_parameters = All_parameters, 
                                 All_model_inputs = All_model_inputs,
                                 All_user_inputs = All_user_inputs, 
                                 Traj_name = "Traj_Preclinical_SCLC_extensive", 
                                 Pathway_no = 10),
           
           # 11. Clinical NSCLC IA1
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_NSCLC_IA1",
             Pathway = "clinical_NSCLC_IA1",
             Route_detection = "Not_screening"),
           
           # 12. Clinical NSCLC IA2
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_NSCLC_IA2",
             Pathway = "Clinical_NSCLC_IA2",
             Route_detection = "Not_screening"),
           
           # 13. Clinical NSCLC IA3
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_NSCLC_IA3",
             Pathway = "Clinical_NSCLC_IA3",
             Route_detection = "Not_screening"),

           # 14. Clinical NSCLC IB
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_NSCLC_IB",
             Pathway = "Clinical_NSCLC_IB",
             Route_detection = "Not_screening"),
           
           # 15. Clinical NSCLC II
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_NSCLC_II",
             Pathway = "Clinical_NSCLC_II",
             Route_detection = "Not_screening"),
           
           # 16. Clinical NSCLC III
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_NSCLC_III",
             Pathway = "Clinical_NSCLC_III",
             Route_detection = "Not_screening"),
           
           # 17. Clinical NSCLC IV
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_NSCLC_IV",
             Pathway = "Clinical_NSCLC_IV",
             Route_detection = "Not_screening"),
           
           # 18. Clinical SCLC limited
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_SCLC_limited",
             Pathway = "Clinical_SCLC_limited",
             Route_detection = "Not_screening"),
           
           # 19. Clinical SCLC extensive
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_model_inputs = All_model_inputs,
             All_user_inputs = All_user_inputs,
             Traj_name = "Traj_clinical_SCLC_extensive",
             Pathway = "Clinical_SCLC_extensive",
             Route_detection = "Not_screening"),
           
           # 20. Screening
           make_screening_traj(Env = Env,
                               All_parameters = All_parameters, 
                               All_model_inputs = All_model_inputs,
                               All_user_inputs = All_user_inputs),
           
           # 21. Other cause death
           trajectory(name = "Other cause death") %>%
             set_attribute(keys = "State", value = 21)
    )  
  
  return(NoCancer_Undiagnosed_traj)
}








make_progression_traj <- function(Env, All_parameters, All_model_inputs,
                                All_user_inputs, Traj_name, Pathway_no) { #All_costs, All_utilities) {
  
  Progression_traj <- trajectory(name = Traj_name) %>%
    
    set_attribute(key = "State", value = Pathway_no) %>%
    set_attribute(key = "Stage", value = Pathway_no) %>%
    set_attribute(key = "No_cancer_states", value = 1, mod = "+") %>%
    
    set_attribute(keys = c("Next_cancer_state", "Next_cancer_state_time"), values = function() {
      
      PLCO <- get_attribute(.env = Env, key = "PLCO")
      Smoking_status_patient <- get_attribute(.env = Env, key = "Smoking_status")
      Current_state <- names(All_model_inputs$States_list[Pathway_no])
      Current_time <- now(.env = Env)
      Next_event <- NH_sample_next_state(Current_state_patient = Current_state, 
                                         States_list = All_model_inputs$States_list,
                                         Event_table = All_model_inputs$Event_table, 
                                         PLCO = PLCO,
                                         Smoking_status = Smoking_status_patient,
                                         All_parameters = All_parameters,
                                         All_user_inputs = All_user_inputs,
                                         All_model_inputs = All_model_inputs,
                                         i = as.integer(sub("patient", "", get_name(.env = Env))) + 1,
                                         j = get_attribute(.env = Env, key = "No_cancer_states"))
      Next_state <- as.numeric(All_model_inputs$States_list[Next_event$Next_state])
      Next_time <- Next_event$Time + Current_time
      
      as.numeric(c(Next_state, Next_time))
    }) %>%
    
    simmer::rollback(target = "No_or_Undiag_state", times = Inf)
  
  return(Progression_traj)
}
















# make_diagnosis_traj

# Creates a simmer trajectory for diagnosis through a specified pathway. The
# trajectory assigns the diagnosis state, applies diagnosis costs, accrues QALYs
# while waiting for diagnosis, and then branches to early or late stage
# diagnosed states with corresponding treatment costs and utilities
#
# @param Env A simmer simulation environment
# @param All_parameters List containing model parameters, including States_list
#   and Disc_rate
# @param All_utilities List containing utility values
# @param Traj_name Character scalar name of the trajectory
# @param Pathway Character scalar of diagnosis pathway/state name
# @param Diagnosis_rate Numeric scalar of rate of diagnosis
# @param Diagnosis_cost Numeric scalar of cost of diagnosis
# @param Treatment_cost Data frame containing treatment costs by diagnosed state
#
# @return A simmer trajectory object

make_diagnosis_traj <- function(Env, All_parameters, #All_costs, All_utilities, 
                                All_model_inputs, All_user_inputs, Traj_name, Pathway, Route_detection) {
  
  Diag_trajectory <- trajectory(name = Traj_name) %>%
    
    # Set diagnosis state
    set_attribute(key = "State", value = function() {
      State <- All_model_inputs$States_list[Pathway]
      
      as.numeric(State)
    }) %>%
    
    # Add discounted diagnosis costs 
    #set_attribute(key = "Total_cost", value = function() {
      #Current_time <- simmer::now(.env = Env)
      #Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      #Cost <- disc_cost(Cost = Diagnosis_cost, Time_in_days = Current_time, Discount_rate = Disc_rate)
      
      #as.numeric(Cost)
    #}, mod = "+") %>%
    
    # Apply diagnosis route specific utility
    #set_attribute(key = "Current_utility", value = function() {
      #Utility <- All_utilities$Utilities %>% filter(State == Pathway) %>% pull(Utility)
      
      #as.numeric(Utility)
    #}, mod = "*") %>%
    
    # Discounted QALYs over this period
    #set_attribute(keys = "Total_QALYs", values = function() {
      #Q_utility <- get_attribute(.env = Env, keys = "Current_utility") 
      #Q_time <- get_attribute(.env = Env, keys = "Time_in_state") 
      #Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      #QALYs <- disc_QALYs(Utility=Q_utility, Time_total=simmer::now(.env = Env), Current_time=Q_time, Discount_rate = Disc_rate)
      
      #as.numeric(QALYs)
    #}, mod = "+") %>%
    
    # Staging - early or late
    set_attribute(key = "Diagnosed_stage", value = function() {
      Stage <- get_attribute(.env = Env, key = "Stage") #- 1
      
      as.numeric(Stage)
    }) %>%
    
    set_attribute(key = "LCM", value = function() {
      
      Age <- get_attribute(.env = Env, key = 'Age')
      Stage <- names(All_model_inputs$States_list[Pathway])
      
      if (Route_detection == "Screening") {
        Route <- Route_detection
      } else {
        #Next_screening_time <- get_attribute(.Env = Env, key = 'Next_screening_time')
        Last_screening_time <- get_attribute(.env = Env, key = "Last_screening_time")
        Current_time <- now(.env = Env)
        
        if (Last_screening_time == 0) {
          Route <- "Outside"
        } else {
          if (Current_time - Last_screening_time < 1) {
            Route <- "Interval"
          } else {
            Route <- "Outside"
          }
        }
        
      }

      LCM <- lung_cancer_mortality(Age = Age, 
                                   Stage = Stage, 
                                   Route = Route, 
                                   All_parameters = All_parameters,
                                   All_user_inputs = All_user_inputs,
                                   All_model_inputs = All_model_inputs,
                                   i = as.integer(sub("patient", "", get_name(.env = Env))) + 1)
      
      as.numeric(LCM)
    }) %>%
    
    # Check if next state death from cancer or other causes
    set_attribute(keys = c("Next_state", "Time_in_state"), value = function() {
      
      Current_time <- now(.env = Env)
      Next_state <- names(States_list[22])
      Next_state_time <- get_attribute(.env = Env, key = "LCM") + Current_time
      
      Time_of_death <- get_attribute(.env = Env, key = "Time_of_death")
      
      Possible_states <- c(Next_state, "Other_cause_death")
      Time_of_events <- c(Next_state_time, Time_of_death)
      
      Next_state <- Possible_states[which.min(Time_of_events)]
      Next_state <- as.numeric(All_model_inputs$States_list[Next_state])
      
      Time_in_state <- min(Time_of_events) - Current_time
      
      c(Next_state, Time_in_state)
    }) %>%
    
    timeout_from_attribute(key = "Time_in_state") %>%
    
    set_attribute(key = "Age", value = function() {
      Time_age <- get_attribute(.env = Env, key = 'Time_in_state')
      as.numeric(Time_age)
    }, mod = "+") %>%
    
    set_attribute(key = "State", value = function() {
      Next_state <- get_attribute(.env = Env, key = "Next_state")
      
      as.numeric(Next_state)
    })
  
  return(Diag_trajectory)
}


















