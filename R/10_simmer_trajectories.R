#########################################################################################################################################################
#########################################################################################################################################################
#
# 10: Simmer trajectories
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

make_natural_history_traj <- function(Env, All_parameters, All_costs, All_utilities) {
  
  NoCancer_Undiagnosed_traj <- trajectory("NoCancer_Undiagnosed") %>%
    
    # Sample patient characteristics
    set_attribute(keys = c("Age", "Gender", "Age_gender_score", "Smoking_status", "PLCO"), values = function() {
      Pmale <- All_parameters %>% filter(Parameter_label == "p_male") %>% pull(Modelled_value)
      Age_mean <- All_parameters %>% filter(Parameter_label == "pop_age_mean") %>% pull(Modelled_value)
      Age_sd <- All_parameters %>% filter(Parameter_label == "pop_age_sd") %>% pull(Modelled_value)
      Age_LL <- All_parameters %>% filter(Parameter_label == "pop_age_LL") %>% pull(Modelled_value)
      Age_UL <- All_parameters %>% filter(Parameter_label == "pop_age_UL") %>% pull(Modelled_value)
      
      Patient_characteristics <- sample_patient_characteristics(n = 1, 
                                                                Pmale = Pmale, 
                                                                Age_mean = Age_mean, 
                                                                Age_sd = Age_sd, 
                                                                Age_LL = Age_LL, 
                                                                Age_UL = Age_UL,
                                                                Never_switch = 0, 
                                                                All_parameters = All_parameters)
      as.numeric(Patient_characteristics)  
    })
  
    # Age/time to death from other cause
    set_attribute(keys = c("Age_of_death", "Time_to_death"), values = function() {
      Comb_method <- "MAX"
      Age_patient <- get_attribute(.env = Env, key = "Age")
      Gender <- get_attribute(.env = Env, key = "Gender")
      PLCO <- get_attribute(.env = Env, key = "PLCO")
      Age_step <- 0.5
      General_mortality_fits <- General_mortality_fits
      OCM <- other_cause_mortality(Comb_method = Comb_method, 
                            Age_patient = Age_patient,
                            Gender = Gender, 
                            PLCO = PLCO,
                            Age_step = Age_step, 
                            General_mortality_fits = General_mortality_fits)
      
      as.numeric(OCM)
    })
    
    # Initialise cost, utility and QALY attributes
    #set_attribute(keys = c("Total_cost", "Current_utility", "Total_QALYs"), value = function() {
      
    #  Utility_start <- All_utilities$Baseline_utility %>% pull(Utility)
    #  c(0, Utility_start, 0)
    #}) %>%
    
    # Sample initial state 
    set_attribute(keys = c("State","Stage"), value = function() {
      
      PLCO <- get_attribute(.env = Env, key = "PLCO")
      Samp_state <- sample_initial_state(States_list = States_list, 
                                         PLCO = PLCO, 
                                         All_parameters = All_parameters)
      Stage <- Samp_state
      
      c(Samp_state, Stage)
    }) %>%
    
    # Initialise first screening time
    #set_attribute(key = 'Next_screening_time', value = All_parameters$Screening_times) %>%
      # WILL DEPEND ON SCREENING STRATEGY
    
    # Mark return point for natural history loop
    log_(message = "", tag = "No_or_Undiag_state", level = 2) %>%
    
    # Sample next state and time to state
    set_attribute(keys = c("Time_in_state", "Next_state"), value = function() {
      
      PLCO <- get_attribute(.env = Env, key = "PLCO")
      Smoking_status_patient <- get_attribute(.env = Env, key = "PLCO")
      Current_state <- get_attribute(.env = Env, key = "Smoking_status")
      Next_event <- NH_sample_next_state(Current_state_patient = Current_state, 
                                         States_list = States_list,
                                         Event_table = Event_table, 
                                         PLCO = PLCO,
                                         Smoking_status = Smoking_status_patient,
                                         All_parameters = All_parameters)
      Age_current <- get_attribute(.env = Env, key = "Age")
      Age_of_death <- get_attribute(.env = Env, key = "Age_of_death")
      
      Possible_states <- c(Next_event$Next_state, "Other_cause_death")
      Age_of_event <- c((Age_current + Next_event$Time), Age_of_death)
      
      Next_state <- Possible_states[which.min(Age_of_event)]
      Time_in_state <- min(Age_of_event) - Age_current
      
      # REMOVE SCREENING FOR NOW
      # If screening occurs before the next natural-history transition,
      # move to the screening state instead
      #if (Next_event$Time + Current_time > get_attribute(.env = Env, key = 'Next_screening_time')) {
      #  Next_state <- 6
      #  Time_in_state <- get_attribute(.env = Env, key = 'Next_screening_time') - Current_time
      #} else {
      #  Next_state <- as.numeric(All_parameters$States_list[Next_event$Next_state])
      #  Time_in_state <- Next_event$Time
      #}
      
      c(Time_in_state, Next_state)
    }) %>%
    
    timeout_from_attribute(key = "Time_in_state") %>%
    
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
           
           # 2. Undiagnosed early stage
           trajectory(name = "Undiagnosed_early") %>%
             
             set_attribute(keys = "State", value = 2) %>%
             set_attribute(keys = "Stage", value = 2) %>%
             simmer::rollback(target = "No_or_Undiag_state", times = Inf),
           
           # 3. Undiagnosed late stage
           trajectory(name = "Undiagnosed_late") %>%
             
             set_attribute(keys = "State", value = 3) %>%
             set_attribute(keys = "Stage", value = 3) %>%
             simmer::rollback(target = "No_or_Undiag_state", times = Inf),
           
           # 4. Diagnosed clinical presentation
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_utilities = All_utilities,
             Traj_name = "traj_diagnosed_clinical",
             Pathway = "Diag_clin",
             Diagnosis_rate = All_parameters$Diagnosis_rates %>% filter(State == "Diag_clin") %>% pull(Rate),
             Diagnosis_cost = All_costs$Diagnosis_costs %>% filter(State == "Diag_clin") %>% pull(Cost),
             Treatment_cost = All_costs$Treatment_costs),
           
           # 5. Diagnosed incidental
           make_diagnosis_traj(
             Env = Env,
             All_parameters = All_parameters,
             All_utilities = All_utilities,
             Traj_name = "traj_diagnosed_incidental",
             Pathway = "Diag_inc",
             Diagnosis_rate = All_parameters$Diagnosis_rates %>% filter(State == "Diag_inc") %>% pull(Rate),
             Diagnosis_cost = All_costs$Diagnosis_costs %>% filter(State == "Diag_inc") %>% pull(Cost),
             Treatment_cost = All_costs$Treatment_costs),
           
           # 6. Screening
           make_screening_traj(Env = Env,
                               All_parameters = All_parameters, 
                               All_costs = All_costs,
                               All_utilities = All_utilities),
           
           # 7. Death
           trajectory(name = "Death") %>%
             set_attribute(keys = "State", value = 7)
    )  
  
  return(NoCancer_Undiagnosed_traj)
}














# make_stage_traj

# Creates a simmer trajectory for a diagnosed cancer stage. The trajectory
# sets the current state, adds discounted treatment cost, applies the relevant
# utility value, waits for a fixed period, and accrues discounted QALYs
#
# @param Env A simmer simulation environment
# @param All_parameters List containing model parameters, including States_list
#   and Disc_rate
# @param All_utilities List containing utility values
# @param Traj_name Character scalar name of the trajectory
# @param Pathway Character scalar of diagnosis pathway/state name
# @param Treatment_cost Treatment cost for early/late stage cancer
#
# @return A `simmer` trajectory object.

make_stage_traj <- function(Env, All_parameters, All_utilities, Traj_name, Pathway, Treatment_cost) {
  
  Stage_trajectory <- trajectory(name = Traj_name) %>%
    
    # Set current diagnosed state
    set_attribute(key = "State", value = function() {
      Stage_state <- All_parameters$States_list[Pathway]
      
      as.numeric(Stage_state)
    }) %>%
    
    # Add discounted treatment cost at current model time
    set_attribute(key = "Total_cost", value = function() {
      Current_time <- simmer::now(.env = Env)
      Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      Cost <- disc_cost(Cost = Treatment_cost, Time_in_days = Current_time, Discount_rate = Disc_rate)
      
      as.numeric(Cost)
    }, mod = "+") %>%
    
    # Apply stage-specific utility
    set_attribute(key = "Current_utility", value = function() {
      Utility <- All_utilities$Utilities %>% filter(State == Pathway) %>% pull(Utility)
      
      as.numeric(Utility)
    }, mod = "*") %>%
    
    # Set time in treatment
    timeout(task = 20) %>%
    
    # Discounted QALYs over this period
    set_attribute(keys = "Total_QALYs", values = function() {
      Q_utility <- get_attribute(.env = Env, keys = "Current_utility") 
      Q_time <- 20
      Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      QALYs <- disc_QALYs(Utility=Q_utility, Time_total=simmer::now(.env = Env), Current_time=Q_time, Discount_rate = Disc_rate)
      
      as.numeric(QALYs)
    }, mod = "+")
  
  return(Stage_trajectory)
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

make_diagnosis_traj <- function(Env, All_parameters, All_utilities, Traj_name, Pathway, Diagnosis_rate, Diagnosis_cost, Treatment_cost) {
  
  Diag_trajectory <- trajectory(name = Traj_name) %>%
    
    # Set diagnosis state
    set_attribute(key = "State", value = function() {
      State <- All_parameters$States_list[Pathway]
      
      as.numeric(State)
    }) %>%
    
    # Add discounted diagnosis costs 
    set_attribute(key = "Total_cost", value = function() {
      Current_time <- simmer::now(.env = Env)
      Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      Cost <- disc_cost(Cost = Diagnosis_cost, Time_in_days = Current_time, Discount_rate = Disc_rate)
      
      as.numeric(Cost)
    }, mod = "+") %>%
    
    # Apply diagnosis route specific utility
    set_attribute(key = "Current_utility", value = function() {
      Utility <- All_utilities$Utilities %>% filter(State == Pathway) %>% pull(Utility)
      
      as.numeric(Utility)
    }, mod = "*") %>%
    
    # Time spent in diagnosis testing
    set_attribute(key = "Time_in_state", value = function() {
      Time_to_next_state <- rexp(1, rate = Diagnosis_rate)
      
      as.numeric(Time_to_next_state)
    }) %>%
    
    timeout_from_attribute(key = "Time_in_state") %>%
    
    # Discounted QALYs over this period
    set_attribute(keys = "Total_QALYs", values = function() {
      Q_utility <- get_attribute(.env = Env, keys = "Current_utility") 
      Q_time <- get_attribute(.env = Env, keys = "Time_in_state") 
      Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      QALYs <- disc_QALYs(Utility=Q_utility, Time_total=simmer::now(.env = Env), Current_time=Q_time, Discount_rate = Disc_rate)
      
      as.numeric(QALYs)
    }, mod = "+") %>%
    
    # Staging - early or late
    set_attribute(key = "Diagnosed_stage", value = function() {
      Stage <- get_attribute(.env = Env, key = "Stage") - 1
      
      as.numeric(Stage)
    }) %>%
    
    # BRANCH
    
    branch(option = function() get_attribute(.env = Env, "Diagnosed_stage"), continue = TRUE,
           
           # 1. Early stage
           make_stage_traj(Env = Env, 
                           All_parameters = All_parameters, 
                           All_utilities = All_utilities, 
                           Traj_name = "traj_early_stage", 
                           Pathway = "Diag_early", 
                           Treatment_cost = Treatment_cost %>% filter(State == "Diag_early") %>% pull(Cost)),
           
           # 2. Late stage
           make_stage_traj(Env = Env, 
                           All_parameters = All_parameters, 
                           All_utilities = All_utilities, 
                           Traj_name = "traj_late_stage", 
                           Pathway = "Diag_late", 
                           Treatment_cost = Treatment_cost %>% filter(State == "Diag_late") %>% pull(Cost))
    )
  
  return(Diag_trajectory)
}





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

make_screening_traj <- function(Env, All_parameters, All_costs, All_utilities) {
  
  traj_screening <- trajectory("screening") %>%
    
    # Set current state to screening
    set_attribute(key = "State", value = 6) %>% 
    
    # Schedule next screening time
    set_attribute(key = 'Next_screening_time', value = All_parameters$Screening_times, mod = "+") %>%
    
    # Add discounted screening test costs
    set_attribute(key = "Total_cost", value = function() {
      
      Current_time <- simmer::now(.env = Env)
      Cost <- All_costs$Diagnosis_costs %>% filter(State == "Screening", Step == 1) %>% pull(Cost)
      Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      Cost <- disc_cost(Cost = Cost, Time_in_days = Current_time, Discount_rate = Disc_rate)
      
      as.numeric(Cost)
    }, mod = "+") %>%
    
    # Apply temporary screening disutility
    set_attribute(key = "Current_utility", value = function() {
      Utility <- -(All_utilities$Utilities %>% filter(State == "Screening") %>% pull(Utility))
      
      as.numeric(Utility)
    }, mod = "+") %>%
    
    # Decide whether screening detects existing cancer
    set_attribute("screen_detected", function() {
      
      Stage <- names(All_parameters$Staging[get_attribute(.env = Env, key = "Stage")])
      S_result <- screening_result(Stage = Stage, 
                                   Early_sensitivity = All_parameters$Screening_parameters %>% filter(Parameter == "Early_sensitivity") %>% pull(Value), 
                                   Late_sensitivity = All_parameters$Screening_parameters %>% filter(Parameter == "Late_sensitivity") %>% pull(Value), 
                                   Specificity = All_parameters$Screening_parameters %>% filter(Parameter == "Specificity") %>% pull(Value))
      
      as.numeric(All_parameters$Screening[S_result])
    }) %>%
    
    # Wait for screening result
    timeout(task = 7) %>% 
    
    # QALYs during screening-result waiting period
    set_attribute(keys = "Total_QALYs", values = function() {
      
      Q_utility <- get_attribute(.env = Env, keys = "Current_utility") 
      Q_time <- 7
      Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
      QALYs <- disc_QALYs(Utility = Q_utility, Time_total = simmer::now(.env = Env), Current_time = Q_time, Discount_rate = Disc_rate)
      
      as.numeric(QALYs)
    }, mod = "+") %>%
    
    # BRANCH
    
    branch(option = function() get_attribute(.env = Env, "screen_detected"), continue = TRUE,
           
           # 1. Screen negative: return to natural history (includes false negative)
           trajectory("screen_negative") %>%
             
             # Restore state to underlying cancer state
             set_attribute(key = "State", value = function() {
               
               Stage <- get_attribute(.env = Env, key = "Stage")
               State <- Stage # 1 - no cancer, 2 - early, 3 - late
               
               as.numeric(State)
             }) %>%
             # Remove temporary screening disutility
             set_attribute(key = "Current_utility", value = function() {
               
               Utility <- (All_utilities$Utilities %>% filter(State == "Screening") %>% pull(Utility))
               
               as.numeric(Utility)
             }, mod = "+") %>%
             simmer::rollback(target = "No_or_Undiag_state", times = Inf),
           
           # 2. Screen positive: proceed to further testing/diagnosis
           trajectory("screen_positive") %>%
             
             # Add diagnostic assessment cost following positive screening
             set_attribute(key = "Total_cost", value = function() {
               
               Current_time <- simmer::now(.env = Env)
               Cost <- All_costs$Diagnosis_costs %>% filter(State == "Screening", Step == 2) %>% pull(Cost)
               Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
               Cost <- disc_cost(Cost = Cost, Time_in_days = Current_time, Discount_rate = Disc_rate)
               
               as.numeric(Cost)
             }, mod = "+") %>%
             
             # Sample time from positive screen to diagnosis outcome
             set_attribute(key = "Time_in_state", value = function() {
               
               Time_to_next_state <- rexp(1, rate = All_parameters$Screening_parameters %>% filter(Parameter == "Rate1") %>% pull(Value))
               
               as.numeric(Time_to_next_state)
             }) %>%
             
             # Determine whether this is a true positive or a false positive
             set_attribute(key = "Diagnosed", value = function() {
               
               Stage <- names(All_parameters$Staging[get_attribute(.env = Env, key = "Stage")])
               if (Stage == 'None') {
                 Diagnosed <- 1
               } else {
                 Diagnosed <- 2
               }
               
               as.numeric(Diagnosed)
             }) %>%
             
             timeout_from_attribute(key = "Time_in_state") %>%
             
             # QALYs between positive screening and further testing
             set_attribute(keys = "Total_QALYs", values = function() {
               
               Q_utility <- get_attribute(.env = Env, keys = "Current_utility") 
               Q_time <- get_attribute(.env = Env, keys = "Time_in_state") 
               Disc_rate <- All_parameters$Disc_rate %>% pull(Rate)
               QALYs <- disc_QALYs(Utility=Q_utility, Time_total=simmer::now(.env = Env), Current_time=Q_time, Discount_rate = Disc_rate)
               
               as.numeric(QALYs)
             }, mod = "+") %>%
             
             # Remove temporary screening disutility
             set_attribute(key = "Current_utility", value = function() {
               Utility <- (All_utilities$Utilities %>% filter(State == "Screening") %>% pull(Utility))
               
               as.numeric(Utility)
             }, mod = "+") %>%
             
             # BRANCH
             
             branch(option = function() get_attribute(.env = Env, "Diagnosed"), continue = TRUE,
                    
                    # 1. False positive: return to no cancer state
                    trajectory("false_positive") %>%
                      set_attribute(key = "State", value = 1) %>%
                      simmer::rollback(target = "No_or_Undiag_state", times = Inf),
                    
                    # 2. True positive: follow screen-detected diagnosis pathway
                    make_diagnosis_traj(
                      Env = Env,
                      All_parameters = All_parameters,
                      All_utilities = All_utilities,
                      Traj_name = "traj_diagnosed_screening",
                      Pathway = "Diag_screen",
                      Diagnosis_rate = All_parameters$Screening_parameters %>% filter(Parameter == "Rate2") %>% pull(Value),
                      Diagnosis_cost = All_costs$Diagnosis_costs %>% filter(State == "Screening", Step == 2) %>% pull(Cost),
                      Treatment_cost = All_costs$Treatment_costs)
                    
             )
    )
  
  return(traj_screening)
}















