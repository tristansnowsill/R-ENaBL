#########################################################################################################################################################
#########################################################################################################################################################
#
# 7: Model run
#
#########################################################################################################################################################
#########################################################################################################################################################


########################
# 1. LOAD R PACKAGES
########################

library(simmer)
library(simmer.plot)
library(dplyr)
library(survival)
library(ggplot2)
library(MASS)
library(msm)
library(nnet)
library(SHELF)
library(tidyr)


# Load in all functions, parameters, inputs

Patient_df <- read.csv("Patient_characteristics_df.csv", header=TRUE)

source("06_simmer_trajectories.R")
source("05_screening.R")
source("02_parameters.R")
source("03_population.R")
source("04_natural_history.R")


Model <- run_model(Screening_strategy = "Biennial_screen",
                   Patient_number = 10,
                   Respond_or_join = TRUE,
                   Model_method = "Fixed",
                   Model_route = "Fixed",
                   Screen_design = "Design_1")

Model <- run_model(Screening_strategy = "Annual_screen",
                   Patient_number = 10,
                   Respond_or_join = TRUE,
                   Model_method = "Random",
                   Model_route = "Risk_in",
                   Screen_design = "None")

  
Model$patient0 %>% print(n=200)
Model$patient1 %>% print(n=200)
Model$patient2 %>% print(n=200)
Model$patient3 %>% print(n=200)
Model$patient4 %>% print(n=200)
Model$patient5 %>% print(n=200)
Model$patient6 %>% print(n=200)
Model$patient7 %>% print(n=200)
Model$patient8 %>% print(n=200)
Model$patient9 %>% print(n=200)
















All_user_inputs <- list(Mortality_comb_method = "SUM",
                        Patient_chars = NULL,
                        Screening_strategy = "Annual_screen",
                        Initial_time_to_screen = 0.08,
                        Patient_number = 10,
                        Respond_or_join = TRUE,
                        #Model_method = "Fixed",
                        Model_method = "Random",
                        #Model_route = "Fixed",
                        Model_route = "Risk_in",
                        Screen_design = "None")



source("10_simmer_trajectories.R")
source("05_screening.R")
source("02_parameters.R")
source("03_population.R")
source("04_natural_history.R")


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
#Atts


#Atts %>% filter(name == "patient0")
#Atts %>% filter(name == "patient1")
#Atts %>% filter(name == "patient2")
#Atts %>% filter(name == "patient3")
#Atts %>% filter(name == "patient4")
#Atts %>% filter(name == "patient5")
#Atts %>% filter(name == "patient6")
#Atts %>% filter(name == "patient7")
#Atts %>% filter(name == "patient8")
#Atts %>% filter(name == "patient9")

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

Atts_results %>% print(n=200)

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


Atts_labelled %>% filter(name == "patient0") %>% print(n=200)
Atts_labelled %>% filter(name == "patient1") %>% print(n=200)
Atts_labelled %>% filter(name == "patient2") %>% print(n=200)
Atts_labelled %>% filter(name == "patient3") %>% print(n=200)
Atts_labelled %>% filter(name == "patient4") %>% print(n=200)
Atts_labelled %>% filter(name == "patient5") %>% print(n=200)
Atts_labelled %>% filter(name == "patient6") %>% print(n=200)
Atts_labelled %>% filter(name == "patient7") %>% print(n=200)
Atts_labelled %>% filter(name == "patient8") %>% print(n=200)
Atts_labelled %>% filter(name == "patient9") %>% print(n=200)

