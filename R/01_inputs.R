#########################################################################################################################################################
#########################################################################################################################################################
#
# 1: Patient initial characteristics + non-intervention effected inputs
#
#########################################################################################################################################################
#########################################################################################################################################################


# Try to keep this file separate
# So that its not being run with the model and patient inputs are selected previously


source("03_population.R")
  

Patient_df <- as.data.frame(t(sapply(1:5, function(i) {
  patient_characteristics(
      All_parameters = All_parameters,
      All_model_inputs = All_model_inputs,
      All_user_inputs = All_user_inputs,
      i_patient1 = i)})))


Patient_df
colnames(Patient_df) <- c("Age", "Gender", "Age_gender_score", "Smoking_status", 
  "PLCO", "IMD", "OtherEth", "Attend_LHC", "Attend_LDCT", "Persist",
  "Age_of_death", "Time_of_death", "State","Stage", "No_cancer_states",
  "Eligibility", "Next_cancer_state", "Next_cancer_state_time")


write.csv(Patient_df, "Patient_characteristics_df.csv", row.names = FALSE)







  