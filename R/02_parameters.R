#########################################################################################################################################################
#########################################################################################################################################################
#
# 2: Parameters
#
#########################################################################################################################################################
#########################################################################################################################################################



All_parameters <- read.csv("Parameters.csv", header=TRUE)





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
    "nh_prog_lambda_7", "nh_prog_xi_7", NA, NA,
    "nh_prog_mu_1", "nh_prog_phi_1", NA, NA,
    "nh_prog_mu_2", "nh_prog_phi_2", NA, NA,

    rep(NA, 37)
  )
)
