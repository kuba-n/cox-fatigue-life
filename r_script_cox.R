#SURVIVAL FUNCTION
surv_function <- Surv(r_dataframe$fatigue_life, r_dataframe$fatigue_survival)

#COX FUNCTION
#create a temporary formula, so a different number of variables can be parsed in coxph() function
temp_formula <- as.formula(paste("surv_function~", paste(dataframe_keys[3:length(dataframe_keys)], collapse = "+") ))

cox_function <- coxph(temp_formula, data = r_dataframe)

# test for proportional hazards assumption
cox_zph <- cox.zph(cox_function)