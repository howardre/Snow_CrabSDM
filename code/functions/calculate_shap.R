calculate_shap <- function(brt_abun, brt_pres, data, vars){
  explain <- data[vars] # contains just variables to calculate SHAP for
  
  # Baseline/expected values
  abun_baseline = mean(pred_fun(brt_abun$model,
                                newdata = explain))
  pres_baseline = mean(pred_fun(brt_pres$model, 
                                newdata = explain))
  
  # Calculate individual model SHAP values
  set.seed(1993)
  abun_shap <- fastshap::explain(brt_abun$model, 
                                 X = explain,
                                 pred_wrapper = pred_fun,
                                 nsim = 50,
                                 shap_only = FALSE,
                                 adjust = TRUE,
                                 parallel = TRUE)
  pres_shap <- fastshap::explain(brt_pres$model,
                                 X = explain,
                                 pred_wrapper = pred_fun,
                                 nsim = 50,
                                 shap_only = FALSE,
                                 adjust = TRUE,
                                 parallel = TRUE)
  
  # Combine SHAP values
  mshap <- mshap(shap_1 = abun_shap$shapley_values,
                 shap_2 = pres_shap$shapley_values,
                 ex_1 = abun_baseline,
                 ex_2 = pres_baseline)
  
  # Return list of all three types of SHAP values
  return(list(abun_shap, pres_shap, mshap))
}
