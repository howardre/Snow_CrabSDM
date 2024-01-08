## Legal Males ----
### TreeSHAP ----
library(treeshap)
library(waterfalls)
leg_male_data <- crab_trans %>% # can use full data set
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station) %>%
  tidyr::drop_na(lncount_leg_male, ice_mean)
leg_male_explain <- leg_male_data[vars] %>% tidyr::drop_na() # cannot use year factor unless hot coded
leg_male_gbm <- gbm.unify(brt_leg_male_abun$model, leg_male_explain)
leg_male_gbm_abun <- treeshap(leg_male_gbm, leg_male_explain)
plot_contribution(leg_male_gbm_abun)
leg_male_sv <- shapviz(leg_male_gbm_abun)

sv_importance(leg_male_sv)
sv_importance(leg_male_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
sv_waterfall(leg_male_sv, 1) # one observation
sv_force(leg_male_sv, 3) # one observation
sv_dependence(leg_male_sv, 
              v = "temperature", 
              color_var = "ice_mean") # specific variable relationships
sv_dependence(leg_male_sv, 
              v = leg_male_names, 
              color_var = NULL,
              color = "#6666CC30")

# Try waterfall split out by year
leg_male_shap <- leg_male_gbm_temp$shaps
leg_male_df <- cbind(leg_male_shap, year = leg_male_data[, 14])
leg_male_final <- leg_male_df[, c(11, 1:10)]
leg_male_waterfall <- leg_male_final %>% 
  group_by(year) %>%
  summarise(across(everything(), list(sum)))

waterfall(leg_male_waterfall)
ggplot(leg_male_df, aes(x = temperature, y = temp_shap)) +
  geom_point(aes(color = ice_mean))

### mSHAP ----
library(mshap)
library(shapper)
library(DALEX)
library(r2pmml)

# May need to specify using anaconda if on Windows machine before attempting to use shap package
# reticulate::use_condaenv("C:/Users/howar/anaconda3")
leg_male_data <- crab_train %>% # can use full data set
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian, legal_male_loading,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year, legal_male_loading_station) %>%
  tidyr::drop_na(lncount_leg_male, ice_mean)
leg_male_explain <- leg_male_data[vars]
leg_male_pres <- leg_male_data[, 13]
leg_male_abun <- leg_male_data[, 11]

# Convert to model usable in treeshap
leg_male_gbm_abun <- gbm.unify(brt_leg_male_abun$model, leg_male_explain)
leg_male_gbm_pres <- gbm.unify(brt_leg_male_base$model, leg_male_explain)

# Calculate SHAP values for each model
leg_male_tshap_abun <- treeshap(leg_male_gbm_abun, leg_male_explain)
leg_male_tshap_pres <- treeshap(leg_male_gbm_pres, leg_male_explain,)

leg_male_abun_sv <- shapviz(leg_male_tshap_abun)
leg_male_pres_sv <- shapviz(leg_male_tshap_pres)
leg_male_msv <- mshapviz(c(leg_male_abun_sv, leg_male_pres_sv))
sv_waterfall(leg_male_abun_sv, 50)
sv_force(leg_male_pres_sv, 50)
sv_dependence(leg_male_abun_sv, 
              v = "temperature", 
              color_var = "ice_mean")

test_baseline <- get_baseline(leg_male_tshap_abun)
test <- get_shap_values(leg_male_abun_sv)

# fastshap
library(fastshap)
pred_fun <- function(X_model, newdata){
  predict.gbm(X_model, newdata)
}

leg_male_abun_baseline = mean(pred_fun(brt_leg_male_abun$model, 
                                       newdata = leg_male_explain))

leg_male_abun_shap <- fastshap::explain(brt_leg_male_abun$model, 
                                        X = leg_male_explain, 
                                        pred_wrapper = pred_fun, 
                                        nsim = 50,
                                        shap_only = FALSE,
                                        adjust = TRUE)
shapviz_ob <- shapviz(shap)
sv_importance(shapviz_ob, kind = "bee")
sv_importance(leg_male_abun_sv, kind = "bee")
sv_dependence(shapviz_ob, 
              v = "temperature", 
              color_var = "ice_mean")
sv_force(shapviz_ob, 30)
sv_force(leg_male_abun_sv, 30)

leg_male_abun_shap <- leg_male_tshap_abun$shaps
leg_male_pres_shap <- leg_male_tshap_pres$shaps

leg_male_mshap <- mshap(shap_1 = leg_male_abun_shap,
                        shap_2 = leg_male_pres_shap,
                        ex_1 = leg_male_abun_baseline, 
                        ex_2 = 0)

mshap::summary_plot(variable_values = leg_male_explain,
                    shap_values = leg_male_mshap$shap_vals)

leg_male_shaps <- as.matrix(leg_male_mshap$shap_vals)
leg_male_mshap_sv <- shapviz(leg_male_shaps, X = leg_male_explain)
sv_importance(leg_male_mshap_sv)
sv_importance(leg_male_mshap_sv, kind = "bee")

sv_dependence(leg_male_mshap_sv, 
              v = "temperature", 
              color_var = "ice_mean")
sv_dependence(leg_male_mshap_sv, v = leg_male_names)

### KernelSHAP ----
leg_male_explain <- leg_male_train[vars] # only use columns in model
leg_male_x <- leg_male_explain[sample(nrow(leg_male_explain), 500), ]
leg_male_shap_abun <- kernelshap(brt_leg_male_abun$model, 
                                 leg_male_explain, 
                                 bg_X = leg_male_x)
saveRDS(leg_male_shap_abun, file = here('data', 'leg_male_shap_full.rds'))

# Visualize
leg_male_sv <- shapviz(leg_male_shap_abun)
# Gives the global effect of variables (absolute value, not directional)
sv_importance(leg_male_sv)
sv_importance(leg_male_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
sv_waterfall(leg_male_sv, 1) # one observation
sv_waterfall(leg_male_sv, leg_male_sv$X$year == "1996") # observations in one year
# Force plots 
# Yellow means variable pushes prediction higher, purple means variable pushes prediction lower
# Scores close to the f(x) value have more of an impact (indicated by SHAP magnitude too)
sv_force(leg_male_sv, 2) # one observation
sv_force(leg_male_sv, leg_male_sv$X$year == "1995") # observations in one year
sv_dependence(leg_male_sv, 
              v = "temperature", 
              color_var = "ice_mean") # specific variable relationships
sv_dependence(leg_male_sv, v = leg_male_names)

# Could I make a heatmap for years by variable? Show change in SHAP over time
# Could just do heatmap by stage/sex