# Title: Snow Crab SHAP Values
# Purpose: Calculate and plot SHAP values from BRTs
# Data created: 07/20/2022

# Libraries and functions ----
library(fastshap) # calculate SHAP values quickly
library(mshap) # combine SHAP values for two-part models
library(shapviz) # visualize SHAP values
library(doParallel)
library(here)
library(dplyr)
library(ggplot2)
source(here('code/functions', 'sv_dependence2D2.R'))
source(here('code/functions', 'sv_dependence2.R'))
source(here('code/functions', 'calculate_shap.R'))

# Read in BRTs ----
brt_mat_female_abun <- readRDS(file = here('data', 'brt_mat_female_abun.rds'))
brt_mat_female_base <- readRDS(file = here('data', 'brt_mat_female_base.rds'))
brt_imm_female_abun <- readRDS(file = here('data', 'brt_imm_female_abun.rds'))
brt_imm_female_base <- readRDS(file = here('data', 'brt_imm_female_base.rds'))
brt_leg_male_abun <- readRDS(file = here('data', 'brt_leg_male_abun.rds'))
brt_leg_male_base <- readRDS(file = here('data', 'brt_leg_male_base.rds'))
brt_sub_male_abun <- readRDS(file = here('data', 'brt_sub_male_abun.rds'))
brt_sub_male_base <- readRDS(file = here('data', 'brt_sub_male_base.rds'))

# Load and prepare data ----
crab_summary <- readRDS(here('data', 'crab_pca.rds')) %>%
  dplyr::select(-geometry)

# Transform female and male data
crab_trans <- mutate(crab_summary,
                     lncount_mat_female = log(mature_female + 1),
                     lncount_imm_female = log(immature_female + 1),
                     lncount_leg_male = log(legal_male + 1),
                     lncount_sub_male = log(sublegal_male + 1),
                     pres_imm_female = ifelse(immature_female > 0, 1, 0),
                     pres_mat_female = ifelse(mature_female > 0, 1, 0),
                     pres_leg_male = ifelse(legal_male > 0, 1, 0),
                     pres_sub_male = ifelse(sublegal_male > 0, 1, 0),
                     year_f = as.factor(year),
                     log_pcod_cpue = log(pcod_cpue + 1)) %>%
  filter(!is.na(temperature),
         !is.na(julian),
         !is.na(depth),
         !is.na(ice_mean)) 

# Customize legend labels - must be done for each plot individually
change_legend_breaks <- function(the_plot, aesthetic, breaks, labels){
  plot_sc <- as.list(the_plot$scales)$scales
  plot_color <- sapply(plot_sc, function(x) x[["aesthetics"]][1])
  plot_idx <- which(aesthetic == plot_color)
  the_plot$scales$scales[[plot_idx]][["breaks"]] <- breaks
  the_plot$scales$scales[[plot_idx]][["labels"]] <- labels
  return(the_plot)
}

leg_male_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                    "latitude", "julian", 
                    "bcs_legal_male", "log_pcod_cpue")
sub_male_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                    "latitude", "julian", 
                    "bcs_sublegal_male", "log_pcod_cpue")
mat_female_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                      "latitude", "julian", 
                      "bcs_mature_female", "log_pcod_cpue")
imm_female_names <- c("depth", "temperature", "phi", "ice_mean", "longitude",
                      "latitude", "julian", 
                      "bcs_immature_female", "log_pcod_cpue")
pred_fun <- function(X_model, newdata){
  gbm::predict.gbm(X_model, 
                   newdata,
                   n.trees = X_model$model$gbm.call$best.trees,
                   type = "response")
}

num_cores <- detectCores() - 2
vars <- c(1:8, 10, 16)

# Legal Males ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
leg_male_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_legal_male,
                longitude, latitude, julian,
                log_pcod_cpue, lncount_leg_male, legal_male, 
                pres_leg_male, year_f, year) 

leg_male_shaps <- calculate_shap(brt_leg_male_abun, 
                                 brt_leg_male_base, 
                                 leg_male_data,
                                 leg_male_names)
stopCluster(cl)

saveRDS(leg_male_shaps, file = here('data', 'leg_male_shaps.rds'))

# Visualize
# Gives the global effect of variables (absolute value, not directional)
# leg_male_shaps <- readRDS(file = here('data', 'leg_male_shaps.rds'))
leg_male_mshap <- as.matrix(leg_male_shaps[[3]]$shap_vals)
leg_male_mshap_sv <- shapviz(leg_male_mshap, X = leg_male_data)
leg_male_pshap <- as.matrix(leg_male_shaps[[2]]$shapley_values)
leg_male_pshap_sv <- shapviz(leg_male_pshap, X = leg_male_data)
leg_male_ashap <- as.matrix(leg_male_shaps[[1]]$shapley_values)
leg_male_ashap_sv <- shapviz(leg_male_ashap, X = leg_male_data)

# Variable importance
sv_importance(leg_male_mshap_sv)
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(leg_male_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Figure for presentation
sv_dependence(leg_male_mshap_sv,
              v = "depth",
              color_var = NULL,
              color = "#00868B35") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Legal Male Crab") +
  theme_classic() +
  ggplot2::theme(legend.box.spacing = grid::unit(0, "pt"),
                 axis.ticks = element_blank(),
                 plot.title = element_text(size = 22, family = "serif", face = "bold"),
                 axis.text = element_text(family = "serif", size = 18),
                 axis.title = element_text(family = "serif", size = 22),
                 axis.text.x = element_text(angle = 0, vjust = 0.7),
                 strip.text = element_text(family = "serif", size = 21),
                 legend.title = element_text(family = "serif", size = 18),
                 legend.text = element_text(family = "serif", size = 17))
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_depth_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Temperature
sv_dependence2(leg_male_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Legal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(leg_male_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Legal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Temperature (phi)
sv_dependence2(leg_male_mshap_sv, 
               v = "temperature",
               color_var = "phi") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Legal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_temperature_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(leg_male_mshap_sv, 
              v = leg_male_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Change to make limits the same for each species
options(shapviz.brewer_args = list(low = "darkslateblue",
                                   mid = "gainsboro",
                                   high = "darkred",
                                   midpoint = 0,
                                   limits = c(-2.2, 3.75))) # use to change color in sv_dependence2D2 function

# Spatial dependence
sv_dependence2D2(leg_male_mshap_sv, 
                 x = "longitude", 
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5) +  
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Legal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Added environmental effects
sv_dependence2D2(leg_male_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("phi", "temperature", "ice_mean", "depth")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Legal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_spatial_shap_env.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Test to see if effects show with bio variables
sv_dependence2D2(leg_male_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("bcs_legal_male", "log_pcod_cpue")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Legal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'leg_male_spatial_shap_bio.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()


# Force plots 
# Yellow means variable pushes prediction higher, purple means variable pushes prediction lower
# Scores close to the f(x) value have more of an impact (indicated by SHAP magnitude too)
# Would like to figure out how to separate by year without including in the model
sv_force(leg_male_mshap_sv, 2) # one observation

## Sublegal Males ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
sub_male_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_sublegal_male,
                longitude, latitude, julian, 
                log_pcod_cpue, lncount_sub_male, sublegal_male, 
                pres_sub_male, year_f, year) %>%
  tidyr::drop_na(lncount_sub_male, ice_mean) 

sub_male_shaps <- calculate_shap(brt_sub_male_abun, 
                                 brt_sub_male_base, 
                                 sub_male_data,
                                 sub_male_names)
stopCluster(cl)

saveRDS(sub_male_shaps, file = here('data', 'sub_male_shaps.rds'))

# Visualize
sub_male_shaps <- readRDS(file = here('data', 'sub_male_shaps.rds'))

sub_male_mshap <- as.matrix(sub_male_shaps[[3]]$shap_vals)
sub_male_mshap_sv <- shapviz(sub_male_mshap, X = sub_male_data)
sub_male_pshap <- as.matrix(sub_male_shaps[[2]]$shapley_values)
sub_male_pshap_sv <- shapviz(sub_male_pshap, X = sub_male_data)
sub_male_ashap <- as.matrix(sub_male_shaps[[1]]$shapley_values)
sub_male_ashap_sv <- shapviz(sub_male_ashap, X = sub_male_data)

# Variable importance
sv_importance(sub_male_mshap_sv)
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(sub_male_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Temperature
sv_dependence2(sub_male_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Sublegal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(sub_male_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Sublegal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Temperature (phi)
sv_dependence2(sub_male_mshap_sv, 
               v = "temperature",
               color_var = "phi") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Sublegal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_temperature_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(sub_male_mshap_sv, 
              v = sub_male_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Spatial dependence
# Just location effect
# Change to make limits the same for each species
options(shapviz.brewer_args = list(low = "darkslateblue",
                                   mid = "gainsboro",
                                   high = "darkred",
                                   midpoint = 0,
                                   limits = c(-2.8, 5.9))) # use to change color in sv_dependence2D2 function

# Spatial dependence
sv_dependence2D2(sub_male_mshap_sv, 
                 x = "longitude", 
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5) +  
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Sublegal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Added environmental effects
sv_dependence2D2(sub_male_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("phi", "temperature", "ice_mean", "depth")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Sublegal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_spatial_shap_env.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Test to see if effects show with bio variables
sv_dependence2D2(sub_male_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("bcs_sublegal_male", "log_pcod_cpue")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Sublegal Male Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'sub_male_spatial_shap_bio.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## Mature Females ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
mat_female_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_mature_female,
                longitude, latitude, julian,
                log_pcod_cpue, lncount_mat_female, mature_female, 
                pres_mat_female, year_f, year) %>%
  tidyr::drop_na(lncount_mat_female, ice_mean) 

mat_female_shaps <- calculate_shap(brt_mat_female_abun, 
                                   brt_mat_female_base, 
                                   mat_female_data,
                                   mat_female_names)
stopCluster(cl)

saveRDS(mat_female_shaps, file = here('data', 'mat_female_shaps.rds'))

# Visualize
mat_female_shaps <- readRDS(file = here('data', 'mat_female_shaps.rds'))

mat_female_mshap <- as.matrix(mat_female_shaps[[3]]$shap_vals)
mat_female_mshap_sv <- shapviz(mat_female_mshap, X = mat_female_data)
mat_female_pshap <- as.matrix(mat_female_shaps[[2]]$shapley_values)
mat_female_pshap_sv <- shapviz(mat_female_pshap, X = mat_female_data)
mat_female_ashap <- as.matrix(mat_female_shaps[[1]]$shapley_values)
mat_female_ashap_sv <- shapviz(mat_female_ashap, X = mat_female_data)

# Variable importance
sv_importance(mat_female_mshap_sv)
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(mat_female_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Temperature
sv_dependence2(mat_female_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Mature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(mat_female_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Mature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Temperature (phi)
sv_dependence2(mat_female_mshap_sv, 
               v = "temperature",
               color_var = "phi") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Mature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_temperature_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(mat_female_mshap_sv, 
              v = mat_female_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Spatial dependence
# Just location effects
# Change to make limits the same for each species
options(shapviz.brewer_args = list(low = "darkslateblue",
                                   mid = "gainsboro",
                                   high = "darkred",
                                   midpoint = 0,
                                   limits = c(-2, 5.6))) # use to change color in sv_dependence2D2 function

# Spatial dependence
sv_dependence2D2(mat_female_mshap_sv, 
                 x = "longitude", 
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5) +  
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Mature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Added environmental effects
sv_dependence2D2(mat_female_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("phi", "temperature", "ice_mean", "depth")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Mature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_spatial_shap_env.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Test to see if effects show with bio variables
sv_dependence2D2(mat_female_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("bcs_mature_female", "log_pcod_cpue")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Mature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'mat_female_spatial_shap_bio.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

## Immature Females ----
cl <- makeCluster(num_cores)
registerDoParallel(cl)
imm_female_data <- crab_trans %>% 
  dplyr::select(depth, temperature, phi, ice_mean, bcs_immature_female,
                longitude, latitude, julian,
                log_pcod_cpue, lncount_imm_female, immature_female, 
                pres_imm_female, year_f, year) %>%
  tidyr::drop_na(lncount_imm_female, ice_mean) 

imm_female_shaps <- calculate_shap(brt_imm_female_abun, 
                                   brt_imm_female_base, 
                                   imm_female_data,
                                   imm_female_names)
stopCluster(cl)

saveRDS(imm_female_shaps, file = here('data', 'imm_female_shaps.rds'))

# Visualize
imm_female_shaps <- readRDS(file = here('data', 'imm_female_shaps.rds'))

imm_female_mshap <- as.matrix(imm_female_shaps[[3]]$shap_vals)
imm_female_mshap_sv <- shapviz(imm_female_mshap, X = imm_female_data)
imm_female_pshap <- as.matrix(imm_female_shaps[[2]]$shapley_values)
imm_female_pshap_sv <- shapviz(imm_female_pshap, X = imm_female_data)
imm_female_ashap <- as.matrix(imm_female_shaps[[1]]$shapley_values)
imm_female_ashap_sv <- shapviz(imm_female_ashap, X = imm_female_data)

# Variable importance
sv_importance(imm_female_mshap_sv)
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_importance_shap.jpg'),
         height = 8,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Swarm importance
sv_importance(imm_female_mshap_sv, kind = "bee") # Use for explaining SHAP values, overall not as useful
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_bee_shap.jpg'),
         height = 8,
         width = 12,
         res = 200,
         units = 'in')
dev.off()

# Dependence plots for individual variables
# Temperature
sv_dependence2(imm_female_mshap_sv, 
               v = "temperature",
               color_var = "ice_mean") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Immature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_temperature_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Phi
sv_dependence2(imm_female_mshap_sv, 
               v = "phi",
               color_var = "depth") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Immature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Temperature (phi)
sv_dependence2(imm_female_mshap_sv, 
               v = "temperature",
               color_var = "phi") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             color = "black") +
  labs(title = "Immature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_temperature_phi_shap.jpg'),
         height = 6,
         width = 8,
         res = 200,
         units = 'in')
dev.off()

# Full set of dependence plots
sv_dependence(imm_female_mshap_sv, 
              v = imm_female_names,
              color_var = NULL,
              color = "aquamarine3",
              alpha = 0.3)
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_mshap.jpg'),
         height = 6,
         width = 10,
         res = 200,
         units = 'in')
dev.off()

# Spatial dependence
# Just location effects
# Change to make limits the same for each species
options(shapviz.brewer_args = list(low = "darkslateblue",
                                   mid = "gainsboro",
                                   high = "darkred",
                                   midpoint = 0,
                                   limits = c(-1.5, 5.8))) # use to change color in sv_dependence2D2 function

# Spatial dependence
sv_dependence2D2(imm_female_mshap_sv, 
                 x = "longitude", 
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5) +  
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Immature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_spatial_shap.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Added environmental effects
sv_dependence2D2(imm_female_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("phi", "temperature", "ice_mean", "depth")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Immature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_spatial_shap_env.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

# Test to see if effects show with bio variables
sv_dependence2D2(imm_female_mshap_sv, 
                 x = "longitude",
                 y = "latitude",
                 size = 2.5,
                 jitter_width = 0.5,
                 jitter_height = 0.5,
                 add_vars = c("bcs_immature_female", "log_pcod_cpue")) +
  labs(x = "Longitude \u00B0W", 
       y = "Latitude \u00B0N",
       title = "Immature Female Crab")
dev.copy(jpeg,
         here('results/SHAP',
              'imm_female_spatial_shap_bio.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
