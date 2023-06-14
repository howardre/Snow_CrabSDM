brt_grid_preds <- function(spatial_grid, abun_brt, base_brt){
  spatial_grid$pred_abun <- predict.gbm(abun_brt$model,
                                        spatial_grid,
                                        n.trees = abun_brt$model$gbm.call$best.trees,
                                        type =  "response")
  spatial_grid$pred_base <- predict.gbm(base_brt$model,
                                        spatial_grid,
                                        n.trees = base_brt$model$gbm.call$best.trees,
                                        type = "response")
  spatial_grid$pred_brt <- spatial_grid$pred_base * spatial_grid$pred_abun
  
  spatial_grid$pred_brt[spatial_grid$dist > 28000] <- NA # Remove predictions too far from samples
  return(spatial_grid)
}

brt_deviance <- function(brt){
  deviance <- (brt$model$self.statistics$mean.null - brt$model$cv.statistics$deviance.mean) /
    brt$model$self.statistics$mean.null
  return(deviance)
}