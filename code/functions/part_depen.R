part_depen <- function(abun_brt){
  windows()
  gbm.plot(abun_brt$model,
           plot.layout = c(3, 4),
           write.title = F,
           smooth = T,
           common.scale = T,
           cex.axis = 1.7,
           cex.lab = 1.7,
           lwd = 1.5,
           show.contrib = FALSE,
           family = "serif")
}