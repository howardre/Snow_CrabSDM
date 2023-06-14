plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  plot(gam,
       pages = 0,
       select = covariate, 
       shade = T,
       shade.col = "#eb8055b2",
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       seWithMean = T,
       scale = 0,
       cex.axis = 6.5,
       cex.lab = 6.5,
       family = "serif",
       lwd = 2.5)
}