variable_figure <- function(gam, range){
  par(mfrow = c(3, 3),
      mar = c(17, 21, .5, 0.6) + 0.1,
      oma = c(3, 1, 1, 1),
      mgp =  c(9, 3, 0))
  plot_variable(gam,
                covariate = 2,
                bounds = range,
                "Day of Year",
                "Species Abundance Anomalies",
                "s")
  plot_variable(gam,
                covariate = 3,
                bounds = range,
                "Depth",
                " ",
                "n")
  plot_variable(gam,
                covariate = 4,
                bounds = range,
                "Phi",
                " ",
                "n")
  plot_variable(gam,
                covariate = 5,
                bounds = range,
                "Temperature",
                "Species Abundance Anomalies",
                "s")
  plot_variable(gam,
                covariate = 6,
                bounds = range,
                "Ice Concentration",
                " ",
                "n")
  plot_variable(gam,
                covariate = 7,
                bounds = range,
                "PCA Loading",
                " ",
                "n")
  plot_variable(gam,
                covariate = 8,
                bounds = range,
                "Cod Abundance",
                "Species Abundance Anomalies",
                "s")
  plot_variable(gam,
                covariate = 9,
                bounds = range,
                "BCS Prevalence",
                " ",
                "n")
}
