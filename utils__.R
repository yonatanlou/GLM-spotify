rescale_0_to_10 <- function(x) {
  scaled <- scale(x, center = min(x), scale = max(x) - min(x))
  rescaled <- as.numeric(scaled * 10)
  
  return(rescaled)
}

calculate_dispersion <- function(model){
  sum(residuals(model, type="pearson")^2)/df.residual(model)
  
}

calculate_comparison_df <- function(models, model_name) {
  
  deviance <- sapply(models, deviance)
  dispersion_param <- sapply(models, calculate_dispersion)
  aic <- sapply(models, AIC)
  bic <- sapply(models, BIC)
  df_res <- sapply(models, function(model)
    sum(model$df.residual))
  model_names <- sprintf("%s%d", model_name, 1:length(models))
  comparison_df <-
    data.frame(model_names,deviance, dispersion_param, aic, bic, df_res)
  comparison_df$scaled_deviance <-
    comparison_df$deviance / comparison_df$dispersion_param
  
  return(comparison_df)
}


extractMeaningfullCoeff <- function(model, q) {
  coefficients <- coef(summary(poisson_model))[, "Estimate"]
  p_values <- coef(summary(poisson_model))[, "Pr(>|z|)"]
  standard_errors <- coef(summary(poisson_model))[, "Std. Error"]
  res <- data.frame(p_values, standard_errors, coefficients)
  threshold <- quantile(abs(coefficients), q)
  significant_coefficients <- coefficients[p_values < 0.05 & abs(coefficients) > threshold]
  
  # Print the significant coefficients
  return(data.frame(coef=sort(abs(significant_coefficients)))  )
}

extractFromCoefTest <- function(model) {
  res <- coeftest(model, vcov = sandwich::vcovHC)[,] %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "coef")
  return(res)
}

library(ggplot2)

generateCoefPlot <- function(models) {
  i <- 1
  plots_list <- list()
  
  for (model in models) {
    coef_summary <- extractFromCoefTest(model)
    
    # Remove the intercept (if present)
    coef_summary <- coef_summary[coef_summary$coef != "(Intercept)", ]
    
    coef_exp <- exp(coef_summary[, c("Estimate")])
    conf_intervals <- exp(confint(model))
    conf_intervals <- conf_intervals[rownames(conf_intervals) != "(Intercept)", ]
    
    result_df <- data.frame(
      Feature = coef_summary$coef,
      coefExp = coef_exp,
      LowerCI = na.omit(conf_intervals[, 1]),
      UpperCI = na.omit(conf_intervals[, 2])
    )
        resStats <- t(data.frame(aic=round(model$aic,1), deviance=round(model$aic,1)))
    g <- ggplot(result_df, aes(x = coefExp, y = Feature)) +
      geom_point(size = 1) +
      geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0) +
      labs(
        title = paste("M", i,"dev:",round(model$deviance,0) ,",aic:",round(model$aic,0) ),
        x = "Relative risk exp(beta)",
        y = "Feature"
      )
      # annotation_custom(tableGrob(resStats), xmin=max(coef_exp) * 0.989, xmax=max(coef_exp)*0.997, ymin=max(result_df[result_df$Feature!=max(result_df$Feature),"Feature"]), ymax=max(result_df$Feature))
    plots_list <- c(plots_list, list(g))
    i <- 1+1
  }
  
  return(plots_list)
}

createCoefPlot <- function(model) {
  coef_summary <- extractFromCoefTest(model)
  
  # Remove the intercept (if present)
  coef_summary <- coef_summary[coef_summary$coef != "(Intercept)", ]
  
  coef_exp <- exp(coef_summary[, c("Estimate")])
  conf_intervals <- exp(confint(model,method="Wald"))
  conf_intervals <- conf_intervals[rownames(conf_intervals) != "(Intercept)", ]
  
  result_df <- data.frame(
    Feature = coef_summary$coef,
    coefExp = coef_exp,
    LowerCI = na.omit(conf_intervals[, 1]),
    UpperCI = na.omit(conf_intervals[, 2])
  )
  
  resStats <- t(data.frame(aic = round(model$aic, 1), deviance = round(model$aic, 1)))
  
  g <- ggplot(result_df, aes(x = coefExp, y = Feature)) +
    geom_point(size = 1) +
    geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0) +
    labs(
      title = paste("Dev:", round(model$deviance, 0), ", AIC:", round(model$aic, 0)),
      x = "Relative risk exp(beta)",
      y = "Feature"
    )
  
  return(g)
}

CreatePlotMultipleModels <- function(models) {
  plots_list <- list()
  
  for (model in models) {
    coef_plot <- createCoefPlot(model)
    plots_list <- c(plots_list, list(coef_plot))
  }
  
  return(plots_list)
}
