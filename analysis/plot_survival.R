library(survminer)
library(survival)

#' plot_survival
#' plot survival 
#' @param data
#' @param variables
#' @param cutoff cutoff methods. When the variables is numerical, we divide the variables into two groups if setting the cutoff methods. cutoff methods is either "median", "mean", "best" or the certain values.
#' @param minprop Setting minprop values if cutoff methods is "best"
#' @param time columns names specifying time.
#' @param status columns names specifying status.
#' @param palette colors for different groups.
#' @param return.data whether return the data frame, Default is FALSE.
#' @param ... parameters from ggsurvplot to control the styles of figures


plot_survival = function(data, 
                         variables,
                         cutoff,
                         time = "time", 
                         status = "status",
                         palette = NULL,
                         minprop = 0.4,
                         return.data = FALSE,
                         ...
                         ){
  if(is.null(palette)){
    palette = c("#E41A1C", "#E7B800", "blue", "green")
  }
  
  
  #add data
  data$time <- data[[time]]
  data$status <- data[[status]]
  data <- dplyr::filter(data, !is.na(.data$time), !is.na(.data$status))
  
  
  #cutoffs: setting the cutoffs for the continuous variables
  
  #see data types.
  types = lapply( variables, function(x) is.numeric(data[[x]]) ) %>%
    purrr::reduce(c)
  
  if(!is.null(cutoff)){
    
    if(sum(types) > 0 ){
      
      message("Set cutoffs for numerical variables by the methods: ", cutoff )
      
      if(cutoff == "mean"){
        cutoff_values = apply(data[, variables[types]] , 2, mean, na.rm = T)
      }else if(cutoff == "median"){
        cutoff_values = apply(data[, variables[types]] , 2, median, na.rm = T)
      }else if(cutoff == "best"){
        cutoff_values = survminer::surv_cutpoint(
          data,
          time = time,
          event = status,
          variables = variables[types],
          minprop = minprop,
          progressbar = TRUE
        )
        cutoff_values = setNames(cutoff_values$cutpoint[, 1], nm = variables[types] )
      }else if(is.numeric(cutoff_values)){
        cutoff_values = cutoff
      }
      
      for(i in variables[types]){
        data[[paste0(i, ".g")]] = ifelse(data[[i]] > cutoff_values[i], "High","Low")
      }
    }
    
  }
  
  variables[types] = str_c(variables[types], ".g")
  
  
  fm = lapply(variables, function(x) as.formula(sprintf("Surv(time, status) ~ %s", x) )
              )

  #fm = as.formula(sprintf("Surv(time, status) ~ %s",  str_c(variables, collapse = " + ") ))
  #model <- coxph(fm , data = data) 
  #broom::tidy(model, exp = T) %>% print()
  
  fit = lapply(fm, function(x) surv_fit(x, data = data ) )
  
  plt = lapply(fit, function(x) 
    ggsurvplot(x ,
             data = data,
             censor = TRUE,
             risk.table = T, 
             risk.table.height = 0.3, 
             linetype = 1, 
             pval = TRUE,
             palette = palette,
             #legend.labs = c("High","Low"),
             ylab = "Overall Survival (OS)",
             #conf.int = TRUE, # Add confidence interval
             ncensor.plot = FALSE,
             ...
  )
  )
  
  names(plt) = variables
  
  if(return.data){
    list(
      data = data,
      plot.list = plt
    )
  }else{
    plt
  }
  
}



