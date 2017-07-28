#' Function to generate data transformation functions for the T1 methods
#' 
#' \code{generates_transformation_functions_T1} is the argument for the 
#'   parameter \code{generates_transform_functions} in \code{genera_T}, which 
#'   is used in the T1 method. In addtion, the Ta method also uses this function 
#'   for the argument. 
#' 
#' @param unit_space_data Matrix with n rows (samples) and (p + 1) columns 
#'                          (variables). Data to generate the unit space. All 
#'                          data should be continuous values and should not have 
#'                          missing values. 
#' 
#' @return \code{generates_transformation_functions_T1} returns a list 
#'           containing three functions. For the first component, the data 
#'           transformation function for independent variables is a function 
#'           that subtracts the mean of each independent variable. For the 
#'           second component, the data transformation function for a dependent 
#'           variable is a function that subtracts the mean of a dependent 
#'           variable. For the third component, the inverse function of the data 
#'           transformation function for a dependent variable is a function that 
#'           adds the mean of a dependent variable. The mean used is the mean of 
#'           the \code{unit_space_data}.  
#' 
#' @seealso \code{\link{T1}} and \code{\link{Ta}}
#' 
#' @examples    
#' # The value of the dependent variable of the following samples mediates  
#' # in the stackloss dataset.
#' stackloss_center <- stackloss[c(9, 10, 11, 20, 21), ] 
#'       
#' tmp <- generates_transformation_functions_T1(stackloss_center)
#' mean_subtraction_function <- tmp[[1]]
#' subtracts_M_0 <- tmp[[2]]
#' adds_M_0 <- tmp[[3]] 
#' 
#' is.function(mean_subtraction_function) # TRUE
#' is.function(subtracts_M_0) # TRUE
#' is.function(adds_M_0) # TRUE
#' 
#' @export
generates_transformation_functions_T1 <- function(unit_space_data) {
  
  center <- apply(unit_space_data, 2, mean)
  unit_space_center <- center[-length(center)]  
  M_0 <- center[length(center)]
  
  subtracts_mean <- 
         generates_normalization_function(unit_space_center = unit_space_center,
                                          is_scaled = FALSE)
  
  subtracts_M_0 <- function(x) x - M_0
  
  adds_M_0 <- function(x) x + M_0
  
  return(list(subtracts_mean, subtracts_M_0, adds_M_0))
  
}


#' Function to generate data transformation functions for the Tb methods
#' 
#' \code{generates_transformation_functions_Tb} is the argument for the 
#'   parameter \code{generates_transform_functions} in \code{genera_T}, which 
#'   is used in the Tb method.
#' 
#' @param sample_data Matrix with n rows (samples) and (p + 1) columns 
#'                      (variables). The Tb method uses all data to generate the 
#'                      unit space. All data should be continuous values and 
#'                      should not have missing values. 
#' 
#' @return \code{generates_transformation_functions_Tb} returns a list 
#'           containing three functions. For the first component, the data 
#'           transformation function for independent variables is a function 
#'           that subtracts the center of each independent variable. For the 
#'           second component, the data transformation function for a dependent 
#'           variable is a function that subtracts the weighted mean of a 
#'           dependent variable. For the third component, the inverse function 
#'           of the data transformation function for a dependent variable is a 
#'           function that adds the weighted mean of a dependent variable. The 
#'           center is determined in a specific manner for the Tb method. The 
#'           center consists of each sample value which maximizes the 
#'           signal-to-noise ratio (S/N) per variable. The values are determined 
#'           independently so that different samples may be selected for 
#'           different variables. The weighted mean is calculated as the 
#'           frequency of being selected in independent variables, similar to 
#'           the center of the dependent variable.
#' 
#' @references
#'   Inou, A., Nagata, Y., Horita, K., & Mori, A. (2012). Prediciton Accuracies 
#'     of Improved Taguchi's T Methods Compared to those of Multiple Regresssion 
#'     Analysis. \emph{Journal of the Japanese Society for Quality Control, 
#'     42}(2), 103-115. (In Japanese) 
#' 
#'   Kawada, H., & Nagata, Y. (2015). An application of a generalized inverse 
#'     regression estimator to Taguchi's T-Method. \emph{Total Quality Science, 
#'     1}(1), 12-21.
#' 
#' @seealso \code{\link{Tb}}
#' 
#' @examples
#' # The value of the dependent variable of the following samples mediates  
#' # in the stackloss dataset.
#' stackloss_center <- stackloss[c(9, 10, 11, 20, 21), ] 
#'     
#' tmp <- generates_transformation_functions_Tb(stackloss_center)
#' center_subtraction_function <- tmp[[1]]
#' subtracts_M_0 <- tmp[[2]]
#' adds_M_0 <- tmp[[3]] 
#' 
#' is.function(center_subtraction_function) # TRUE
#' is.function(subtracts_M_0) # TRUE
#' is.function(adds_M_0) # TRUE
#' 
#' @export
generates_transformation_functions_Tb <- function(sample_data) {
  
  get_eta <- function(one_sample_data, all_sample_data) {
    model <- general_T(unit_space_data = one_sample_data, 
                       signal_space_data = all_sample_data,
                       generates_transform_functions = 
                                          generates_transformation_functions_T1, 
                       includes_transformed_data = FALSE)
      
    return(model$eta_hat)
      
  }    
    
  etas <- apply(sample_data, 1, get_eta, sample_data)
  # apply per row (=1) because the rows and columns were transposed in above. 
  max_eta_index <- apply(etas, 1, which.max)

  unit_space_center <- 
                 diag(as.matrix(sample_data[max_eta_index, -ncol(sample_data)]))
  
  #M_0 <- sum(etas / sum(etas) * sample_data[max_eta_index, ncol(sample_data)])
  max_eta <- diag(as.matrix(etas[, max_eta_index]))
  M_0 <- 
     sum(max_eta / sum(max_eta) * sample_data[max_eta_index, ncol(sample_data)])
  
  subtracts_center <- 
         generates_normalization_function(unit_space_center = unit_space_center,
                                          is_scaled = FALSE)
  
  subtracts_M_0 <- function(x) x - M_0
  
  adds_M_0 <- function(x) x + M_0
  
  return(list(subtracts_center, subtracts_M_0, adds_M_0))
  
}
