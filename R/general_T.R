#' General function to generate a prediction expression for a family of Taguchi 
#'   (T) methods
#' 
#' \code{general_T} is a (higher-order) general function that generates a 
#'   prediction expression for a family of Taguchi (T) methods. Each T method 
#'   can be implemented by setting the parameters of this function appropriately.
#'
#' @param unit_space_data Matrix with n rows (samples) and (p + 1) columns 
#'                          (variables). The 1 ~ p th columns are independent 
#'                          variables and the (p + 1) th column is a dependent 
#'                          variable. Underlying data to obtain a representative 
#'                          point for the normalization of the 
#'                          \code{signal_space_data}. All data should be 
#'                          continuous values and should not have missing values. 
#' @param signal_space_data Matrix with m rows (samples) and (p + 1) columns 
#'                            (variables). The 1 ~ p th columns are independent 
#'                            variables and the (p + 1) th column is a dependent 
#'                            variable. Underlying data to generate a prediction 
#'                            expression. All data should be continuous values 
#'                            and should not have missing values. 
#' @param generates_transform_functions A function that takes the 
#'                                        \code{unit_space_data} as an (only) 
#'                                        argument and returns a list containing 
#'                                        three functions. A data transformation 
#'                                        function for independent variables is 
#'                                        the first component, a data 
#'                                        transformation function for a 
#'                                        dependent variable is the second 
#'                                        component, and an inverse function of 
#'                                        the data transformation function for a 
#'                                        dependent variable is the third 
#'                                        component. The data transformation 
#'                                        function for independent variables 
#'                                        takes independent variable data (a 
#'                                        matrix of p columns) as an (only) 
#'                                        argument and returns the transformed 
#'                                        independent variable data. The data 
#'                                        transformation function for a 
#'                                        dependent variable takes dependent 
#'                                        variable data (a vector) as an (only) 
#'                                        argument and returns the transformed 
#'                                        dependent variable data. The inverse 
#'                                        function of the data transformation 
#'                                        for a dependent variable takes the 
#'                                        transformed dependent variable data (a 
#'                                        vector) as an (only) argument and 
#'                                        returns the untransformed dependent 
#'                                        variable data.
#' @param includes_transformed_data If \code{TRUE}, then the transformed data 
#'                                    are included in a return object.
#' 
#' @return A list containing the following components is returned.
#' 
#'  \item{beta_hat}{Vector with length q. Estimated proportionality constants 
#'                   between each independent variable and the dependent 
#'                   variable.}
#'  \item{eta_hat}{Vector with length q. Estimated squared signal-to-noise 
#'                  ratios (S/N) coresponding to \code{beta_hat}.}
#'  \item{M_hat}{Vector with length n. The estimated values of the dependent 
#'                variable after the data transformation for 
#'                \code{signal_space_data}.} 
#'  \item{overall_prediction_eta}{Numeric. The overall squared signal-to-noise 
#'                                 ratio (S/N).}
#'  \item{transforms_independent_data}{Data transformation function generated 
#'                                      from \code{generates_transform_functions} 
#'                                      based on \code{unit_space_data}. The 
#'                                      function for independent variables takes 
#'                                      independent variable data (a matrix of p 
#'                                      columns) as an (only) argument and 
#'                                      returns the transformed independent 
#'                                      variable data.}
#'  \item{transforms_dependent_data}{Data transformation function generated in 
#'                                    \code{generates_transform_functions} based 
#'                                    on the \code{unit_space_data}. The 
#'                                    function for a dependent variable takes 
#'                                    dependent variable data (a vector) as an 
#'                                    (only) argument and returns the 
#'                                    transformed dependent variable data.}
#'  \item{inverses_transformed_dependent_data}{Inverse function generated in the 
#'                                              \code{generates_transform_functions} 
#'                                              based on \code{unit_space_data}. 
#'                                              The function of the takes the 
#'                                              transformed dependent variable 
#'                                              data (a vector) as an (only) 
#'                                              argument and returns the 
#'                                              dependent variable data inversed 
#'                                              from the transformed dependent 
#'                                              variable data.}
#'  \item{m}{The number of samples for \code{signal_space_data}.}
#'  \item{q}{The number of independent variables after the data transformation. 
#'            According to the data transoformation function, q may be equal to 
#'            p.}
#'  \item{X}{If \code{includes_transformed_data} is \code{TRUE}, then the 
#'            independent variable data after the data transformation for the 
#'            \code{signal_space_data} are included.}
#'  \item{M}{If \code{includes_transformed_data} is \code{TRUE}, then the (true) 
#'            value of the dependent variable after the data transformation for 
#'            the \code{signal_space_data} are included.}
#'    
#' @seealso \code{\link{T1}}, \code{\link{Ta}}, and \code{\link{Tb}}
#' 
#' @examples 
#' # The value of the dependent variable of the following samples mediates  
#' # in the stackloss dataset.
#' stackloss_center <- stackloss[c(9, 10, 11, 20, 21), ] 
#' 
#' # The following samples are data other than the unit space data and the test 
#' # data.   
#' stackloss_signal <- stackloss[-c(2, 9, 10, 11, 12, 19, 20, 21), ] 
#' 
#' # The following settings are same as the T1 method.             
#' model <- general_T(unit_space_data = stackloss_center, 
#'                    signal_space_data = stackloss_signal,
#'                    generates_transform_functions = 
#'                                        generates_transformation_functions_T1, 
#'                    includes_transformed_data = TRUE)
#'                          
#' (model$M_hat)
#' 
#' @export
general_T <- function(unit_space_data, 
                      signal_space_data, 
                      generates_transform_functions,
                      includes_transformed_data = FALSE) {
  
  check_data(unit_space_data)
  check_data(signal_space_data)
  
  if (is.vector(unit_space_data)) {
    unit_space_data <- matrix(unit_space_data, ncol = length(unit_space_data))
  }
  
  if (is.vector(signal_space_data)) {
    signal_space_data <- 
                     matrix(signal_space_data, ncol = length(signal_space_data))
  }
  
  functions <- generates_transform_functions(unit_space_data)
  transforms_independent_data <- functions[[1]]
  transforms_dependent_data <- functions[[2]]
  inverses_transformed_dependent_data <- functions[[3]]
  
  p <- ncol(signal_space_data) - 1
  
  signal_space_data_x <- signal_space_data[, -(p + 1)]
  signal_space_data_y <- signal_space_data[, (p + 1)]
  
  X <- transforms_independent_data(signal_space_data_x)
  M <- transforms_dependent_data(signal_space_data_y)
  
  r <- sum(M ^ 2)
  S_T <- apply(X ^ 2, 2, sum)
  S_beta <- apply(X * M, 2, sum) ^ 2 / r
  S_e <- S_T - S_beta
  V_e <- S_e / (nrow(X) - 1)
  
  beta_hat <- apply(X * M, 2, sum) / r
  
  eta_hat <- r ^ -1 * (S_beta - V_e) / V_e
  eta_hat[eta_hat < 0] <- 0
  
  M_hat <- calc_M_hat(X, beta_hat, eta_hat)
  
  #y_hat <- M_hat + M_0
  
  overall_prediction_eta <- calc_overall_predicton_eta(M, M_hat)
  
  model <- list(beta_hat = beta_hat, eta_hat = eta_hat, 
                M_hat = M_hat, #M_0 = M_0, #y_hat = y_hat, 
                overall_prediction_eta = overall_prediction_eta,
                transforms_independent_data = transforms_independent_data,
                transforms_dependent_data = transforms_dependent_data,
                inverses_transformed_dependent_data = 
                                            inverses_transformed_dependent_data)
  
  model$m <- nrow(X)
  model$q <- ncol(X)
  
  if (includes_transformed_data) {
    model$X <- X
    model$M <- M
  }  
  
  return(model)
  
}

#' General function to implement a forecasting method for a family of Taguchi (T) 
#'   methods
#' 
#' \code{general_forecasting.T} is the general function that implements a 
#'   forecasting method for a family of Taguchi (T) methods. Each forecasting 
#'   method of a family of T methods can be implemented by setting the 
#'   parameters of this function appropriately.
#' 
#' @param model Object generated as a model. 
#' @param newdata Matrix with n rows (samples) and p columns (variables). The 
#'                  data are used to calculate the desired distances from the 
#'                  unit space. All data should be continuous values and should 
#'                  not have missing values.  
#' @param includes_transformed_newdata If \code{TRUE}, then the transformed data 
#'                                       for \code{newdata} are included in a 
#'                                       return object.
#' 
#' @return A list containing the following components is returned.
#' 
#'  \item{M_hat}{Vector with length n. The estimated values of the dependent 
#'                variable after the data trasformation.}
#'  \item{y_hat}{Vector with length n. The estimated values after the inverse 
#'                transformation from \code{M_hat}.}
#'  \item{model}{Object passed by \code{model}.}
#'  \item{n}{The number of samples for \code{newdata}.}
#'  \item{q}{The number of variables after the data transformation.}
#'  \item{X}{If \code{includes_transformed_newdata} is \code{TRUE}, then the 
#'            transformed data for \code{newdata} are included.}
#' 
#' @seealso \code{\link{forecasting.T1}}, \code{\link{forecasting.Ta}}, and 
#'            \code{\link{forecasting.Tb}}
#' 
#' @examples
#' # The value of the dependent variable of the following samples mediates  
#' # in the stackloss dataset.
#' stackloss_center <- stackloss[c(9, 10, 11, 20, 21), ] 
#' 
#' # The following samples are data other than the unit space data and the test 
#' # data.   
#' stackloss_signal <- stackloss[-c(2, 9, 10, 11, 12, 19, 20, 21), ] 
#' 
#' # The following settings are same as the T1 method. 
#' model <- general_T(unit_space_data = stackloss_center, 
#'                    signal_space_data = stackloss_signal,
#'                    generates_transform_functions = 
#'                                        generates_transformation_functions_T1, 
#'                    includes_transformed_data = TRUE)
#'                    
#' # The following test samples are chosen casually. 
#' stackloss_test <- stackloss[c(2, 12, 19), -4] 
#' 
#' forecasting <- general_forecasting.T(model = model, 
#'                                      newdata = stackloss_test, 
#'                                      includes_transformed_newdata = TRUE)
#'                               
#' (forecasting$y_hat) # Estimated values
#' (stackloss[c(2, 12, 19), 4]) # True values
#' 
#' @export
general_forecasting.T <- function(model, 
                                  newdata, 
                                  includes_transformed_newdata = FALSE) {
  
  if (!missing(newdata)) {
    
    check_data(newdata)

    # For the case newdata is given as vector when the number of sample = 1.
    if (is.vector(newdata)) {
      newdata <- matrix(newdata, ncol = length(newdata))
    }
    
    transformed_newdata <- model$transforms_independent_data(newdata)
    
    M_hat <- calc_M_hat(transformed_newdata, model$beta_hat, model$eta_hat)    

  } else {
    
    M_hat <- model$M_hat     
      
  }
  
  y_hat <- model$inverses_transformed_dependent_data(M_hat)

  forecasting <- list(M_hat = M_hat, y_hat = y_hat)
  
  forecasting$model <- model
  
  if (!missing(newdata)) {
    forecasting$n <- nrow(transformed_newdata)
    forecasting$p <- ncol(transformed_newdata)
  
    if (includes_transformed_newdata) {
      forecasting$X <- transformed_newdata
    }
    
  } else {
    forecasting$n <- model$n
    forecasting$p <- model$p
    
    if (includes_transformed_newdata) {
      if (!is.null(model$X)){
        forecasting$X <- model$X
      } else {
        warning(paste("Data values will not be included in the return object,",
                      "because data values are not included in the model object."))
      }
    }
  }
  
  return(forecasting)
  
}

#' Function to estimate M value (M hat) for a family of T methods.
#' 
#' \code{calc_M_hat} estimates M values (M hat) for the T method.
#' 
#' @param X Matrix with n rows (samples) and q columns (variables). The 
#'            independent variable data after the data transformation. All data 
#'            should be continuous values and should not have missing values.
#' @param beta_hat Vector with length q. Estimated proportionality constants 
#'                   between each independent variable and the dependent 
#'                   variable.
#' @param eta_hat Vector with length q. Estimated squared signal-to-noise ratios 
#'                  (S/N) coresponding to \code{beta_hat}.
#' 
#' @return Vector with length n. Estimated M values (M hat).
#' 
#' @seealso \code{\link{general_T}} and \code{\link{general_forecasting.T}}
#' 
#' @examples
#' # The value of the dependent variable of the following samples mediates  
#' # in the stackloss dataset.
#' stackloss_center <- stackloss[c(9, 10, 11, 20, 21), ] 
#' 
#' # The following samples are data other than the unit space data and the test 
#' # data.   
#' stackloss_signal <- stackloss[-c(2, 9, 10, 11, 12, 19, 20, 21), ] 
#' 
#' # The following settings are same as the T1 method.             
#' model <- general_T(unit_space_data = stackloss_center, 
#'                    signal_space_data = stackloss_signal,
#'                    generates_transform_functions = 
#'                                        generates_transformation_functions_T1, 
#'                    includes_transformed_data = TRUE)
#'         
#' modified_eta_hat <- model$eta_hat
#' modified_eta_hat[3] <- 0
#'         
#' (modified_M_hat <- calc_M_hat(model$X, model$beta_hat, modified_eta_hat))
#' 
#' @export
calc_M_hat <- function(X, beta_hat, eta_hat) {
  
  if (ncol(X) != length(beta_hat)) {
    stop("X and beta_hat should have the same length.")
  }
  
  if (length(beta_hat) != length(eta_hat)) {
    stop("beta_hat and eta_hat should have the same length.")
  }
  
  M_hat <- apply(eta_hat * t(X) / beta_hat, 2, sum) / sum(eta_hat)
  
  return(M_hat)
  
}

#' Function to calculate overall prediction eta for the T method
#' 
#' \code{calc_M_hat} calculates the overall prediction eta for the T method.
#' 
#' @param M Vector with length n. The (true) value of the dependent 
#'            variable after the data trasformation.
#' @param M_hat Vector with length n. The estimated values of the dependent 
#'                variable after the data trasformation.
#' 
#' @return Numeric. Overall prediction eta which is used to measure the 
#'           estimation accuracy.
#' 
#' @seealso \code{\link{general_T}} and \code{\link{general_forecasting.T}}
#' 
#' @examples
#' 
#' # The value of the dependent variable of the following samples mediates  
#' # in the stackloss dataset.
#' stackloss_center <- stackloss[c(9, 10, 11, 20, 21), ] 
#' 
#' # The following samples are data other than the unit space data and the test 
#' # data.   
#' stackloss_signal <- stackloss[-c(2, 9, 10, 11, 12, 19, 20, 21), ] 
#' 
#' # The following settings are same as the T1 method.             
#' model <- general_T(unit_space_data = stackloss_center, 
#'                    signal_space_data = stackloss_signal,
#'                    generates_transform_functions = 
#'                                        generates_transformation_functions_T1, 
#'                    includes_transformed_data = TRUE)
#'                
#' modified_eta_hat <- model$eta_hat
#' modified_eta_hat[3] <- 0
#'         
#' modified_M_hat <- calc_M_hat(model$X, model$beta_hat, modified_eta_hat)
#' 
#' (modified_overall_predicton_eta <- calc_overall_predicton_eta(model$M, 
#'                                                              modified_M_hat))
#' 
#' @export
calc_overall_predicton_eta <- function(M, M_hat) {
  
  if (length(M) != length(M_hat)) {
    stop("M and M_hat should be the same length.")
  }
  
  L <- sum(M * M_hat) 
  r <- sum(M ^ 2)
  S_T <- sum(M_hat ^ 2)
  S_beta <- L ^ 2 / r 
  S_e <- S_T - S_beta
  V_e <- S_e / (length(M) - 1)
  
  eta_hat <- r ^ -1 * (S_beta - V_e) / V_e
  #log10ed_eta_hat <- 10 * log(eta_hat, base = 10)
  
  return(eta_hat = eta_hat)
  
}
