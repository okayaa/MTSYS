#' Function to generate a prediction expression for the Tb method
#' 
#' \code{Tb} generates a prediction expression for the Tb method. In 
#'   \code{\link{general_T}}, the data are normalized by subtracting the center
#'   and without scaling based on \code{sample_data}. The center is determined 
#'   by the specific way for the Tb method. For details, please see 
#'   \code{\link{generates_transformation_functions_Tb}}. All the sample data 
#'   are used for both unit space and signal space. 
#'
#' @param sample_data Matrix with n rows (samples) and (p + 1) columns 
#'                      (variables). The 1 ~ p th columns are independent 
#'                      variables and the (p + 1) th column is a dependent 
#'                      variable. All data should be continuous values and 
#'                      should not have missing values. 
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
#'                variable after the data transformation for \code{sample_data}.}
#'  \item{overall_prediction_eta}{Numeric. The overall squared signal-to-noise 
#'                                 ratio (S/N).}
#'  \item{transforms_independent_data}{Data transformation function generated 
#'                                      from \code{generates_transform_functions} 
#'                                      based on the \code{unit_space_data}. The 
#'                                      function for independent variables takes 
#'                                      independent variable data (a matrix of p 
#'                                      columns) as an (only) argument and 
#'                                      returns the transformed independent 
#'                                      variable data.}
#'  \item{transforms_dependent_data}{Data transformation function generated from 
#'                                    \code{generates_transform_functions} based 
#'                                    on the \code{unit_space_data}. The 
#'                                    function for a dependent variable takes 
#'                                    dependent variable data (a vector) as an 
#'                                    (only) argument and returns the 
#'                                    transformed dependent variable data.}
#'  \item{inverses_dependent_data}{Data transformation function generated 
#'                                  from \code{generates_transform_functions} 
#'                                  based on the \code{unit_space_data}. The 
#'                                  function of the takes the transformed 
#'                                  dependent variable data (a vector) as an 
#'                                  (only) argument and returns the dependent 
#'                                  variable data inversed from the transformed 
#'                                  dependent variable data.}
#'  \item{m}{The number of samples for \code{sample_data}.}
#'  \item{q}{The number of independent variables after the data transformation. 
#'            q equals p.}
#'  \item{X}{If \code{includes_transformed_data} is \code{TRUE}, then the 
#'            independent variable data after the data transformation for the 
#'            \code{sample_data} are included.}
#'  \item{M}{If \code{includes_transformed_data} is \code{TRUE}, then the (true) 
#'            value of the dependent variable after the data transformation for 
#'            the \code{sample_data} are included.}
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
#' @seealso \code{\link{general_T}}, 
#'            \code{\link{generates_transformation_functions_Tb}}, and 
#'            \code{\link{forecasting.Tb}}    
#' 
#' @examples  
#' model_Tb <- Tb(sample_data = stackloss[-c(2, 12, 19), ], 
#'                includes_transformed_data = TRUE)
#'                          
#' (model_Tb$M_hat)
#' 
#' @export
Tb <- function(sample_data, includes_transformed_data = FALSE) {
  
  model_Tb <- general_T(unit_space_data = sample_data, 
                        signal_space_data = sample_data,
                        generates_transform_functions = 
                                          generates_transformation_functions_Tb, 
                        includes_transformed_data = includes_transformed_data)
  
  class(model_Tb) <- "Tb"
  
  return(model_Tb)
  
}

#' Forecasting method for the Tb method
#' 
#' \code{forecasting.Tb} (via \code{\link{forecasting}}) estimates the dependent 
#'   values based on the Tb model.
#' 
#' @param model Object of class "Tb" generated by \code{\link{Tb}} or 
#'                \code{\link{generates_model}}(..., method = "Tb").
#' @param newdata Matrix with n rows (samples) and p columns (variables). The 
#'                  Data to be estimated. All data should be continuous values 
#'                  and should not have missing values.  
#' @param includes_transformed_newdata If \code{TRUE}, then the transformed data 
#'                                       for \code{newdata} are included in a 
#'                                       return object.
#' 
#' @return A list containing the following components is returned.
#' 
#'  \item{M_hat}{Vector with length n. The estimated values of the dependent 
#'                variable after the data transformation.}
#'  \item{y_hat}{Vector with length n. The estimated values after the inverse 
#'                transformation from \code{M_hat}.}
#'  \item{model}{Object of class "Tb" passed by \code{model}.}
#'  \item{n}{The number of samples for \code{newdata}.}
#'  \item{q}{The number of variables after the data transformation. q equals p.} 
#'  \item{X}{If \code{includes_transformed_newdata} is \code{TRUE}, then the 
#'            transformed data for \code{newdata} are included.}
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
#' @seealso \code{\link{general_forecasting.T}} and \code{\link{Tb}} 
#'     
#' @examples     
#' model_Tb <- Tb(sample_data = stackloss[-c(2, 12, 19), ], 
#'                includes_transformed_data = TRUE)
#'                               
#' forecasting_Tb <- forecasting(model = model_Tb, 
#'                               newdata = stackloss[c(2, 12, 19), -4], 
#'                               includes_transformed_newdata = TRUE)
#'                               
#' (forecasting_Tb$y_hat) # Estimated values
#' (stackloss[c(2, 12, 19), 4]) # True values                        
#' 
#' @export
forecasting.Tb <- function(model, 
                           newdata, 
                           includes_transformed_newdata = FALSE) {
  
  if (!inherits(model, "Tb")) {
    warning("calling forecasting.Tb(<fake-Tb-object>) ...")
  }
  
  general_forecasting.T(model = model, 
                        newdata = newdata, 
                        includes_transformed_newdata = 
                                                   includes_transformed_newdata)
  
}
