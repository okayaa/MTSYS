#' General function to generate a unit space for a family of 
#'   Mahalanobis-Taguchi (MT) methods
#' 
#' \code{general_MT} is a (higher-order) general function that generates a unit 
#'   space for a family of Mahalanobis-Taguchi (MT) methods. Each MT method can 
#'   be implemented by setting the parameters of this function appropriately.
#'
#' @param unit_space_data Matrix with n rows (samples) and p columns (variables).
#'                          Data to generate the unit space. All data should be 
#'                          continuous values and should not have missing values. 
#' @param generates_transform_function Function that takes \code{unit_space_data} 
#'                                       as an (only) argument and returns a 
#'                                       data transformation function. The data 
#'                                       transformation function takes data as 
#'                                       an (only) argument and returns the 
#'                                       transformed data.          
#' @param calc_A Function that returns A in a quadratic form x'Ax. \code{calc_A} 
#'                 takes the transformed data as an (only) argument.   
#' @param includes_transformed_data If \code{TRUE}, then the transformed data 
#'                                    are included in a return object. 
#' 
#' @return A list containing the following components is returned.
#' 
#' \item{A}{q x q matrix calculated by \code{calc_A}.}
#' \item{calc_A}{Function passed by \code{calc_A}.}
#' \item{transforms_data}{Data transformation function generated from 
#'                         \code{generates_transform_function} based on 
#'                         \code{unit_space_data}.}
#' \item{distance}{Vector with length n. Distances from the unit space to each 
#'                  sample.}
#' \item{n}{The number of samples.}
#' \item{q}{The number of independent variables after the data transformation. 
#'           According to the data transoformation function, q may be equal to p.} 
#' \item{x}{If \code{includes_transformed_data} is \code{TRUE}, then the 
#'           transformed data are included.}
#'   
#' @seealso \code{\link{MT}}, \code{\link{MTA}} and \code{\link{RT}}
#' 
#' @examples 
#' # 40 data for versicolor in the iris dataset                            
#' iris_versicolor <- iris[61:100, -5] 
#' 
#' # The following settings are same as the MT method.                          
#' unit_space <- general_MT(unit_space_data = iris_versicolor, 
#'                          generates_transform_function = 
#'                                             generates_normalization_function,
#'                          calc_A = function(x) solve(cor(x)),  
#'                          includes_transformed_data = TRUE)
#'                          
#' (unit_space$distance)
#' 
#' @export
general_MT <- function(unit_space_data, 
                       calc_A,
                       generates_transform_function,
                       includes_transformed_data = FALSE) {
  
  check_data(unit_space_data)
  
  transforms_data <- generates_transform_function(unit_space_data)
  
  transformed_unit_space_data <- transforms_data(unit_space_data)
  
  A <- calc_A(transformed_unit_space_data)

  x <- as.matrix(transformed_unit_space_data)
  distance <- sqrt(rowSums((x %*% A) * x) / ncol(x))
  
  unit_space <- list(A = A, calc_A = calc_A, transforms_data = transforms_data,
                     distance = distance)
  
  unit_space$n <- nrow(transformed_unit_space_data)
  unit_space$q <- ncol(transformed_unit_space_data)

  if (includes_transformed_data) {
    unit_space$x <- transformed_unit_space_data
  }
  
  return(unit_space)

}

#' General function to implement a diagnosis method for a family of 
#'   Mahalanobis-Taguchi (MT) methods
#' 
#' \code{general_diagnosis.MT} is the general function that implements a 
#'   diagnosis method for a family of Mahalanobis-Taguchi (MT) methods. Each 
#'   diagnosis method of a family of MT methods can be implemented by setting 
#'   the parameters of this function appropriately.
#' 
#' @param unit_space Object generated as a unit space. 
#' @param newdata Matrix with n rows (samples) and p columns (variables). The 
#'                  data are used to calculate the desired distances from the 
#'                  unit space. All data should be continuous values and should 
#'                  not have missing values.  
#' @param threshold Numeric specifying the threshold value to classify each 
#'                    sample into positive (\code{TRUE}) or negative 
#'                    (\code{FALSE}).
#' @param includes_transformed_newdata If \code{TRUE}, then the transformed data 
#'                                       for \code{newdata} are included in a 
#'                                       return object.
#' 
#' @return A list containing the following components is returned.
#' 
#' \item{distance}{Vector with length n. Distances from the unit space to each 
#'                  sample.}
#' \item{le_threshold}{Vector with length n. Logical values indicating the 
#'                      distance of each sample is less than or equal to the 
#'                      threhold value (\code{TRUE}) or not (\code{FALSE}).}
#' \item{threshold}{Numeric value to classify the sample into positive or 
#'                   negative. }
#' \item{unit_space}{Object passed by \code{unit_space}.}
#' \item{n}{The number of samples for \code{newdata}.}
#' \item{q}{The number of independent variables after the data transformation. 
#'           According to the data transoformation function, q may be equal to p.}  
#' \item{x}{If \code{includes_transformed_newdata} is \code{TRUE}, then the 
#'           transformed data for \code{newdata} are included.}
#'   
#' @seealso \code{\link{diagnosis.MT}}, \code{\link{diagnosis.MTA}}, and 
#'            \code{\link{diagnosis.RT}}
#' 
#' @examples                           
#' # 40 data for versicolor in the iris dataset                            
#' iris_versicolor <- iris[61:100, -5] 
#' 
#' # The following settings are same as the MT method.                          
#' unit_space <- general_MT(unit_space_data = iris_versicolor, 
#'                          generates_transform_function = 
#'                                             generates_normalization_function,
#'                          calc_A = function(x) solve(cor(x)),  
#'                          includes_transformed_data = TRUE)
#'
#' # 10 data for each kind (setosa, versicolor, virginica) in the iris dataset                         
#' iris_test <- iris[c(1:10, 51:60, 101:111), -5]
#'                          
#' diagnosis <- general_diagnosis.MT(unit_space = unit_space, 
#'                                   newdata = iris_test, 
#'                                   threshold = 4,
#'                                   includes_transformed_newdata = TRUE)
#'                               
#' (diagnosis$distance)
#' (diagnosis$le_threshold)                          
#' 
#' @export
general_diagnosis.MT <- function(unit_space, 
                               newdata, 
                               threshold,
                               includes_transformed_newdata = FALSE) {

  if (missing(threshold)) {
    warning("Parameter \"threshold\" is missing. The threshold value will be NA.")
  }
  
  if (!missing(newdata)) {
    
    check_data(newdata)

    # For the case newdata is given as vector when the number of sample = 1.
    if (is.vector(newdata)) {
      newdata <- matrix(newdata, ncol = length(newdata))
    }
    
    transformed_newdata <- unit_space$transforms_data(newdata)
    
    x <- as.matrix(transformed_newdata)
    distance <- sqrt(rowSums((x %*% unit_space$A) * x) / ncol(x))

  } else {
    
    distance <- unit_space$distance
      
  }

  if (missing(threshold)) {
    threshold <- NA
  }
  
  le_threshold <- (distance <= threshold)

  diagnosis <- list(distance = distance, le_threshold = le_threshold, 
                     threshold = threshold)
  
  if (!missing(newdata)) {
    diagnosis$n <- nrow(transformed_newdata)
    diagnosis$q <- ncol(transformed_newdata)
  
    if (includes_transformed_newdata) {
      diagnosis$x <- transformed_newdata
    }
  } else {
    diagnosis$n <- unit_space$n
    diagnosis$q <- unit_space$q
    
    if (includes_transformed_newdata) {
      if (!is.null(unit_space$x)){
        diagnosis$x <- unit_space$x
      } else {
        warning(paste("Data values will not be included in the return object,",
                      "because data values are not included in the unit_space",
                      "object."))
      }
    }
  }
  
  diagnosis$unit_space <- unit_space
  
  return(diagnosis)
  
}
