#' Function to generate the data normalization function
#'
#' \code{generates_normalization_function} returns the data normalization
#'   function. The data normalization function is generated based on
#'   \code{unit_space_data}.
#'
#' @param unit_space_data Matrix with n rows (samples) and p columns (variables).
#'                          Data to generate the unit space. All data should be
#'                          continuous values and should not have missing values.
#' @param unit_space_center Vector with length p. The values are subtrahends in
#'                            normalization. If missing, the mean for each
#'                            column of \code{unit_space_data} is used for
#'                            normalization.
#' @param unit_space_scale Vector with length p. The values are divisors in
#'                           normalization. If missing and \code{is_scaled} is
#'                           \code{TRUE}, then the unbiased standard deviation
#'                           for each column of \code{unit_space_data} is used
#'                           for normalization.
#' @param is_scaled Logical. If \code{TRUE} (default value), normalization is
#'                    conducted by subtracting \code{unit_space_center} and
#'                    dividing by \code{unit_space_scale}. If \code{FALSE},
#'                    normalization is conducted by subtracting
#'                    \code{unit_space_center} only.
#' @return Function is returned which takes an n x p matrix as an (only)
#'           argument and returns a normalized n x p matrix. The normalization
#'           is conducted based on \code{unit_space_data}.
#'
#' @seealso \code{\link{MT}} and \code{\link{MTA}}
#'
#' @examples
#' # 40 data for versicolor in the iris dataset
#' iris_versicolor <- iris[61:100, -5]
#'
#' normalizes_data <- generates_normalization_function(iris_versicolor)
#'
#' is.function(normalizes_data) # TRUE
#'
#' @importFrom stats sd
#' @export
generates_normalization_function <- function(unit_space_data,
                                             unit_space_center,
                                             unit_space_scale,
                                             is_scaled = TRUE) {

  if (!missing(unit_space_data) && missing(unit_space_center)) {
    unit_space_center <- apply(unit_space_data, 2, mean)
  }

  if (!missing(unit_space_data) && missing(unit_space_scale) && is_scaled) {
    unit_space_scale <- apply(unit_space_data, 2, sd)
  }

  normalizes_data <- function(data) {

    if (is.vector(data)) {
      data <- matrix(data, ncol = length(data))
    }

    if (is_scaled) {
      normalized_data <- scale(data, unit_space_center, unit_space_scale)
    } else {
      normalized_data <- scale(data, unit_space_center, FALSE)
    }

    return(normalized_data)
  }

  return(normalizes_data)

}

#' Function to generate a data transformation function for the
#'   Recognition-Taguchi (RT) method
#'
#' \code{generates_dimensionality_reduction_function} returns the data
#'   transformation function for the Recognition-Taguchi (RT) method based on
#'   the \code{unit_space_data}. The function reduces the dimensionality of data
#'   into 2 synthetic variables.
#'
#' @param unit_space_data Matrix with n rows (samples) and p columns (variables).
#'                          Data to generate the unit space. All data should be
#'                          continuous values and should not have missing values.
#'
#' @return Function is returned which takes an n x p matrix as an (only)
#'           argument and returns a dimensionality-reduced n x 2 data frame with
#'           named columns; Y_1 and Y_2.
#'
#' @references
#'   Taguchi, G. (2006). Objective Function and Generic Function (11).
#'     \emph{Journal of Quality Engineering Society, 14}(2), 5-9. (In Japanese)
#'
#'   Huda, F., Kajiwara, I., Hosoya, N., & Kawamura, S. (2013). Bolt loosening
#'     analysis and diagnosis by non-contact laser excitation vibration tests.
#'     \emph{Mechanical systems and signal processing, 40}(2), 589-604.
#'
#' @seealso \code{\link{RT}}
#'
#' @examples
#' # 40 data for versicolor in the iris dataset
#' iris_versicolor <- iris[61:100, -5]
#'
#' reduces_dimensionality <-
#'                  generates_dimensionality_reduction_function(iris_versicolor)
#'
#' is.function(reduces_dimensionality) # TRUE
#'
#' @export
generates_dimensionality_reduction_function <- function(unit_space_data) {

  generates_dimension_reduction_data <- function(data, center) {

    r <- sum(center ^ 2)
    L <- data %*% center
    S_T <- apply(data ^ 2, 1, sum)
    S_beta <- L ^ 2 / r
    S_e <- S_T - S_beta
    V_e <- S_e / (ncol(data) - 1)

    return(data.frame(Y_1 = L / r, Y_2 = sqrt(V_e)))

  }

  unit_space_data <- as.matrix(unit_space_data)

  unit_space_center <- apply(unit_space_data, 2, mean)
  unit_space_dimension_reduction_data <-
          generates_dimension_reduction_data(unit_space_data, unit_space_center)

  unit_space_dimension_reduction_center <-
                             apply(unit_space_dimension_reduction_data, 2, mean)

  reduces_dimensionality <- function(data) {

    if (is.vector(data)) {
      data <- matrix(data, ncol = length(data))
    }

    data <- as.matrix(data)

    dimension_reduction_data <-
                     generates_dimension_reduction_data(data, unit_space_center)

    return(scale(x = dimension_reduction_data,
                 center = unit_space_dimension_reduction_center,
                 scale = FALSE))

  }

  return(reduces_dimensionality)

}
