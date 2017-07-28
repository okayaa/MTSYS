#' Function to predict a diagnosis for a family of Mahalanobis-Taguchi (MT)
#'   methods
#'
#' \code{diagnosis} is a generic function. For details, see
#'   \code{\link{diagnosis.MT}}, \code{\link{diagnosis.MTA}},
#'   \code{\link{diagnosis.RT}} or \code{\link{general_diagnosis.MT}}.
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
#'  \item{distance}{Vector with length n. Distances from the unit space to each
#'                   sample.}
#'  \item{le_threshold}{Vector with length n. Logical values indicating the
#'                       distance of each sample is less than or equal to the
#'                       threhold value (\code{TRUE}) or not (\code{FALSE}).}
#'  \item{threshold}{Numeric value to classify the sample into positive or
#'                    negative.}
#'  \item{unit_space}{Object passed by \code{unit_space}.}
#'  \item{n}{The number of samples for \code{newdata}.}
#'  \item{q}{The number of variables after the data transformation.}
#'  \item{x}{If \code{includes_transformed_newdata} is \code{TRUE}, then the
#'            transformed data for \code{newdata} are included.}
#'
#' @seealso \code{\link{diagnosis.MT}}, \code{\link{diagnosis.MTA}}, and
#'            \code{\link{diagnosis.RT}}
#'
#' @export
diagnosis <- function(unit_space,
                      newdata,
                      threshold,
                      includes_transformed_newdata) {

  UseMethod("diagnosis", unit_space)

}

#' Function to predict a forecasting for a family of Taguchi (T) methods
#'
#' \code{forecasting} is a generic function. For details, see
#'   \code{\link{forecasting.T1}}, \code{\link{forecasting.Ta}},
#'   \code{\link{forecasting.Tb}} or \code{\link{general_forecasting.T}}.
#'
#' @param model Object generated as a model.
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
#' @export
forecasting <- function(model, newdata, includes_transformed_newdata) {
  UseMethod("forecasting", model)
}

#' @importFrom stats complete.cases
check_data <- function(data) {
  if (!all(complete.cases(data))) {
    stop(paste("Some missing values are including in data.",
               "Data without missing are needed."))
  }

  if (sum(!is.numeric(as.matrix(data))) != 0) {
    warning("Data include some values which are not numeric.")
  }

}

# Used at a test of generates_dimensionality_reduction_function
restores_Ys <- function(object) {

  center <- attr(object, "scaled:center")

  attr(object, "dimname") <- NULL
  attr(object, "scaled:center") <- NULL

  return(as.data.frame(t(t(object) + center)))
}
