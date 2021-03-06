% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{diagnosis}
\alias{diagnosis}
\title{Function to predict a diagnosis for a family of Mahalanobis-Taguchi (MT)
  methods}
\usage{
diagnosis(unit_space, newdata, threshold, includes_transformed_newdata)
}
\arguments{
\item{unit_space}{Object generated as a unit space.}

\item{newdata}{Matrix with n rows (samples) and p columns (variables). The
data are used to calculate the desired distances from the
unit space. All data should be continuous values and should
not have missing values.}

\item{threshold}{Numeric specifying the threshold value to classify each
sample into positive (\code{TRUE}) or negative
(\code{FALSE}).}

\item{includes_transformed_newdata}{If \code{TRUE}, then the transformed data
for \code{newdata} are included in a
return object.}
}
\value{
A list containing the following components is returned.

 \item{distance}{Vector with length n. Distances from the unit space to each
                  sample.}
 \item{le_threshold}{Vector with length n. Logical values indicating the
                      distance of each sample is less than or equal to the
                      threhold value (\code{TRUE}) or not (\code{FALSE}).}
 \item{threshold}{Numeric value to classify the sample into positive or
                   negative.}
 \item{unit_space}{Object passed by \code{unit_space}.}
 \item{n}{The number of samples for \code{newdata}.}
 \item{q}{The number of variables after the data transformation.}
 \item{x}{If \code{includes_transformed_newdata} is \code{TRUE}, then the
           transformed data for \code{newdata} are included.}
}
\description{
\code{diagnosis} is a generic function. For details, see
  \code{\link{diagnosis.MT}}, \code{\link{diagnosis.MTA}},
  \code{\link{diagnosis.RT}} or \code{\link{general_diagnosis.MT}}.
}
\seealso{
\code{\link{diagnosis.MT}}, \code{\link{diagnosis.MTA}}, and
           \code{\link{diagnosis.RT}}
}
