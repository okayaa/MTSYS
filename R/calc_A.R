#' Function to calculate a cofactor matrix
#' 
#' \code{calc_cofactor} calculates a cofactor matrix.
#' 
#' @param data Matrix with n rows (samples) and p columns (variables). All data 
#'               should be continuous values and should not have missing values.  
#' 
#' @return \code{calc_cofactor} returns a cofactor matrix of size p x p.
#'           
#' @seealso \code{\link{MTA}}
#' 
#' @examples
#' # 40 data for versicolor in the iris dataset                            
#' iris_versicolor <- iris[61:100, -5] 
#'                             
#' calc_cofactor(cov(iris_versicolor))
#'  
#' @export
calc_cofactor <- function(data) {  
  p <- ncol(data)
  
  #Function "det (determinant)" does not correspond to a 1 x 1 matrix.
  #Function "det1" returns the element itself as the determinant in the case.
  det1 <- function(x, ...) {
    if (length(x) == 1) {
      D <- x
    } else {
      D <- det(x)
    }
    return(D)
  }
  
  cofactor <- 
        outer(seq_len(p), seq_len(p), 
              Vectorize(function(i, j) (-1) ^ (i + j) * det1(data[-i, -j])))
  
  return(cofactor)
}
