#' Loss functions
#'
#' @param x,y Vectors whose loss to calculate
#' @param tol A tolerance below which things are considered neighbors and above
#'  which is considered error
#' @returns
#' Most of these functions return the loss between the elements of \code{x} and \code{y}.
#'   Note that the misclassification loss is simply \code{!=}.
#'
#' \code{loss_tolerance} is actually a function factory, which returns a loss
#'   function.
#' @name nca_loss
NULL

#' @rdname nca_loss
#' @export
loss_sq_error <- function(x, y) {
  (x - y)^2
}

#' @rdname nca_loss
#' @export
loss_abs_error <- function(x, y) {
  abs(x - y)
}

#' @rdname nca_loss
#' @export
loss_tolerance <- function(tol) {
  function(x, y) {
    abs(x - y) > tol
  }
}


