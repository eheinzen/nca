#' Loss functions
#'
#' @param x,y Vectors whose loss to calculate
#' @param tol A tolerance below which things are considered neighbors and above
#'  which is considered error
#' @param ... Other arguments (not in use)
#' @name nca_loss
NULL


#' @rdname nca_loss
#' @export
loss_misclassification <- function(x, y, ...) {
  # this function is only needed because it has dots
  x != y
}

#' @rdname nca_loss
#' @export
loss_sq_error <- function(x, y, ...) {
  (x - y)^2
}

#' @rdname nca_loss
#' @export
loss_abs_error <- function(x, y, ...) {
  abs(x - y)
}

#' @rdname nca_loss
#' @export
loss_tolerance <- function(x, y, tol, ...) {
  abs(x - y) > tol
}


