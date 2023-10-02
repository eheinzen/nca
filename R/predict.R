

#' @rdname nca
#' @export
predict.nca <- function(object, newdata, ..., na.action = stats::na.pass) {
  if(missing(newdata) || is.null(newdata)) {
    return(object$fitted)
  }

  if(is.data.frame(newdata)) {
    if(is.null(object$terms)) stop("The 'nca' object was not originally built with a formula and data.frame. Make sure 'newdata=' is a matrix.")
    Terms <- stats::delete.response(object$terms)
    mf <- stats::model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
    if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, mf)
    X <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)

    stopifnot(
      ncol(X) == ncol(object$coefficients)
    )

  } else if(!is.matrix(X)) {
    stop("'X' must be a data.frame or a matrix")
  } else if(!is.null(object$terms)) {
    warning("It is recommended that 'newdata' be a data.frame if the original 'nca' object was fit with one")
  } else {
    X <- newdata
    stopifnot(
      "The size (number of columns) of 'newdata' must be compatible with the data used to train the 'nca' object" = ncol(X) == ncol(object$coefficients)
    )
  }

  N <- nrow(X)
  AX <- object$coefficients %*% t(X)
  pij <- calculate_pij(AX = AX, ref = object$projected_X)
  if(object$classification) {
    classes <- object$classes
    names(classes) <- classes
    pred <- outer(1:N, classes, FUN = Vectorize(function(i, yy) {
      sum(pij[i, object$y == yy])
    }))
  } else {
    pred <- rowSums(pij * matrix(object$y, nrow = N, ncol = N, byrow = TRUE))
  }

  pred
}


