nca_fitted_regression <- function(pp, yy) {
  stopifnot(ncol(pp) == length(yy))
  rowSums(pp * matrix(yy, nrow = nrow(pp), ncol = ncol(pp), byrow = TRUE))
}

nca_fitted_classification <- function(pp, yy, classes) {
  stopifnot(ncol(pp) == length(yy))
  vapply(classes, function(yyy) {
    rowSums(pp[, yy == yyy])
  }, numeric(nrow(pp)))
}

#' @rdname nca
#' @export
predict.nca <- function(object, newdata, ..., neighborhood = NULL, na.action = stats::na.pass) {
  Call <- match.call()
  if(missing(newdata) || is.null(newdata)) {
    return(object$fitted)
  }

  if(is.data.frame(newdata)) {
    if(is.null(object$terms)) stop("The 'nca' object was not originally built with a formula and data.frame. Make sure 'newdata=' is a matrix.")
    Terms <- stats::delete.response(object$terms)

    idx <- match(c("object", "newdata", "neighborhood", "na.action"), names(Call), nomatch = 0)
    tmp.call <- Call[c(1, idx)]
    names(tmp.call)[2:3] <- c("formula", "data")
    tmp.call[[1L]] <- quote(stats::model.frame)
    tmp.call$formula <- Terms
    tmp.call$xlev <- object$xlevels
    mf <- eval(tmp.call, parent.frame())
    neighborhood <- mf[["(neighborhood)"]]
    if(is.null(neighborhood) && length(object$neighborhood_names) == 1) {
      neighborhood <- rep.int(object$neighborhood_names, nrow(mf))
    } else if(is.null(neighborhood)) stop("Please specify which neighborhood the new data belongs to")

    if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, mf)
    X <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)

    stopifnot(
      ncol(X) == ncol(object$coefficients)
    )

  } else if(!is.matrix(newdata)) {
    stop("'newdata' must be a data.frame or a matrix")
  } else if(!is.null(object$terms)) {
    warning("It is recommended that 'newdata' be a data.frame if the original 'nca' object was fit with one")
  } else {
    X <- newdata
    if(ncol(X) != ncol(object$coefficients)) {
      stop("The size (number of columns) of 'newdata' must be compatible with the data used to train the 'nca' object")
    }
    if(is.null(neighborhood) && length(object$neighborhood_names) == 1) {
      neighborhood <- rep.int(object$neighborhood_names, nrow(X))
    } else if(is.null(neighborhood)) stop("Please specify which neighborhood the new data belongs to")
  }

  if(!all(neighborhood %in% object$neighborhood_names)) {
    stop("Detected new neighborhoods: ", paste0(setdiff(neighborhood, object$neighborhood_names), collapse = ", "))
  }

  idx <- split(seq_len(nrow(X)), f = neighborhood, drop = TRUE)
  Xsplit <- split.data.frame(X, f = neighborhood, drop = TRUE)

  AX <- lapply(Xsplit, function(xx) object$coefficients %*% t(xx))
  pij <- Map(AX, ref = object$projected_X[names(Xsplit)], f = calculate_pij)

  reord <- order(do.call(c, idx))
  if(object$classification) {
    classes <- object$classes
    names(classes) <- classes
    fitted <- Map(pij, object$y[names(Xsplit)], f = nca_fitted_classification, MoreArgs = list(classes = classes))
    fitted <- do.call(rbind, unname(fitted))
    fitted <- fitted[reord, ]
  } else {
    fitted <- Map(pij, object$y[names(Xsplit)], f = nca_fitted_regression)
    fitted <- do.call(c, unname(fitted))
    fitted <- fitted[reord]
  }

  fitted
}


