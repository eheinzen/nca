#' Neighborhood Components Analysis
#'
#' @param formula A formula denoting the response and features.
#' @param data An optional data.frame from which to get the variables in \code{formula}.
#' @param subset An optional vector specifying the subset of rows to use.
#' @param na.action The action to take when \code{NA}s are present.
#' @param ... Other arguments passed to \code{model.matrix} or \code{nca.fit} or \code{loss}.
#' @param y The response vector whose loss is to be measured. This can be a factor
#'   or character vector, for which classification will be performed, or a numeric,
#'   for which regression will be performed.
#' @param X A matrix of N data points (rows) by K features (columns).
#' @param n_components How many components to use for the NCA algorithm.
#' @param neighborhood An indicator of a distinct neighborhood,
#'   for when, a priori, the data is completely separable into distinct groups
#'   (that is, \code{p_ij} can be rearranged to be block diagonal); having this
#'   prior knowledge of which points structurally can and can't be neighbors speeds
#'   up computation
#' @param init How to initialize the transformation matrix
#' @param loss A vectorized function fed to \code{\link{outer}} for determining
#'   the loss between two elements of \code{y}. It is assumed (but not checked)
#'   that the loss is symmetric. For regression, this defaults to \code{\link{loss_sq_error}},
#'   and for classification it defaults to \code{\link{loss_inaccuracy}}.
#' @param lambda A penalty parameter to penalize the transformation matrix back to 0. The penalty applied
#'   is \code{1/2 * lambda * sum(transformation^2)}.
#' @param optim_method The method passed to \code{\link{optim}}.
#' @param optim_control The control passed to \code{\link{optim}}. It can be
#'   useful for, e.g., increasing verbosity of the optimization.
#' @param debug A logical, for debugging
#' @param newdata New data to predict
#' @param object An object of class \code{'nca'}
#' @details
#' This differs from the NCA publication (Goldberger et al.) in a few ways.
#'  First, it uses a vectorized
#' gradient; second, it minimizes loss instead of maximizing accuracy (hence it
#' can support regression); third, it supports a penalty parameter (\code{lambda},
#' as in Yang et al.).
#'
#' @examples
#' library(datasets)
#' data(iris)
#'
#' nca.iris <- nca(Species ~ ., data = iris, n_components = 1)
#' pred <- fitted(nca.iris)
#' table(colnames(pred)[apply(pred, 1, which.max)] == iris$Species) # 2% Leave-One-Out error
#'
#' @references J. Goldberger, G. Hinton, S. Roweis, R. Salakhutdinov.
#' "Neighbourhood Components Analysis". Advances in Neural Information
#' Processing Systems. 17, 513-520, 2005.
#'
#' Yang, W., K. Wang, W. Zuo. "Neighborhood Component Feature Selection for
#' High-Dimensional Data." Journal of Computers. Vol. 7, Number 1, January, 2012.
#' @name nca
NULL

#' @rdname nca
#' @export
nca <- function(formula, data, subset, na.action, ...) {
  Call <- match.call()
  idx <- match(c("formula", "data", "subset", "na.action"), names(Call), nomatch = 0)
  if(idx[1] == 0) stop("A formula argument is required")
  temp.call <- Call[c(1, idx)]
  temp.call[[1L]] <- quote(stats::model.frame)
  temp.call$drop.unused.levels <- TRUE
  mf <- eval(temp.call, parent.frame())
  if(nrow(mf) == 0) stop("No (non-missing) observations")
  Terms <- stats::terms(mf)
  attr(Terms, "intercept") <- FALSE
  naaction <- stats::na.action(mf)

  if(attr(Terms, "response") != 1) stop("You must provide a left-hand side to the formula")
  y <- mf[[1]]
  X <- stats::model.matrix(Terms, data = mf, ...)

  out <- nca.fit(y = y, X = X, ...)
  out$terms <- Terms
  out$na.action <- stats::na.action(mf)
  out$contrasts <- attr(X, "contrasts")
  out$xlevels <- .getXlevels(Terms, mf)
  out
}

#' @rdname nca
#' @export
nca.fit <- function(y, X, n_components, init = c("pca", "identity"), loss = NULL, ...,
                     neighborhood = NULL,
                     lambda = 0, optim_method = "L-BFGS-B", optim_control = list(), debug = FALSE) {
  # set.seed(20230920)
  # X <- matrix(rnorm(1000), 100, 10)
  # y <- drop(X %*% rnorm(10))
  # n_components <- 2
  # lambda = 0
  # loss <- NULL
  # optim_method <- "L-BFGS-B"
  # optim_control <- list(REPORT = 1L, trace = 1)
  # init <- "pca"
  # nca.fit(y = y, X = X, n_components = n_components, optim_control = optim_control)
  stopifnot(
    n_components >= 1,
    is.matrix(X),
    is.numeric(X),
    n_components <= ncol(X),
    length(y) == nrow(X),
    lambda >= 0,
    !anyNA(X),
    !anyNA(y),
    is.null(neighborhood) || length(neighborhood) == nrow(X),
    is.list(optim_control)
  )
  if(debug) {
    if(is.null(optim_control$REPORT)) optim_control$REPORT <- 1
    if(is.null(optim_control$trace)) optim_control$trace <- 1
  }

  if(is.null(loss) && (is.factor(y) || is.character(y))) {
    classification <- TRUE
    loss <- "loss_inaccuracy"
  } else {
    stopifnot(is.numeric(y))
    classification <- FALSE
    loss <- "loss_sq_error"
  }
  loss <- match.fun(loss)

  N <- nrow(X)
  k <- ncol(X)

  has_nbh <- !is.null(neighborhood)
  if(!has_nbh) {
    neighborhood <- rep.int(1, N)
  }
  idx <- split(seq_len(N), f = neighborhood, drop = TRUE)
  ysplit <- split(y, f = neighborhood, drop = TRUE)
  yiyj <- lapply(ysplit, function(yy) outer(yy, yy, FUN = loss, ...))

  init <- match.arg(init)
  if(init == "pca") {
    X.pca <- stats::prcomp(X)
    A.init <- unname(t(X.pca$rotation[, paste0("PC", seq_len(n_components)), drop = FALSE]))
  } else if(init == "identity") {
    A.init <- matrix(0, nrow = n_components, ncol = k)
    idx <- seq_len(min(k, n_components))
    A.init[cbind(idx, idx)] <- 1
  }

  Xsplit <- split.data.frame(X, f = neighborhood, drop = TRUE)

  env <- new.env()
  calculate_once <- function(A, verbose = debug) {
    if(identical(A, env$A)) {
      if(verbose) cat("Already have it\n")
      return(NULL)
    }
    if(verbose) cat("Recalculating\n")
    env$A <- A
    env$AX <- lapply(Xsplit, function(xx) A %*% t(xx))
    env$pij <- lapply(env$AX, calculate_pij)
    return(NULL)
  }

  calculate_objective <- function(A) {
    dim(A) <- dim(A.init)
    calculate_once(A)
    # stopifnot(env$A == A)
    do.call(sum, Map(env$pij, yiyj, f = function(pp, yy) sum(pp * yy)))/N + 0.5*lambda*sum(A^2) # since we're minimizing, we want to *add* the penalty
  }

  calculate_gradient <- function(A) {
    dim(A) <- dim(A.init)
    calculate_once(A)
    # stopifnot(env$A == A)

    AX <- env$AX
    pij <- env$pij
    tmp <- Map(AX, pij, yiyj, Xsplit, f = function(ax, pp, yy, xx) {
      py <- pp * yy
      pi <- rowSums(py)
      W <- py - pp * matrix(pi, nrow = nrow(pp), ncol = ncol(pp), byrow = FALSE)
      stopifnot(all.equal(rowSums(W), rep.int(0, nrow(pp))))
      W2 <- W + t(W)
      diag(W2) <- -colSums(W)
      ax %*% W2 %*% xx
    })
    (2/N)*Reduce(`+`, tmp) + lambda*A
  }

  out <- stats::optim(
    par = A.init,
    fn = calculate_objective,
    gr = calculate_gradient,
    method = optim_method,
    control = optim_control
  )

  A <- out$par
  out$par <- NULL # don't store it twice
  calculate_once(A)

  if(classification) {
    classes <- sort(unique(y))
    names(classes) <- classes
    fitted <- do.call(rbind, Map(env$pij, ysplit, f = function(pp, yy) {
      outer(seq_len(nrow(pp)), classes, FUN = Vectorize(function(i, yyy) {
        sum(pij[i, yy == yyy])
      }))
    }))
    fitted <- fitted[order(do.call(c, idx)), ]

  } else {
    fitted <- Map(env$pij, ysplit, f = function(pp, yy) rowSums(pp * matrix(y, nrow = nrow(pp), ncol = ncol(pp), byrow = TRUE)))
    fitted <- fitted[order(do.call(c, idx))]
    classes <- NULL
  }

  structure(list(
    optim = out,
    coefficients = A,
    projected_X = AX,
    y = y,
    neighborhood = neighborhood,
    neighborhood_names = names(ysplit),
    fitted = fitted,
    classification = classification,
    classes = classes,
    lambda = lambda
  ), class = "nca")
}


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


