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
#' @param init How to intialize the transformation matrix
#' @param loss A vectorized function fed to \code{\link{outer}} for determining
#'   the loss between two elements of \code{y}. It is assumed (but not checked)
#'   that the loss is symmetric. For regression, this defaults to \code{\link{loss_sq_error}},
#'   and for classification it defaults to \code{\link{loss_inaccuracy}}.
#' @param lambda A penalty parameter to penalize the transformation matrix back to 0.
#' @param optim.method The method passed to \code{\link{optim}}.
#' @param optim.control The control passed to \code{\link{optim}}. It can be
#'  useful for, e.g., increasing verbosity of the optimization.
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
  out
}


#' @rdname nca
#' @export
nca.fit <- function(y, X, n_components, init = c("pca", "identity"), loss = NULL,
                    lambda = 0, optim.method = "L-BFGS-B", optim.control = list(), ...) {
  # set.seed(20230920)
  # X <- matrix(rnorm(1000), 100, 10)
  # y <- drop(X %*% rnorm(10))
  # n_components <- 2
  # lambda = 0
  # loss <- NULL
  # optim.method <- "L-BFGS-B"
  # optim.control <- list(REPORT = 1L, trace = 1)
  # init <- "pca"
  # nca.fit(y = y, X = X, n_components = n_components, optim.control = optim.control)
  stopifnot(
    n_components >= 1,
    is.matrix(X),
    n_components <= ncol(X),
    length(y) == nrow(X),
    lambda >= 0,
    !anyNA(X),
    !anyNA(y)
  )
  if(is.null(loss) && (is.factor(y) || is.character(y))) {
    classification <- TRUE
    loss <- "loss_inaccuracy"
  } else {
    stopifnot(is.numeric(y))
    classification <- FALSE
    loss <- "loss_sq_error"
  }
  loss <- match.fun(loss)
  yiyj <- outer(y, y, FUN = loss, ...)
  N <- nrow(X)
  k <- ncol(X)

  init <- match.arg(init)
  if(init == "pca") {
    X.pca <- stats::prcomp(X)
    A.init <- t(X.pca$rotation[, paste0("PC", seq_len(n_components)), drop = FALSE])
  } else if(init == "identity") {
    A.init <- matrix(0, nrow = n_components, ncol = k)
    idx <- seq_len(min(k, n_components))
    A.init[cbind(idx, idx)] <- 1
  }

  calculate_objective <- function(A) {
    dim(A) <- dim(A.init)
    AX <- A %*% t(X)
    pij <- calculate_pij(AX)
    mean(colSums(t(pij) * yiyj)) + lambda*sum(A^2) # since we're minimizing, we want to *add* the penalty
  }

  calculate_gradient <- function(A) {
    dim(A) <- dim(A.init)
    AX <- A %*% t(X)
    pij <- calculate_pij(AX)
    py <- pij * yiyj
    pi <- rowSums(py)
    W <- py - pij * matrix(pi, nrow = N, ncol = N)
    stopifnot(all.equal(rowSums(W), rep.int(0, N)))
    W2 <- W + t(W)
    diag(W2) <- -colSums(W)
    (2/N)*AX %*% W2 %*% X + 2*lambda*A
  }

  out <- stats::optim(
    par = A.init,
    fn = calculate_objective,
    gr = calculate_gradient,
    method = optim.method,
    control = optim.control
  )

  A <- out$par
  out$par <- NULL # don't store it twice
  AX <- A %*% t(X)
  pij <- calculate_pij(AX)
  if(classification) {
    classes <- sort(unique(y))
    names(classes) <- classes
    fitted <- outer(1:N, classes, FUN = Vectorize(function(i, yy) {
      sum(pij[i, y == yy])
    }))
  } else {
    fitted <- rowSums(pij * matrix(y, nrow = N, ncol = N, byrow = TRUE))
    classes <- NULL
  }

  structure(list(
    optim = out,
    coefficients = A,
    projected_X = AX,
    y = y,
    fitted = fitted,
    classification = classification,
    classes = classes,
    lambda = lambda
  ), class = "nca")
}
