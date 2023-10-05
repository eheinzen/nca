
#' @rdname nca
#' @export
nca3.fit <- function(y, X, n_components, init = c("pca", "identity"), loss = NULL, ...,
                    neighborhood = NULL,
                    lambda = 0, optim_method = "L-BFGS-B", optim_control = list(),
                    sparse = NA, debug = FALSE) {
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
    is.matrix(X),
    is.numeric(X),
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

  if(!missing(init) && is.matrix(init)) {
    stopifnot(
      is.numeric(init),
      ncol(init) == k,
      nrow(init) <= k,
      "The inits can't be all zero for any row" = rowSums(init != 0) >= 1,
      "The rows of the inits should all be distinct" = nrow(unique(init)) == nrow(init)
    )
    A.init <- init
  } else {
    stopifnot(
      n_components >= 1,
      n_components <= k
    )
    init <- match.arg(init)
    if(init == "pca") {
      X.pca <- stats::prcomp(X)
      A.init <- unname(t(X.pca$rotation[, paste0("PC", seq_len(n_components)), drop = FALSE]))
    } else if(init == "identity") {
      A.init <- matrix(0, nrow = n_components, ncol = k)
      idx <- seq_len(min(k, n_components))
      A.init[cbind(idx, idx)] <- 1
    }
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
    sp <- if(is.na(sparse)) {
      if(is.null(env$pij)) FALSE else is_nca_sparse(env$pij[[1]]) || mean(env$pij[[1]] == 0) > 0.7
    } else sparse
    env$pij <- lapply(env$AX, calculate_pij2, sparse = sp)
    return(NULL)
  }

  calculate_objective <- function(A) {
    dim(A) <- dim(A.init)
    calculate_once(A)
    # stopifnot(env$A == A)
    do.call(sum, Map(env$pij, yiyj, f = function(pp, yy) nca_sum(nca_mult(pp, yy))))/N + 0.5*lambda*sum(A^2) # since we're minimizing, we want to *add* the penalty
  }

  calculate_gradient <- function(A) {
    dim(A) <- dim(A.init)
    calculate_once(A)
    # stopifnot(env$A == A)

    AX <- env$AX
    pij <- env$pij
    tmp <- Map(AX, pij, yiyj, Xsplit, f = function(ax, pp, yy, xx) {
      py <- nca_mult(pp, yy)
      pi <- nca_rowSums(py)
      W <- nca_minus(py, nca_mult(pp, matrix(pi, nrow = nrow(yy), ncol = ncol(yy), byrow = FALSE)))
      stopifnot(all.equal(nca_rowSums(W), rep.int(0, nrow(yy))))
      if(is_nca_sparse(W)) {
        do_the_W(AX = ax, X = xx, W = W)
      } else {
        # W <- as.matrix.nca_sparse(W)
        W2 <- W + spam::t(W)
        diag(W2) <- -spam::colSums(W)
        ax %*% W2 %*% xx
      }
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

  if(out$convergence > 0) warning("The algorithm didn't converge")
  cat("Made it")
  A <- out$par
  out$par <- NULL # don't store it twice
  calculate_once(A)

  reord <- order(do.call(c, idx))
  if(classification) {
    classes <- sort(unique(y))
    names(classes) <- classes
    fitted <- Map(env$pij, ysplit, f = nca_fitted_classification2, MoreArgs = list(classes = classes))
    fitted <- do.call(rbind, unname(fitted))
    fitted <- fitted[reord, ]

  } else {
    fitted <- Map(env$pij, ysplit, f = nca_fitted_regression2)
    fitted <- do.call(c, unname(fitted))
    fitted <- fitted[reord]
    classes <- NULL
  }

  structure(list(
    optim = out,
    coefficients = A,
    projected_X = env$AX,
    y = ysplit,
    neighborhood_membership = idx,
    neighborhood_names = names(ysplit),
    fitted = fitted,
    classification = classification,
    classes = classes,
    lambda = lambda
  ), class = "nca")
}

#' @rdname nca
#' @export
print.nca <- function(x, ...) {
  cat(
    "\nAn object of class 'nca': ",
    if(x$classification) paste("classification on", length(x$classes), "classes") else "regression",
    " using ", nrow(x$coefficients), " component(s).\n",
    sep = ""
  )
}

