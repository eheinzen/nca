
calculate_pij2 <- function(AX, ref = NULL, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  N <- ncol(AX)
  Nj <- ncol(ref)
  base <- t(AX) %*% ref
  if(nr) {
    idx <- 1:N
    idx <- cbind(idx, idx)
    # stopifnot(isSymmetric(base))
    # mj <- diag(base)
    mj <- base[idx] ## should be the fastest
    # mj <- colSums(ref^2)

  } else {
    mj <- colSums(ref^2)
  }

  # (AX_i - AX_j) %*% (Ax_i - AX_j) = AX_i %*% AX_i - AX_j %*% AX_i - AX_j %*% AX_i + AX_j %*% AX_j
  # mi <- colSums(AX^2) ## we're subtracting the max from rows anyway; don't subtract twice
  tmp <-  2*base - matrix(mj, nrow = N, ncol = Nj, byrow = TRUE) #+ matrix(mi, N, Nj)
  if(nr) {
    # diag(tmp) <- -Inf # can't be a neighbor to oneself
    tmp[idx] <- -Inf ## should be faster
  }
  mx <- apply(tmp, 1, max)
  tmp <- tmp - matrix(mx, nrow = N, ncol = Nj)
  if(unnormalized) return(tmp)
  tmp <- exp(tmp)
  tmp <- tmp / matrix(rowSums(tmp), nrow = N, ncol = Nj)
  stopifnot(!anyNA(tmp))
  tmp
}


calculate_pij3 <- function(AX, ref = NULL, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX

  # here we're doing j, i until the very end
  N <- ncol(AX)
  Nj <- ncol(ref)
  base <- t(ref) %*% AX
  if(nr) {
    idx <- 1:N
    idx <- cbind(idx, idx)
    # stopifnot(isSymmetric(base))
    # mj <- diag(base)
    mj <- base[idx] ## should be the fastest
    # mj <- colSums(ref^2)

  } else {
    mj <- colSums(ref^2)
  }

  # mi <- colSums(AX^2) ## we're subtracting the max from columns anyway; don't subtract twice
  tmp <-  2*base - mj#matrix(mj, nrow = N, ncol = Nj) #+ matrix(mi, N, Nj, byrow = T)
  if(nr) {
    # diag(tmp) <- -Inf # can't be a neighbor to oneself
    tmp[idx] <- -Inf ## should be faster
  }
  tmp <- apply(tmp, 2, function(xx) {
    # xx <- 2*xx - mj - mi[i]
    xxx <- exp(xx - max(xx))
    xxx / sum(xxx)
  })
  stopifnot(!anyNA(tmp))
  t(tmp)
}


calculate_pij4 <- function(AX, ref = NULL, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  N <- ncol(AX)
  Nj <- ncol(ref)
  base <- t(ref) %*% AX
  if(nr) {
    idx <- 1:N
    idx <- cbind(idx, idx)
    # stopifnot(isSymmetric(base))
    # mj <- diag(base)
    mj <- base[idx] ## should be the fastest
    # mj <- colSums(ref^2)

  } else {
    mj <- colSums(ref^2)
  }

  # mi <- colSums(AX^2) ## we're subtracting the max from rows anyway; don't subtract twice
  tmp <-  2*base - mj#matrix(mj, nrow = N, ncol = Nj) #+ matrix(mi, N, Nj, byrow = TRUE)
  if(nr) {
    # diag(tmp) <- -Inf # can't be a neighbor to oneself
    tmp[idx] <- -Inf ## should be faster
  }
  mx <- apply(tmp, 2, max)
  tmp <- tmp - matrix(mx, nrow = N, ncol = Nj, byrow = TRUE)
  if(unnormalized) return(tmp)
  tmp <- exp(tmp)
  tmp <- tmp / matrix(colSums(tmp), nrow = N, ncol = Nj, byrow = TRUE)
  stopifnot(!anyNA(tmp))
  t(tmp)
}


calculate_pij5 <- function(AX, ref = NULL, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  N <- ncol(AX)
  Nj <- ncol(ref)
  base <- t(AX) %*% ref
  if(nr) {
    idx <- 1:N
    idx <- cbind(idx, idx)
    # stopifnot(isSymmetric(base))
    # mj <- diag(base)
    mj <- base[idx] ## should be the fastest
    # mj <- colSums(ref^2)

  } else {
    mj <- colSums(ref^2)
  }

  # mi <- colSums(AX^2) ## we're subtracting the max from rows anyway; don't subtract twice
  tmp <-  2*base - matrix(mj, nrow = N, ncol = Nj, byrow = TRUE) #+ matrix(mi, N, Nj)
  if(nr) {
    # diag(tmp) <- -Inf # can't be a neighbor to oneself
    tmp[idx] <- -Inf ## should be faster
  }
  mx <- apply(tmp, 1, max)
  tmp <- tmp - mx#matrix(mx, nrow = N, ncol = Nj)
  if(unnormalized) return(tmp)
  tmp <- exp(tmp)
  tmp <- tmp / rowSums(tmp)#matrix(rowSums(tmp), nrow = N, ncol = Nj)
  stopifnot(!anyNA(tmp))
  tmp
}




calculate_pij6 <- function(AX, ref = NULL, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  if(nrow(ref) != nrow(AX)) stop("'AX' and 'ref' must have the same number of rows")
  # here we're doing j, i until the very end
  N <- ncol(AX)
  Nj <- ncol(ref)

  diff <- array(AX, dim = c(nrow(AX), N, Nj)) -
    aperm(array(ref, dim = c(nrow(ref), Nj, N)), c(1, 3, 2))
  tmp <- -colSums(diff^2)
  if(nr) {
    # diag(tmp) <- -Inf # can't be a neighbor to oneself
    idx <- 1:N
    idx <- cbind(idx, idx)
    tmp[idx] <- -Inf ## should be faster
  }
  tmp <- apply(tmp, 1, function(xx) {
    # xx <- 2*xx - mj - mi[i]
    xxx <- exp(xx - max(xx))
    xxx / sum(xxx)
  })
  stopifnot(!anyNA(tmp))
  t(tmp)
}


calculate_pij7 <- function(AX, ref = NULL, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  N <- ncol(AX)
  Nj <- ncol(ref)
  base <- t(AX) %*% ref
  if(nr) {
    idx <- 1:N
    idx <- cbind(idx, idx)
    # stopifnot(isSymmetric(base))
    # mj <- diag(base)
    mj <- base[idx] ## should be the fastest
    # mj <- colSums(ref^2)

  } else {
    mj <- colSums(ref^2)
  }

  # mi <- colSums(AX^2) ## we're subtracting the max from columns anyway; don't subtract twice
  #tmp <-  2*base - matrix(mj, nrow = N, ncol = Nj) #+ matrix(mi, N, Nj, byrow = T)
  if(nr) {
    # diag(tmp) <- -Inf # can't be a neighbor to oneself
    base[idx] <- -Inf ## should be faster
  }
  tmp <- apply(2*base, 1, function(xx) {
    # xx <- 2*xx - mj - mi[i]
    xx <- xx - mj
    xxx <- exp(xx - max(xx))
    xxx / sum(xxx)
  })
  stopifnot(!anyNA(tmp))
  t(tmp)
}
