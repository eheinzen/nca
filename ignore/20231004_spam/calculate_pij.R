calculate_pij <- function(AX, ref = NULL, sparse = FALSE, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  N <- ncol(AX)
  FUN <- function(i) {
    tmp <- -colSums((ref - AX[, i])^2)
    if(nr) {
      tmp[i] <- -Inf # can't be a neighbor to oneself
    }
    tmp <- tmp - max(tmp) # this is okay because we're multiplying both top and bottom
    if(unnormalized) return(tmp)
    tmp <- exp(tmp)
    tmp <- tmp / sum(tmp)
    if(!sparse) return(tmp)
    idx <- tmp > 0
    cbind(i = i, j = which(idx), v = tmp[idx])
  }
  if(!sparse) {
    pij <- t(vapply(1:N, FUN, numeric(ncol(ref))))
    stopifnot(!anyNA(pij))
  } else {
    pij <- do.call(rbind, lapply(1:N, FUN))
    pij <- spam::spam(list(i = pij[, "i"], j = pij[, "j"], value = pij[, "v"]), nrow = N, ncol = ncol(ref))
    stopifnot(!anyNA(pij@entries))
  }
  pij
}

calculate_pij2 <- function(AX, ref = NULL, sparse = FALSE, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  N <- ncol(AX)
  FUN <- function(i) {
    tmp <- -colSums((ref - AX[, i])^2)
    if(nr) {
      tmp[i] <- -Inf # can't be a neighbor to oneself
    }
    tmp <- tmp - max(tmp) # this is okay because we're multiplying both top and bottom
    if(unnormalized) return(tmp)
    tmp <- exp(tmp)
    tmp <- tmp / sum(tmp)
    if(!sparse) return(tmp)
    idx <- tmp > 0
    cbind(i = i, j = which(idx), v = tmp[idx])
  }
  if(!sparse) {
    pij <- t(vapply(1:N, FUN, numeric(ncol(ref))))
    stopifnot(!anyNA(pij))
  } else {
    pij <- do.call(rbind, lapply(1:N, FUN))
    pij <- nca_sparse(entries = pij, nrow = N, ncol = ncol(ref))
  }
  pij
}
