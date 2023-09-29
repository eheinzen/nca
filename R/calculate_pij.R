calculate_pij <- function(AX, unnormalized = FALSE) {
  N <- ncol(AX)
  pij <- t(vapply(1:N, function(i) {
    tmp <- -colSums((AX - AX[, i])^2)
    tmp[i] <- -Inf
    tmp <- tmp - max(tmp) # this is okay because we're multiplying both top and bottom
    if(unnormalized) return(tmp)
    tmp <- exp(tmp)
    tmp / sum(tmp)
  }, numeric(N)))
  stopifnot(!anyNA(pij))
  pij
}
