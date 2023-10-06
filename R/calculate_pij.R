calculate_pij <- function(AX, ref = NULL, unnormalized = FALSE) {
  if(nr <- is.null(ref)) ref <- AX
  N <- ncol(AX)
  pij <- t(vapply(1:N, function(i) {
    tmp <- -colSums((ref - AX[, i])^2)
    if(nr) {
      tmp[i] <- -Inf # can't be a neighbor to oneself
    }
    if(unnormalized) return(tmp)
    tmp <- tmp - max(tmp) # this is okay because we're multiplying both top and bottom
    tmp <- exp(tmp)
    tmp / sum(tmp)
  }, numeric(ncol(ref))))
  stopifnot(!anyNA(pij))
  pij
}
