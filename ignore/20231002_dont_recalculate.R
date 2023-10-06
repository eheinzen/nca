dat <- do.call(rbind, lapply(1:20, function(x) iris))
y <- dat$Species
X <- as.matrix(dat[names(iris) != "Species"])
microbenchmark::microbenchmark(
  nca.fit(y = y, X = X, n_components = 3),
  nca2.fit(y = y, X = X, n_components = 3),
  times = 50L,
  check = "equal"
)


nca.iris <- nca.fit(y = y, X = X, n_components = 3)

microbenchmark::microbenchmark(
  calculate_pij(nca.iris$projected_X[[1]]),
  # calculate_pij6(nca.iris$projected_X[[1]]),
  calculate_pij3(nca.iris$projected_X[[1]]),
  # calculate_pij7(nca.iris$projected_X[[1]]),
  times = 100,
  check = "equal"
)
