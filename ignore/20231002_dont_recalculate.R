dat <- do.call(rbind, lapply(1:5, function(x) cbind(iris, iris[names(iris) != "Species"],
                                                    iris[names(iris) != "Species"], iris[names(iris) != "Species"],
                                                    iris[names(iris) != "Species"], iris[names(iris) != "Species"])))
y <- dat$Species
X <- as.matrix(dat[names(dat) != "Species"])
microbenchmark::microbenchmark(
  diag(nca.fit(y = y, X = X, n_components = "diagonal")$coefficients),
  nca.fit2(y = y, X = X, n_components = "diagonal")$coefficients,
  times = 10L,
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
