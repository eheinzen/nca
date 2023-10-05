dat <- do.call(rbind, lapply(1:10, function(x) iris))
y <- dat$Species
X <- as.matrix(dat[names(iris) != "Species"])
microbenchmark::microbenchmark(
  nca.fit(y = y, X = X, n_components = 3),
  # nca2.fit(y = y, X = X, n_components = 3, sparse = FALSE),
  # nca2.fit(y = y, X = X, n_components = 1, sparse = TRUE),
  nca2.fit(y = y, X = X, n_components = 3, sparse = NA),
  # nca3.fit(y = y, X = X, n_components = 1, sparse = FALSE),
  # nca3.fit(y = y, X = X, n_components = 1, sparse = TRUE),
  # nca3.fit(y = y, X = X, n_components = 1, sparse = NA),
  times = 5L
)


