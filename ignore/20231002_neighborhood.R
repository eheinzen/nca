
y <- iris$Species
X <- as.matrix(iris[1:4])
nca.fit(y = y, X = X, n_components = 1)$coefficients
nca.fit2(y = y, X = X, n_components = 1)$coefficients
