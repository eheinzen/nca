data(iris)

test_that("matrix and formula interfaces give the same answer", {
  y <- iris$Species
  X <- as.matrix(iris[names(iris) != "Species"])
  nca.iris <- nca.fit(y = y, X = X, n_components = 1)
  nca.iris3 <- nca(Species ~ ., data = iris, n_components = 1)
  expect_equal(nca.iris$coefficients, nca.iris3$coefficients)
})

# -------------------------------------------------------------------------


test_that("two neighborhoods work", {
  y <- iris$Species
  X <- as.matrix(iris[names(iris) != "Species"])
  nca.iris <- nca.fit(y = y, X = X, n_components = 1)
  nca.iris2 <- nca.fit(y = c(y, y), X = rbind(X, X), n_components = 1,
                       neighborhood = rep(0:1, each = nrow(iris)))

  expect_equal(nca.iris$coefficients, nca.iris2$coefficients, tolerance = 1e-6)

  nca.iris3 <- nca(Species ~ ., data = rbind(iris, iris), n_components = 1, neighborhood = rep(0:1, each = nrow(iris)))
  expect_error(predict(nca.iris2, newdata = X), "Please specify which neighborhood")
  expect_error(predict(nca.iris3, newdata = iris), "Please specify which neighborhood")

  expect_error(predict(nca.iris3, newdata = iris, neighborhood = rep.int(2, nrow(iris))), "Detected new neighborhoods")

  pred <- predict(nca.iris3, newdata = iris, neighborhood = rep.int(0, nrow(iris)))
  expect_identical(dim(pred), c(nrow(iris), 3L))

  pred <- predict(nca.iris2, newdata = X, neighborhood = rep.int(0, nrow(iris)))
  expect_identical(dim(pred), c(nrow(iris), 3L))
})



# -------------------------------------------------------------------------



test_that("prediction works for data.frame", {
  nca.iris <- nca(Species ~ ., data = iris, n_components = 1)
  expect_false(all(fitted(nca.iris) == predict(nca.iris, newdata = iris)))

  nca.iris <- nca(Sepal.Length ~ ., data = iris, n_components = 1)
  expect_error(predict(nca.iris, newdata = iris), NA)

  d <- data.frame(Sepal.Width = 3.5, Petal.Length = 1.4, Petal.Width = 0.2)
  expect_error(predict(nca.iris, newdata = d), "'Species' not found")
  lvls <- sort(unique(iris$Species))
  d$Species <- factor("Not a species", levels = c("Not a species", lvls))
  expect_error(predict(nca.iris, newdata = d), "factor Species has new level Not a species")

  d$Species <- factor("setosa", levels = lvls)
  d$Sepal.Width <- "1"
  expect_error(predict(nca.iris, newdata = d), "variable 'Sepal.Width' was fitted with type \"numeric\" but type \"character\" was supplied")
})

test_that("prediction works for matrix", {
  nca.iris <- nca.fit(y = iris$Species, X = as.matrix(iris[names(iris) != "Species"]), n_components = 1)

  expect_error(predict(nca.iris, newdata = as.matrix(iris[names(iris) != "Species"])), NA)
  expect_error(predict(nca.iris, newdata = matrix(1, 1, 1)), "must be compatible")

})



# -------------------------------------------------------------------------


test_that("custom inits work", {
  expect_error(nca(Species ~ ., data = iris, init = matrix(0, 1, 4)), "can't be all zero")
  expect_error(nca(Species ~ ., data = iris, init = matrix(1, 2, 4)), "distinct")
  expect_error(nca(Species ~ ., data = iris, init = cbind(1:0, 0:1, 0, 0)), NA)
})



# -------------------------------------------------------------------------


test_that("A warning is thrown if it doesn't converge", {
  expect_warning(nca(Species ~ ., data = iris, n_components = 1, optim_control = list(maxit = 1L)), "didn't converge")
})



# -------------------------------------------------------------------------


test_that("custom loss function works", {
  expect_message(nca(Species ~ ., data = iris, n_components = 1, loss = `!=`), "You passed a custom")
  expect_message(nca(Species ~ ., data = iris, n_components = 1, loss = `!=`, mode = "classification"), NA)
})


