data(iris)
test_that("prediction works for data.frame", {
  nca.iris <- nca(Species ~ ., data = iris, n_components = 1)
  expect_false(all(fitted(nca.iris) == predict(nca.iris, newdata = iris)))

  nca.iris <- nca(Sepal.Length ~ ., data = iris, n_components = 1)
  expect_error(predict(nca.iris, newdata = iris), NA)

  d <- data.frame(Sepal.Width = 3.5, Petal.Length = 1.4, Petal.Width = 0.2)
  expect_error(predict(nca.iris, newdata = d), "'Species' not found")
  d$Species <- "Not a species"
  expect_error(predict(nca.iris, newdata = d), "factor Species has new level Not a species")

  d$Species <- "setosa"
  d$Sepal.Width <- "1"
  expect_error(predict(nca.iris, newdata = d), "variable 'Sepal.Width' was fitted with type \"numeric\" but type \"character\" was supplied")
})

test_that("prediction works for matrix", {
  nca.iris <- nca.fit(y = iris$Species, X = as.matrix(iris[names(iris) != "Species"]), n_components = 1)

  expect_error(predict(nca.iris, newdata = as.matrix(iris[names(iris) != "Species"])), NA)
  expect_error(predict(nca.iris, newdata = matrix(1, 1, 1)), "must be compatible")

})
