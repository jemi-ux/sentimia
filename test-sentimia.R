test_that("dichotomie trouve la racine de x^3 - 2x - 5", {
  f <- function(x) x^3 - 2 * x - 5
  expect_equal(dichotomie(f, 2, 3), 2.0945515, tolerance = 1e-4)
})

test_that("point_fixe trouve la solution de x = cos(x)", {
  expect_equal(point_fixe(cos, 1), 0.7390851, tolerance = 1e-4)
})

test_that("newton_raphson calcule la racine carrée de 2", {
  f  <- function(x) x^2 - 2
  df <- function(x) 2 * x
  expect_equal(newton_raphson(f, df, 2), sqrt(2), tolerance = 1e-6)
})

test_that("newton_raphson signale une dérivée nulle", {
  f  <- function(x) x^2 + 1
  df <- function(x) 2 * x
  expect_error(newton_raphson(f, df, 0))
})

test_that("colebrook donne le bon coefficient de frottement", {
  lambda <- colebrook(eps = 0.00026, D = 0.1, Re = 1e5)
  expect_equal(lambda, 0.02657, tolerance = 1e-3)
})
