f  <- function(x) x^2 - 2
df <- function(x) 2 * x

test_that("la bissection trouve racine de 2", {
  r <- bissection(f, 1, 2)
  expect_equal(r[["racine"]], sqrt(2), tolerance = 1e-8)
})

test_that("la bissection refuse deux bornes de meme signe", {
  expect_error(bissection(function(x) x^2 + 1, 0, 1))
})

test_that("Newton trouve racine de 2 en peu d'iterations", {
  r <- newton(f, df, 1)
  expect_equal(r[["racine"]], sqrt(2), tolerance = 1e-8)
  expect_lte(r[["iterations"]], 6)
})

test_that("la secante trouve racine de 2", {
  r <- secante(f, 1, 2)
  expect_equal(r[["racine"]], sqrt(2), tolerance = 1e-8)
})

test_that("le point fixe trouve la solution de cos x = x", {
  r <- point_fixe(cos, 1)
  expect_equal(r[["racine"]], 0.7390851332, tolerance = 1e-7)
})
