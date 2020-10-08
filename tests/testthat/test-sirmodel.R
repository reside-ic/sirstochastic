test_that("test correct default parameters returned from get_parameters()", {

  pars <- get_parameters()
  expect_identical(pars$beta, 0.5)
  expect_identical(pars$nu, 0.3)
  expect_identical(pars$mu, 0.001)
  expect_identical(pars$prop_immune, 0)
  expect_identical(pars$I0_at_steady_state, 0)
  expect_identical(pars$N, 10000)
  expect_identical(pars$I0, 5)
  expect_identical(pars$dt, 0.01)

})

test_that("test that func sirmodel when given no data runs with default data", {

  pars <- list()
  res <- sirmodel(pars)
  expect_true(!is.null(res))

})

test_that("test that a single SIR model can be run and give a valid plot with default data, N = 10000", {

  pars <- get_parameters()
  res <- sirmodel(pars)
  check <- displaythemodel(res)
  expect_equal(check$check, NULL)

})

test_that("test that a single SIR model can be run and give a valid plot for N = 1000", {

  pars <- get_parameters()
  res <- sirmodel(pars)
  pars[["N"]] <- 1000
  check <- displaythemodel(res)
  expect_equal(check$check, NULL)

})

test_that("test that a single SIR run gives increasing R and decreasing S", {

  pars <- get_parameters()
  res <- sirmodel(pars)

  for( i in length(res$S))
  {
    if(i > 0){
      expect_gte(res$S[i], res$S[i-1])
    }
  }

  for( i in length(res$R))
  {
    if(i > 0){
      expect_gte(res$R[i-1], res$S[i])
    }
  }

})

test_that("test for non-valid numbers and number of data points", {

  pars <- get_parameters()
  res <- sirmodel(pars)
  NT <- pars$N + 1
  expect_equal(length(res$time), NT)
  expect_equal(length(res$S), NT)
  expect_equal(length(res$I), NT)
  expect_equal(length(res$R), NT)
  expect_false(any(is.na(res)))

})

test_that("there are no infections when beta is 0", {

  pars <- get_parameters()
  pars[["beta"]] <- 0
  res <- sirmodel(pars)
  NT <- pars$N + 1
  expect_true(all(res$I[1:NT] <= 5))

})

test_that("test that everyone is infected when beta is very high", {

  pars <- get_parameters()
  pars[["beta"]] <- 1e100
  res <- sirmodel(pars)
  NT <- pars$N + 1
  expect_true(all(res$S[2:NT] == 0))

})

test_that("test that no one is infected if I is 0 at t = 0", {

  pars <- get_parameters()
  pars[["I0"]] <- 0
  res <- sirmodel(pars)
  NT <- pars$N + 1
  expect_true(all(res$I[1:NT] == 0))

})
