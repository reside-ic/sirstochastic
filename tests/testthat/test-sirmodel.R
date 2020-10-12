test_that("test correct default parameters returned from get_parameters()", {

  pars <- get_parameters()
  expect_identical(pars$beta, 0.5)
  expect_identical(pars$nu, 0.3)
  expect_identical(pars$mu, 0.001)
  expect_identical(pars$prop_immune, 0)
  expect_identical(pars$N, 10000)
  expect_identical(pars$I0, 5)
  expect_identical(pars$dt, 0.01)

})

test_that("can set parameters", {
  p <- get_parameters(list(beta = 0))
  expect_equal(p$beta, 0)
})

test_that("test that func sirmodel when given no data runs with default data", {

  pars <- list()
  end_time <- 100
  res <- sirmodel(end_time, pars)
  expect_true(!is.null(res))

})

test_that("that when no parameters are set parameter list defaults to sir_model_parameters_defaults", {
  end_time <- 100
  pars <- get_parameters(NULL)
  expect_equal(pars, sir_model_parameters_defaults())

})

test_that("test overrides works for parameter list", {

  dfs <- run_with_repetitions(
    100,
    1
  )

  expect_true(is.data.frame(dfs) == TRUE)

})

test_that("test run_simulation_with_repetitions works with multiple runs", {

  dfs <- run_with_repetitions(
    100,
    10
  )

  expect_true(is.data.frame(dfs) == TRUE)

})

test_that("test that a single SIR model can be run and give a valid plot for N = 1000", {

  pars <- get_parameters()
  end_time <- 100
  pars[["N"]] <- 1000
  res <- sirmodel(end_time, pars)
  expect_true(all(res$S + res$I + res$R == 1000))
  expect_equal(max(res$time),100)

})

test_that("test that a single SIR run gives increasing R and decreasing S", {

  pars <- get_parameters()
  pars[["mu"]] <- 0
  end_time <- 100
  res <- sirmodel(end_time, pars)

  expect_true(all(diff(res$S) <= 0))
  expect_true(all(diff(res$R) >= 0))

})

test_that("test for non-valid numbers and number of data points", {

  pars <- get_parameters()
  end_time <- 100
  res <- sirmodel(end_time, pars)

  expect_equal(length(res$S), length(res$time))
  expect_equal(length(res$I), length(res$time))
  expect_equal(length(res$R), length(res$time))
  expect_false(any(is.na(res)))

})

test_that("there are no infections when beta is 0", {

  pars <- get_parameters()
  pars[["beta"]] <- 0
  end_time <- 100
  res <- sirmodel(end_time, pars)
  expect_true(all(res$I[1:length(res$time)] <= 5))

})

test_that("test that everyone is infected when beta is very high", {

  pars <- get_parameters()
  pars[["beta"]] <- 1e100
  pars[["mu"]] <- 0
  end_time <- 100
  res <- sirmodel(end_time, pars)
  expect_true(all(res$S[2:length(res$time)] == 0))

})

test_that("test that no one is infected if I is 0 at t = 0", {

  pars <- get_parameters()
  pars[["I0"]] <- 0
  end_time <- 100
  res <- sirmodel(end_time, pars)
  expect_true(all(res$I[1:length(res$time)] == 0))

})
