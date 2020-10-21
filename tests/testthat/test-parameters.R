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

test_that("test that parameters can be set", {

  p <- get_parameters(list(beta = 0))
  expect_equal(p$beta, 0)
})

test_that("test that setting an incorrect parameter gives an error", {

  expect_error(
    p <- get_parameters(list(beta = -1.0)),
    '*'
  )
})

test_that("that when no parameters are set, parameter list defaults to sir_model_parameters_defaults", {

  pars <- get_parameters(NULL)
  expect_equal(pars, sir_model_parameters_defaults())

})
