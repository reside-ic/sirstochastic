test_that("test that a single SIR run gives a valid plot", {

  pars <- get_parameters()
  res <- sirmodel(pars)
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

  expect_equal(length(res$time), 10001)
  expect_equal(length(res$S), 10001)
  expect_equal(length(res$I), 10001)
  expect_equal(length(res$R), 10001)
  expect_false(any(is.na(res)))

})

test_that("there are no infections when beta is 0", {

  pars <- get_parameters()
  pars[["beta"]] <- 0
  res <- sirmodel(pars)
  expect_true(all(res$I[1:10001] == 0))

})

test_that("test that everyone is infected when beta is very high", {

  pars <- get_parameters()
  pars[["beta"]] <- 1e100
  res <- sirmodel(pars)
  expect_true(all(res$S[1:10001] == 0))
})

test_that("test that no one is infected if I is 0 at t = 0", {

  pars <- get_parameters()
  pars[["I0"]] <- 0
  res <- sirmodel(pars)
  ## Susceptible population is never drawn down:
  ##expect_equal(s, array(s[, , 1], c(17, 1, 101)))
  expect_true(all(res$I[1:10001] == 0))

})
