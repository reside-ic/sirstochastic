test_that("test that func compartmental_sirmodel when given no data runs with default data", {

  pars <- list()
  end_time <- 100
  res <- compartmental_sirmodel(end_time, pars)
  expect_true(!is.null(res))

})

test_that("test overrides works for parameter list", {

  dfs <- run_with_repetitions(
    100,
    1,
    get_parameters()
  )

  expect_true(is.data.frame(dfs) == TRUE)

})

test_that("test run_simulation_with_repetitions works with multiple runs", {

  dfs <- run_with_repetitions(
    100,
    10,
    get_parameters()
  )

  expect_true(is.data.frame(dfs) == TRUE)

})

test_that("test that a single SIR model can be run and give a valid plot for N = 1000", {

  pars <- get_parameters()
  end_time <- 100
  pars[["N"]] <- 1000
  res <- compartmental_sirmodel(end_time, pars)
  expect_true(all(res$S + res$I + res$R == 1000))
  expect_equal(max(res$time),100)

})

test_that("test that a single SIR run gives increasing R and decreasing S", {

  pars <- get_parameters()
  pars[["mu"]] <- 0
  end_time <- 100
  res <- compartmental_sirmodel(end_time, pars)

  expect_true(all(diff(res$S) <= 0))
  expect_true(all(diff(res$R) >= 0))

})

test_that("test for non-valid numbers and number of data points", {

  pars <- get_parameters()
  end_time <- 100
  res <- compartmental_sirmodel(end_time, pars)

  expect_equal(length(res$S), length(res$time))
  expect_equal(length(res$I), length(res$time))
  expect_equal(length(res$R), length(res$time))
  expect_false(any(is.na(res)))

})

test_that("there are no infections when beta is 0", {

  pars <- get_parameters()
  pars[["beta"]] <- 0
  end_time <- 100
  res <- compartmental_sirmodel(end_time, pars)
  expect_true(all(res$I[1:length(res$time)] <= 5))

})

test_that("test that everyone is infected when beta is very high", {

  pars <- get_parameters()
  pars[["beta"]] <- 1e100
  pars[["mu"]] <- 0
  end_time <- 100
  res <- compartmental_sirmodel(end_time, pars)
  expect_true(all(res$S[2:length(res$time)] == 0))

})

test_that("test that no one is infected if I is 0 at t = 0", {

  pars <- get_parameters()
  pars[["I0"]] <- 0
  end_time <- 100
  res <- compartmental_sirmodel(end_time, pars)
  expect_true(all(res$I[1:length(res$time)] == 0))

})

test_that("test that an empty simulation exits ok for individual model", {
  library(individual)
  library(reshape2)

  population <- 10000

  S <- State$new('S', population)
  I <- State$new('I', 0)
  R <- State$new('R', 0)

  human <- Individual$new('human', list(S, I, R))

  output <- simulate(human, list(), 1)

  true_output <- data.frame(timestep = c(1))

  expect_equal(true_output, output)

  expect_error(
    simulate(human, list(), 0),
    '*'
  )

})

test_that("test individual model with 10000 humans", {
  library(individual)
  library(reshape2)

  pars <- get_parameters()

  population <- pars$N
  NI <- pars$I0
  NR <- 2
  pops <- population - NI - NR
  timestep <- pars$num/pars$dt

  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', NR)

  immunity <- Variable$new('immunity',  rep(0, pars$N))
  age  <- Variable$new('age', rep(0, pars$N))
  location <- Variable$new('location', rep(0, pars$N))
  human <- Individual$new('human', list(S, I, R), variables = list(immunity, age, location))

  processes <- list(
    sirstochastic::individual_S_to_I(S, I, human, immunity, age, location, pars),
    sirstochastic::individual_I_to_R(I, R, human, immunity, age, location, pars),
    sirstochastic::individual_R_to_S(S, R, human, immunity, age, location, pars),
    sirstochastic::render_state_sizes(S, I, R, human)
  )

  output <- simulate(human, processes, timestep)

  df <-   data.frame(S = output$susceptable_counts, I = output$infected_counts, R = output$recovered_counts, time = output$time, type = "Individual",  legend = "Individual", stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))

})


test_that("test individual model with 10000 humans with immunity", {
  library(individual)
  library(reshape2)

  pars <- get_parameters()
  pars[["indludeimmune"]] <- TRUE

  population <- pars$N
  NI <- pars$I0
  NR <- 2
  pops <- population - NI - NR
  timestep <- pars$num/pars$dt

  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', NR)

  immunity <- Variable$new('immunity', runif(population, 0, .5))
  age  <- Variable$new('age', rep(0, pars$N))
  location <- Variable$new('location', rep(0, pars$N))
  human <- Individual$new('human', list(S, I, R), variables = list(immunity, age, location))

  processes <- list(
    sirstochastic::individual_S_to_I(S, I, human, immunity, age, location, pars),
    sirstochastic::individual_I_to_R(I, R, human, immunity, age, location, pars),
    sirstochastic::individual_R_to_S(S, R, human, immunity, age, location, pars),
    sirstochastic::render_state_sizes(S, I, R, human)
  )

  output <- simulate(human, processes, timestep, parameters = list(immunity_level = .6))

  df <-   data.frame(S = output$susceptable_counts, I = output$infected_counts, R = output$recovered_counts, time = output$time, type = "Individual",  legend = "Individual", stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))

})

test_that("test individual model with 10000 humans with immunity and age effects", {
  library(individual)
  library(reshape2)

  pars <- get_parameters()
  pars[["indludeimmune"]] <- TRUE
  pars[["includeage"]] <- TRUE

  population <- pars$N
  NI <- pars$I0
  NR <- 2
  pops <- population - NI - NR
  timestep <- pars$num/pars$dt

  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', NR)

  immunity <- Variable$new('immunity', runif(population, 0, .1))
  rate=1/pars$average_age
  age <- Variable$new('age', rexp(pars$N, rate))
  location <- Variable$new('location', rep(0, pars$N))
  human <- Individual$new('human', list(S, I, R), variables = list(immunity, age, location))

  processes <- list(
    sirstochastic::individual_S_to_I(S, I, human, immunity, age, location, pars),
    sirstochastic::individual_I_to_R(I, R, human, immunity, age, location, pars),
    sirstochastic::individual_R_to_S(S, R, human, immunity, age, location, pars),
    sirstochastic::render_state_sizes(S, I, R, human)
  )

  output <- simulate(human, processes, timestep, parameters = list(immunity_level = .2, age_level=0.3))

  df <-   data.frame(S = output$susceptable_counts, I = output$infected_counts, R = output$recovered_counts, time = output$time, type = "Individual",  legend = "Individual", stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))

})

test_that("test individual model with 10000 humans with age effects", {
  library(individual)
  library(reshape2)

  pars <- get_parameters()
  pars[["includeage"]] <- TRUE

  population <- pars$N
  NI <- pars$I0
  NR <- 2
  pops <- population - NI - NR
  timestep <- pars$num/pars$dt

  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', NR)

  immunity <- Variable$new('immunity', rep(0, pars$N))
  rate=1/pars$average_age
  age <- Variable$new('age', rexp(pars$N, rate))
  location <- Variable$new('location', rep(0, pars$N))
  human <- Individual$new('human', list(S, I, R), variables = list(immunity, age, location))

  processes <- list(
    sirstochastic::individual_S_to_I(S, I, human, immunity, age, location, pars),
    sirstochastic::individual_I_to_R(I, R, human, immunity, age, location, pars),
    sirstochastic::individual_R_to_S(S, R, human, immunity, age, location, pars),
    sirstochastic::render_state_sizes(S, I, R, human)
  )

  output <- simulate(human, processes, timestep, parameters = list(age_level=0.3))

  df <-   data.frame(S = output$susceptable_counts, I = output$infected_counts, R = output$recovered_counts, time = output$time, type = "Individual",  legend = "Individual", stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))

})

test_that("test individual model with 10000 humans with effects due to location", {
  library(individual)
  library(reshape2)

  pars <- get_parameters()
  pars[["includelocation"]] <- TRUE

  population <- pars$N
  NI <- pars$I0
  NR <- 2
  pops <- population - NI - NR
  timestep <- pars$num/pars$dt

  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', NR)

  immunity <- Variable$new('immunity', rep(0, pars$N))
  age <- Variable$new('age', rep(0, pars$N))
  location <- Variable$new('location', runif(population, 0, .2))
  human <- Individual$new('human', list(S, I, R), variables = list(immunity, age, location))

  processes <- list(
    sirstochastic::individual_S_to_I(S, I, human, immunity, age, location, pars),
    sirstochastic::individual_I_to_R(I, R, human, immunity, age, location, pars),
    sirstochastic::individual_R_to_S(S, R, human, immunity, age, location, pars),
    sirstochastic::render_state_sizes(S, I, R, human)
  )

  output <- simulate(human, processes, timestep, parameters = list(location_level=0.2))

  df <-   data.frame(S = output$susceptable_counts, I = output$infected_counts, R = output$recovered_counts, time = output$time, type = "Individual",  legend = "Individual", stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))

})

test_that("test individual model with 10000 humans with immunity, age and location effects", {
  library(individual)
  library(reshape2)

  pars <- get_parameters()
  pars[["indludeimmune"]] <- TRUE
  pars[["includeage"]] <- TRUE
  pars[["includelocation"]] <- TRUE

  population <- pars$N
  NI <- pars$I0
  NR <- 2
  pops <- population - NI - NR
  timestep <- pars$num/pars$dt

  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', NR)

  immunity <- Variable$new('immunity', runif(population, 0, .1))
  rate=1/pars$average_age
  age <- Variable$new('age', rexp(pars$N, rate))
  location <- Variable$new('location', runif(population, 0, .2))
  human <- Individual$new('human', list(S, I, R), variables = list(immunity, age, location))

  processes <- list(
    sirstochastic::individual_S_to_I(S, I, human, immunity, age, location, pars),
    sirstochastic::individual_I_to_R(I, R, human, immunity, age, location, pars),
    sirstochastic::individual_R_to_S(S, R, human, immunity, age, location, pars),
    sirstochastic::render_state_sizes(S, I, R, human)
  )

  output <- simulate(human, processes, timestep, parameters = list(immunity_level = .2, age_level=0.3, location_level = 0.4))

  df <-   data.frame(S = output$susceptable_counts, I = output$infected_counts, R = output$recovered_counts, time = output$time, type = "Individual",  legend = "Individual", stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))

})
