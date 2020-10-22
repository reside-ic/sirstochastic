#' @title Stochastic SIR model, compartmental
#'
#' This takes N iterations for time step dt
#'
#' @param end_time en time of data
#' @param pars parameter list
#'
#' @return dataframe
#' @export
#'
#' @examples
#' compartmental_sirmodel(100, list())
#' run_with_repetitions(100, 1, list(), FALSE)
compartmental_sirmodel <- function(end_time, pars = NULL) {

  pars <- get_parameters(pars)

  time <- seq(0, end_time, by = pars$dt)
  NT <- length(time)
  S <- rep(NA_integer_, NT)
  I <- rep(NA_integer_, NT)
  R <- rep(NA_integer_, NT)

  ## Initialise S, I, R
  initsir <- initialisesir(pars)

  S[1] <- initsir$S
  I[1] <- initsir$I
  R[1] <- initsir$R

  snew <- 0
  inew <- 0
  rnew <- 0

  for (j in seq_along(time))
  {
    if( j > 1){
      S[j] <- news
      I[j] <- newi
      R[j] <- newr
    }

    # calculate information for infections, recoveries and births
    inf <- infections(I[j], S[j], pars)
    pars$n_deaths_S <- inf$n_deaths_S
    pars$n_infections_S <- inf$n_infections_S

    rec <- recoveries(I[j], pars)
    pars$n_deaths_I <- rec$n_deaths_I
    pars$n_recoveries_I <- rec$n_recoveries_I

    bir <- births(R[j], pars)
    pars$n_deaths_R <- bir$n_deaths_R
    pars$n_births <- bir$n_births

    # update for next time step
    # Find next values of S, I and R
    updated <- update(S[j], I[j], R[j], pars)

    news <- updated$news
    newi <- updated$newi
    newr <- updated$newr
  }

  data.frame(S = S, I = I, R = R, time = time, type = "Compartmental", legend = "Compartment", stringsAsFactors = FALSE)
}

initialisesir <- function(pars){
  ## Stochastic solution
  R0 <- pars$beta/pars$nu + pars$mu
  # SIR: new definition of number of infected people at endemic equilibrium state
  I_star <- pars$N*pars$mu*(pars$beta - pars$nu - pars$mu)/(pars$beta*(pars$mu + pars$nu))
  # SIR: new definition of number of susceptible people at endemic equilibrium state
  S_star <- pars$N/R0

  # initial susceptible people
  S0 <- (pars$N-pars$I0)*(1 - pars$prop_immune)
  R0 <- pars$N - pars$I0 - S0

  S <- 0
  I <- 0
  R <- 0

  # initialise S, I and R
  if (pars$I0_at_steady_state){
    S <- round(S_star)
    I <- round(I_star)
    R <- pars$N - round(I_star) - round(S_star)
  }
  else
  {
    S <- S0
    I <- pars$I0
    R <- R0
  }

  list(S = S, I = I, R = R)
}

#' @importFrom stats rbinom
infections <- function(I, S, pars){

  # SIR: two types of events for S, so competing hazards. A fraction of
  # S events are deaths and the rest are infections.

  FOI <- (pars$beta * I)/pars$N
  fmu <- FOI + pars$mu
  fmudt <- fmu * pars$dt

  prob1 <- 1.0 - exp(-1.0 * fmudt)
  prob2 <- pars$mu/fmu

  n_events_S <- rbinom(1, S, prob1)

  if(n_events_S > 0){
    n_deaths_S <- rbinom(1, n_events_S, prob2)
    n_infections_S <- n_events_S - n_deaths_S
  }
  else
  {
    n_deaths_S <- 0
    n_infections_S <- 0
  }

  list(n_deaths_S = n_deaths_S, n_infections_S = n_infections_S)
}

recoveries <- function(I, pars){
  # SIR: two types of events for I, so competing hazards and a fraction of
  # I events are deaths and the rest are recoveries
  coeff <- pars$nu + pars$mu
  coeffdt <- coeff * pars$dt

  prob1 <- 1.0 - exp(-1.0 * coeffdt)
  prob2 <- 1.0 - exp(-1.0 * pars$mu/coeff)

  n_events_I <- rbinom(1, I, prob1)

  n_deaths_I <- rbinom(1, n_events_I, prob2)
  n_recoveries_I <- n_events_I - n_deaths_I

  list(n_deaths_I = n_deaths_I, n_recoveries_I = n_recoveries_I)
}

#' @importFrom stats rbinom
births <- function(R, pars){

  coeffdt <- pars$mu * pars$dt
  prob <- 1.0 - exp(-1.0 * coeffdt)
  n_deaths_R <- rbinom(1, R, prob)
  n_births <- pars$n_deaths_S  + pars$n_deaths_I + n_deaths_R

  list(n_deaths_R = n_deaths_R, n_births = n_births)
}

update <- function(S, I, R, pars){

    news <- S - pars$n_deaths_S - pars$n_infections_S + pars$n_births
    newi <- I + pars$n_infections_S  - pars$n_recoveries_I  - pars$n_deaths_I
    newr <- R + pars$n_recoveries_I - pars$n_deaths_R
    newN <- news + newi + newr

    list(news = news, newi = newi, newr = newr)
}

#' @title Run the simulation with repetitions
#'
#' @param end_time end time for run
#' @param repetitions n times to run the simulation
#' @param pars parameter list
#' @param parallel execute runs in parallel, TRUE or FALSE
#' @return dataframe
#' @export
run_with_repetitions <- function(
  end_time,
  repetitions,
  pars,
  parallel = FALSE
) {
  if (parallel) {
    fapply <- parallel::mclapply
  } else {
    fapply <- lapply
  }
  dfs <- fapply(
    seq(repetitions),
    function(repetition) {
      df <- compartmental_sirmodel(end_time, pars)
      df$repetition <- repetition
      df
    }
  )

  do.call("rbind", dfs)
}

displaythemodel <- function(df) {

  # This function displays data in a list. df must be in the form of a list.

  # Check if df is a dataframe. If yes then turn it into a list
  numruns <- 0
  datapoints <- 0
  subtitle <- ""

  if(is.data.frame(df)){
    numdatapoints <- paste(length(df$time))
    numruns <- 1
    df <- list(df)
    subtitle <- paste('Simulation for', numruns, 'run and', numdatapoints, 'data points')
  }
  else{
    numruns <- length(df)
    numdatapoints <- length(df[[1]][[1]])-1
    subtitle <- paste('Simulation for', numruns, 'runs,', numdatapoints, 'data points per run')
  }

  # Create group id for data
  df <- dplyr::bind_rows(df, .id = "group")

  # Convert to long format
  df <- tidyr::pivot_longer(tibble::as_tibble(df), c("S", "I", "R"))

  strname <- paste("SIR", df$type, "Model Simulation")

  ggplot2::ggplot(df, ggplot2::aes(x=df$time, y=df$value, group=interaction(df$group, df$name), colour=df$name ) ) +
    ggplot2::geom_line(size=0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = strname, subtitle = subtitle, color=df$legend) +
    ggplot2::labs(y ="S, I, & R", x="time") +
    ggplot2::theme(
      legend.justification = c("right", "top"),
      legend.box = c("horizontal", "vertical")
    ) +
    ggplot2::scale_colour_manual(values=c("blue", "red", "green")) +
    ggplot2::theme(text = ggplot2::element_text(color = "#444444", family = 'Lucida Bright'),
    plot.title = ggplot2::element_text(size = 26, color = '#333333'),
    plot.subtitle = ggplot2::element_text(size = 13),
    axis.title.x = ggplot2::element_text(size = 16, color = '#333333'),
    axis.title.y = ggplot2::element_text(angle = 0, vjust = .5))
}

#' @title Infection_process -> S to I
#'
#' @param S susceptible
#' @param I infected
#' @param human human
#' @param immunity immunity
#' @param age age
#' @param location location
#' @param pars parameter list
#'
#' @export
#'
#' @examples
#' individual_S_to_I(S, I, immunity, age, location, pars)
#' @importFrom stats runif
individual_S_to_I <- function(S, I, human, immunity, age, location, pars = NULL) {
  function(api) {

    pars <- get_parameters(pars)

    # calculate information for infections, recoveries and births
    inf <- infections(length(api$get_state(human, I)), length(api$get_state(human, S)), pars)

    n_to_infect <- inf$n_infections_S
    susceptible <- api$get_state(human, S)

    if(!pars$variations){
        infected <- susceptible[sample.int(length(susceptible), n_to_infect)]
        api$queue_state_update(human, I, infected)
    }
    if(pars$indludeimmune){
      # Get the immunity for susceptible humans and use the complement to modify the
      # infection rate
      rate_modifier <- 1 - api$get_variable(human, immunity, susceptible)
      infected <- susceptible[runif(length(susceptible)) < (pars$infection_rate * rate_modifier)]
      api$queue_state_update(human, I, infected)
    }
    if(pars$includeage){
      # Get the age for susceptible humans and use the complement to modify the
      # infection rate
      rate_modifier <- 1 - api$get_variable(human, age, susceptible)
      infected <- susceptible[runif(length(susceptible)) < (pars$location_rate * rate_modifier)]
      api$queue_state_update(human, I, infected)
    }
    if(pars$includelocation){
      # Get the location for susceptible humans and use the complement to modify the
      # infection rate
      rate_modifier <- 1 - api$get_variable(human, location, susceptible)
      infected <- susceptible[runif(length(susceptible)) < (pars$location_rate * rate_modifier)]
      api$queue_state_update(human, I, infected)
    }
  }
}

#' @title Recovery_process -> I to R
#'
#' @param I infected
#' @param R recovered
#' @param human human
#' @param immunity immunity
#' @param age age
#' @param location location
#' @param pars parameter list
#'
#' @export
#'
#' @examples
#' individual_I_to_R(I, R, human, immunity, age, location, pars)
individual_I_to_R <- function(I, R, human, immunity, age, location, pars = NULL) {
  function(api) {

    pars <- get_parameters(pars)
    rec <- recoveries(length(api$get_state(human, I)), pars)
    n_to_recover <- rec$n_recoveries_I
    infected <- api$get_state(human, I)

    if(!pars$variations){
      recovered <- infected[sample.int(length(infected), n_to_recover)]
      api$queue_state_update(human, R, recovered)
    }
    if(pars$includeage){
      rate_modifier <- 1 - api$get_variable(human, age, infected)
      recovered <- infected[runif(length(infected)) < (pars$age_rate * rate_modifier)]
      api$queue_state_update(human, I, recovered)
    }
    if(pars$includelocation){
      rate_modifier <- 1 - api$get_variable(human, location, infected)
      recovered <- infected[runif(length(infected)) < (pars$location_rate * rate_modifier)]
      api$queue_state_update(human, I, recovered)
    }
  }
}

#' @title Recovered to susceptible -> R to S
#'
#' @param S susceptible
#' @param R recovered
#' @param human human
#' @param immunity immunity
#' @param age age
#' @param location location
#' @param pars parameter list
#'
#' @export
#'
#' @examples
#' individual_R_to_S(S, R, human, immunity, age, location, pars)
#' @importFrom stats runif
individual_R_to_S <- function(S, R, human, immunity, age, location, pars = NULL) {
  function(api) {

    pars <- get_parameters(pars)
    bir <- births(length(api$get_state(human, R)), pars)
    n_to_susceptible <- bir$n_births
    from_state <- api$get_state(human, R)

    if(!pars$variations){
      if(length(from_state) != 0 && length(from_state) > n_to_susceptible)
      {
        thenewsusceptible <- from_state[sample.int(length(from_state), n_to_susceptible)]
        api$queue_state_update(human, S, thenewsusceptible)
      }
    }
    if(pars$indludeimmune){
      recovered <- from_state[runif(length(from_state)) < pars$recovery_rate]
      api$queue_state_update(human, R, recovered)
      api$queue_variable_update(human, immunity, api$get_parameters()$immunity_level, recovered)
    }
    if(pars$includeage){
      recovered <- from_state[runif(length(from_state)) < pars$age_rate]
      api$queue_state_update(human, R, recovered)
      api$queue_variable_update(human, age, api$get_parameters()$age_level, recovered)
    }
    if(pars$includelocation){
      recovered <- from_state[runif(length(from_state)) < pars$recovery_rate]
      api$queue_state_update(human, R, recovered)
      api$queue_variable_update(human, age, api$get_parameters()$location_level, recovered)
    }
  }
}

#' @title Renders the sizes for S, I, R
#' @param S S
#' @param I I
#' @param R R
#' @param human human
#' @export
#' @examples
#' render_state_sizes(S, I, R, human)
#'
render_state_sizes <- function(S, I, R, human) {
  function(api) {
    api$render('susceptable_counts', length(api$get_state(human, S)))
    api$render('infected_counts', length(api$get_state(human, I)))
    api$render('recovered_counts', length(api$get_state(human, R)))
  }
}




















