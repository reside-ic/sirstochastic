#' Stochastic SIR model, compartmental
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
#' individual_sirmodel(100,list())
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
    inf <- infections(pars, I[j], S[j])
    rec <- recoveries(pars, I[j])
    bir <- births(pars, R[j], inf, rec)

    # update for next time step
    # Find next values of S, I and R
    updated <- update(pars, inf, rec, bir, S[j], I[j], R[j])

    news <- updated$news
    newi <- updated$newi
    newr <- updated$newr
  }

  data.frame(S = S, I = I, R = R, time = time, stringsAsFactors = FALSE)
}

initialisesir <- function(pars){
  ## Stochastic solution
  R0 <- pars$beta/( pars$nu + pars$mu )
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
    R <-pars$N - round(I_star) - round(S_star)
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
infections <- function(pars, IJ, SJ){
  # SIR: two types of events for S, so competing hazards. A fraction of
  # S events are deaths and the rest are infections.
  FOI <- pars$beta * IJ / pars$N
  fmu <- FOI + pars$mu
  fmudt <- fmu *pars$dt

  prob1 <- 1.0 - exp(-1.0 * fmudt)
  prob2 <- pars$mu/fmu

  n_events_S <- rbinom(1, SJ, prob1)
  if(n_events_S >0){
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

recoveries <- function(pars, IJ){
  # SIR: two types of events for I, so competing hazards and a fraction of
  # I events are deaths and the rest are recoveries
  coeff <- pars$nu + pars$mu
  coeffdt <- coeff * pars$dt

  prob1 <- 1.0 - exp(-1.0 * coeffdt)
  prob2 <- 1.0 - exp(-1.0 * pars$mu/coeff)

  n_events_I <- rbinom(1, IJ, prob1)
  n_deaths_I <- rbinom(1, n_events_I, prob2)
  n_recoveries_I <- n_events_I - n_deaths_I

  list(n_deaths_I = n_deaths_I, n_recoveries_I = n_recoveries_I)
}

births <- function(pars, RJ, inf, rec){

  coeffdt <- pars$mu*pars$dt
  prob <- 1.0 - exp(-1.0 * coeffdt)

  n_deaths_R <- rbinom(1, RJ, prob)
  n_births <- inf$n_deaths_S + rec$n_deaths_I + n_deaths_R

  list(n_deaths_R = n_deaths_R, n_births = n_births)
}

update <- function(pars, inf, rec, bir, SJ, IJ, RJ){

    news <- SJ - inf$n_deaths_S - inf$n_infections_S + bir$n_births
    newi <- IJ + inf$n_infections_S - rec$n_recoveries_I - rec$n_deaths_I
    newr <- RJ + rec$n_recoveries_I - bir$n_deaths_R
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
    datapoints <- paste(nrow(df) - 1)
    numruns <- 1
    df <- list(df)
    subtitle <- paste('Simulation for', numruns, 'run and', datapoints, 'data points')
  }
  else{
    numruns <- length(df)
    datapoints <- length(df[[1]][[1]])-1
    subtitle <- paste('Simulation for', numruns, 'runs,', datapoints, 'data points per run')
  }

  # Create group id for data
  df <- dplyr::bind_rows(df, .id = "group")

  # Convert to long format
  df <- tidyr::pivot_longer(tibble::as_tibble(df), c("S", "I", "R"))

  ggplot2::ggplot(df, ggplot2::aes(x=df$time, y=df$value, group=interaction(df$group, df$name), colour=df$name ) ) +
    ggplot2::geom_line(size=0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "SIR model simulation", subtitle = subtitle, color="Compartment") +
    ggplot2::labs(y ="S, I, & R", x="time") +
    ggplot2::theme(
      legend.justification = c("right", "top"),
      legend.box = c("horizontal", "vertical")
    ) +
    ggplot2::scale_colour_manual(values=c("blue", "red", "green")) +
    ggplot2::theme(text = ggplot2::element_text(color = "#444444", family = 'Helvetica Neue'),
    plot.title = ggplot2::element_text(size = 26, color = '#333333'),
    plot.subtitle = ggplot2::element_text(size = 13),
    axis.title.x = ggplot2::element_text(size = 16, color = '#333333'),
    axis.title.y = ggplot2::element_text(angle = 0, vjust = .5))
}

plot_states <- function(output) {
  ggplot(
    melt(output, 'timestep'),
    aes(x = timestep, y = value, group = variable)) +
    geom_line(aes(color = variable), size=0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "SIR model simulation (individual)", color="Individual") +
    ggplot2::labs(y ="S, I, & R", x="time") +
    ggplot2::scale_colour_manual(values=c("blue", "red", "green")) +
    ggplot2::theme(text = ggplot2::element_text(color = "#444444", family = 'Helvetica Neue'),
                   plot.title = ggplot2::element_text(size = 26, color = '#333333'),
                   plot.subtitle = ggplot2::element_text(size = 13),
                   axis.title.x = ggplot2::element_text(size = 16, color = '#333333'),
                   axis.title.y = ggplot2::element_text(angle = 0, vjust = .5))
}

#' Stochastic SIR model, individual
#'
#' Takes timestep steps over time for population N
#'
#' @param end_time end time for data
#' @param pars parameter list
#'
#' @return dataframe
#'
#' @export
#'
#' @examples
#' individual_sirmodel(100, list())
individual_sirmodel <- function(end_time, pars = NULL) {

  pars <- get_parameters(pars)
  population <- pars$N

  process <- function(api) {
    # calculate information for infections, recoveries and births
    inf <- infections(pars, length(api$get_state(human, I)), length(api$get_state(human, S)))
    rec <- recoveries(pars, length(api$get_state(human, I)))
    bir <- births(pars, length(api$get_state(human, R)), inf, rec)
    #infection_process -> S to I
    n_to_infect <- inf$n_infections_S
    susceptible <- api$get_state(human, S)
    infected <- susceptible[sample.int(length(susceptible), n_to_infect)]
    api$queue_state_update(human, I, infected)
    #recovery_process -> I to R
    n_to_recover <- rec$n_recoveries_I
    infected <- api$get_state(human, I)
    recovered <- infected[sample.int(length(infected), n_to_recover)]
    api$queue_state_update(human, R, recovered)
    #Recovered to susceptible -> R to S
    # n_to_susceptible <- n_to_recover - bir$n_deaths_R
    # from_state <- api$get_state(human, R)
    #
    # if(length(from_state) != 0 && length(from_state) > n_to_susceptible){
    #   thenewsusceptible <- from_state[sample.int(length(from_state), n_to_susceptible)]
    #   api$queue_state_update(human, S, thenewsusceptible)
    # }
  }

  render_state_sizes <- function(api) {
    api$render('susceptable_counts', length(api$get_state(human, S)))
    api$render('infected_counts', length(api$get_state(human, I)))
    api$render('recovered_counts', length(api$get_state(human, R)))
  }

  # Set no of infections, NI, no of susceptibles, pops,
  # no of recovered at the start and no of time points
  NI <- pars$I0
  pops <- population - NI
  timestep <- pars$num/pars$dt
  S <- State$new('S', pops)
  I <- State$new('I', NI)
  R <- State$new('R', 0)

  human <- Individual$new('human', list(S, I, R))

  output <- simulate(human,
                     list(process,
                     render_state_sizes),
                     timestep)

  plot_states <- function(output) {
    ggplot(
      melt(output, 'timestep'),
      aes(x = timestep, y = value, group = variable)
    ) + geom_line(aes(color = variable)) +
      ggplot2::labs(title = "SIR model simulation (individual)", color="Individual") +
      ggplot2::labs(y ="S, I, & R", x="time")
  }
  plot_states(output)

  output

}


















