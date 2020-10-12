
#' Stochastic SIR model
#'
#' This takes N iterations for time step dt
#'
#' @return
#' @export
#'
#' @examples
#' sirmodel(pars)
#' initialisesir(pars)
#' infections(pars, IJ, SJ)
#' recoveries(pars, IJ)
#' births(pars, RJ, n_deaths_S, n_deaths_I)
#' update(pars, inf, rec, bir, SJ, IJ, RJ, j)
#' displaythemodel(df2)
sirmodel <- function(end_time, pars) {
  warnings()

  if(length(pars) == 0){
    pars <- get_parameters()
  }

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

  data.frame(S = S, I = I, R = R, time = time, plotsave = pars$plotsave, stringsAsFactors = FALSE)
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

  # initialise S, I and R
  S <- if (pars$I0_at_steady_state) round(S_star) else S0
  I <- if (pars$I0_at_steady_state) round(I_star) else pars$I0
  R <- if (pars$I0_at_steady_state) pars$N - round(I_star) - round(S_star) else pars$N - pars$I0 - S0

  list(S = S, I = I, R = R)
}

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

  ggplot2::ggplot(df, ggplot2::aes(x=time, y=value, group=interaction(group, name), colour=name ) ) +
    ggplot2::geom_line(size=1) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "SIR against time", subtitle = subtitle, color="Category") +
    ggplot2::labs(y ="S, I, & R") +
    ggplot2::theme(
      legend.justification = c("right", "top"),
      legend.box = c("horizontal", "vertical")
    ) +
    ggplot2::scale_colour_manual(values=c("blue", "red", "purple")) +
    ggplot2::theme( plot.title = ggplot2::element_text(color="black", size=14,
      face="bold.italic"), axis.title.x = ggplot2::element_text(color="blue", size=14, face="bold"),
      axis.title.y = ggplot2::element_text(color="blue", size=14, face="bold.italic"))

  if(df$plotsave[[1]]){
    savedfilename <- paste(numruns, 'run(s)', datapoints, 'pts.png')
    ggplot2::ggsave(
      savedfilename, h = 9/3, w = 16/3, type = "cairo-png")
  }

}



