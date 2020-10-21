#' @title Default parameters for SIR model
#'
#' @param overrides use a named parameter list instead of defaults
#' Parameters defined below
#'
#' * `pars` list of parameters
#' Compartmental
#' * `beta` contact rate
#' * `nu rate` of recovery
#' * `mu rate` of death
#' * `prop_immune` proportion of people who are immune
#' * `N number` of people being investigated
#' * `num` used for countlim/num time points
#' * `I0` initial number of infected people
#' * `dt` time step
#' * `I0_at_steady_state` boolean value
#' * `n_deaths_S` number of deaths at susceptible stage
#' * `n_infections_S` number of infections at susceptible stage
#' * `n_deaths_I` number of deaths at infected stage
#' * `n_recoveries_I` number of recoveries at infected stage
#' * `n_deaths_R` number of deaths at recovered stage
#' * `n_births` number of births
#' Individual only
#' * `average_age` average age for population
#' * `includeage` TRUE if age used
#' * `includebirth` TRUE if immunity used
#' * `indludeimmune` TRUE if immunity used
#' * `includelocation` TRUE if location used
#'
#' @return list
#' @export
get_parameters <- function(overrides = list()) {

  pars <- sir_model_parameters_defaults()

  # Override pars with any client specified ones
  if (!is.list(overrides) && !is.null(overrides)) {
    stop('overrides must be a list')
  }

  for (name in names(overrides)) {
    if (!(name %in% names(pars))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    pars[[name]] <- overrides[[name]]
  }

  if(pars$beta < 0){
    stop("'beta' must be positive")
  }

  if(pars$nu < 0){
    stop("'nu' must be positive")
  }

  if(pars$mu < 0){
    stop("'mu' must be positive")
  }

  if(pars$prop_immune > 0 || pars$prop_immune < 0){
    stop("'prop_immune' must be between 0 and 1 (inclusive)")
  }

  if(pars$N <= 0){
    stop("'N' must be positive")
  }

  if(pars$num <= 0){
    stop("'num' must be positive")
  }

  if(pars$I0 > pars$N || pars$I0 < 0){
    stop("'I0' must be positive and never greater than N")
  }

  if(pars$dt <= 0){
    stop("'dt' must be positive and greater than 0")
  }

  if(pars$n_deaths_S < 0){
    stop("'n_deaths_S' must be positive and greater than or equal to 0")
  }

  if(pars$n_infections_S < 0){
    stop("'n_infections_S' must be positive and greater than or equal to 0")
  }

  if(pars$n_deaths_I < 0){
    stop("'n_deaths_I' must be positive and greater than or equal to 0")
  }

  if(pars$n_recoveries_I < 0){
    stop("'n_recoveries_I' must be positive and greater than or equal to 0")
  }

  if(pars$n_deaths_R < 0){
    stop("'n_deaths_R' must be positive and greater than or equal to 0")
  }

  if(pars$n_births < 0){
    stop("'n_births' must be positive and greater than or equal to 0")
  }

  if(pars$average_age <= 0){
    stop("'average_age' must be positive and greater than 0")
  }

  pars

}

sir_model_parameters_defaults <- function() {

  pars <- list(
    # Compartmental
    beta = 0.5,
    nu = 0.3,
    mu = 0.001,
    prop_immune = 0,
    I0_at_steady_state = FALSE,
    N = 10000,
    num = 100,
    I0 = 5,
    dt = 0.01,
    n_deaths_S = 0,
    n_infections_S = 0,
    n_deaths_I = 0,
    n_recoveries_I = 0,
    n_deaths_R = 0,
    n_births = 0,
    # individual only
    average_age = 30,
    includeage = FALSE,
    includebirth = FALSE,
    indludeimmune = FALSE,
    includelocation = FALSE
    )

  pars
}

#' Events list
#'
#' @param overrides use a named events list instead of defaults
#'
#' @return events
#' @export
get_events <- function(overrides = list()) {

  events <- sir_individual_model_events_defaults()

  # Override pars with any client specified ones
  if (!is.list(overrides) && !is.null(overrides)) {
    stop('overrides must be a list')
  }

  for (name in names(overrides)) {
    if (!(name %in% names(events))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    events[[name]] <- overrides[[name]]
  }
}

sir_individual_model_events_defaults <- function() {

  events <- list(
    # individual only
    name = '',
    listener = NULL)

  events
}
