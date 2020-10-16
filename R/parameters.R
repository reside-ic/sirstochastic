#' Default parameters for SIR model
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
#'
#' Individual - note tis uses
#' * `average_age` average age for population
#' * `includeage` if age is to be included, is TRUE
#' * `includeimmunity` if immunity is to be included, is TRUE
#' * `includebirth` if birth is to be included, is TRUE
#'
#' @return parameter list
#' @export
#'

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

  if(pars$N == 0 || pars$N < 0){
    stop("'N' must be positive")
  }

  if(pars$I0 > pars$N || pars$I0 < 0){
    stop("'I0' must be positive and never greater than N")
  }

  if(pars$dt <= 0){
    stop("'dt' must be positive and greater than 0")
  }

  if(pars$average_age <= 0){
    stop("'age' must be positive and greater than 0")
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
    # individual only
    includeage = FALSE,
    average_age = 30,
    includebirth = FALSE,
    includeimmunity = FALSE)

  pars
}

#' Variable list: immune, age, birth or location
#'
#' @param overrides use a named variable list instead of defaults
#'
#' @return variable list
#' @export
get_variables <- function(overrides = list()) {

  vars <- sir_individual_model_variables_defaults()

  # Override pars with any client specified ones
  if (!is.list(overrides) && !is.null(overrides)) {
    stop('overrides must be a list')
  }

  for (name in names(overrides)) {
    if (!(name %in% names(vars))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    vars[[name]] <- overrides[[name]]
  }
}

sir_individual_model_variables_defaults <- function() {

  vars <- list(
   # individual only
   immune = NULL,
   age = NULL,
   birth = NULL,
   location = NULL)

  vars
}

#' Events list
#'
#' @param overrides use a named events list instead of defaults
#'
#' @return events list
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
