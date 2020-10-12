#' Default parameters for SIR model
#'
#' @param overrides use a named parameter list instead of defaults
#' Parameters defined below
#'
#' * `pars` list of parameters
#' * `beta` contact rate
#' * `nu rate` of recovery
#' * `mu rate` of death
#' * `prop_immune` proportion of people who are immune
#' * `N number` of people being investigated
#' * `num` used for countlim/num time points
#' * `I0` initial number of infected people
#' * `dt` time step
#' * `I0_at_steady_state` boolean value
#' * `plotsave` TRUE if plotted data is to be saved to file
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

  pars

}

sir_model_parameters_defaults <- function() {

  pars <- list(
    beta = 0.5,
    nu = 0.3,
    mu = 0.001,
    prop_immune = 0,
    I0_at_steady_state = FALSE,
    N = 10000,
    num = 100,
    I0 = 5,
    dt = 0.01,
    plotsave = TRUE)

  pars
}

