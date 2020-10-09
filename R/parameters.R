#' Default parameters for SIR model
#'
#' @param overrides
#' @param beta contact rate
#' @param nu rate of recovery
#' @param mu rate of death
#' @param prop_immune proportion of people who are immune
#' @param N number of people being investigated
#' @param num used for countlim/num time points
#' @param I0 initial number of infected people
#' @param dt time step
#' @param I0_at_steady_state boolean value
#'
#' @return
#' @export
#'
#' @examples
#'
#'


get_parameters <- function(overrides = list()) {

  parameters <- list(
    beta = 0.5,
    nu = 0.3,
    mu = 0.001,
    prop_immune = 0,
    I0_at_steady_state = FALSE,
    N = 10000,
    num = 100,
    I0 = 5,
    dt = 0.01)


  if(parameters$beta < 0){
    stop("'beta' must be positive")
  }

  if(parameters$nu < 0){
    stop("'nu' must be positive")
  }

  if(parameters$mu < 0){
    stop("'mu' must be positive")
  }

  if(parameters$prop_immune > 0 || parameters$prop_immune < 0){
    stop("'prop_immune' must be between 0 and 1 (inclusive)")
  }

  if(parameters$N == 0 || parameters$N < 0){
    stop("'N' must be positive")
  }

  if(parameters$I0 > parameters$N || parameters$I0 < 0){
    stop("'I0' must be positive and never greater than N")
  }

  if(parameters$dt <= 0){
    stop("'dt' must be positive and greater than 0")
  }

  parameters

}
