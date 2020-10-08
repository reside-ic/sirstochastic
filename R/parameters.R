#' Default parameters for SIR model
#'
#' @param overrides
#' @param beta contact rate
#' @param nu rate of recovery
#' @param mu rate of death
#' @param prop_immune proportion of people who are immune
#' @param I0_at_steady_state initial number of infected at steady state
#' @param N number of people being investigated
#' @param num used for N/num time points
#' @param I0 initial number of infected people
#' @param dt time step
#' @param problim limit for exponent in probabilities
#' @param countlim limit for number of points to be plotted
#'
#' @return
#' @export
#'
#' @examples
get_parameters <- function(overrides = list()) {

  parameters <- list(
    beta = 0.5,
    nu = 0.3,
    mu = 0.001,
    prop_immune = 0,
    I0_at_steady_state = 0,
    N = 10000,
    num = 100,
    I0 = 5,
    dt = 0.01,
    problim = 1000,
    countlim = 10000)

  parameters

}
