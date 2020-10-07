#' Default parameters for SIR model
#'
#' @param overrides
#' @param beta contact rate
#' @param nu rate of recovery
#' @param mu rate of death
#' @param prop_immune proportion of people who are immune
#' @param I0_at_steady_state initial number of infected at steady state
#' @param N number of people being investigated
#' @param I0 initial number of infected people
#' @param dt time step
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
    I0 = 5,
    dt = 0.01,
    n_events_S = NULL,
    n_deaths_S = NULL,
    n_infections_S = NULL,
    n_events_S = NULL,
    n_deaths_S = NULL,
    n_recoveries_S = NULL,
    n_deaths_R = NULL,
    n_deaths_R = NULL)

  parameters

}
