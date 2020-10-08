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
#' @param n_events_S number of susceptible events
#' @param n_deaths_S number of susceptible deaths
#' @param n_infections_S number of susceptible infections
#' @param n_events_I number of infected events
#' @param n_deaths_I number of infected deaths
#' @param n_recoveries_I number of infected recoveries
#' @param n_births_R number of recovered births
#' @param n_deaths_R number of recovered deaths
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
    n_events_S = NULL,
    n_deaths_S = NULL,
    n_infections_S = NULL,
    n_events_I = NULL,
    n_deaths_I = NULL,
    n_recoveries_I = NULL,
    n_births_R = NULL,
    n_deaths_R = NULL)

  parameters

}
