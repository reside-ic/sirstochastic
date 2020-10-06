get_parameters <- function(overrides = list()) {

  parameters <- list(
    beta = 0.5,
    nu = 0.3,
    mu = 0.001,
    prop_immune = 0,
    I0_at_steady_state = 0,
    N = 10000,
    I0 = 5,
    dt = 0.01)

  parameters

}
