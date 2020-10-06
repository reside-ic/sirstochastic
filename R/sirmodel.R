
#' Stochastic SIR model
#'
#' This takes 10000 iterations for time step dt = 0.01
#'
#' @return
#' @export
#'
#' @examples
#' sirmodel(pars)
#' initialisesir()
#' infections(IJ, SJ)
#' recoveries(IJ)
#' births(RJ, n_deaths_S, n_deaths_I)
#' displaythemodel(results)
#  sirmodel <- function(pars) {
sirmodel <- function(pars) {
  warnings()
  NT <- pars$N + 1
  S <- rep(NA_integer_, NT)
  I <- rep(NA_integer_, NT)
  R <- rep(NA_integer_, NT)

  ## Initialise S, I, R
  initsir <- initialisesir()

  S[1] <- initsir$S
  I[1] <- initsir$I
  R[1] <- initsir$R

  count <- seq(1, 10000, by=1)
  for (j in count)
  {
    inf <- infections(I[j], S[j])
    rec <- recoveries(I[j])
    bir <- births(R[j], inf$n_deaths_S, rec$n_deaths_I)

    # update for next time step
    if (j == NT){
      break
    }
    else{
      S[j+1] <- S[j] - inf$n_deaths_S - inf$n_infections_S + bir$n_births
      I[j+1] <- I[j] + inf$n_infections_S - rec$n_recoveries_I - rec$n_deaths_I
      R[j+1] <- R[j] + rec$n_recoveries_I - bir$n_deaths_R
      newN <- S[j+1] + I[j+1] + R[j+1]
      if(newN != pars$N){
        thetime <- pars$dt * (j + 1)
        stop(sprintf("newN is not equal to N for time = %f",
                     thetime))
      }
    }
  }

  time <- seq(0, 100, by = pars$dt)
  list(S = S, I = I, R = R, time = time)
}

initialisesir <- function(){
  pars <- get_parameters()
  ## Stochastic solution
  R0 <- pars$beta/( pars$nu + pars$mu )
  # SIR: new definition of number of infected people at endemic equilibrium state
  I_star <- pars$N*pars$mu*(pars$beta - pars$nu - pars$mu)/(pars$beta*(pars$mu + pars$nu))
  # SIR: new definition of number of susceptible people at endemic equilibrium state
  S_star <- pars$N/R0

  S0 <- (pars$N-pars$I0)*(1 - pars$prop_immune)  # initial susceptible people
  R0 <- pars$N - pars$I0 - S0

  S <- if (pars$I0_at_steady_state > 0) round(S_star) else S0
  I <- if (pars$I0_at_steady_state > 0) round(I_star) else pars$I0
  R <- if (pars$I0_at_steady_state > 0) pars$N - round(I_star) - round(S_star) else pars$N - pars$I0 - S0

  list(S = S, I = I, R = R)
}

infections <- function(IJ, SJ){

  pars <- get_parameters()
  FOI <- pars$beta * IJ / pars$N
  fmu <- FOI + pars$mu
  prob1 <- fmu *pars$dt
  prob2 <- (pars$mu/fmu)
  n_events_S <- rbinom(1, SJ, prob1) # SIR: two types of events for S, so competing hazards.
  n_deaths_S <- rbinom(1, n_events_S, prob2) # SIR: a fraction of S events are deaths.
  n_infections_S <- n_events_S - n_deaths_S	# SIR: ...the rest are infections.

  list(n_deaths_S = n_deaths_S, n_infections_S = n_infections_S)
}

recoveries <- function(IJ){
  pars <- get_parameters()
  prob3 <- (pars$nu + pars$mu)* pars$dt
  prob4 <- pars$mu/(pars$mu + pars$nu)
  n_events_I <- rbinom(1, IJ, prob3) # SIR: two types of events for I, so competing hazards.
  n_deaths_I <- rbinom(1, n_events_I, prob4) # SIR: a fraction of I events are deaths.
  n_recoveries_I <- n_events_I - n_deaths_I # SIR: ...the rest are recoveries.

  list(n_deaths_I = n_deaths_I, n_recoveries_I = n_recoveries_I)
}

births <- function(RJ, n_deaths_S, n_deaths_I){
  pars <- get_parameters()
  prob5 <- pars$mu*pars$dt
  n_deaths_R <- rbinom(1, RJ, prob5)
  n_births <- n_deaths_S + n_deaths_I + n_deaths_R

  list(n_deaths_R = n_deaths_R, n_births = n_births)
}

displaythemodel <- function(results) {

  sir_col <- c("#8c8cd9", "#cc0044", "#999966")
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  check = matplot(results$time[1:10000], cbind(results$S[1:10000], results$I[1:10000], results$R[1:10000]), type = "l", col = sir_col, lty = 1, xlab="time", ylab="S, I, R")
  legend("topright", inset = 0.05, lwd = 1, col = sir_col, legend = c("S", "I", "R"))
  title("SIR model, single run, N = 10,000, dt = 0.01", line = -2, adj = 0.2, col.main = "gray", font.main = 4)

  list(check)
}

