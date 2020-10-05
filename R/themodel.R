
#' Stochastic SIR model
#'
#' This takes 10000 iterations for time step dt = 0.01
#'
#' @return
#' @export
#'
#' @examples
#' themodel(dt)
#' displaythemodel(dt)
themodel <- function(dt) {
  warnings()
  beta <- 0.5 # contact rate
  nu <- 0.3 # recovery
  mu <- 0.001 # death rate
  prop_immune <- 0 # proportion of population initially immune
  I0_at_steady_state <- 0
  N <- 10000  # total population.
  NT <- 10001
  xS <- rep(NA_integer_, NT)
  yI <- rep(NA_integer_, NT)
  zR <- rep(NA_integer_, NT)

  ## Stochastic solution
  R0 <- beta/( nu + mu )
  I_star <- N*mu*(beta - nu - mu)/(beta*(mu + nu)) # SIR: new definition of number of infecteds at endemic equilibrium state
  S_star <- N/R0 # SIR: new definition of number of susceptibles at endemic equilibrium state
  I0 <- 5 # initial infecteds
  S0 <- (N-I0)*(1 - prop_immune) # initial susceptibles
  R0 <- N - I0 - S0

  xS[1] <- if (I0_at_steady_state > 0) round(S_star) else S0
  yI[1] <- if (I0_at_steady_state > 0) round(I_star) else I0
  zR[1] <- if (I0_at_steady_state > 0) N - round(I_star) - round(S_star) else N - I0 - S0

  ##dt <- 0.01

  count <- seq(1, 10000, by=1)
  for (j in count)
  {
    FOI <- beta * yI[j] / N
    fmu <- FOI + mu
    prob1 <- 1.0 - exp(-1.0* fmu *dt)
    prob2 <- 1.0 - exp(-1.0*(mu/fmu))
    prob11 <- fmu *dt
    prob21 <- (mu/fmu)
    n_events_S <- rbinom(1, xS[j], prob11) # SIR: two types of events for S, so competing hazards.
    n_deaths_S <- rbinom(1, n_events_S, prob21) # SIR: a fraction of S events are deaths.
    n_infections_S <- n_events_S - n_deaths_S	# SIR: ...the rest are infections.

    prob3 <- 1.0 - exp(-1.0*(nu + mu)* dt)
    prob4 <- 1.0 - exp(-1.0*(mu/(mu + nu)))
    prob31 <- (nu + mu)* dt
    prob41 <- mu/(mu + nu)
    n_events_I <- rbinom(1, yI[j], prob31) # SIR: two types of events for I, so competing hazards.
    n_deaths_I <- rbinom(1, n_events_I, prob41) # SIR: a fraction of I events are deaths.
    n_recoveries_I <- n_events_I - n_deaths_I # SIR: ...the rest are recoveries.
    prob5 <- mu*dt
    n_deaths_R <- rbinom(1, zR[j], prob5)
    n_births <- n_deaths_S + n_deaths_I + n_deaths_R

    # update for next time step
    if (j == NT){
      break
    }
    else{
      xS[j+1] <- xS[j] - n_deaths_S - n_infections_S + n_births
      yI[j+1] <- yI[j] + n_infections_S - n_recoveries_I - n_deaths_I
      zR[j+1] <- zR[j] + n_recoveries_I - n_deaths_R
      newN <- xS[j+1] + yI[j+1] + zR[j+1]
      if(newN != N){
        stop(sprintf("newN is not equal to N for time = %f", thetime))
      }
    }
  }
  time <- seq(0, 100, by = dt)
  list(S = xS, I = yI, R = zR, time = time)
}

displaythemodel <- function(results) {
  warnings()

  sir_col <- c("#8c8cd9", "#cc0044", "#999966")
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  matplot(results$time[1:10000], cbind(results$S[1:10000], results$I[1:10000], results$R[1:10000]), type = "l", col = sir_col, lty = 1, xlab="time", pch=c(1,16,17))
  legend(1, lwd = 1, col = sir_col, legend = c("S", "I", "R"))

}

