
#' Stochastic SIR model
#'
#' This takes N iterations for time step dt
#'
#' @return
#' @export
#'
#' @examples
#' sirmodel(pars)
#' initialisesir(pars)
#' infections(pars, IJ, SJ)
#' recoveries(pars, IJ)
#' births(pars, RJ, n_deaths_S, n_deaths_I)
#' update(pars, inf, rec, bir, SJ, IJ, RJ, j)
#' displaythemodel(results)
sirmodel <- function(pars) {
  warnings()

  # Take default parameters if pars is empty
  if (length(pars) == 0 ){
    pars <- get_parameters()
  }

  NT <- pars$N + 1
  S <- rep(NA_integer_, NT)
  I <- rep(NA_integer_, NT)
  R <- rep(NA_integer_, NT)

  ## Initialise S, I, R
  initsir <- initialisesir(pars)

  S[1] <- initsir$S
  I[1] <- initsir$I
  R[1] <- initsir$R

  count <- seq(1, pars$N, by=1)
  for (j in count)
  {
    # calculate information for infections, recoveries and births
    inf <- infections(pars, I[j], S[j])
    rec <- recoveries(pars, I[j])
    bir <- births(pars, R[j], inf, rec)

    # update for next time step
    if (j == NT){
      # break if we have reached the end
      break
    }
    else{
      # Find next values of S, I and R
      updated <- update(pars, inf, rec, bir, S[j], I[j], R[j], j)
    }

    S[j + 1] <- updated$news
    I[j + 1] <- updated$newi
    R[j + 1] <- updated$newr
  }

  tlimit <- pars$N/pars$num
  time <- seq(0, tlimit, by = pars$dt)
  list(S = S, I = I, R = R, time = time, countlim = pars$countlim, dt = pars$dt)
}

initialisesir <- function(pars){

  ## Stochastic solution
  R0 <- pars$beta/( pars$nu + pars$mu )
  # SIR: new definition of number of infected people at endemic equilibrium state
  I_star <- pars$N*pars$mu*(pars$beta - pars$nu - pars$mu)/(pars$beta*(pars$mu + pars$nu))
  # SIR: new definition of number of susceptible people at endemic equilibrium state
  S_star <- pars$N/R0

  # initial susceptible people
  S0 <- (pars$N-pars$I0)*(1 - pars$prop_immune)
  R0 <- pars$N - pars$I0 - S0

  # initialise S, I and R
  S <- if (pars$I0_at_steady_state > 0) round(S_star) else S0
  I <- if (pars$I0_at_steady_state > 0) round(I_star) else pars$I0
  R <- if (pars$I0_at_steady_state > 0) pars$N - round(I_star) - round(S_star) else pars$N - pars$I0 - S0

  list(S = S, I = I, R = R)
}

infections <- function(pars, IJ, SJ){
  # SIR: two types of events for S, so competing hazards. A fraction of
  # S events are deaths and the rest are infections.
  FOI <- pars$beta * IJ / pars$N
  fmu <- FOI + pars$mu
  fmudt <- fmu *pars$dt

  prob1 <- 1.0 - exp(-1.0 * fmudt)
  #prob2 <- 1.0 - exp(-1.0 * pars$mu/fmu)
  prob2 <- 1.0 - exp(-1.0 * pars$mu/fmu)

  # print("n_events_S")
  # print(prob1)
  print("n_deaths_S")
  print(prob2)

  n_events_S <- rbinom(1, SJ, prob1)
  n_deaths_S <- rbinom(1, n_events_S, prob2)
  n_infections_S <- n_events_S - n_deaths_S

  list(n_deaths_S = n_deaths_S, n_infections_S = n_infections_S)
}

recoveries <- function(pars, IJ){
  # SIR: two types of events for I, so competing hazards and a fraction of
  # I events are deaths and the rest are recoveries
  coeff <- pars$nu + pars$mu
  coeffdt <- coeff * pars$dt

  prob1 <- 1.0 - exp(-1.0 * coeffdt)
  prob2 <- 1.0 - exp(-1.0 * pars$mu/coeff)
  # print("n_events_I")
  # print(prob1)
  # print("n_deaths_I")
  # print(prob2)
  n_events_I <- rbinom(1, IJ, prob1)
  n_deaths_I <- rbinom(1, n_events_I, prob2)
  n_recoveries_I <- n_events_I - n_deaths_I

  list(n_deaths_I = n_deaths_I, n_recoveries_I = n_recoveries_I)
}

births <- function(pars, RJ, inf, rec){

  coeffdt <- pars$mu*pars$dt
  prob <- 1.0 - exp(-1.0 * coeffdt)
  # print("n_deaths_R")
  # print(prob)

  n_deaths_R <- rbinom(1, RJ, prob)
  n_births <- inf$n_deaths_S + rec$n_deaths_I + n_deaths_R

  list(n_deaths_R = n_deaths_R, n_births = n_births)
}

update <- function(pars, inf, rec, bir, SJ, IJ, RJ, j){

    news <- SJ - inf$n_deaths_S - inf$n_infections_S + bir$n_births
    newi <- IJ + inf$n_infections_S - rec$n_recoveries_I - rec$n_deaths_I
    newr <- RJ + rec$n_recoveries_I - bir$n_deaths_R
    newN <- news + newi + newr

    if(newN != pars$N){
      thetime <- pars$dt * (j + 1)
      stop(sprintf("newN is not equal to N for time = %f",
                   thetime))
    }

    list(news = news, newi = newi, newr = newr)
}

displaythemodel <- function(res) {

  sir_col <- c("#8c8cd9", "#cc0044", "#999966")
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  check = matplot(res$time[1:res$countlim], cbind(res$S[1:res$countlim], res$I[1:res$countlim], res$R[1:res$countlim]), type = "l", col = sir_col, lty = 1, xlab="time", ylab="S, I, R")
  legend("topright", inset = 0.05, lwd = 1, col = sir_col, legend = c("S", "I", "R"))
  a <- paste("N = ", res$N)
  b <- paste("dt = ", res$dt)
  c <- paste(a, b)
  d <- paste("SIR model,",  c )
  title(d, line = -2, adj = 0.2, col.main = "gray", font.main = 4)

  list(check)
}

