
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sirstochastic

<!-- badges: start -->

<!-- badges: end -->

The goal of package sirstochastic is to run simulations of single and
multiple iterations of the SIR stochastic model.

## Installation

You can install the released version of sirstochastic from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sirstochastic")
install.packages("png")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("reside-ic/sirstochastic")
```

## Equations

The equations defining the model are:

\#\#dS/dt = - beta \* S \* I / N \#\#dI/dt = beta \* S \* I / N - sigma
\* I \#\#dR/dt = sigma \* I

and these are converted to the discrete stochastic model by imagining
that in a small period of time `dt` the chance of an event happening
follows a Bernoulli trial. In a discrete time step of size `dt` we model
the number of infections (i.e., the number of individuals who move from
`S` to `I`) as a binomial draw with `n = S` and `p = beta * I / N * dt`)
and the number of recoveries (movements from `I` to `R`) as a binomial
draw with `n = I` and `p = sigma * dt`. Note, N = S + I + R.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sirstochastic)
library(png)

model <- function() {

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

  dt <- 0.01

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
    n_infections_S <- n_events_S - n_deaths_S   # SIR: ...the rest are infections.

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
}
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

An example of an SIR plot is:

<img src="man/figures/README-stochastic-1.png" width="100%" /> \#\#
License [MIT](https://choosealicense.com/licenses/mit/)

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!