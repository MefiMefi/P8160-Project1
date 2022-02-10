library(survival)
library(tidyverse)

exp_base <- function(t, lambda = 0.5) lambda * 1 * t^0
weibull_base <- function(t, lambda = 0.5, gamma = 1.5) lambda * gamma * t^(gamma - 1)
# cox no base

#curve(exp_base, from = 0, to = 5, lty = 1, ylim = c(0, 2), ylab = expression(h[0](t)), xlab = "Follow-up time t")
#curve(weibull_base, from = 0, to = 5, lty = 2, add = TRUE)
#legend(x = "topleft", lty = 1:2, legend = c("Exponential baseline hazard", "Weibull baseline hazard"), bty = "n")


#simulate data
survival_data <- function(n, beta, lambda, gamma, alpha, distribution) {
  exponential = function(u, x, lambda, beta){
    time = -log(u)/(exp(x*beta)*lambda)
  }
  
  weibull = function(u, x, lambda, gamma, beta){
    time = ( -log(u) / (exp(x * beta) * lambda) ) ^ (1 / gamma)
  }
  
  #cox
  
  #treatment & times
  x = rbinom(n = n, size = 1, prob = 0.5)
  u = runif(n)
  
  if (distribution == "exponential") {
    time = exponential(u, x, lambda, beta) } 
  else if (distribution == "weibull") {
    time = weibull(u, x, lambda, gamma, beta) }
  
  #censor
  C = rexp(n, lambda)
  
  # follow-up times
  observed_time = pmin(time, C)
  status = 1 * (time <= C)
  
  survival_data = data.frame(id = 1:n,
                             time = time,
                             observed_time = observed_time,
                             status = status,
                             x = x)
}

simulate = function(sim, n, beta, dist = "exponential") {

  exp_beta = rep(NA, sim)
  weibull_beta = rep(NA, sim)
  cox_beta = rep(NA, sim)

  for (i in 1:sim) {
      data = survival_data(n, beta, distribution = dist, 
                           lambda = 0.1, gamma = 4, alpha = 4)
      # fit
      fit_exp = survreg(Surv(data$observed_time, data$status) ~ data$x, dist = "exponential") 
      fit_weibull = survreg(Surv(data$observed_time, data$status) ~ data$x, dist = "weibull")
      fit_cox = coxph(Surv(data$observed_time, data$status) ~ data$x)
      
      # Save coefficients 
      exp_beta[i] = -fit_exp$coefficients[-1]
      weibull_beta[i] = -fit_weibull$coefficients[-1] / fit_weibull$scale
      cox_beta[i] = fit_cox$coefficients[1]
      }
  
  coef = tibble(exp = exp_beta,
                weibull = weibull_beta,
                cox = cox_beta)
  
  results = list(coef)
  names(results) = c("coefficients")
  return(results)
}

tibble(sample_size = c(300, 350, 400, 450, 500)) %>% 
  mutate(
    output = map(.x = sample_size, 
                          ~simulate(sim = 1000, n = .x, beta = 4))) 

