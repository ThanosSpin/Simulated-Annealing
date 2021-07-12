# Clear R environment
rm(list = ls())


# Schaffer function
schaffer <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- (sin(x1^2 - x2^2))^2 - 0.5
  fact2 <- (1 + 0.001*(x1^2 + x2^2))^2
  
  y <- 0.5 + fact1/fact2
  y
}

# SA funnction 
simulated_annealing <- function(func, s0, niter, step = 0.1) {
  
  # Initialize
  ## s stands for state
  ## f stands for function value
  ## b stands for best
  ## c stands for current
  ## n stands for neighbor
  s_b <- s_c <- s_n <- s0
  f_b <- f_c <- f_n <- func(s_n)
  message("It\tBest\tCurrent\tNeigh\tTemp")
  message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", 0L, f_b, f_c, f_n, 1))
  
  for (k in 1:niter) {     
    Temp <- (1 - step)^k
    # consider a random neighbor
    s_n <- rnorm(2, s_c, 1)
    f_n <- func(s_n)
      # update current state
    if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
      s_c <- s_n
      f_c <- f_n
    }
    # update best state
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n         
    }
    message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", k, f_b, f_c, f_n, Temp))
  }
  return(list(iterations = niter, best_value = f_b, best_state = s_b))
}

set.seed(12345)
sol <- simulated_annealing(schaffer, c(-10, 10), 1000)
sol 
#global optimum for schaffer function is a global minimum with f(xx) = 0, at (0,0)
# SA algorithm approaches the true global optimum
