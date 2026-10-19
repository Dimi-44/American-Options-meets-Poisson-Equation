if (!requireNamespace("RQuantLib", quietly = TRUE)) install.packages("RQuantLib")
library(RQuantLib)
price_american_option_sor <- function(
    S0 = 100.0, K = 100.0, T = 1.0, r = 0.05, sigma = 0.2,
    option_type = "put", S_max = 300.0, M = 100, N = 1000,
    omega = 1.2, tol = 1e-6, align_strike = FALSE
) {
  if (align_strike) {
    k_node <- round((K / S_max) * M)
    dS <- K / k_node
    S_max <- dS * M
    omega <- 2 / (1 + sin(pi / M))
  } else {
    dS <- S_max / M
  }
  dt <- T / N
  S <- seq(0, S_max, length.out = M + 1)
  payoff <- if (tolower(option_type) == "put") pmax(K - S, 0) else pmax(S - K, 0)
  V <- payoff
  i <- 1:(M - 1)
  a <- -0.5 * dt * (sigma^2 * i^2 - r * i)
  b <- 1 + dt * (sigma^2 * i^2 + r)
  c <- -0.5 * dt * (sigma^2 * i^2 + r * i)
  for (n in 1:N) {
    V_old <- V
    tau <- n * dt
    if (tolower(option_type) == "put") {
      V[1] <- K * exp(-r * tau)
      V[M + 1] <- 0.0
    } else {
      V[1] <- 0.0
      V[M + 1] <- S_max - K * exp(-r * tau)
    }
    for (iter in 1:10000) {
      max_diff <- 0.0
      for (k in 2:M) {
        idx <- k - 1
        gs_val <- (V_old[k] - a[idx] * V[k - 1] - c[idx] * V[k + 1]) / b[idx]
        V_new  <- (1 - omega) * V[k] + omega * gs_val
        V_new  <- max(V_new, payoff[k])
        diff <- abs(V_new - V[k])
        if (diff > max_diff) max_diff <- diff
        V[k] <- V_new
      }
      if (max_diff < tol) break
    }
  }
  return(approx(S, V, xout = S0)$y)
}
run_r_benchmark <- function() {
  cat("\n=============================================================================\n")
  cat(" R AMERICAN OPTION BENCHMARK: BASE SOR vs. OPTIMIZED SOR vs. RQUANTLIB\n")
  cat("=============================================================================\n\n")
  base_put  <- price_american_option_sor(M = 100, N = 1000, omega = 1.2, option_type = "put")
  opt_put   <- price_american_option_sor(M = 300, N = 3000, align_strike = TRUE, option_type = "put")
  ql_put    <- AmericanOption("put", 100, 100, 0, 0.05, 1.0, 0.2, timeSteps = 1000)$value
  base_call <- price_american_option_sor(M = 100, N = 1000, omega = 1.2, option_type = "call")
  opt_call  <- price_american_option_sor(M = 300, N = 3000, align_strike = TRUE, option_type = "call")
  ql_call   <- AmericanOption("call", 100, 100, 0, 0.05, 1.0, 0.2, timeSteps = 1000)$value
  cat(sprintf("%-6s | %-12s | %-12s | %-12s | %-12s | %-12s\n", 
              "Type", "Base SOR", "Opt SOR", "RQuantLib", "Err (Base)", "Err (Opt)"))
  cat(paste(rep("-", 78), collapse = ""), "\n")
  cat(sprintf("%-6s | %-12.4f | %-12.4f | %-12.4f | %-12.6f | %-12.6f\n", 
              "Put", base_put, opt_put, ql_put, abs(base_put - ql_put), abs(opt_put - ql_put)))
  cat(sprintf("%-6s | %-12.4f | %-12.4f | %-12.4f | %-12.6f | %-12.6f\n", 
              "Call", base_call, opt_call, ql_call, abs(base_call - ql_call), abs(opt_call - ql_call)))
  cat("\n")
}
run_r_benchmark()
