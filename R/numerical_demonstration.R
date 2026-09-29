#!/usr/bin/env Rscript

# Numerical illustration for "Approximation by mixtures of multivariate
# Erlang distributions".  The script uses base R only.  It regenerates the
# manuscript figure and the compact CSV summaries from a fixed random seed.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_file))
results_dir <- file.path(repo_root, "results")
figures_dir <- file.path(repo_root, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
seed <- 20260926L
set.seed(seed)

if (quick) {
  B_density <- 20L
  B_resolution <- 20L
  B_functional <- 100L
  grid_size <- 20L
} else {
  B_density <- 150L
  B_resolution <- 120L
  B_functional <- 1000L
  grid_size <- 100L
}

# Target law: a nonsingular bivariate lognormal model from the same family as
# the insurance example discussed in the manuscript.  Its density is known;
# aggregate-loss benchmarks are evaluated by one-dimensional quadrature.
mu <- c(0, 0.25)
sigma <- c(0.60, 0.80)
rho <- 0.55
conditional_sigma <- sigma[2L] * sqrt(1 - rho^2)

r_target <- function(n) {
  z1 <- rnorm(n)
  z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  cbind(exp(mu[1L] + sigma[1L] * z1),
        exp(mu[2L] + sigma[2L] * z2))
}

d_target <- function(x1, x2) {
  z1 <- (log(x1) - mu[1L]) / sigma[1L]
  z2 <- (log(x2) - mu[2L]) / sigma[2L]
  exponent <- -(z1^2 - 2 * rho * z1 * z2 + z2^2) /
    (2 * (1 - rho^2))
  exp(exponent) /
    (2 * pi * sigma[1L] * sigma[2L] * sqrt(1 - rho^2) * x1 * x2)
}

target_tail <- function(u) {
  integrand <- function(z1) {
    x1 <- exp(mu[1L] + sigma[1L] * z1)
    conditional_mu <- mu[2L] + rho * sigma[2L] * z1
    conditional_tail <- rep(1, length(z1))
    below <- x1 < u
    conditional_tail[below] <- plnorm(
      u - x1[below], meanlog = conditional_mu[below],
      sdlog = conditional_sigma, lower.tail = FALSE
    )
    dnorm(z1) * conditional_tail
  }
  integrate(integrand, -10, 10, rel.tol = 1e-10,
            subdivisions = 2000L)$value
}

conditional_lognormal_stop_loss <- function(threshold, x1, meanlog, sdlog) {
  cutoff <- threshold - x1
  mean_x2 <- exp(meanlog + 0.5 * sdlog^2)
  answer <- x1 - threshold + mean_x2
  positive <- cutoff > 0
  answer[positive] <- (x1[positive] - threshold) *
    plnorm(cutoff[positive], meanlog = meanlog[positive], sdlog = sdlog,
           lower.tail = FALSE) +
    mean_x2[positive] * pnorm(
      (meanlog[positive] + sdlog^2 - log(cutoff[positive])) / sdlog
    )
  answer
}

target_layer <- function(u, v) {
  integrand <- function(z1) {
    x1 <- exp(mu[1L] + sigma[1L] * z1)
    conditional_mu <- mu[2L] + rho * sigma[2L] * z1
    layer <- conditional_lognormal_stop_loss(
      u, x1, conditional_mu, conditional_sigma
    ) - conditional_lognormal_stop_loss(
      u + v, x1, conditional_mu, conditional_sigma
    )
    dnorm(z1) * layer
  }
  integrate(integrand, -10, 10, rel.tol = 1e-10,
            subdivisions = 2000L)$value
}

# Posterior-mean density on a rectangular midpoint grid.  Observations sharing
# a grid cell are aggregated before the product-Erlang kernels are evaluated.
erlang_density_grid <- function(x, n, grid) {
  m1 <- 1L + floor(n * x[, 1L])
  m2 <- 1L + floor(n * x[, 2L])
  key <- paste(m1, m2, sep = ":")
  cell_count <- table(key)
  index <- do.call(rbind, strsplit(names(cell_count), ":", fixed = TRUE))
  index <- matrix(as.integer(index), ncol = 2L)
  weight <- as.numeric(cell_count) / nrow(x)

  d1 <- vapply(index[, 1L], function(k) dgamma(grid, shape = k, rate = n),
               numeric(length(grid)))
  d2 <- vapply(index[, 2L], function(k) dgamma(grid, shape = k, rate = n),
               numeric(length(grid)))
  if (is.null(dim(d1))) d1 <- matrix(d1, ncol = 1L)
  if (is.null(dim(d2))) d2 <- matrix(d2, ncol = 1L)
  estimate <- tcrossprod(sweep(d1, 2L, weight, `*`), d2)
  grid_step <- grid[2L] - grid[1L]
  grid_limit <- max(grid) + 0.5 * grid_step
  squared_kernel_norm <- function(m) {
    log_constant <- log(n) + lgamma(2 * m - 1) -
      (2 * m - 1) * log(2) - 2 * lgamma(m)
    exp(log_constant) * pgamma(grid_limit, shape = 2 * m - 1,
                               rate = 2 * n)
  }
  mean_component_norm2 <- sum(
    weight * squared_kernel_norm(index[, 1L]) *
      squared_kernel_norm(index[, 2L])
  )
  attr(estimate, "occupied_cells") <- length(cell_count)
  attr(estimate, "mean_component_norm2") <- mean_component_norm2
  estimate
}

functional_inputs <- function(x, n, exceedance_u, layer_u, layer_v) {
  k <- 2L + floor(n * x[, 1L]) + floor(n * x[, 2L])
  tail_value <- pgamma(exceedance_u, shape = k, rate = n,
                       lower.tail = FALSE)
  stop_loss <- function(u) {
    (k / n) * pgamma(u, shape = k + 1, rate = n, lower.tail = FALSE) -
      u * pgamma(u, shape = k, rate = n, lower.tail = FALSE)
  }
  layer_value <- stop_loss(layer_u) - stop_loss(layer_u + layer_v)
  cbind(exceedance = tail_value, layer = layer_value)
}

posterior_functional_summary <- function(b) {
  n_obs <- nrow(b)
  estimate <- colMeans(b)
  centered <- sweep(b, 2L, estimate, `-`)
  covariance <- crossprod(centered) / (n_obs * (n_obs + 1))
  list(estimate = estimate, covariance = covariance)
}

M <- 6
dx <- M / grid_size
grid <- (seq_len(grid_size) - 0.5) * dx
true_density <- outer(grid, grid, d_target)
area <- dx^2

density_experiment <- function(N, n, B) {
  sum_estimate <- matrix(0, grid_size, grid_size)
  sum_ise <- 0
  sum_integrated_posterior_loss <- 0
  sum_occupied <- 0
  for (b in seq_len(B)) {
    estimate <- erlang_density_grid(r_target(N), n, grid)
    sum_occupied <- sum_occupied + attr(estimate, "occupied_cells")
    sum_estimate <- sum_estimate + estimate
    ise <- sum((estimate - true_density)^2) * area
    posterior_spread <- (attr(estimate, "mean_component_norm2") -
                           sum(estimate^2) * area) / (N + 1)
    sum_ise <- sum_ise + ise
    sum_integrated_posterior_loss <- sum_integrated_posterior_loss +
      ise + max(0, posterior_spread)
  }
  mean_estimate <- sum_estimate / B
  bias2 <- sum((mean_estimate - true_density)^2) * area
  mise <- sum_ise / B
  c(bias2 = bias2, variance = max(0, mise - bias2), mise = mise,
    integrated_posterior_loss = sum_integrated_posterior_loss / B,
    mean_occupied_cells = sum_occupied / B)
}

# At each sample size the density resolution follows Corollary 5.5 exactly:
# n_D = ceiling(sqrt(N)) for d = 2 and alpha = 1.  A separate diagnostic
# sweep includes this value and round comparison resolutions at N = 1000.
density_resolution <- function(N) as.integer(ceiling(sqrt(N)))
resolution_N <- 1000L
resolution_n <- sort(unique(c(20L, density_resolution(resolution_N),
                              50L, 100L, 200L)))
resolution_rows <- lapply(resolution_n, function(n) {
  out <- density_experiment(resolution_N, n, B_resolution)
  data.frame(N = resolution_N, n = n, bias2 = out[["bias2"]],
             variance = out[["variance"]], mise = out[["mise"]],
             integrated_posterior_loss =
               out[["integrated_posterior_loss"]],
             mean_occupied_cells = out[["mean_occupied_cells"]])
})
resolution_summary <- do.call(rbind, resolution_rows)

# Use round sample sizes with the exact integer square-root density rule.
# Bounded functionals use the sufficient undersmoothing choice n = N.
N_values <- c(200L, 500L, 1000L, 2000L)
exceedance_u <- 4
layer_u <- 2.5
layer_v <- 2.5
truth <- c(exceedance = target_tail(exceedance_u),
           layer = target_layer(layer_u, layer_v))

rate_rows <- vector("list", length(N_values))
for (j in seq_along(N_values)) {
  N <- N_values[j]
  n_density <- density_resolution(N)
  density_out <- density_experiment(N, n_density, B_density)

  estimates <- matrix(NA_real_, B_functional, 2L,
                      dimnames = list(NULL, names(truth)))
  posterior_sd <- estimates
  posterior_corr <- numeric(B_functional)
  posterior_covariance <- numeric(B_functional)
  for (b in seq_len(B_functional)) {
    x <- r_target(N)
    inputs <- functional_inputs(x, N, exceedance_u, layer_u, layer_v)
    post <- posterior_functional_summary(inputs)
    estimates[b, ] <- post$estimate
    posterior_sd[b, ] <- sqrt(diag(post$covariance))
    posterior_covariance[b] <- post$covariance[1L, 2L]
    posterior_corr[b] <- post$covariance[1L, 2L] /
      sqrt(post$covariance[1L, 1L] * post$covariance[2L, 2L])
  }

  bias <- colMeans(estimates) - truth
  mse <- colMeans(sweep(estimates, 2L, truth, `-`)^2)
  rate_rows[[j]] <- data.frame(
    N = N,
    n_density = n_density,
    density_bias2 = density_out[["bias2"]],
    density_variance = density_out[["variance"]],
    density_mise = density_out[["mise"]],
    density_integrated_posterior_loss =
      density_out[["integrated_posterior_loss"]],
    density_mean_occupied_cells = density_out[["mean_occupied_cells"]],
    n_functional = N,
    exceedance_truth = truth[["exceedance"]],
    exceedance_bias = bias[["exceedance"]],
    exceedance_rmse = sqrt(mse[["exceedance"]]),
    exceedance_mean_posterior_sd = mean(posterior_sd[, "exceedance"]),
    exceedance_integrated_posterior_rmse = sqrt(mean(
      (estimates[, "exceedance"] - truth[["exceedance"]])^2 +
        posterior_sd[, "exceedance"]^2
    )),
    layer_truth = truth[["layer"]],
    layer_bias = bias[["layer"]],
    layer_rmse = sqrt(mse[["layer"]]),
    layer_mean_posterior_sd = mean(posterior_sd[, "layer"]),
    layer_integrated_posterior_rmse = sqrt(mean(
      (estimates[, "layer"] - truth[["layer"]])^2 +
        posterior_sd[, "layer"]^2
    )),
    mean_posterior_covariance = mean(posterior_covariance),
    empirical_estimator_covariance = cov(estimates[, 1L], estimates[, 2L]),
    mean_posterior_correlation = mean(posterior_corr),
    empirical_estimator_correlation = cor(estimates[, 1L], estimates[, 2L])
  )
}
rate_summary <- do.call(rbind, rate_rows)

density_slope <- unname(coef(lm(log(density_mise) ~ log(N),
                                data = rate_summary))[2L])
exceedance_mse_slope <- unname(coef(lm(log(exceedance_rmse^2) ~ log(N),
                                       data = rate_summary))[2L])
layer_mse_slope <- unname(coef(lm(log(layer_rmse^2) ~ log(N),
                                  data = rate_summary))[2L])

write.csv(resolution_summary,
          file.path(results_dir, "density_resolution.csv"), row.names = FALSE)
write.csv(rate_summary, file.path(results_dir, "rate_summary.csv"),
          row.names = FALSE)

metadata <- c(
  sprintf("seed: %d", seed),
  sprintf("mode: %s", if (quick) "quick" else "full"),
  sprintf("density replications per design: %d", B_density),
  sprintf("resolution replications per design: %d", B_resolution),
  sprintf("functional replications per design: %d", B_functional),
  sprintf("density midpoint grid: %d x %d on [0, %g]^2", grid_size,
          grid_size, M),
  sprintf("sample sizes: %s", paste(N_values, collapse = ", ")),
  "density resolution rule: n_D = ceiling(sqrt(N))",
  sprintf("selected density resolutions: %s",
          paste(rate_summary$n_density, collapse = ", ")),
  sprintf("resolution sweep sample size: %d", resolution_N),
  sprintf("resolution sweep values: %s", paste(resolution_n, collapse = ", ")),
  "functional resolution rule: n = N",
  sprintf("empirical log-log slope, density MISE: %.4f", density_slope),
  sprintf("empirical log-log slope, exceedance MSE: %.4f",
          exceedance_mse_slope),
  sprintf("empirical log-log slope, layer MSE: %.4f", layer_mse_slope),
  "",
  capture.output(sessionInfo())
)
writeLines(metadata, file.path(results_dir, "run_metadata.txt"))

draw_figure <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.1, 5.0, 1.8, 0.8),
      mgp = c(3.2, 0.75, 0), tcl = -0.25, las = 1,
      cex = 0.84)

  density_components <- as.matrix(
    resolution_summary[, c("bias2", "variance", "mise")]
  )
  matplot(resolution_summary$n, density_components,
          type = "b", log = "xy", pch = c(1, 2, 16), lty = c(2, 3, 1),
          col = c("#0072B2", "#D55E00", "#000000"),
          ylim = c(min(density_components) * 0.8,
                   max(density_components) * 3),
          xlab = expression("resolution " * n),
          ylab = "integrated squared error")
  abline(v = density_resolution(resolution_N), lty = 3, col = "grey45")
  legend("topright", c("squared bias", "sampling variance", "MISE"),
         pch = c(1, 2, 16), lty = c(2, 3, 1),
         col = c("#0072B2", "#D55E00", "#000000"), bty = "n",
         cex = 0.78)
  mtext(sprintf("(a) Bias-variance, N = %d", resolution_N),
        side = 3, line = 0.4,
        font = 2, cex = 0.86)

  y <- cbind(rate_summary$density_mise,
             rate_summary$exceedance_rmse^2,
             rate_summary$layer_rmse^2)
  matplot(rate_summary$N, y, type = "b", log = "xy",
          pch = c(16, 17, 15), lty = 1,
          col = c("#000000", "#009E73", "#CC79A7"),
          ylim = c(min(y) * 0.8, max(y) * 3),
          xlab = expression("sample size " * N), ylab = "mean squared error")
  ref_density <- y[1L, 1L] * (rate_summary$N / rate_summary$N[1L])^(-0.5)
  ref_functional <- y[1L, 2L] * (rate_summary$N / rate_summary$N[1L])^(-1)
  lines(rate_summary$N, ref_density, lty = 2, col = "grey35")
  lines(rate_summary$N, ref_functional, lty = 3, col = "grey35")
  legend("topright",
         c("density", "exceedance", "capped layer",
           expression(N^{-1/2}), expression(N^{-1})),
         pch = c(16, 17, 15, NA, NA), lty = c(1, 1, 1, 2, 3),
         col = c("#000000", "#009E73", "#CC79A7", "grey35", "grey35"),
         bty = "n", cex = 0.76)
  mtext("(b) Risk at prescribed resolutions", side = 3, line = 0.4,
        font = 2, cex = 0.86)
}

pdf(file.path(figures_dir, "numerical_demonstration.pdf"),
    width = 7.2, height = 3.25, useDingbats = FALSE)
draw_figure()
dev.off()

png(file.path(figures_dir, "numerical_demonstration.png"),
    width = 1800, height = 810, res = 250)
draw_figure()
dev.off()

message("Wrote results to ", results_dir)
message("Wrote figure to ", figures_dir)
