library(genlasso)

# ==== SETUP ====
set.seed(42)
n_data <- 200
x <- seq(1, n_data, 1)
y <- sin(2 * pi * x / n_data) + rnorm(n_data, 0, 0.4)


# ==== FINITE DIFFERENCE MATRIX ====
build_diff_matrix <- function(n, order = 1) {
  # For trend filtering, we use a difference operator of order (order + 1)
  diff(diag(n), differences = order + 1)
}

D <- build_diff_matrix(n_data, order = 3)

# ==== SMOOTHED QUANTILE LOSS ====
smooth_check_loss <- function(u, tau, eps = 0.05) {
  tau * u + eps * log(exp(-u / eps))
}

grad_smooth_check <- function(y, theta, tau, eps = 0.05) {
  u <- y - theta
  sigmoid <- 1 / (1 + exp(u / eps))
  -(tau - sigmoid)
}


# === NEGATIVE LOG POSTERIOR ===
neg_log_posterior <- function(theta, y, tau, lambda) {
  resid <- y - theta
  f_theta <- sum(smooth_check_loss(resid, tau, eps = 0.05))
  g_theta <- lambda * sum(abs(as.vector(D %*% theta)))
  return(f_theta + g_theta)
}


# ==== PROX USING PATH INDEX ====
prox_g <- function(z, D, gamma, lambda) {
  eff <- gamma * lambda
  fit <- genlasso(z, D = D)
  idx <- which.min(abs(fit$lambda - eff))
  as.numeric(coef(fit, i = idx)$beta[,1])
}


# ==== P-ULA SAMPLER ====
pula_sampler <- function(y, tau = 0.5, lambda = 50,
                         D, delta = 0.15, gamma = 0.1, n_iter = 500,
                         eps = 0.05, burn_in = 200) {
  n <- length(y)
  theta_chain <- matrix(NA, nrow = n_iter, ncol = n)
  theta <- rep(0, n)

  for (i in 1:n_iter) {
    grad <- grad_smooth_check(y, theta, tau, eps)
    z <- theta - delta * grad + sqrt(2 * delta) * rnorm(n)
    theta <- prox_g(z, D, gamma, lambda)
    theta_chain[i, ] <- theta
  }

  theta_post <- theta_chain[(burn_in + 1):n_iter, , drop = FALSE]
  list(samples = theta_post,
       estimate = colMeans(theta_post))
}

# === P-ULA WARM-UP ===

theta_current <- pula_sampler(y, tau = 0.5, D = D)

theta_current <- theta_current$estimate


# ==== GRADIENT OF (f + g) ====
grad_f_plus_g <- function(y, theta, tau, eps, gamma, D, lambda) {
  # Step 1: Gradient of smooth likelihood f
  grad_f_theta <- grad_smooth_check(y, theta, tau, eps)
  # Step 2: Forward gradient step
  x_tilde <- theta - gamma * grad_f_theta

  # Step 3: Proximal mapping for nonsmooth prior g
  prox_step <- prox_g(x_tilde, D, gamma, lambda)

  # Step 4: Moreau gradient approximation
  grad_total <- (theta - prox_step) / gamma

  return(grad_total)
}

acceptance_count <- 0

pmala_sampler <- function(y, tau = 0.5, lambda = 50, gamma = 0.1,
                         D, delta = 0.15, n_iter = 2000,
                         eps = 0.05, burn_in = 500) {
  n <- length(y)
  theta_chain <- matrix(NA, nrow = n_iter, ncol = n)

  for (i in 1:n_iter) {
    grad <- grad_smooth_check(y, theta_current, tau, eps)
    z <- theta_current - delta * grad + sqrt(2 * delta) * rnorm(n)
    theta_proposal <- prox_g(z, D, gamma, lambda)

    log_posterior_current <- -neg_log_posterior(theta_current, y, tau, lambda)
    log_posterior_proposal <- -neg_log_posterior(theta_proposal, y, tau, lambda)
    log_posterior_ratio <- log_posterior_proposal - log_posterior_current

    # Compute mu(θ) for transition density using fast approximation
    grad_psi_current <- grad_f_plus_g(y, theta_current, tau, eps, gamma, D, lambda)
    mu_current <- theta_current - delta * grad_psi_current

    grad_psi_proposal <- grad_f_plus_g(y, theta_proposal, tau, eps, gamma, D, lambda)
    mu_proposal <- theta_proposal - delta * grad_psi_proposal

    # Compute log proposal densities
    log_q_proposal_given_current <- -sum((theta_proposal - mu_current)^2) / (4*delta)
    log_q_current_given_proposal <- -sum((theta_current - mu_proposal)^2) / (4*delta)
    log_proposal_ratio <- log_q_current_given_proposal - log_q_proposal_given_current

    # Compute log acceptance probability
    log_alpha <- min(0, log_posterior_ratio + log_proposal_ratio)
    #print(exp(log_posterior_ratio + log_proposal_ratio))

    # Step 4: Accept/reject
    if(log(runif(1)) < log_alpha) {
      theta_current <- theta_proposal
      if(i > burn_in) acceptance_count <- acceptance_count + 1
    }

    # Store post burn-in samples
    if(i > burn_in) {
      samples[i - burn_in, ] <- theta_current
    }
    # Print progress
    if(i %% 100 == 0) {
      if(i <= burn_in) {
        cat("Burn-in iation:", i, "\n")
      } else {
        cat("Sampling iation:", i - burn_in,
            "Acceptance rate:", acceptance_count/(i - burn_in), "\n")
      }
    }

    theta_chain[i, ] <- theta_current

  }
  theta_post <- theta_chain[(burn_in + 1):n_iter, , drop = FALSE]
  list(samples = theta_post,
       estimate = colMeans(theta_post))

}

median_result <- pmala_sampler(y, tau = 0.5, lambda = 50, D = D, n_iter = 500, burn_in = 200) # Just to see if the code works, iter will be increased
theta_median <- median_result$estimate
