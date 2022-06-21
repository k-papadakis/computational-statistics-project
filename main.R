library(ggplot2)

seed <- 7

# PART 1

# 1.A
# f/g maximized at x = +-1 with value sqrt(2*pi/e)
# M := sqrt(2*pi / e)
# 1/M probability of acceptance
# M expected number of trials per accepted sample


rejection_sample_normal <- function(N, seed = NULL) {
  # Rejection sample from the standard normal
  # using the standard Cauchy scaled optimally
  set.seed(seed)
  M <- sqrt(2 * pi / exp(1))
  n_sampled <- 0
  n_trials <- 1
  samples <- rep(NA, N)
  while (n_sampled < N) {
    n_trials <- n_trials + 1
    x <- rcauchy(1)
    f <- dnorm(x)
    g <- dcauchy(x)
    p_accept <- f / (M * g)
    u <- runif(1)
    if (u < p_accept) {
      n_sampled <- n_sampled + 1
      samples[n_sampled] <- x
    }
  }
  return(list(
    samples = samples,
    n_trials = n_trials,
    acceptance_rate = N / n_trials,
    M = M,
    N = N
  ))
}


rs <- rejection_sample_normal(N = 1000, seed = seed)

cat(sprintf(
  '
True Acceptance Rate: %.4f
Estimated Acceptance Rate: %.4f

True Expected #Trials/Sample: %.4f:
Estimated Expected #Trials/Sample: %.4f

True Mean: %.4f
Estimated Mean: %.4f

True SD: %.4f
Estimated SD: %.4f
',
  1 / rs$M, rs$acceptance_rate,
  rs$M, 1 / rs$acceptance_rate,
  0, mean(rs$samples), 1, sd(rs$samples)
))

ggplot(data.frame(rs$samples), aes(x = rs$samples, y = ..density..)) +
  geom_histogram(binwidth = 3.491 * 1 * rs$N^(-1 / 3)) +  # slides 1, slide 11
  geom_density(aes(color = 'Epanechnikov KDE'),
               kernel = 'epanechnikov',
               key_glyph = draw_key_path) +
  stat_function(fun = dnorm, inherit.aes = FALSE,
                aes(color = 'Normal Distribution'),
                key_glyph = draw_key_path) +
  labs(color = '')


# 1.B
raoblackwell_nbinom_mean <- function(n, r, p, n_iter = 1, seed = NULL) {
  # p is the probability of failure
  set.seed(seed)

  original <- array(dim = c(n, n_iter))
  improved <- array(dim = c(n, n_iter))
  ns <- 1 : n

  for (iter in 1:n_iter) {
    lambda <- rgamma(n, shape = r, scale = (1 - p) / p)
    x <- rpois(n, lambda)
    original[, iter] <- cumsum(x) / ns
    improved[, iter] <- cumsum(lambda) / ns
  }

  return(list(
    ns = ns,
    original = original,
    improved = improved,
    est_mean_original = apply(original, 1, mean),
    est_mean_improved = apply(improved, 1, mean),
    true_mean = r * (1 - p) / p,
    est_var_original = apply(original, 1, var),
    est_var_improved = apply(improved, 1, var),
    true_var_original = function(n) {r * p/(1-p^2) / n},
    true_var_improved = function(n) {r * p^2/(1-p^2) / n}
  ))
}


rb <- raoblackwell_nbinom_mean(n = 5000, r = 1, p = 0.5, n_iter = 100, seed = seed)
# var_ratio <- rb$est_var_improved / rb$est_var_original

n_skip <- 100
realization_id <- 2

# Plot a single realization for delta_n (all n)
ggplot(
  data.frame(
    n = rb$ns[-(1:n_skip)],
    estimator_value = c(
      rb$original[-(1:n_skip), realization_id],
      rb$improved[-(1:n_skip), realization_id]
    ),
    estimator_type = rep(
      c('Original', 'Rao-Blackwellized'),
      each = length(rb$est_var_original) - n_skip
    )
  ),
  aes(x = n, y = estimator_value, color = estimator_type),
) +
  geom_line() +
  geom_hline(
    aes(
      yintercept = rb$true_mean,
      color = 'Estimator Expected Value')
  ) +
  labs(
    x = 'n', y = 'Estimator Value',
    title = 'Estimation of the Mean',
    color = 'Estimator Type',
    linetype = 'Value'
  )

# Plot the variance of the realizations of the delta_n's (all n)
ggplot(
  data.frame(
    n = rb$ns[-(1:n_skip)],
    estimated_estimator_var = c(
      rb$est_var_original[-(1:n_skip)],
      rb$est_var_improved[-(1:n_skip)]
    ),
    estimator_type = rep(
      c('Original', 'Rao-Blackwellized'),
      each = length(rb$est_var_original) - n_skip
    ),
    true_estimator_var = c(
      rb$true_var_original(rb$ns[-(1:n_skip)]),
      rb$true_var_improved(rb$ns[-(1:n_skip)])
    )
  ),
  aes(x = n, color = estimator_type)
) +
  geom_line(
    aes(
      y = estimated_estimator_var,
      linetype = 'Estimated'
    )
  ) +
  geom_line(aes(y = true_estimator_var, linetype = 'Real')) +
  labs(
    x = 'n', y = 'Estimator Variance',
    title = 'Rao-Blackwellization',
    color = 'Estimator Type',
    linetype = 'Estimator Variance'
  )





# 1.C
# E[sum(X_i)^2] = V(sum(X_i)) + (E[sum(X_i)])^2 = n/12 + n^2/4
# E[sum(X_i)^2 / n] = 1/12 + n/4
# V[sum(X_i)^2 / n] = 1/n^2 (E[sum(X_i)^4] - E[sum(X_i)^2])
# E[sum(X_i)^4] = sum(E[X_i X_j X_k X_l]) =
#   = (n)_4 E[X]^4 +  # all indices i,j,k,l different
#      + choose(4, 2) (n)_3 E[X]^2 E[X^2] +  # two same, two different
#      + choose(4, 2)/2 (n)_2 E[X^2] E[X^2] +  # two same, two same
#      + choose(4, 3) (n)_2 E[X^3] E[X^3] E[X] +  # three same
#      + n E[X^4]  # all same
# Using the fact that E[X^k] = 1 / (k+1) when X ~ U(0,1)
# we finally get
# V[sum(X_i)^2 / n] = (5 - 3/n + 30*n) / 360
#
# Now, considering the distribution of T:
# sum of n=80 uniforms is essentially Normal(n/2, n/12) = sqrt(n/12) Normal(n/2, 1)
# https://en.wikipedia.org/wiki/Irwin%E2%80%93Hall_distribution
# The square of Normal(n/2, 1) is non-central chi squared with k=1, lambda = (n/2)^2.
# For large lambda, this becomes approximately N(k+lambda, 2(k + 2 lambda))
# https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Properties
# Thus, approximately T ~ 1/n * sqrt(n/12) * Normal(1 + (n/2)^2, 2(1 + 2 (n/2)^2)
# So, roughly we have sigma^2 ~= 2(1 + 2*(n/2)^2) * (1/n^2) * n/12 = n/12 + 1/(6n)
# For n = 80, this evaluates to 2.582392.

generate_ts <- function(num_samples, n, seed=NULL) {
  set.seed(seed)
  ts <- array(dim=num_samples)
  for (i in 1:num_samples) {
    xs <- runif(n)
    t <- sum(xs)^2 / n
    ts[i] <- t
  }
  return(ts)
}
get_t_mean <- function(n) {n/4 + 1/12}
get_t_sd <- function(n) {sqrt((5 - 3/n + 30*n) / 360)}  # Exact value


ts <- generate_ts(num_samples = 10000, n = 80, seed=seed)

cat(sprintf(
  '
Estimated mean: %.4f
True mean: %.4f

Estimated SD: %.4f
True SD: %.4f
  ',
  mean(ts), get_t_mean(80),
  sd(ts), get_t_sd(80)
))

ggplot(data.frame(ts), aes(x = ts, y = ..density..)) +
  geom_histogram(binwidth = 3.491 * min(IQR(ts)/1.345, sd(ts)) * length(ts)^(-1/3)) +  # slides 1, slide 11
  geom_density(aes(color = 'Gaussian KDE'),
               kernel = 'epanechnikov',
               key_glyph = draw_key_path) +
  stat_function(fun = dnorm, args = list(mean = get_t_mean(80), sd = get_t_sd(80)),
                inherit.aes = FALSE,
                aes(color = 'Normal Distribution'),
                key_glyph = draw_key_path) +
  labs(x = 'T', color = '')


# 1.D
bootstrap <- function(data, statistic, R, seed = NULL) {
  set.seed(seed)
  t <- numeric(R)
  for (i in 1:R) {
    indices <- sample(seq_along(data), length(data), replace = TRUE)
    t[i] <- statistic(data, indices)
  }
  return(list(
    t = t,
    data = data,
    statistic = statistic
  ))
}


jackknife <- function(data, statistic) {
  t_leaveoneout <- numeric(length(data))
  for (i in seq_along(data)) {
    t_leaveoneout[i] <- statistic(data[-i])
  }
  t_entire <- statistic(data)
  n <- length(data)

  t_leaveoneout_mean <- mean(t_leaveoneout)
  jkval <- n * t_entire - (n - 1) * t_leaveoneout_mean
  # Estimates about properties of the original estimator
  bias <- (n - 1) * (t_leaveoneout_mean - t_entire)
  se <- sqrt((n-1) * mean((t_leaveoneout - t_leaveoneout_mean)^2))
  # values <- t_leaveoneout
  return(list(
    val = jkval,
    vals = t_leaveoneout,
    bias = bias,
    se = se,
    statistic = statistic
  ))
}


t_statistic <- function(data, indices = NULL) {
  if (is.null(indices)) {
    indices <- seq_along(data)
  }
  sample <- data[indices]
  t <- sum(sample)^2 / length(sample)
  return(t)
}


xs <- readRDS('data1.rds')

# Bootstrap
bs <- bootstrap(xs, t_statistic, 10000, seed = seed)
ts_boot <- bs$t

# Jackknife
jk <- jackknife(xs, t_statistic)

cat(sprintf(
  '
Bootstrap Estimate of SE: %.4f
Jackknife Estimate of SE: %.4f
  ',
  sd(bs$t), jk$se
))

ggplot(
  data.frame(
    T = c(ts, bs$t),
    Method = rep(c('Regular', 'Bootstrap'), each=length(ts))
  ),
  aes(T, ..density.., fill = Method)
) +
  geom_histogram(
    binwidth = 3.491 * min(IQR(ts)/1.345, sd(ts)) * length(ts)^(-1/3),
    alpha = 0.3,
    position = 'identity'
  )  +
  stat_function(
    fun = dnorm, args = list(mean = get_t_mean(80), sd = get_t_sd(80)),
    inherit.aes = FALSE,
    color = 'black',
    alpha = 0.9,
    key_glyph = draw_key_path
  )


# PART 2


