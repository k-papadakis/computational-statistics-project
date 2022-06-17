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

raoblackwell_nbinom_mean <- function(n, r, p, seed = NULL) {
  # p is the probability of failure
  set.seed(seed)
  x <- rnbinom(n, size = r, prob = p)
  lambda <- rgamma(n, shape = r, scale = (1 - p) / p)
  ns <- 1 : n
  original <- cumsum(x) / ns
  raoblackwell <- cumsum(lambda) / ns
  expected_value <- r * (1 - p) / p
  return(list(
    original = original,
    raoblackwell = raoblackwell,
    expected_value = expected_value
  ))
}


rb <- raoblackwell_nbinom_mean(n = 5000, r = 1, p = 0.5, seed = seed)

ggplot(data = data.frame(ns = seq_along(rb$original),
                         original=rb$original,
                         raoblackwell=rb$raoblackwell),
       aes(ns))+
  geom_line(aes(y = original, color = 'Original Estimator')) +
  geom_line(aes(y = raoblackwell, color = 'Rao-Blackwellized Estimator')) +
  geom_hline(aes(yintercept = rb$expected_value, color = 'Estimator Expected Value')) +
  labs(title = 'Estimation of the Mean of a Negative Binomial distribution',
       x = 'Number of Samples', y = 'Estimator Value', color = '')


# 1.C

