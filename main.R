library(ggplot2)

# PART 1

# (A)
# f/g maximized at x = +-1 with value sqrt(2*pi/e)
# M := sqrt(2*pi / e)
# 1/M probability of acceptance
# M expected number of trials per accepted sample

# RANDOM_STATE <- 41
# set.seed(RANDOM_STATE)

N <- 1000
M <- sqrt(2*pi / exp(1))
n_sampled <- 0
n_trials <- 1
samples <- rep(NA, N)
while (n_sampled < N) {
  n_trials <- n_trials + 1
  x <- rcauchy(1)
  f <- dnorm(x)
  g <- dcauchy(x)
  p_accept <- f / (M*g)
  u <- runif(1)
  if (u < p_accept) {
    n_sampled <- n_sampled + 1
    samples[n_sampled] <- x
  }
}
acceptance_rate <- n_sampled / n_trials

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
  1/M, acceptance_rate,
  M, 1/acceptance_rate,
  0, mean(samples), 1, sd(samples)
))


ggplot(data.frame(samples), aes(x=samples, y=..density..)) +
  geom_histogram(binwidth = 3.491 * 1 * N^(-1/3)) +  # slides 1, slide 11
  geom_density(aes(color='Epanechnikov KDE'),
               kernel='epanechnikov',
               key_glyph=draw_key_path) +
  stat_function(fun = dnorm, inherit.aes = FALSE,
                aes(color='Normal Distribution'),
                key_glyph=draw_key_path) +
  labs(color="")

cat(sprintf(
  '\n\nTrue Mean: %.4f\nEstimated Mean: %.4f',
  0, mean(samples)
))


# B



