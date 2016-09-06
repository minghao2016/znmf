kl_divergence <- function(a, b, lambda=sqrt(.Machine$double.eps)) {
  sum(a * log(a / (b + lambda) + lambda) - a + b)
}

znmf <- function(v, k, niter=100, tol=1e-6, peek_interval=10, lambda=sqrt(.Machine$double.eps)) {
  n <- nrow(v)
  m <- ncol(v)

  alpha <- sqrt(2 * mean(v) / k)

  w <- matrix(runif(n * k, 0, alpha), nrow=n, ncol=k)
  h <- matrix(runif(k * m, 0, alpha), nrow=k, ncol=m)

  errors <- numeric()
  last_error <- Inf
  for (i in 1:niter) {
    h <- h * (t(w) %*% (v / (w %*% h + lambda))) / colSums(w)
    w <- w * ((v / (w %*% h + lambda)) %*% t(h)) / rowSums(h)
    if (i %% peek_interval == 0) {
      error <- kl_divergence(v, w %*% h)
      message(sprintf('(%d) error: %f; diff: %f', i, error, last_error - error))
      errors <- c(errors, error)
      if (last_error - error < tol) break
      last_error <- error
    }
  }

  list(
    h = h,
    w = w,
    prod = w %*% h,
    errors = errors
  )
}

n <- 5000
m <- 200

v <- matrix(c(rep(c(rep(1, n), rep(0, n)), m),
        rep(c(rep(0, n), rep(1, n)), m)), n * 2)

res <- znmf(v, 7, 300)
saveRDS(res, 'res.RDS')
res <- readRDS('res.RDS')

res$prod %>% zh
res$w %>% zh
res$h %>% zh
rwl %>% zh


se <- readRDS('out-datasets/batchSE3.RDS')
rwl <- log(se@assays$data$counts+1)
dim(rwl)

res <- znmf(rwl, 3, 3000)

errors <- res$errors %>%
  enframe('step', 'error') %>%
  mutate(last_error=lag(error)) %>%
  mutate(diff=last_error - error) %>%
  na.omit

ggplot() +
  geom_line(aes(step, diff), errors) +
  scale_y_log10()
