stop()
q('no')

install.packages('Rcpp')
install.packages('RcppArmadillo')
install.packages('inline')
library(Rcpp)
library(RcppArmadillo)
library(inline)
source('~/.z.R')
source('~/.h.R')


n <- 5000
m <- 200

v <- matrix(c(rep(c(rep(1, n), rep(0, n)), m),
              rep(c(rep(0, n), rep(1, n)), m)), n * 2)

rw <- readRDS('~/rna2/a01-NMFEM_paper/d1-data_processing/x13-process_zeisel_data/temp.RDS')

res <- znmf(rw, 7, 2000, 1e-6, 10, sqrt(.Machine$double.eps))
res$H %>% zh
res$W %>% zh


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
    w <- sweep(w * ((v / (w %*% h + lambda)) %*% t(h)), 2, rowSums(h), `/`)
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

res <- znmf(v, 2, 3000)


saveRDS(res, 'res.RDS')
res <- readRDS('res.RDS')

res$prod %>% zh
res$w %>% zh
res$h %>% zh
rwl %>% zh


se <- readRDS('~/rna2/a01-NMFEM_paper/d1-data_processing/out-datasets/batchSE3.RDS')
rwl <- log(se@assays$data$counts+1)
dim(rwl)

sourceCpp('../src/znmf.cpp')
znmf2 <- function(v, k, niter=100, tol=1e-6, peek_interval=10, lambda=sqrt(.Machine$double.eps)) {
  c_znmf(v, k, niter, tol, peek_interval, lambda)
}

microbenchmark::microbenchmark(
 znmf(rwl, 2, 30),
 znmf2(rwl, 2, 30),
 times = 3
)
 
res <- znmf(rwl, 2, 30)
res <- znmf2(rwl, 2, 30)

res$w %>% zh
res$h %>% zh
cbind(rbind(matrix(rep(0, 4), 2), res$w),
      rbind(res$h, rwl)) %>% zh
print(proc.time() - time)

errors <- res$errors %>%
  enframe('step', 'error') %>%
  mutate(last_error=lag(error)) %>%
  mutate(diff=last_error - error) %>%
  na.omit

ggplot() +
  geom_line(aes(step, diff), errors) +
  scale_y_log10()
