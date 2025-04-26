test_that("Identical case returns 0", {
  mu1 <- c(0, 0)
  sigma1 <- matrix(c(1, 0, 0, 1), 2, 2)
  mu2 <- c(0, 0) 
  sigma2 <- matrix(c(1, 0, 0, 1), 2, 2)
  
  kl_div <- kl_divergence_gaussian(mu1, sigma1, mu2, sigma2)
  expect_equal(kl_div, -Inf, tolerance = 1e-10)
})

test_that("Invalid inputs throw an error", {
  mu1 <- c(0, 0)
  sigma1 <- matrix(c(1, 0), 2, 1)  # Not a square matrix
  mu2 <- c(0, 0, 0) 
  sigma2 <- matrix(c(1, 0, 0, 1), 2, 2)
  
  expect_error(kl_divergence_gaussian(mu1, sigma1, mu2, sigma2))
})

test_that("Known good test case", {
  mu1 <- c(0, 0) 
  sigma1 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  mu2 <- c(1, 1)
  sigma2 <- matrix(c(2, 0, 0, 2), 2, 2)
  
  kl_div <- kl_divergence_gaussian(mu1, sigma1, mu2, sigma2)
  expect_equal(kl_div, log(0.8369882), tolerance = 1e-6)
})


test_that("Univariate setting for KL Gaussian approx", {
  mu1 <- c(0)
  sigma1 <- matrix(1, 1, 1)

  mu2 <- c(1)
  sigma2 <- matrix(2, 1, 1)
  kl_div <- kl_divergence_gaussian(mu1, sigma1, mu2, sigma2)
  expect_equal(kl_div, -1.05966, tolerance = 1e-6)
})
