################################################################
#. Gaussian  ARMA
################################################################


test_that("Epps' statistics  works for GPs", {
  set.seed(169721)
  y = arima.sim(n = 600, model = list(ar = 0.3, ma = 0.2))
  ht = nortsTest::epps.statistic(y = y,lambda = c(1,2))
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(ht, 2.919187,tolerance = 1e-3)
})

test_that("Epps' statistics  works for GPs, and custom lambda", {
  set.seed(169721)
  y = arima.sim(n = 600, model = list(ar = 0.3, ma = 0.2))
  ht = nortsTest::epps.statistic(y = y,lambda = c(1,2,3,4,5))
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(ht, 7.205453,tolerance = 1e-3)
})


test_that("Epps' statistics  works for non GPs", {
  set.seed(169721)
  y = arima.sim(n = 600, model = list(ar = 0.3, ma = 0.2),
                rand.gen = rgamma, shape = 0.4)
  ht = nortsTest::epps.statistic(y = y,lambda = c(1,2))
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(ht, 91.6286,tolerance = 1e-3)
})

test_that("Epps' statistics  works for non GPs, and custom lambda", {
  set.seed(169721)
  y = arima.sim(n = 600, model = list(ar = 0.3, ma = 0.2),
                rand.gen = rgamma, shape = 0.4)
  ht = nortsTest::epps.statistic(y = y,lambda = c(1,2,3,4,5))
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(ht, 312.4624,tolerance = 1e-3)
})
