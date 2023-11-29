set.seed(169721)
y = arima.sim(n = 1000,model = list(ar = 0.3,ma = 0.2),)

test_that("Vavra's statistic", {
  ht = nortsTest::vavra.sample(y)
  expect_equal(mean(ht), 0.4242066, tolerance = 0.3)
})

test_that("Lobato's statistic", {
  ht = nortsTest::lobato.statistic(y)
  expect_equal(ht, 1.3353, tolerance = 0.3)
})

test_that("Epps' statistic", {
  ht = nortsTest::epps.statistic(y)
  expect_equal(ht, 2.72, tolerance = 0.3)
})

test_that("Random Projections' statistics", {
  ht = lapply(nortsTest::rp.sample(y),mean)
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(ht$lobato, 1.146889, tolerance = 0.3)
  #check computations are less than 2.5s
  expect_equal(ht$epps, 2.219923, tolerance = 0.3)
})
