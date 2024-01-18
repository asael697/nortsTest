set.seed(169721)
y = arima.sim(n = 1000, model = list(ar = 0.3, ma = 0.2))

test_that("Vavra's statistic", {
  ht = nortsTest::vavra.sample(y)
  expect_equal(mean(ht), 0.4242066, tolerance = 0.3)
})

test_that("sieve_bootstrap",{
  ht = nortsTest::sieve.bootstrap(y)
  expect_equal(ncol(ht), length(y))
  expect_equal(nrow(ht), 1000)
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
  expect_equal(ht$lobato, 2.47123, tolerance = 0.3)
  #check computations are less than 2.5s
  expect_equal(ht$epps, 1.399825, tolerance = 0.3)
})

test_that("Random Projections' samples", {
  k = sample(1:10, 1)
  ht = nortsTest::rp.sample(y, k = k)
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(length(ht$lobato), k)
  #check computations are less than 2.5s
  expect_equal(length(ht$epps), k)
})

set.seed(169721)
n = 3000
y = rnorm(n)
x = rnorm(n)

test_that("2-D El Bouch's statistics", {
  ht = nortsTest::elbouch.statistic(y, x)
  expect_equal(ht[1], 8 * (n - 1) / (n + 1), tolerance = 0.5)
  expect_equal(ht[2], 64/n, tolerance = 0.5)
  expect_equal(ht[3] < 1, TRUE)
})

test_that("1-D El Bouch's statistics", {
  ht = nortsTest::elbouch.statistic(y)
  expect_equal(ht[1], 3 * (n - 1) / (n + 1), tolerance = 0.5)
  expect_equal(ht[2], 24/n, tolerance = 0.5)
  expect_equal(ht[3] < 1, TRUE)
})
