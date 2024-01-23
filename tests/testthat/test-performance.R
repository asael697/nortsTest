################################################################
#. Time performance
################################################################

set.seed(1697231)
y = arima.sim(n = 3000, model = list(ar = 0.3, ma = 0.2))

test_that("Vavra's test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::vavra.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Bootstrap Jarque Beras' test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::jb_bootstrap.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Bootstrap Shapiro's test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::shapiro_bootstrap.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Bootstrap Cramer Von Mises' test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::cvm_bootstrap.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Lobato's test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::lobato.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Bootstrap Lobato's test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::lobato_bootstrap.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Epps' test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::epps.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Bootstrap Epps' test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::epps_bootstrap.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})

test_that("Random Projections' test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::rp.test(y)
  end = difftime(Sys.time(), start, units="secs")
  expect_equal(end < 30, TRUE)
})
