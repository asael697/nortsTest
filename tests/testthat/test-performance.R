################################################################
#. Time performance
################################################################

set.seed(1697231)
y = arima.sim(n = 3000,model = list(ar = 0.3,ma = 0.2))

test_that("Vavra's test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::vavra.test(y)
  end = Sys.time() - start
  expect_equal(end < 10, TRUE)
})

test_that("Lobato's test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::lobato.test(y)
  end = Sys.time() - start
  expect_equal(end < 10, TRUE)
})

test_that("Epps' test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::epps.test(y)
  end = Sys.time() - start
  expect_equal(end < 10, TRUE)
})

test_that("Random Projections' test performance lower than 10s", {
  start = Sys.time()
  ht = nortsTest::rp.test(y)
  end = Sys.time() - start
  expect_equal(end < 10, TRUE)
})
