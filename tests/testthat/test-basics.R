################################################################
#. Gaussian  ARMA
################################################################

set.seed(169721)
y = arima.sim(n = 600,model = list(ar = 0.3,ma = 0.2))

test_that("Vavra's test works for GPs", {
  start = Sys.time()
  ht = nortsTest::vavra.test(y)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})

test_that("Lobato's test works for GPs", {
  start = Sys.time()
  ht = nortsTest::lobato.test(y)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})

test_that("Epps' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::epps.test(y)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})

test_that("Random Projections' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::rp.test(y,k = 16)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})

################################################################
#. NON Gaussian  ARMA
################################################################

set.seed(169721)
y = arima.sim(n = 600,model = list(ar = 0.3,ma = 0.2),rand.gen = rgamma,shape = 0.4)

test_that("Vavra's test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::vavra.test(y)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})

test_that("Lobato's test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::lobato.test(y)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})

test_that("Epps' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::epps.test(y)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})

test_that("Random Projections' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::rp.test(y,k = 16)
  end = Sys.time() - start
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 2.5, TRUE)
})
