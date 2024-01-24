################################################################
#. Gaussian  ARMA
################################################################

set.seed(169721)
y = arima.sim(n = 600, model = list(ar = 0.3, ma = 0.2))

test_that("Vavra's test works for GPs", {
  start = Sys.time()
  ht = nortsTest::vavra.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Jarque Beras' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::jb_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Shapiro's test works for GPs", {
  start = Sys.time()
  ht = nortsTest::shapiro_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Cramer Von Mises' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::cvm_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Lobato's test works for GPs", {
  start = Sys.time()
  ht = nortsTest::lobato.test(y)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Lobato's test works for GPs", {
  start = Sys.time()
  ht = nortsTest::lobato_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Epps' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::epps.test(y)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Epps' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::epps_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Random Projections' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::rp.test(y,k = 16)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("El Bouch' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::elbouch.test(y)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

set.seed(169721)
y = rnorm(200)
x = rnorm(200)

test_that("El Bouch' test works for GPs", {
  start = Sys.time()
  ht = nortsTest::elbouch.test(y, x)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value >= 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

################################################################
#. NON Gaussian  ARMA
################################################################

set.seed(169721)
y = arima.sim(n = 600, model = list(ar = 0.3, ma = 0.2),
              rand.gen = rgamma, shape = 0.4)

test_that("Vavra's test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::vavra.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Jarque Beras' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::jb_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Shapiro test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::shapiro_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.1), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Cramer Von Mises' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::cvm_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Lobato's test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::lobato.test(y)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Lobato's test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::lobato_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Epps' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::epps.test(y)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Bootstrap Epps' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::epps_bootstrap.test(y,reps = 200)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("Random Projections' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::rp.test(y,k = 16)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

test_that("El Bouch' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::elbouch.test(y)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a gamma ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})

set.seed(169721)
y = rgamma(200, shape = 2)
x = rgamma(200, shape = 3)

test_that("El Bouch' test works for non GPs", {
  start = Sys.time()
  ht = nortsTest::elbouch.test(y, x)
  end = difftime(Sys.time(), start, units="secs")
  # check the test choose the right hypothesis when using a Gaussian ARMA
  expect_equal(unname(ht$p.value < 0.05), TRUE)
  #check computations are less than 2.5s
  expect_equal(end < 10, TRUE)
})
