test_that("parameters are coerced to integer", {
  res <- run_twolife_simulation(steps = 5.9, n = 3.1)
  expect_equal(res$steps, 5L) # C++ receives as integer
  expect_equal(res$n, 3L)
})
