test_that("output has expected fields", {
  res <- run_twolife_simulation(steps = 1, n = 2)
  expect_true(all(c("population","habitat","steps","n") %in% names(res)))
})
