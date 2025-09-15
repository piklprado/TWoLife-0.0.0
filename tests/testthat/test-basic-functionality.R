test_that("Basic simulation runs without error", {
  # Create simple test landscape
  habitat <- matrix(rep(1, 25), nrow = 5, ncol = 5)
  
  landscape_params <- list(
    cells_per_side = 5,
    cell_size = 1.0,
    boundary_condition = 1
  )
  
  individual_params <- list(
    initial_population_size = 10,
    neighbor_radius = 2.0,
    step_length = 1.0,
    base_birth_rate = 0.5,
    base_mortality_rate = 0.1
  )
  
  simulation_params <- list(
    max_events = 100,
    neutral_mode = FALSE
  )
  
  # Test that simulation runs
  expect_no_error({
    result <- twolife_simulation(
      landscape_params, 
      individual_params, 
      simulation_params, 
      habitat
    )
  })
})

test_that("Simulation returns expected structure", {
  habitat <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  
  landscape_params <- list(cells_per_side = 3, cell_size = 1.0)
  individual_params <- list(initial_population_size = 5)
  simulation_params <- list(max_events = 50)
  
  result <- twolife_simulation(landscape_params, individual_params, simulation_params, habitat)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("final_population_size" %in% names(result))
  expect_true("survivor_x" %in% names(result))
  expect_true("survivor_y" %in% names(result))
  expect_true("event_times" %in% names(result))
  expect_true("event_types" %in% names(result))
  
  # Check types
  expect_true(is.numeric(result$final_population_size))
  expect_true(is.numeric(result$survivor_x))
  expect_true(is.numeric(result$event_times))
})

test_that("Parameter validation works", {
  habitat <- matrix(rep(1, 4), nrow = 2, ncol = 2)
  
  # Test invalid parameters
  expect_error({
    twolife_simulation("not_a_list", list(), list(), habitat)
  })
  
  expect_error({
    twolife_simulation(
      list(cells_per_side = -1), 
      list(initial_population_size = 10), 
      list(max_events = 100), 
      habitat
    )
  })
})