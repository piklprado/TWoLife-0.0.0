# functions.R
# Purpose: Complete TWoLife R interface for individual-based simulations with landscape
# generation, simulation execution, analysis, and visualization.

#' Run TWoLife Individual-Based Simulation
#' 
#' @description 
#' Execute a spatially-explicit individual-based population simulation where organisms
#' move, reproduce, and die across heterogeneous landscapes. Supports genetic evolution,
#' habitat selection behavior, and density-dependent demography.
#' 
#' @param landscape_params List containing landscape parameters:
#'   \describe{
#'     \item{\code{habitat}}{Square matrix of habitat quality values (required)}
#'     \item{\code{cell_size}}{Size of each cell in world units (default: 1.0)}
#'     \item{\code{boundary_condition}}{Edge behavior: 0=absorbing, 1=periodic, 2=reflective (default: 1)}
#'     \item{\code{density_type}}{Density dependence: 0=global, 1=local (default: 1)}
#'     \item{\code{matrix_mortality_multiplier}}{Death rate multiplier in matrix habitat (default: 2.0)}
#'     \item{\code{matrix_dispersal_multiplier}}{Movement rate multiplier in matrix habitat (default: 0.5)}
#'   }
#' @param individual_params List containing individual parameters:
#'   \describe{
#'     \item{\code{initial_population_size}}{Starting number of individuals (default: 200)}
#'     \item{\code{neighbor_radius}}{Radius for density interactions (default: 2.0)}
#'     \item{\code{vision_angle}}{Angular range sampled during random walk in radians (default: pi)}
#'     \item{\code{step_length}}{Maximum distance moved per dispersal event (default: 1.0)}
#'     \item{\code{base_dispersal_rate}}{Base rate of movement events (default: 0.1)}
#'     \item{\code{base_birth_rate}}{Base rate of reproduction (default: 0.3)}
#'     \item{\code{base_mortality_rate}}{Base rate of death (default: 0.2)}
#'     \item{\code{birth_density_slope}}{Density effect on birth rate (default: 0.02)}
#'     \item{\code{mortality_density_slope}}{Density effect on death rate (default: 0.02)}
#'   }
#' @param genetic_params List containing genetic parameters:
#'   \describe{
#'     \item{\code{genotype_means}}{Optimal habitat values for each individual (default: all 1.0)}
#'     \item{\code{genotype_sds}}{Niche width (habitat tolerance) for each individual (default: all 0.0)}
#'     \item{\code{mutation_rates}}{Mutation rate per generation for each individual (default: all 0.0)}
#'     \item{\code{plasticities}}{Phenotypic flexibility for each individual (default: all 0.0)}
#'     \item{\code{sampling_points}}{Number of locations sampled for habitat selection (default: all 0 = random walk)}
#'   }
#' @param simulation_params List containing simulation control parameters:
#'   \describe{
#'     \item{\code{max_events}}{Maximum number of demographic events (auto-calculated if NULL)}
#'     \item{\code{neutral_mode}}{All individuals act with fitness of the average initial population (default: FALSE)}
#'   }
#' @param ... Additional named arguments (deprecated - use parameter lists instead)
#' @param output_file Optional file path for detailed event output (advanced usage)
#' 
#' @return A list of class 'twolife_result' containing:
#'   \describe{
#'     \item{\code{summary}}{List with final_population_size, total_events, duration, status}
#'     \item{\code{survivors}}{Data frame with id, x, y, genotype for surviving individuals}
#'     \item{\code{spatial}}{List with world_size and num_patches}
#'     \item{\code{events}}{List with complete event history (times, types, coordinates, genotypes)}
#'     \item{\code{parameters}}{All parameters used in the simulation}
#'   }
#' 
#' @details
#' TWoLife simulates individual-based population dynamics where each organism has:
#' \itemize{
#'   \item Spatial location with movement behavior
#'   \item Genetic traits determining habitat preferences  
#'   \item Demographic rates affected by local density
#'   \item Optional behavioral habitat selection
#' }
#' 
#' The simulation proceeds as a continuous-time Markov process where individuals
#' experience birth, death, or dispersal events with rates determined by their
#' genetic traits, local environment, and nearby population density.
#' 
#' @section Parameter Guidelines:
#' \itemize{
#'   \item Ensure \code{base_birth_rate > base_mortality_rate} for population viability
#'   \item Use \code{step_length <= neighbor_radius} for realistic spatial interactions  
#'   \item Set \code{sampling_points > 0} to enable habitat selection behavior
#'   \item Match genetic parameter vector lengths to \code{initial_population_size}
#' }
#' 
#' @examples
#' # Basic simulation with default parameters
#' habitat <- create_fractal_landscape(8, fractality = 0.6, habitat_proportion = 0.4)
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = habitat),
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 50)
#' )
#' print(result)
#' 
#' \donttest{
#' # Advanced simulation with genetic evolution
#' result_genetic <- twolife_simulation(
#'   landscape_params = list(
#'     habitat = habitat,
#'     boundary_condition = 1,  # Periodic boundaries
#'     density_type = 1         # Local density dependence
#'   ),
#'   individual_params = list(
#'     initial_population_size = 20,
#'     base_birth_rate = 0.35,
#'     base_mortality_rate = 0.25
#'   ),
#'   genetic_params = list(
#'     genotype_means = runif(20, 0.3, 0.7),    # Diverse starting genotypes
#'     genotype_sds = rep(0.15, 20),            # Moderate specialists
#'     mutation_rates = rep(0.03, 20),          # Moderate evolution
#'     plasticities = rep(0.02, 20),            # Low plasticity
#'     sampling_points = rep(5, 20)             # Habitat selection
#'   ),
#'   simulation_params = list(max_events = 200)
#' )
#' 
#' # Visualize results
#' plot_simulation_on_landscape(result_genetic)
#' }
#' 
#' @seealso 
#' \code{\link{create_fractal_landscape}} for landscape generation,
#' \code{\link{plot_simulation_on_landscape}} for result visualization,
#' \code{\link{compute_population_size}} for trajectory analysis
#' 
#' @export
twolife_simulation <- function(landscape_params = list(), 
                               individual_params = list(), 
                               genetic_params = list(),
                               simulation_params = list(), 
                               ..., 
                               output_file = NULL) {
  
  # Check for deprecated habitat_grid argument
  dots <- list(...)
  if ("habitat_grid" %in% names(dots)) {
    stop("habitat_grid argument is deprecated. Please pass habitat inside landscape_params$habitat", call. = FALSE)
  }
  
  # Validate required landscape_params
  if (!is.list(landscape_params)) {
    stop("landscape_params must be a list", call. = FALSE)
  }
  
  if (is.null(landscape_params$habitat)) {
    stop("habitat matrix is required in landscape_params$habitat", call. = FALSE)
  }
  
  # Extract and validate habitat matrix
  habitat_grid <- landscape_params$habitat
  
  if (!is.matrix(habitat_grid)) {
    stop("landscape_params$habitat must be a matrix", call. = FALSE)
  }
  
  if (nrow(habitat_grid) != ncol(habitat_grid)) {
    stop("landscape_params$habitat must be square", call. = FALSE)
  }
  
  if (nrow(habitat_grid) <= 0) {
    stop("landscape_params$habitat must have positive dimensions", call. = FALSE)
  }
  
  if (any(!is.finite(habitat_grid))) {
    stop("landscape_params$habitat contains non-finite values", call. = FALSE)
  }
  
  # TWoLife canonical defaults
  defaults <- list(
    # Spatial parameters
    cell_size = 1.0,
    boundary_condition = 1,                    # 0=absorbing, 1=periodic, 2=reflective
    density_type = 1,                         # 0=global, 1=local
    matrix_mortality_multiplier = 2.0,
    matrix_dispersal_multiplier = 0.5,
    
    # Movement parameters
    neighbor_radius = 2.0,
    vision_angle = pi,                        # radians
    step_length = 1.0,
    base_dispersal_rate = 0.1,
    
    # Demography parameters
    base_birth_rate = 0.3,
    base_mortality_rate = 0.20,
    birth_density_slope = 0.02,
    mortality_density_slope = 0.02,
    
    # Population parameters
    initial_population_size = 200,
    initial_placement_mode = 1,
    
    # Simulation parameters
    neutral_mode = FALSE
  )
  
  # Combine parameters with defaults
  all_params <- modifyList(defaults, c(
    landscape_params[names(landscape_params) != "habitat"], 
    individual_params, 
    simulation_params
  ))
  
  # Calculate dynamic max_events based on population size
  pop_size <- all_params$initial_population_size
  if (is.null(simulation_params$max_events) && !("max_events" %in% names(dots))) {
    all_params$max_events <- (50 * pop_size)
  } else if (is.null(all_params$max_events)) {
    all_params$max_events <- (50 * pop_size)
  }
  
  # Parameter constraints validation - use stop() for critical errors
  if (all_params$base_birth_rate <= all_params$base_mortality_rate) {
    stop("Invalid parameters: base_birth_rate (", all_params$base_birth_rate, 
         ") must be > base_mortality_rate (", all_params$base_mortality_rate, ")", 
         call. = FALSE)
  }
  
  if (all_params$step_length > all_params$neighbor_radius) {
    warning("Constraint violation: step_length (", all_params$step_length, 
            ") should be <= neighbor_radius (", all_params$neighbor_radius, ")", 
            call. = FALSE)
  }
  
  if (all_params$matrix_dispersal_multiplier <= 0) {
    stop("Constraint violation: matrix_dispersal_multiplier must be > 0, got: ", 
         all_params$matrix_dispersal_multiplier, call. = FALSE)
  }
  
  if (all_params$vision_angle <= 0 || all_params$vision_angle > 2 * pi) {
    warning("Constraint violation: vision_angle (", all_params$vision_angle, 
            ") should be in (0, 2*pi]", call. = FALSE)
  }
  
  # Process genetic parameters with improved handling
  process_genetic_parameter <- function(param_value, param_name, pop_size, default_value) {
    if (is.null(param_value)) {
      return(rep(default_value, pop_size))
    } else if (length(param_value) == 1) {
      return(rep(param_value, pop_size))
    } else if (length(param_value) == pop_size) {
      return(param_value)
    } else {
      stop("Parameter '", param_name, "' must have length 1 or ", pop_size, 
           ", got length ", length(param_value), call. = FALSE)
    }
  }
  
  genotype_means <- process_genetic_parameter(
    genetic_params$genotype_means, "genotype_means", pop_size, 1
  )
  
  genotype_sds <- process_genetic_parameter(
    genetic_params$genotype_sds, "genotype_sds", pop_size, 0
  )
  
  mutation_rates <- process_genetic_parameter(
    genetic_params$mutation_rates, "mutation_rates", pop_size, 0
  )
  
  plasticities <- process_genetic_parameter(
    genetic_params$plasticities, "plasticities", pop_size, 0
  )
  
  # Special handling for integer parameter
  sampling_points_raw <- genetic_params$sampling_points
  if (is.null(sampling_points_raw)) {
    sampling_points <- rep(0L, pop_size)
  } else if (length(sampling_points_raw) == 1) {
    sampling_points <- rep(as.integer(sampling_points_raw), pop_size)
  } else if (length(sampling_points_raw) == pop_size) {
    sampling_points <- as.integer(sampling_points_raw)
  } else {
    stop("Parameter 'sampling_points' must have length 1 or ", pop_size, 
         ", got length ", length(sampling_points_raw), call. = FALSE)
  }
  
  # Process initial coordinates
  if (all_params$initial_placement_mode == 3) {
    if (is.null(individual_params$initial_x_coordinates) || 
        is.null(individual_params$initial_y_coordinates)) {
      stop("initial_x_coordinates and initial_y_coordinates required for placement mode 3", call. = FALSE)
    }
    initial_x <- individual_params$initial_x_coordinates
    initial_y <- individual_params$initial_y_coordinates
  } else {
    initial_x <- rep(0, pop_size)
    initial_y <- rep(0, pop_size)
  }
  
  # Execute C++ simulation via internal function
  result <- run_twolife_simulation(
    neighbor_radius = all_params$neighbor_radius,
    initial_population_size = as.integer(all_params$initial_population_size),
    vision_angle = all_params$vision_angle,
    step_length = all_params$step_length,
    base_dispersal_rate = all_params$base_dispersal_rate,
    base_birth_rate = all_params$base_birth_rate,
    base_mortality_rate = all_params$base_mortality_rate,
    birth_density_slope = all_params$birth_density_slope,
    mortality_density_slope = all_params$mortality_density_slope,
    habitat = habitat_grid,
    cell_size = all_params$cell_size,
    density_type = as.integer(all_params$density_type),
    matrix_mortality_multiplier = all_params$matrix_mortality_multiplier,
    matrix_dispersal_multiplier = all_params$matrix_dispersal_multiplier,
    initial_placement_mode = as.integer(all_params$initial_placement_mode),
    boundary_condition = as.integer(all_params$boundary_condition),
    max_events = all_params$max_events,
    initial_x_coordinates = initial_x,
    initial_y_coordinates = initial_y,
    genotype_means = genotype_means,
    genotype_sds = genotype_sds,
    mutation_rates = mutation_rates,
    plasticities = plasticities,
    sampling_points = sampling_points,
    neutral_mode = all_params$neutral_mode,
    output_file = output_file
  )
  
  # Compute essential summary stats
  final_pop <- as.integer(length(result$survivor_x))
  total_events <- length(result$event_times)
  duration <- if(total_events > 0) max(result$event_times) else 0
  
  # Lean result structure
  lean_result <- list(
    # Essential summary
    summary = list(
      final_population_size = final_pop,
      total_events = total_events,
      duration = duration,
      status = if(final_pop > 0) "surviving" else "extinct"
    ),
    
    # Survivors
    survivors = if(final_pop > 0) {
      data.frame(
        id = result$survivor_ids,
        x = result$survivor_x,
        y = result$survivor_y,
        genotype = result$survivor_genotypes,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(id = integer(0), x = numeric(0), y = numeric(0), 
                 genotype = numeric(0), stringsAsFactors = FALSE)
    },
    
    # Spatial info
    spatial = list(
      world_size = result$world_size,
      num_patches = result$num_patches
    ),
    
    # Raw events
    events = list(
      times = result$event_times,
      types = result$event_types,
      individual_ids = result$individual_ids,
      patch_ids = result$patch_ids,
      x_coordinates = result$x_coordinates,
      y_coordinates = result$y_coordinates,
      genotypes = result$genotypes
    ),
    
    # Parameters used
    parameters = list(
      landscape = list(
        habitat = habitat_grid,
        cell_size = all_params$cell_size,
        boundary_condition = all_params$boundary_condition,
        density_type = all_params$density_type,
        matrix_mortality_multiplier = all_params$matrix_mortality_multiplier,
        matrix_dispersal_multiplier = all_params$matrix_dispersal_multiplier
      ),
      individual = list(
        initial_population_size = all_params$initial_population_size,
        neighbor_radius = all_params$neighbor_radius,
        vision_angle = all_params$vision_angle,
        step_length = all_params$step_length,
        base_dispersal_rate = all_params$base_dispersal_rate,
        base_birth_rate = all_params$base_birth_rate,
        base_mortality_rate = all_params$base_mortality_rate,
        birth_density_slope = all_params$birth_density_slope,
        mortality_density_slope = all_params$mortality_density_slope,
        initial_placement_mode = all_params$initial_placement_mode,
        initial_x_coordinates = if(all_params$initial_placement_mode == 3) initial_x else NULL,
        initial_y_coordinates = if(all_params$initial_placement_mode == 3) initial_y else NULL
      ),
      genetic = list(
        genotype_means = genotype_means,
        genotype_sds = genotype_sds,
        mutation_rates = mutation_rates,
        plasticities = plasticities,
        sampling_points = sampling_points
      ),
      simulation = list(
        max_events = all_params$max_events,
        neutral_mode = all_params$neutral_mode
      )
    )
  )
  
  # S3 class for methods
  class(lean_result) <- c("twolife_result", "list")
  
  return(lean_result)
}

#' Calculate Population Trajectory Over Time
#' 
#' @description
#' Extract and compute population size changes throughout a TWoLife simulation.
#' This function processes the event history to create a time series of 
#' population size, useful for analyzing population dynamics and persistence.
#' 
#' @param result Result object from \code{\link{twolife_simulation}}
#' 
#' @return Data frame with columns:
#'   \describe{
#'     \item{\code{time}}{Simulation time points}
#'     \item{\code{population_size}}{Number of individuals at each time point}
#'   }
#' 
#' @details
#' This function reconstructs the population trajectory by processing all
#' demographic events (births, deaths, dispersal, emigration) in chronological
#' order. The resulting trajectory shows how population size changed over
#' the course of the simulation.
#' 
#' Event types processed:
#' \describe{
#'   \item{-1}{Initial placement (+1 individual)}
#'   \item{0}{Death (-1 individual)} 
#'   \item{1}{Birth (+1 individual)}
#'   \item{2}{Dispersal (no change in population)}
#'   \item{3}{Emigration (-1 individual)}
#' }
#' 
#' @examples
#' # Run simulation
#' habitat <- create_fractal_landscape(8, fractality = 0.5, habitat_proportion = 0.3)
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = habitat),
#'   individual_params = list(initial_population_size = 15),
#'   simulation_params = list(max_events = 100)
#' )
#' 
#' # Extract trajectory
#' trajectory <- compute_population_size(result)
#' head(trajectory)
#' 
#' # Plot population over time
#' plot(trajectory$time, trajectory$population_size, 
#'      type = "l", main = "Population Trajectory",
#'      xlab = "Time", ylab = "Population Size")
#' 
#' \donttest{
#' # Analyze trajectory characteristics
#' max_pop <- max(trajectory$population_size)
#' final_pop <- tail(trajectory$population_size, 1)
#' mean_pop <- mean(trajectory$population_size)
#' 
#' cat("Peak population:", max_pop, "\n")
#' cat("Final population:", final_pop, "\n")
#' cat("Mean population:", round(mean_pop, 1), "\n")
#' 
#' # Check for population crashes
#' min_pop <- min(trajectory$population_size)
#' if (min_pop == 0) {
#'   extinction_time <- trajectory$time[trajectory$population_size == 0][1]
#'   cat("Population went extinct at time:", extinction_time, "\n")
#' }
#' }
#' 
#' @seealso 
#' \code{\link{twolife_simulation}} for generating simulation results,
#' \code{\link{plot_simulation_on_landscape}} for spatial visualization
#' 
#' @export
compute_population_size <- function(result) {
  
  # Extract event data
  events_df <- data.frame(
    time = result$events$times,
    event_type = result$events$types,
    stringsAsFactors = FALSE
  )
  
  # Ensure chronological ordering
  events_df <- events_df[order(events_df$time), ]
  
  # Map event types to population changes
  events_df$pop_change <- ifelse(events_df$event_type == -1, 1,    # Initial placement
                                 ifelse(events_df$event_type == 0, -1,     # Death
                                        ifelse(events_df$event_type == 1, 1,      # Birth  
                                               ifelse(events_df$event_type == 2, 0,      # Dispersal (no net change)
                                                      ifelse(events_df$event_type == 3, -1, 0))))) # Emigration
  
  # Calculate running population total
  events_df$population_size <- cumsum(events_df$pop_change)
  
  return(events_df[, c("time", "population_size")])
}

#' Create Fractal Landscape
#' 
#' @description
#' Generate realistic heterogeneous landscapes using fractal algorithms.
#' Creates spatial patterns with controllable fragmentation and habitat proportion,
#' suitable for testing spatial population models.
#' 
#' @param cells_per_side Integer. Number of cells per side of the square landscape (must be >= 2)
#' @param fractality Numeric. Fractal parameter controlling spatial correlation (0-1).
#'   Higher values (0.7-0.9) create smoother, more connected patterns.
#'   Lower values (0.1-0.4) create more fragmented, random-like patterns.
#' @param min_value Numeric. Minimum habitat value for continuous landscapes (default: 0.0)
#' @param max_value Numeric. Maximum habitat value for continuous landscapes (default: 1.0)
#' @param habitat_proportion Numeric. Proportion of cells classified as habitat (0-1).
#'   If specified, creates binary landscape (0/1). If NULL, creates continuous landscape.
#' @param return_as_landscape_params Logical. If TRUE, returns list(habitat = matrix)
#'   suitable for direct use in \code{\link{twolife_simulation}}. If FALSE, returns matrix only (default: FALSE).
#' 
#' @return Square matrix representing habitat grid, or list with habitat component if 
#'   \code{return_as_landscape_params = TRUE}
#' 
#' @details
#' This function uses the diamond-square fractal algorithm to generate spatially
#' autocorrelated landscapes. The fractality parameter controls the degree of
#' spatial correlation:
#' 
#' \describe{
#'   \item{High fractality (0.7-0.9)}{Smooth gradients, large connected patches}
#'   \item{Medium fractality (0.4-0.7)}{Moderate fragmentation, intermediate patches}  
#'   \item{Low fractality (0.1-0.4)}{High fragmentation, small scattered patches}
#' }
#' 
#' For binary landscapes, \code{habitat_proportion} determines the fraction of cells
#' with value 1 (habitat) vs 0 (matrix). The threshold is applied to the
#' continuous fractal pattern, preserving spatial structure.
#' 
#' @section Dependencies:
#' Requires the 'rflsgen' package for fractal generation. If not installed, 
#' install with: \code{install.packages("rflsgen")}
#' 
#' @examples
#' # Binary fragmented landscape (agricultural mosaic)
#' fragmented <- create_fractal_landscape(
#'   cells_per_side = 10,
#'   fractality = 0.3,        # High fragmentation
#'   habitat_proportion = 0.4  # 40% habitat
#' )
#' 
#' # Binary connected landscape (protected area)
#' connected <- create_fractal_landscape(
#'   cells_per_side = 10,
#'   fractality = 0.8,        # Low fragmentation  
#'   habitat_proportion = 0.4  # 40% habitat
#' )
#' 
#' \donttest{
#' # Continuous environmental gradient
#' gradient <- create_fractal_landscape(
#'   cells_per_side = 20,
#'   fractality = 0.7,
#'   habitat_proportion = NULL  # Continuous values
#' )
#' 
#' # Ready for simulation use
#' landscape_params <- create_fractal_landscape(
#'   cells_per_side = 15,
#'   fractality = 0.6,
#'   habitat_proportion = 0.3,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Use directly in simulation
#' result <- twolife_simulation(
#'   landscape_params = landscape_params,
#'   individual_params = list(initial_population_size = 20),
#'   simulation_params = list(max_events = 100)
#' )
#' }
#' 
#' @seealso 
#' \code{\link{create_corner_landscape}} for simpler test landscapes,
#' \code{\link{plot_landscape_world_coords}} for visualization,
#' \code{\link{twolife_simulation}} for using landscapes in simulations
#' 
#' @export
create_fractal_landscape <- function(cells_per_side, 
                                     fractality, 
                                     min_value = 0.0, 
                                     max_value = 1.0,
                                     habitat_proportion = NULL,
                                     return_as_landscape_params = FALSE) {
  
  # INPUT VALIDATION
  if (!is.numeric(cells_per_side) || cells_per_side <= 0 || cells_per_side != round(cells_per_side)) {
    stop("cells_per_side must be a positive integer", call. = FALSE)
  }
  
  if (!is.numeric(fractality) || fractality < 0 || fractality > 1) {
    stop("fractality must be between 0 and 1", call. = FALSE)
  }
  
  if (!is.numeric(min_value) || !is.numeric(max_value)) {
    stop("min_value and max_value must be numeric", call. = FALSE)
  }
  
  if (min_value >= max_value) {
    stop("min_value must be less than max_value", call. = FALSE)
  }
  
  if (!is.null(habitat_proportion)) {
    if (!is.numeric(habitat_proportion) || habitat_proportion < 0 || habitat_proportion > 1) {
      stop("habitat_proportion must be between 0 and 1", call. = FALSE)
    }
  }
  
  # DEPENDENCY CHECK
  if (!requireNamespace("rflsgen", quietly = TRUE)) {
    stop("Package 'rflsgen' is required for fractal landscapes. Install with: install.packages('rflsgen')", 
         call. = FALSE)
  }
  
  # FRACTAL GENERATION
  range01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  raw_fractal <- rflsgen::flsgen_terrain(
    cells_per_side,    
    cells_per_side,    
    fractality         
  )
  
  fractal_matrix <- as.matrix(raw_fractal)
  normalized_fractal <- range01(fractal_matrix)
  
  # DETERMINE OUTPUT TYPE
  if (!is.null(habitat_proportion)) {
    # BINARY LANDSCAPE with edge case handling
    
    # Handle edge cases first
    if (habitat_proportion >= 1.0) {
      # All habitat
      binary_landscape <- matrix(1, nrow = cells_per_side, ncol = cells_per_side)
    } else if (habitat_proportion <= 0.0) {
      # No habitat (all matrix)
      binary_landscape <- matrix(0, nrow = cells_per_side, ncol = cells_per_side)
    } else {
      # Normal case: use quantile-based thresholding
      threshold <- quantile(normalized_fractal, 1 - habitat_proportion, na.rm = TRUE)
      binary_landscape <- matrix(0, nrow = cells_per_side, ncol = cells_per_side)
      # Use >= instead of > to be more inclusive at boundaries
      binary_landscape[normalized_fractal >= threshold] <- 1
    }
    
    if (return_as_landscape_params) {
      return(list(habitat = binary_landscape))
    } else {
      return(binary_landscape)
    }
    
  } else {
    # CONTINUOUS LANDSCAPE
    scaled_fractal <- min_value + (normalized_fractal * (max_value - min_value))
    
    habitat_grid <- matrix(nrow = cells_per_side, ncol = cells_per_side)
    count <- 0
    for (i in 1:cells_per_side) {
      for (j in 1:cells_per_side) {
        count <- count + 1
        habitat_grid[i, j] <- scaled_fractal[count]
      }
    }
    
    if (return_as_landscape_params) {
      return(list(habitat = habitat_grid))
    } else {
      return(habitat_grid)
    }
  }
}

#' Create Corner Test Landscapes
#' 
#' @description 
#' Creates test landscapes with habitat (value 1) concentrated in one of the four corners,
#' with the rest of the landscape as matrix (value 0). Useful for testing spatial effects
#' and population dynamics under different habitat configurations.
#' 
#' @param cells_per_side Integer. Number of cells per side of the square landscape (must be >= 4)
#' @param corner Character. Which corner to place habitat: "top-left", "top-right", "bottom-left", "bottom-right"
#' @param corner_size Integer. Size of the habitat corner in cells (default: 1/4 of cells_per_side)
#' @param return_as_landscape_params Logical. If TRUE, returns list(habitat = matrix) suitable for 
#'   \code{\link{twolife_simulation}}. If FALSE, returns matrix only (default: FALSE)
#' 
#' @return Square matrix representing habitat grid with habitat in specified corner, or list with 
#'   habitat component if \code{return_as_landscape_params = TRUE}
#' 
#' @details
#' Corner landscapes provide simple, controlled test cases for understanding
#' how habitat placement affects population dynamics. They are particularly
#' useful for:
#' \itemize{
#'   \item Testing edge effects and boundary conditions
#'   \item Comparing dispersal strategies and movement patterns
#'   \item Understanding source-sink dynamics
#'   \item Validating simulation behavior in extreme habitat configurations
#' }
#' 
#' The habitat corner is always square-shaped and positioned at the specified corner
#' of the landscape. All other cells are set to 0 (matrix habitat).
#' 
#' @examples
#' # Basic corner landscape
#' habitat_tl <- create_corner_landscape(12, corner = "top-left")
#' 
#' # Different corners for comparison
#' corners <- c("top-left", "top-right", "bottom-left", "bottom-right")
#' landscapes <- lapply(corners, function(c) {
#'   create_corner_landscape(10, corner = c, corner_size = 3)
#' })
#' names(landscapes) <- corners
#' 
#' \donttest{
#' # Test population dynamics with corner habitat
#' habitat_corner <- create_corner_landscape(16, corner = "top-left", corner_size = 4)
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = habitat_corner),
#'   individual_params = list(
#'     initial_population_size = 20,
#'     base_dispersal_rate = 0.2  # Higher dispersal for corner landscapes
#'   ),
#'   simulation_params = list(max_events = 200)
#' )
#' 
#' # Compare all four corners
#' corner_results <- list()
#' for(corner in corners) {
#'   habitat <- create_corner_landscape(12, corner = corner)
#'   corner_results[[corner]] <- twolife_simulation(
#'     landscape_params = list(habitat = habitat),
#'     individual_params = list(initial_population_size = 15),
#'     simulation_params = list(max_events = 150)
#'   )
#' }
#' 
#' # Compare final population sizes
#' sapply(corner_results, function(r) r$summary$final_population_size)
#' }
#' 
#' @seealso 
#' \code{\link{create_fractal_landscape}} for realistic landscapes,
#' \code{\link{plot_landscape_world_coords}} for visualization,
#' \code{\link{twolife_simulation}} for running simulations
#' 
#' @export
create_corner_landscape <- function(cells_per_side, 
                                    corner = "top-left",
                                    corner_size = NULL,
                                    return_as_landscape_params = FALSE) {
  
  # INPUT VALIDATION
  if (!is.numeric(cells_per_side) || cells_per_side < 4 || cells_per_side != round(cells_per_side)) {
    stop("cells_per_side must be an integer >= 4", call. = FALSE)
  }
  
  valid_corners <- c("top-left", "top-right", "bottom-left", "bottom-right")
  if (!corner %in% valid_corners) {
    stop("corner must be one of: ", paste(valid_corners, collapse = ", "), call. = FALSE)
  }
  
  # Set default corner size (1/4 of landscape)
  if (is.null(corner_size)) {
    corner_size <- max(2, round(cells_per_side / 4))
  }
  
  if (!is.numeric(corner_size) || corner_size < 1 || corner_size >= cells_per_side) {
    stop("corner_size must be between 1 and cells_per_side-1", call. = FALSE)
  }
  
  # Initialize landscape with all matrix (0)
  landscape <- matrix(0, nrow = cells_per_side, ncol = cells_per_side)
  
  # Define corner coordinates based on corner choice
  corner_coords <- switch(corner,
                          "top-left" = list(
                            rows = (cells_per_side - corner_size + 1):cells_per_side,
                            cols = 1:corner_size
                          ),
                          "top-right" = list(
                            rows = (cells_per_side - corner_size + 1):cells_per_side,
                            cols = (cells_per_side - corner_size + 1):cells_per_side
                          ),
                          "bottom-left" = list(
                            rows = 1:corner_size,
                            cols = 1:corner_size
                          ),
                          "bottom-right" = list(
                            rows = 1:corner_size,
                            cols = (cells_per_side - corner_size + 1):cells_per_side
                          )
  )
  
  # Set habitat in specified corner
  landscape[corner_coords$rows, corner_coords$cols] <- 1
  
  # Return format
  if (return_as_landscape_params) {
    return(list(habitat = landscape))
  } else {
    return(landscape)
  }
}

#' Plot Landscape with World Coordinates
#' 
#' @description 
#' Enhanced landscape plotting that uses world coordinates matching the C++ simulation
#' and applies consistent matrix transformation for proper orientation.
#' 
#' @param landscape_data Matrix or list(habitat = matrix)
#' @param cell_size Numeric. Size of each cell (default: 1.0)
#' @param filename Character. File path for export (optional)
#' @param main Character. Plot title
#' @param colors Character. Color scheme: "terrain", "viridis", "plasma", "habitat", "quality", "binary"
#' @param show_legend Logical. Whether to show color scale legend
#' @param add_grid Logical. Whether to add grid lines and reference axes
#' 
#' @return Invisibly returns filename if saved, NULL otherwise
#' 
#' @examples
#' habitat <- create_fractal_landscape(20, fractality = 0.6, habitat_proportion = 0.4)
#' plot_landscape_world_coords(habitat, cell_size = 1.0, colors = "habitat")
#' 
#' @seealso 
#' \code{\link{plot_simulation_on_landscape}} for overlaying simulation results
#' 
#' @export
plot_landscape_world_coords <- function(landscape_data, 
                                        cell_size = 1.0,
                                        filename = NULL,
                                        main = "Landscape (World Coordinates)",
                                        colors = "habitat",
                                        show_legend = TRUE,
                                        add_grid = TRUE) {
  
  # Extract habitat grid
  if (is.list(landscape_data) && "habitat" %in% names(landscape_data)) {
    habitat_grid <- landscape_data$habitat
  } else if (is.matrix(landscape_data)) {
    habitat_grid <- landscape_data
  } else {
    stop("landscape_data must be a matrix or list with $habitat component", call. = FALSE)
  }
  
  if (!is.matrix(habitat_grid)) {
    stop("habitat_grid must be a matrix", call. = FALSE)
  }
  
  # Geometry - FIXED coordinate dimensions
  n <- nrow(habitat_grid)
  world_side_length <- n * cell_size
  x_coords <- seq(-world_side_length/2, world_side_length/2, length.out = n)
  y_coords <- seq(-world_side_length/2, world_side_length/2, length.out = n)
  
  # Apply same transformation as plot_simulation_on_landscape
  z <- t(apply(habitat_grid, 2, rev))
  
  # FIXED: Handle uniform landscapes and color inversion
  is_binary <- all(habitat_grid %in% c(0, 1))
  z_range <- range(z, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  
  # Color palette with correct mapping
  if (colors == "habitat") {
    if (is_binary) {
      # FIXED: Correct color order - 0=beige (matrix), 1=green (habitat)
      color_palette <- c("#F5F5DC", "#228B22")  # Beige for 0, Green for 1
      
      # FIXED: Handle uniform binary landscapes
      if (is_uniform) {
        if (z_range[1] == 0) {
          # All matrix - show beige
          color_palette <- "#F5F5DC"
        } else {
          # All habitat - show green  
          color_palette <- "#228B22"
        }
      }
    } else {
      # Continuous landscape
      color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
    }
  } else {
    # Other color schemes
    n_colors <- if (is_binary && !is_uniform) 2 else if (is_uniform) 1 else 100
    color_palette <- switch(colors,
                            "terrain" = terrain.colors(n_colors),
                            "viridis" = {
                              if (requireNamespace("viridisLite", quietly = TRUE)) {
                                viridisLite::viridis(n_colors)
                              } else {
                                terrain.colors(n_colors)
                              }
                            },
                            terrain.colors(n_colors)
    )
  }
  
  # FIXED: Handle plotting for uniform vs non-uniform landscapes
  if (is_uniform) {
    # For uniform landscapes, create a simple filled plot
    image(x = x_coords, y = y_coords, z = z,
          col = color_palette,
          main = main,
          xlab = "X (World Coordinates)",
          ylab = "Y (World Coordinates)",
          asp = 1,
          useRaster = TRUE)
  } else {
    # For variable landscapes, use standard image plotting
    image(x = x_coords, y = y_coords, z = z, 
          col = color_palette,
          main = main,
          xlab = "X (World Coordinates)",
          ylab = "Y (World Coordinates)",
          asp = 1,
          useRaster = TRUE)
  }
  
  # Add grid and reference lines
  if (add_grid) {
    grid(col = "white", lty = 1, lwd = 0.5)
    abline(h = 0, v = 0, col = "red", lty = 2, lwd = 2)
  }
  
  invisible(NULL)
}

#' Plot Simulation Results on Landscape
#' 
#' @description 
#' Visualize TWoLife simulation results by overlaying individual positions
#' on the landscape used in the simulation. Shows spatial distribution of
#' surviving individuals and their genetic characteristics.
#' 
#' @param simulation_result A twolife_result object from \code{\link{twolife_simulation}}
#' @param point_size Numeric. Size of points representing individuals (default: 2)
#' @param point_color Character. Color of survivor points (default: "red")
#' @param point_shape Numeric. Shape of points using R's pch values (default: 16)
#' @param show_genotypes Logical. Color points by genotype values (default: FALSE)
#' @param landscape_colors Character. Color scheme for landscape background: 
#'   "habitat", "terrain", "viridis" (default: "habitat")
#' @param main Character. Plot title (auto-generated if NULL)
#' @param filename Character. Optional filename to save plot
#' @param add_stats Logical. Add text with population statistics (default: TRUE)
#' 
#' @return Invisibly returns the simulation result object
#' 
#' @details
#' This function creates a spatial visualization showing:
#' \itemize{
#'   \item Landscape habitat quality as background colors
#'   \item Surviving individual positions as points
#'   \item Optional genetic diversity through point colors
#'   \item Population statistics in plot title or text
#' }
#' 
#' When \code{show_genotypes = TRUE}, points are colored by their genotype values
#' using a heat color scale, allowing visualization of spatial genetic structure.
#' 
#' The landscape is displayed using world coordinates matching the simulation,
#' with habitat quality shown through the selected color scheme.
#' 
#' @examples
#' # Run basic simulation
#' habitat <- create_fractal_landscape(10, fractality = 0.6, habitat_proportion = 0.4)
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = habitat),
#'   individual_params = list(initial_population_size = 20),
#'   simulation_params = list(max_events = 150)
#' )
#' 
#' # Basic visualization
#' plot_simulation_on_landscape(result)
#' 
#' # Customized visualization
#' plot_simulation_on_landscape(result, 
#'                              point_size = 3,
#'                              point_color = "blue",
#'                              landscape_colors = "terrain",
#'                              main = "Custom Simulation Results")
#' 
#' \donttest{
#' # Simulation with genetic diversity
#' result_genetic <- twolife_simulation(
#'   landscape_params = list(habitat = habitat),
#'   individual_params = list(initial_population_size = 25),
#'   genetic_params = list(
#'     genotype_means = runif(25, 0.2, 0.8),
#'     genotype_sds = rep(0.15, 25),
#'     mutation_rates = rep(0.03, 25),
#'     sampling_points = rep(10, 25)
#'   ),
#'   simulation_params = list(max_events = 200)
#' )
#' 
#' # Show genetic structure
#' plot_simulation_on_landscape(result_genetic, 
#'                              show_genotypes = TRUE,
#'                              main = "Spatial Genetic Structure")
#' 
#' # Compare multiple scenarios
#' par(mfrow = c(1, 2))
#' plot_simulation_on_landscape(result, main = "No Genetics")
#' plot_simulation_on_landscape(result_genetic, show_genotypes = TRUE, 
#'                              main = "With Genetics")
#' par(mfrow = c(1, 1))
#' }
#' 
#' @seealso 
#' \code{\link{twolife_simulation}} for generating results,
#' \code{\link{plot_landscape_world_coords}} for plotting landscapes alone,
#' \code{\link{compute_population_size}} for trajectory analysis
#' 
#' @export
plot_simulation_on_landscape <- function(simulation_result,
                                         point_size = 2,
                                         point_color = "red",
                                         point_shape = 16,
                                         show_genotypes = FALSE,
                                         landscape_colors = "habitat",
                                         main = NULL,
                                         filename = NULL,
                                         add_stats = TRUE) {
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  # Extract data
  habitat_grid <- simulation_result$parameters$landscape$habitat
  cell_size <- simulation_result$parameters$landscape$cell_size
  survivors <- simulation_result$survivors
  n_survivors <- if (is.null(survivors)) 0L else nrow(survivors)
  if (is.null(main)) main <- paste("Simulation Results:", n_survivors, "Survivors")
  
  # FIXED: Geometry - same coordinate system as plot_landscape_world_coords
  n <- nrow(habitat_grid)
  world_side_length <- n * cell_size
  x_coords <- seq(-world_side_length/2, world_side_length/2, length.out = n)
  y_coords <- seq(-world_side_length/2, world_side_length/2, length.out = n)
  
  # Apply same transformation as plot_landscape_world_coords
  z <- t(apply(habitat_grid, 2, rev))
  
  # FIXED: Handle uniform landscapes and color inversion
  is_binary <- all(habitat_grid %in% c(0, 1))
  z_range <- range(z, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  
  # FIXED: Color palette with correct mapping using landscape_colors parameter
  if (landscape_colors == "habitat") {
    if (is_binary) {
      # FIXED: Correct color order - 0=beige (matrix), 1=green (habitat)
      color_palette <- c("#F5F5DC", "#228B22")  # Beige for 0, Green for 1
      
      # FIXED: Handle uniform binary landscapes
      if (is_uniform) {
        if (z_range[1] == 0) {
          # All matrix - show beige
          color_palette <- "#F5F5DC"
        } else {
          # All habitat - show green  
          color_palette <- "#228B22"
        }
      }
    } else {
      # Continuous landscape
      color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
    }
  } else {
    # Other color schemes
    n_colors <- if (is_binary && !is_uniform) 2 else if (is_uniform) 1 else 100
    color_palette <- switch(landscape_colors,
                            "terrain" = terrain.colors(n_colors),
                            "viridis" = {
                              if (requireNamespace("viridisLite", quietly = TRUE)) {
                                viridisLite::viridis(n_colors)
                              } else {
                                terrain.colors(n_colors)
                              }
                            },
                            terrain.colors(n_colors)
    )
  }
  
  # FIXED: Handle plotting for uniform vs non-uniform landscapes with proper coordinates
  if (is_uniform) {
    # For uniform landscapes, create a simple filled plot
    image(x = x_coords, y = y_coords, z = z,
          col = color_palette,
          main = main,
          xlab = "X (World Coordinates)",
          ylab = "Y (World Coordinates)",
          asp = 1,
          useRaster = TRUE)
  } else {
    # For variable landscapes, use standard image plotting
    image(x = x_coords, y = y_coords, z = z, 
          col = color_palette,
          main = main,
          xlab = "X (World Coordinates)",
          ylab = "Y (World Coordinates)",
          asp = 1,
          useRaster = TRUE)
  }
  
  # Add grid and reference lines
  grid(col = "white", lty = 1, lwd = 0.5)
  abline(h = 0, v = 0, col = "red", lty = 2, lwd = 1)
  
  # Overlay survivors with proper coordinate system
  if (n_survivors > 0) {
    if (show_genotypes && !is.null(survivors$genotype)) {
      # Color by genotype
      gpal <- heat.colors(100)
      gr <- range(survivors$genotype, na.rm = TRUE)
      if (gr[1] != gr[2]) {
        gscaled <- (survivors$genotype - gr[1])/(gr[2]-gr[1])
        idx <- pmax(1, pmin(100, round(gscaled*99)+1))
        ptcols <- gpal[idx]
      } else {
        ptcols <- rep(gpal[50], n_survivors)
      }
      points(survivors$x, survivors$y, col = ptcols, pch = point_shape, cex = point_size)
    } else {
      points(survivors$x, survivors$y, col = point_color, pch = point_shape, cex = point_size)
    }
  }
  
  invisible(simulation_result)
}

#' Quick Plot of Simulation Results
#' 
#' @description 
#' Simplified function for quickly visualizing TWoLife simulation results.
#' 
#' @param simulation_result A twolife_result object
#' @param ... Additional arguments passed to \code{\link{plot_simulation_on_landscape}}
#' 
#' @examples
#' result <- twolife_simulation(...)  # your simulation
#' quick_plot_result(result)
#' 
#' @seealso 
#' \code{\link{plot_simulation_on_landscape}} for full customization
#' 
#' @export
quick_plot_result <- function(simulation_result, ...) {
  plot_simulation_on_landscape(simulation_result, ...)
}

#' Print Method for TWoLife Results
#' 
#' @description
#' Custom print method for twolife_result objects.
#' 
#' @param x A twolife_result object
#' @param ... Additional arguments (unused)
#' 
#' @return Invisibly returns the input object
#' 
#' @examples
#' result <- twolife_simulation(...)
#' print(result)  # or just: result
#' 
#' @export
print.twolife_result <- function(x, ...) {
  cat("TWoLife Simulation Result\n")
  cat("========================\n")
  cat("Status:", x$summary$status, "\n")
  cat("Final population:", x$summary$final_population_size, "\n") 
  cat("Duration:", round(x$summary$duration, 2), "\n")
  cat("Total events:", x$summary$total_events, "\n")
  cat("World size:", x$spatial$world_size, "×", x$spatial$world_size, "\n")
  if (x$summary$final_population_size > 0) {
    cat("Survivors: Use result$survivors for details\n")
  }
  cat("\nUse compute_population_size() for trajectory analysis.\n")
  cat("Use plot_simulation_on_landscape() for visualization.\n")
  invisible(x)
}