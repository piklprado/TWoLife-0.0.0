# functions.R - Complete TWoLife R interface
# UPDATED: Added all missing @param documentation and proper imports
# UPDATED: Added show_legend parameter to check_habitat_match functions
# FIXED: Replaced non-ASCII characters with ASCII equivalents

#' @importFrom grDevices colorRampPalette heat.colors terrain.colors
#' @importFrom graphics abline grid image legend mtext par points
#' @importFrom stats cor median quantile rnorm runif sd
#' @importFrom utils modifyList
NULL

# ============================================================================
# CORE SIMULATION FUNCTIONS
# ============================================================================

#' Run TWoLife Individual-Based Simulation
#' 
#' Runs a spatially-explicit individual-based simulation with habitat selection,
#' genetic variation, phenotypic plasticity, and demographic processes. Supports 
#' rectangular landscapes and customizable habitat selection parameters.
#' 
#' @param landscape_params List containing landscape parameters. Required component:
#'   \describe{
#'     \item{habitat}{Matrix. Binary (0/1) or continuous habitat values. Required.}
#'     \item{cell_size}{Numeric. Size of each landscape cell. Default: 1.0}
#'     \item{boundary_condition}{Integer. Boundary behavior: 1 = reflective (individuals bounce 
#'       back into landscape), 2 = absorbing (individuals leave and are removed, generating 
#'       emigration events), 3 = periodic (individuals wrap to opposite edge, torus topology). 
#'       Default: 1}
#'     \item{density_type}{Integer. Density calculation: 1 = local (within neighbor_radius), 
#'       2 = global (entire population). Default: 1}
#'     \item{matrix_mortality_multiplier}{Numeric. Mortality multiplier for matrix cells (habitat=0). 
#'       Values > 1 increase mortality in matrix. Default: 2.0}
#'     \item{matrix_dispersal_multiplier}{Numeric. Dispersal multiplier for matrix cells. 
#'       Values < 1 reduce movement in matrix (movement costs). Default: 0.5}
#'   }
#' @param individual_params List containing individual-level parameters:
#'   \describe{
#'     \item{initial_population_size}{Integer. Starting population size. Default: 200}
#'     \item{neighbor_radius}{Numeric. Radius for density calculations (world units). Default: 2.0}
#'     \item{vision_angle}{Numeric. Vision angle in radians for movement decisions. 
#'       pi (180 degrees) = forward hemisphere, 2*pi (360 degrees) = random walk. Default: pi}
#'     \item{step_length}{Numeric. Distance moved per dispersal event (world units). Default: 1.0}
#'     \item{base_dispersal_rate}{Numeric. Base dispersal probability per time unit (0-1). Default: 0.1}
#'     \item{base_birth_rate}{Numeric. Base birth probability per time unit (0-1). 
#'       Should typically be > base_mortality_rate for population persistence. Default: 0.3}
#'     \item{base_mortality_rate}{Numeric. Base mortality probability per time unit (0-1). Default: 0.20}
#'     \item{birth_density_slope}{Numeric. Density-dependence slope for births (negative effect). Default: 0.02}
#'     \item{mortality_density_slope}{Numeric. Density-dependence slope for mortality (positive effect). Default: 0.02}
#'     \item{initial_placement_mode}{Integer. Placement strategy: 1 = random in habitat cells, 
#'       2 = random anywhere with normal distribution around center, 3 = custom coordinates. Default: 1}
#'     \item{initial_x_coordinates}{Numeric vector. Custom x positions (required if initial_placement_mode = 3)}
#'     \item{initial_y_coordinates}{Numeric vector. Custom y positions (required if initial_placement_mode = 3)}
#'   }
#' @param genetic_params List containing genetic parameters (each can be single value or vector matching initial_population_size):
#'   \describe{
#'     \item{genotype_means}{Numeric. Mean genotype values (underlying genetic optimum). Default: 1}
#'     \item{genotype_sds}{Numeric. Standard deviations (niche width - tolerance to habitat mismatch). 
#'       Smaller values = specialist, larger values = generalist. Default: 0}
#'     \item{mutation_rates}{Numeric. Per-birth mutation rate (SD of normal distribution added to 
#'       offspring genotype). Higher values = more genetic variation between generations. Default: 0}
#'     \item{plasticities}{Numeric. Phenotypic plasticity (SD of normal distribution added to 
#'       individual phenotype from genotype). Higher values = more phenotypic variation from environment. 
#'       Works similarly to mutation_rates but affects phenotype expression rather than inheritance. Default: 0}
#'     \item{sampling_points}{Integer. Number of habitat locations sampled for plasticity calculation. 
#'       Higher values = more accurate environmental assessment but slower. Only relevant if plasticities > 0. 
#'       Default: 0}
#'     \item{habitat_selection_temperatures}{Numeric. Temperature parameter for softmax habitat selection 
#'       (must be positive). Controls selection strength during dispersal: Lower values (< 1) = strong 
#'       selection (nearly always choose best habitat), higher values (> 1) = weak selection (more 
#'       exploration/randomness), 1 = balanced (probability of choosing habitat proportional to its 
#'       relative fitness). Default: 1.0}
#'   }
#' @param simulation_params List containing simulation control parameters:
#'   \describe{
#'     \item{max_events}{Integer. Maximum number of events to simulate. Default: 50 * initial_population_size}
#'     \item{neutral_mode}{Logical. If TRUE, disables habitat selection (all locations equally preferred, 
#'       and all individuals are assigned the average genotype of the initial population). Useful for 
#'       null model comparisons to test effects of habitat selection. Default: FALSE}
#'   }
#' @param history_detail Character. Level of event history detail to record.
#'   Options: "minimal" (only time, event type, individual ID - fastest, smallest memory),
#'   "standard" (adds spatial coordinates, patch ID, genotype - enables most analyses),
#'   "full" (adds phenotype and niche width - enables exact historical reconstruction).
#'   Default: "standard"
#' @param master_seed Integer. Seed for reproducible simulations. If NULL, results are stochastic. Default: NULL
#' 
#' @return A list of class 'twolife_result' with components:
#'   \describe{
#'     \item{summary}{List with status ("surviving" or "extinct"), final_population_size (integer), 
#'       total_events (integer), and duration (numeric time units)}
#'     \item{survivors}{Data frame with columns: id (integer), x (numeric), y (numeric), 
#'       genotype (numeric), phenotype (numeric), width (numeric). Empty data frame if extinct.}
#'     \item{spatial}{List containing world_width (numeric), world_height (numeric), 
#'       world_size (numeric, maximum dimension: max(world_width, world_height)), num_patches (integer)}
#'     \item{events}{List with event history. Content depends on history_detail:
#'       \itemize{
#'         \item Always included: times (numeric vector), types (integer vector: -1=initial, 0=death, 
#'           1=birth, 2=movement, 3=emigration), individual_ids (integer vector)
#'         \item If history_detail is "standard" or "full": patch_ids, x_coordinates, y_coordinates, genotypes
#'         \item If history_detail == "full": phenotypes, widths
#'       }}
#'     \item{parameters}{Nested list preserving all input parameters (landscape, individual, genetic, simulation)}
#'   }
#' 
#' @details
#' Simulation Stopping Conditions:
#'   The simulation automatically stops when any of these conditions is met:
#'   \itemize{
#'     \item Population reaches 0 (extinction)
#'     \item Population exceeds 1,000,000 individuals (overflow protection)
#'     \item max_events is reached
#'   }
#'   Users can interrupt long simulations with Ctrl+C or Esc. The simulation checks 
#'   for user interrupts every 1000 events.
#' 
#' Genetic Architecture:
#'   The model implements a quantitative genetics framework:
#'   \itemize{
#'     \item Genotype: Underlying genetic optimum (heritable, passed to offspring with mutation)
#'     \item Mutation: Offspring genotype = parent genotype + N(0, mutation_rate)
#'     \item Phenotype: Expressed optimum = genotype + N(0, plasticity), where N(0, plasticity) 
#'       is random normal variation representing environmental influence
#'     \item Niche width (genotype_sds): Determines fitness specialism via Gaussian fitness function. 
#'       Smaller values create specialists (fitness drops sharply away from optimum), larger values 
#'       create generalists (fitness remains high across wider habitat range).
#'   }
#'   Higher fitness increases birth rate and decreases mortality rate relative to habitat quality.
#'   Mutation affects inheritance (parent to offspring), while plasticity affects individual 
#'   phenotype expression within a generation.
#' 
#' Boundary Conditions:
#'   When individuals move beyond landscape edges:
#'   \itemize{
#'     \item Reflective (1): Individual bounces back, remains in simulation
#'     \item Absorbing (2): Individual exits permanently (emigration event, type = 3)
#'     \item Periodic (3): Individual wraps to opposite edge (torus topology)
#'   }
#' 
#' Event Types:
#'   Event history records these event types:
#'   \itemize{
#'     \item -1: Initial placement
#'     \item 0: Death (natural or density-dependent)
#'     \item 1: Birth (reproduction)
#'     \item 2: Movement (dispersal within landscape)
#'     \item 3: Emigration (boundary crossing with absorbing boundaries)
#'   }
#' 
#' @examples
#' # Quick example for CRAN checks (< 2 seconds)
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 123
#' )
#' 
#' print(result)
#' head(result$survivors)
#' 
#' # With genetic variation
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   genetic_params = list(
#'     genotype_means = rnorm(15, mean = 0.5, sd = 0.15),
#'     genotype_sds = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 456
#' )
#' 
#' \donttest{
#' # Larger examples (not run during CRAN checks)
#' landscape_full <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.7,
#'   habitat_proportion = 0.4,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Longer simulation for detailed analysis
#' result_long <- twolife_simulation(
#'   landscape_params = landscape_full,
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 5000),
#'   history_detail = "full"
#' )
#' 
#' # Analyze population trajectory
#' pop_trajectory <- population_size(result_long)
#' plot(pop_trajectory$time, pop_trajectory$population_size, type = "l")
#' }
#' 
#' @export
twolife_simulation <- function(landscape_params = list(), 
                               individual_params = list(), 
                               genetic_params = list(),
                               simulation_params = list(),
                               history_detail = "standard",
                               master_seed = NULL) {
  
  # Validate history_detail parameter
  valid_levels <- c("minimal", "standard", "full")
  if (!history_detail %in% valid_levels) {
    stop("history_detail must be one of: ", 
         paste(valid_levels, collapse = ", "), 
         call. = FALSE)
  }
  
  if (!is.list(landscape_params)) {
    stop("landscape_params must be a list", call. = FALSE)
  }
  
  if (is.null(landscape_params$habitat)) {
    stop("habitat matrix is required in landscape_params$habitat", call. = FALSE)
  }
  
  habitat_grid <- landscape_params$habitat
  
  if (!is.matrix(habitat_grid)) {
    stop("landscape_params$habitat must be a matrix", call. = FALSE)
  }
  
  if (nrow(habitat_grid) <= 0 || ncol(habitat_grid) <= 0) {
    stop("landscape_params$habitat must have positive dimensions", call. = FALSE)
  }
  
  if (any(!is.finite(habitat_grid))) {
    stop("landscape_params$habitat contains non-finite values", call. = FALSE)
  }
  
  if (!is.null(master_seed)) {
    if (!is.numeric(master_seed) || length(master_seed) != 1 || master_seed != round(master_seed)) {
      stop("master_seed must be a single integer", call. = FALSE)
    }
  }
  
  defaults <- list(
    cell_size = 1.0,
    boundary_condition = 1,
    density_type = 1,
    matrix_mortality_multiplier = 2.0,
    matrix_dispersal_multiplier = 0.5,
    neighbor_radius = 2.0,
    vision_angle = pi,
    step_length = 1.0,
    base_dispersal_rate = 0.1,
    base_birth_rate = 0.3,
    base_mortality_rate = 0.20,
    birth_density_slope = 0.02,
    mortality_density_slope = 0.02,
    initial_population_size = 200,
    initial_placement_mode = 1,
    neutral_mode = FALSE
  )
  
  all_params <- modifyList(defaults, c(
    landscape_params[names(landscape_params) != "habitat"], 
    individual_params, 
    simulation_params
  ))
  
  pop_size <- all_params$initial_population_size
  if (is.null(simulation_params$max_events)) {
    all_params$max_events <- (50 * pop_size)
  } else if (is.null(all_params$max_events)) {
    all_params$max_events <- (50 * pop_size)
  }
  
  if (all_params$base_birth_rate <= all_params$base_mortality_rate) {
    stop("Invalid parameters: base_birth_rate must be > base_mortality_rate", call. = FALSE)
  }
  
  if (all_params$step_length > all_params$neighbor_radius) {
    warning("step_length should be <= neighbor_radius", call. = FALSE)
  }
  
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
  
  habitat_selection_temperatures <- process_genetic_parameter(
    genetic_params$habitat_selection_temperatures, "habitat_selection_temperatures", pop_size, 1.0
  )
  
  if (any(habitat_selection_temperatures <= 0)) {
    stop("All habitat_selection_temperatures must be positive", call. = FALSE)
  }
  
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
    habitat_selection_temperatures = habitat_selection_temperatures,
    neutral_mode = all_params$neutral_mode,
    history_detail = history_detail,
    master_seed = master_seed
  )
  
  final_pop <- as.integer(length(result$survivor_x))
  total_events <- length(result$event_times)
  duration <- if(total_events > 0) max(result$event_times) else 0
  
  # Build events list conditionally based on history_detail
  events_list <- list(
    times = result$event_times,
    types = result$event_types,
    individual_ids = result$individual_ids
  )
  
  if (history_detail != "minimal") {
    events_list$patch_ids <- result$patch_ids
    events_list$x_coordinates <- result$x_coordinates
    events_list$y_coordinates <- result$y_coordinates
    events_list$genotypes <- result$genotypes
  }
  
  if (history_detail == "full") {
    events_list$phenotypes <- result$phenotypes
    events_list$widths <- result$widths
  }
  
  lean_result <- list(
    summary = list(
      final_population_size = final_pop,
      total_events = total_events,
      duration = duration,
      status = if(final_pop > 0) "surviving" else "extinct"
    ),
    survivors = if(final_pop > 0) {
      data.frame(
        id = result$survivor_ids,
        x = result$survivor_x,
        y = result$survivor_y,
        genotype = result$survivor_genotypes,
        phenotype = result$survivor_phenotypes,
        width = result$survivor_widths,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(id = integer(0), x = numeric(0), y = numeric(0), 
                 genotype = numeric(0), phenotype = numeric(0), 
                 width = numeric(0), stringsAsFactors = FALSE)
    },
    spatial = list(
      world_width = result$world_width,
      world_height = result$world_height,
      world_size = max(result$world_width, result$world_height),
      num_patches = result$num_patches
    ),
    events = events_list,
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
        sampling_points = sampling_points,
        habitat_selection_temperatures = habitat_selection_temperatures
      ),
      simulation = list(
        max_events = all_params$max_events,
        neutral_mode = all_params$neutral_mode,
        history_detail = history_detail,
        master_seed = master_seed
      )
    )
  )
  
  class(lean_result) <- c("twolife_result", "list")
  return(lean_result)
}

#' Calculate Population Trajectory Over Time
#' 
#' Computes population size at each event time from simulation results.
#' Works with all history_detail levels (only requires event times and types).
#' Useful for analyzing population dynamics, detecting extinctions, and 
#' identifying critical time points.
#' 
#' @param result A 'twolife_result' object returned by \code{\link{twolife_simulation}}.
#'   Must contain \code{$events$times} and \code{$events$types} components.
#' 
#' @return A data frame with two columns:
#'   \describe{
#'     \item{time}{Numeric. Event times in chronological order}
#'     \item{population_size}{Integer. Population size at each event time (cumulative sum of births and deaths)}
#'   }
#' 
#' @details
#' Event types and their effects on population size:
#' \itemize{
#'   \item -1 (initial): +1 to population
#'   \item 0 (death): -1 from population
#'   \item 1 (birth): +1 to population
#'   \item 2 (movement): no change to population
#'   \item 3 (emigration): -1 from population (only occurs with absorbing boundaries)
#' }
#' 
#' Note: Emigration events (type 3) only occur when boundary_condition = 2 (absorbing).
#' With reflective boundaries (1), individuals bounce back. With periodic boundaries (3),
#' individuals wrap around.
#' 
#' @examples
#' # Quick example
#' set.seed(300)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Run simulation
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 789
#' )
#' 
#' # Calculate trajectory
#' trajectory <- population_size(result)
#' head(trajectory)
#' 
#' # Plot population over time
#' plot(trajectory$time, trajectory$population_size, 
#'      type = "l", 
#'      xlab = "Time", 
#'      ylab = "Population Size",
#'      main = "Population Dynamics")
#' 
#' # Find peak population
#' peak_time <- trajectory$time[which.max(trajectory$population_size)]
#' peak_size <- max(trajectory$population_size)
#' cat("Peak population:", peak_size, "at time", peak_time, "\n")
#' 
#' @seealso \code{\link{twolife_simulation}} for running simulations,
#'   \code{\link{snapshot_at_time}} for reconstructing population state at specific times
#' 
#' @export
population_size <- function(result) {
  events_df <- data.frame(
    time = result$events$times,
    event_type = result$events$types,
    stringsAsFactors = FALSE
  )
  
  events_df <- events_df[order(events_df$time), ]
  
  events_df$pop_change <- ifelse(events_df$event_type == -1, 1,
                                 ifelse(events_df$event_type == 0, -1,
                                        ifelse(events_df$event_type == 1, 1,
                                               ifelse(events_df$event_type == 2, 0,
                                                      ifelse(events_df$event_type == 3, -1, 0)))))
  
  events_df$population_size <- cumsum(events_df$pop_change)
  
  return(events_df[, c("time", "population_size")])
}

#' Reconstruct Population State at Specific Time
#' 
#' Replays the event log from a simulation to reconstruct the exact positions
#' and states of all living individuals at a specified time point. Useful for
#' temporal analysis, creating animations, and studying population dynamics
#' at critical moments.
#' 
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}.
#'   Must have been created with \code{history_detail = "standard"} or \code{"full"}.
#' @param target_time Numeric. Time point to reconstruct population state. 
#'   Must be >= 0 and <= maximum simulation time. If larger than max time,
#'   automatically uses final time with a warning.
#' @param color_by Character. How to color points in visualization:
#'   \itemize{
#'     \item "genotype" - Color by genotype values
#'     \item "phenotype" - Color by phenotype values (requires history_detail = "full")
#'     \item "none" - Single color (red)
#'   }
#'   Default: "genotype"
#' @param show_plot Logical. If TRUE, displays visualization of population at target_time. Default: TRUE
#' @param point_size Numeric. Size of points in plot. Default: 2
#' 
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{time}{Numeric. The target_time used (may differ from input if capped at max time)}
#'     \item{n_alive}{Integer. Number of living individuals at target_time}
#'     \item{population}{Data frame with columns: id (integer), x (numeric), y (numeric), 
#'       genotype (numeric), phenotype (numeric, NA if history_detail != "full"), 
#'       width (numeric, NA if history_detail != "full"), patch_id (integer)}
#'     \item{events_processed}{Integer. Number of events processed up to target_time}
#'     \item{history_level}{Character. The history_detail level from the simulation 
#'       ("minimal", "standard", or "full")}
#'   }
#' 
#' @details
#' This function requires simulations run with \code{history_detail = "standard"} or 
#' \code{"full"}. The "minimal" history level does not record spatial coordinates,
#' making reconstruction impossible. The function will check and provide a clear
#' error message if the result doesn't support snapshots.
#' 
#' The reconstruction works by:
#' \enumerate{
#'   \item Filtering events up to target_time
#'   \item Tracking births and deaths to determine who is alive
#'   \item Recording the most recent position for each living individual
#'   \item Optionally visualizing the spatial distribution
#' }
#' 
#' @examples
#' # Quick example
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.5,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Run with standard history
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 100),
#'   history_detail = "standard",
#'   master_seed = 123
#' )
#' 
#' # Get trajectory to find interesting time
#' traj <- population_size(result)
#' mid_time <- traj$time[nrow(traj) %/% 2]
#' 
#' # Reconstruct population at midpoint
#' state <- snapshot_at_time(result, target_time = mid_time, show_plot = FALSE)
#' 
#' # Examine the snapshot
#' cat("Population at time", state$time, ":", state$n_alive, "individuals\n")
#' head(state$population)
#' 
#' \donttest{
#' # Larger example with visualization
#' landscape_big <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.7,
#'   habitat_proportion = 0.4,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' result_big <- twolife_simulation(
#'   landscape_params = landscape_big,
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 150),
#'   history_detail = "standard",
#'   master_seed = 456
#' )
#' 
#' traj_big <- population_size(result_big)
#' mid_time_big <- traj_big$time[nrow(traj_big) %/% 2]
#' 
#' # With plotting
#' snapshot_at_time(result_big, target_time = mid_time_big)
#' 
#' # Create animation sequence
#' times <- seq(0, max(traj_big$time), length.out = 10)
#' for (t in times) {
#'   snapshot_at_time(result_big, target_time = t, 
#'                   color_by = "genotype", point_size = 3)
#'   Sys.sleep(0.5)  # Pause between frames
#' }
#' }
#' 
#' @seealso \code{\link{population_size}} for population trajectories,
#'   \code{\link{twolife_simulation}} for running simulations with appropriate history_detail
#' 
#' @export
snapshot_at_time <- function(simulation_result, 
                             target_time,
                             color_by = "genotype",
                             show_plot = TRUE,
                             point_size = 2) {
  
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  # Check history detail level
  history_level <- simulation_result$parameters$simulation$history_detail
  
  if (is.null(history_level)) {
    # Backward compatibility - infer from available fields
    if ("phenotypes" %in% names(simulation_result$events)) {
      history_level <- "full"
    } else if ("x_coordinates" %in% names(simulation_result$events)) {
      history_level <- "standard"
    } else {
      history_level <- "minimal"
    }
  }
  
  # Validate capability
  if (history_level == "minimal") {
    stop("Cannot reconstruct spatial positions with history_detail='minimal'.\n",
         "Re-run simulation with history_detail='standard' or 'full'.", 
         call. = FALSE)
  }
  
  # Validate color_by parameter
  if (!color_by %in% c("genotype", "phenotype", "none")) {
    stop("color_by must be one of: 'genotype', 'phenotype', or 'none'", call. = FALSE)
  }
  
  # Extract available fields based on history level
  events <- data.frame(
    time = simulation_result$events$times,
    type = simulation_result$events$types,
    id = simulation_result$events$individual_ids,
    stringsAsFactors = FALSE
  )
  
  if (history_level != "minimal") {
    events$patch_id <- simulation_result$events$patch_ids
    events$x <- simulation_result$events$x_coordinates
    events$y <- simulation_result$events$y_coordinates
    events$genotype <- simulation_result$events$genotypes
  }
  
  use_exact_phenotypes <- FALSE
  if (history_level == "full") {
    events$phenotype <- simulation_result$events$phenotypes
    events$width <- simulation_result$events$widths
    use_exact_phenotypes <- TRUE
  }
  
  # Validate target_time
  if (target_time < 0) {
    stop("target_time must be non-negative", call. = FALSE)
  }
  
  if (target_time > max(events$time)) {
    warning("target_time exceeds simulation duration. Using final time.", call. = FALSE)
    target_time <- max(events$time)
  }
  
  # Filter events up to target time
  events_up_to <- events[events$time <= target_time, ]
  
  if (nrow(events_up_to) == 0) {
    stop("No events occurred before target_time", call. = FALSE)
  }
  
  # Track population state
  alive_ids <- c()
  alive_data <- list()
  
  for (i in seq_len(nrow(events_up_to))) {
    event <- events_up_to[i, ]
    
    if (event$type == -1) {
      # Initial state
      alive_ids <- c(alive_ids, event$id)
      alive_data[[as.character(event$id)]] <- event
      
    } else if (event$type == 1) {
      # Birth
      alive_ids <- c(alive_ids, event$id)
      alive_data[[as.character(event$id)]] <- event
      
    } else if (event$type == 0 || event$type == 3) {
      # Death or emigration
      alive_ids <- setdiff(alive_ids, event$id)
      alive_data[[as.character(event$id)]] <- NULL
      
    } else if (event$type == 2) {
      # Movement - update position
      if (as.character(event$id) %in% names(alive_data)) {
        alive_data[[as.character(event$id)]]$x <- event$x
        alive_data[[as.character(event$id)]]$y <- event$y
        alive_data[[as.character(event$id)]]$patch_id <- event$patch_id
      }
    }
  }
  
  # Convert to data frame
  if (length(alive_ids) == 0) {
    alive_pop <- data.frame(
      id = integer(0),
      x = numeric(0),
      y = numeric(0),
      genotype = numeric(0),
      phenotype = numeric(0),
      width = numeric(0),
      patch_id = integer(0)
    )
  } else {
    alive_pop <- do.call(rbind, lapply(alive_data, function(row) {
      data.frame(
        id = row$id,
        x = row$x,
        y = row$y,
        genotype = row$genotype,
        phenotype = if(use_exact_phenotypes) row$phenotype else NA,
        width = if(use_exact_phenotypes) row$width else NA,
        patch_id = row$patch_id,
        stringsAsFactors = FALSE
      )
    }))
    rownames(alive_pop) <- NULL
    
    # Calculate approximate phenotypes if not using exact values
    if (!use_exact_phenotypes && nrow(alive_pop) > 0) {
      plasticities <- simulation_result$parameters$genetic$plasticities
      widths <- simulation_result$parameters$genetic$genotype_sds
      
      for (i in seq_len(nrow(alive_pop))) {
        ind_id <- alive_pop$id[i]
        
        if (ind_id <= length(plasticities)) {
          plast <- plasticities[ind_id]
        } else {
          plast <- median(plasticities)
        }
        
        alive_pop$phenotype[i] <- alive_pop$genotype[i] + rnorm(1, 0, plast)
        
        if (ind_id <= length(widths)) {
          alive_pop$width[i] <- widths[ind_id]
        } else {
          alive_pop$width[i] <- median(widths)
        }
      }
      
      # Override with exact values for final survivors
      if (!is.null(simulation_result$survivors) && nrow(simulation_result$survivors) > 0) {
        for (i in seq_len(nrow(alive_pop))) {
          survivor_match <- simulation_result$survivors$id == alive_pop$id[i]
          if (any(survivor_match)) {
            alive_pop$phenotype[i] <- simulation_result$survivors$phenotype[survivor_match][1]
            alive_pop$width[i] <- simulation_result$survivors$width[survivor_match][1]
          }
        }
      }
    }
  }
  
  # Visualization
  if (show_plot && nrow(alive_pop) > 0) {
    habitat_grid <- simulation_result$parameters$landscape$habitat
    cell_size <- simulation_result$parameters$landscape$cell_size
    
    n_rows <- nrow(habitat_grid)
    n_cols <- ncol(habitat_grid)
    world_width <- n_cols * cell_size
    world_height <- n_rows * cell_size
    
    x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
    y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
    z <- t(apply(habitat_grid, 2, rev))
    
    # Determine if binary and set colors accordingly
    is_binary <- all(habitat_grid %in% c(0, 1))
    z_range <- range(habitat_grid)
    
    if (is_binary) {
      # Binary: 0 = white (matrix), 1 = green (habitat)
      landscape_cols <- c("white", "#228B22")
    } else {
      # Continuous: terrain colors
      landscape_cols <- terrain.colors(100)
    }
    
    par(mfrow = c(1, 2))
    
    # Left panel: Landscape
    image(x = x_coords, y = y_coords, z = z,
          col = landscape_cols,
          main = "Landscape",
          xlab = "X", ylab = "Y", asp = 1)
    grid(col = "gray80", lty = 1, lwd = 0.5)
    
    # Right panel: Population with color mapping
    trait_values <- if (color_by == "phenotype") alive_pop$phenotype else alive_pop$genotype
    
    if (color_by != "none" && z_range[1] != z_range[2]) {
      # Normalize trait values to [0, 1] based on landscape range
      trait_normalized <- (trait_values - z_range[1]) / (z_range[2] - z_range[1])
      trait_normalized <- pmax(0, pmin(1, trait_normalized))
      
      # Map to the SAME color palette as the landscape
      if (is_binary) {
        # For binary: traits < 0.5 -> white, >= 0.5 -> green
        color_indices <- ifelse(trait_normalized < 0.5, 1, 2)
        point_colors <- landscape_cols[color_indices]
      } else {
        # For continuous: map to terrain.colors consistently
        color_indices <- pmax(1, pmin(100, round(trait_normalized * 99) + 1))
        point_colors <- landscape_cols[color_indices]
      }
    } else {
      point_colors <- rep("red", nrow(alive_pop))
    }
    
    trait_name <- if (color_by == "phenotype") "Phenotype" else if (color_by == "genotype") "Genotype" else "Trait"
    
    plot(alive_pop$x, alive_pop$y,
         col = point_colors, pch = 16, cex = point_size,
         main = paste("Population at t =", round(target_time, 2), 
                      if (color_by != "none") paste("\n(colored by", trait_name, ")") else ""),
         xlab = "X", ylab = "Y", asp = 1,
         xlim = c(-world_width/2, world_width/2),
         ylim = c(-world_height/2, world_height/2))
    grid(col = "gray80", lty = 1, lwd = 0.5)
    
    par(mfrow = c(1, 1))
  }
  
  return(invisible(list(
    time = target_time,
    n_alive = nrow(alive_pop),
    population = alive_pop,
    events_processed = nrow(events_up_to),
    history_level = history_level
  )))
}

# ============================================================================
# LANDSCAPE GENERATION FUNCTIONS
# ============================================================================

#' Create Fractal Landscape
#' 
#' Generates a spatially autocorrelated landscape using fractal algorithms.
#' Supports both square and rectangular landscapes, and can create either 
#' continuous or binary (habitat/matrix) patterns. Useful for testing
#' habitat selection and spatial processes.
#' 
#' @param cells_per_row Integer. Number of cells per row. Must be positive. Typical values: 10-100.
#' @param cells_per_col Integer. Number of cells per column. If NULL, creates square landscape 
#'   (cells_per_col = cells_per_row). Default: NULL
#' @param fractality Numeric between 0 and 1. Controls spatial autocorrelation:
#'   \itemize{
#'     \item 0 = completely random (no spatial structure)
#'     \item 0.3-0.5 = moderate clumping
#'     \item 0.7-0.9 = highly clumped/fragmented
#'     \item 1 = maximum spatial structure
#'   }
#' @param min_value Numeric. Minimum habitat value for continuous landscapes. Default: 0.0
#' @param max_value Numeric. Maximum habitat value for continuous landscapes. Must be > min_value. Default: 1.0
#' @param habitat_proportion Numeric between 0 and 1. If provided, creates binary landscape:
#'   \itemize{
#'     \item Proportion of cells designated as habitat (value = 1)
#'     \item Remaining cells are matrix (value = 0)
#'     \item Uses fractality to determine spatial arrangement
#'     \item If NULL, creates continuous landscape
#'   }
#'   Default: NULL
#' @param return_as_landscape_params Logical. If TRUE, returns \code{list(habitat = matrix)}
#'   suitable for direct use in \code{twolife_simulation()}. If FALSE, returns matrix only. Default: FALSE
#' 
#' @return If \code{return_as_landscape_params = FALSE}: A numeric matrix with dimensions 
#'   \code{cells_per_row} x \code{cells_per_col}. Values are either continuous (between min_value 
#'   and max_value) or binary (0 and 1).
#'   
#'   If \code{return_as_landscape_params = TRUE}: A list with one component:
#'   \describe{
#'     \item{habitat}{The generated landscape matrix}
#'   }
#' 
#' @details
#' Algorithm:
#'   The fractal generation uses iterative neighbor averaging:
#'   \enumerate{
#'     \item Start with random values at each cell
#'     \item For each iteration (controlled by fractality):
#'       \itemize{
#'         \item Average each cell with its 8 neighbors
#'         \item Add small random noise
#'       }
#'     \item More iterations (higher fractality) -> more spatial clumping
#'   }
#'   
#'   Number of smoothing iterations = round(fractality x 10). For example:
#'   \itemize{
#'     \item fractality = 0.0 -> 0 iterations (completely random)
#'     \item fractality = 0.5 -> 5 iterations (moderate structure)
#'     \item fractality = 0.9 -> 9 iterations (highly clumped)
#'   }
#' 
#' Mathematical Formula:
#'   For each smoothing iteration t, the value at cell (i,j) is updated:
#'   
#'   \deqn{value_{t+1}(i,j) = \alpha \times mean(neighbors_t(i,j)) + (1-\alpha) \times \epsilon}
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{\alpha = 0.9} (smoothing weight)
#'     \item \eqn{neighbors_t(i,j)} = values of 8 surrounding cells at iteration t
#'     \item \eqn{\epsilon \sim Uniform(0, 1)} (random noise)
#'     \item Number of iterations = round(fractality x 10)
#'   }
#'   
#'   Edge cells use available neighbors only (fewer than 8).
#' Binary vs Continuous:
#'   \itemize{
#'     \item If habitat_proportion is NULL: Creates continuous landscape with values 
#'       scaled to [min_value, max_value]
#'     \item If habitat_proportion is provided: Applies threshold to make binary 
#'       (0 = matrix, 1 = habitat), with fractality determining spatial clustering
#'   }
#' 
#' Usage in Simulations:
#'   Binary landscapes (habitat_proportion specified) are typical for most simulations,
#'   representing discrete habitat patches in a hostile matrix. Continuous landscapes
#'   represent gradual environmental gradients. The matrix_mortality_multiplier and
#'   matrix_dispersal_multiplier parameters in twolife_simulation() only affect 
#'   cells with value = 0 in binary landscapes.
#' 
#' @examples
#' # Quick examples
#' landscape1 <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.5,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Visualize
#' plot_landscape(landscape1, main = "Binary Fractal Landscape")
#' 
#' # Continuous landscape
#' landscape2 <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0,
#'   max_value = 1
#' )
#' 
#' # Check dimensions and values
#' dim(landscape2)
#' range(landscape2)
#' 
#' \donttest{
#' # Larger examples (not run during CRAN checks)
#' 
#' # Binary landscape
#' landscape3 <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.7,
#'   habitat_proportion = 0.4,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Rectangular landscape
#' landscape4 <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   cells_per_col = 20,
#'   fractality = 0.6,
#'   habitat_proportion = 0.5,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Use in simulation
#' result <- twolife_simulation(
#'   landscape_params = landscape4,
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 100),
#'   master_seed = 123
#' )
#' 
#' # Compare different fractality levels
#' par(mfrow = c(2, 2))
#' for (frac in c(0.2, 0.5, 0.7, 0.9)) {
#'   land <- create_fractal_landscape(
#'     cells_per_row = 5,
#'     fractality = frac,
#'     habitat_proportion = 0.4
#'   )
#'   plot_landscape(
#'     list(habitat = land),
#'     main = paste("Fractality =", frac)
#'   )
#' }
#' par(mfrow = c(1, 1))
#' }
#' 
#' @seealso \code{\link{plot_landscape}} for visualization
#' 
#' @export
create_fractal_landscape <- function(cells_per_row, 
                                     cells_per_col = NULL,
                                     fractality, 
                                     min_value = 0.0, 
                                     max_value = 1.0,
                                     habitat_proportion = NULL,
                                     return_as_landscape_params = FALSE) {
  
  if (!is.numeric(cells_per_row) || cells_per_row <= 0 || cells_per_row != round(cells_per_row)) {
    stop("cells_per_row must be a positive integer", call. = FALSE)
  }
  
  if (is.null(cells_per_col)) {
    cells_per_col <- cells_per_row
  }
  
  if (!is.numeric(cells_per_col) || cells_per_col <= 0 || cells_per_col != round(cells_per_col)) {
    stop("cells_per_col must be a positive integer", call. = FALSE)
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
  
  range01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  # Generate fractal pattern using iterative smoothing
  fractal_matrix <- matrix(runif(cells_per_row * cells_per_col), nrow = cells_per_row, ncol = cells_per_col)
  n_iterations <- max(1, round(fractality * 10))
  
  for (iter in 1:n_iterations) {
    smoothed <- matrix(0, nrow = cells_per_row, ncol = cells_per_col)
    
    for (i in 1:cells_per_row) {
      for (j in 1:cells_per_col) {
        neighbors <- c()
        if (i > 1) neighbors <- c(neighbors, fractal_matrix[i-1, j])
        if (i < cells_per_row) neighbors <- c(neighbors, fractal_matrix[i+1, j])
        if (j > 1) neighbors <- c(neighbors, fractal_matrix[i, j-1])
        if (j < cells_per_col) neighbors <- c(neighbors, fractal_matrix[i, j+1])
        if (i > 1 && j > 1) neighbors <- c(neighbors, fractal_matrix[i-1, j-1])
        if (i > 1 && j < cells_per_col) neighbors <- c(neighbors, fractal_matrix[i-1, j+1])
        if (i < cells_per_row && j > 1) neighbors <- c(neighbors, fractal_matrix[i+1, j-1])
        if (i < cells_per_row && j < cells_per_col) neighbors <- c(neighbors, fractal_matrix[i+1, j+1])
        
        weight_original <- 1 - fractality
        weight_neighbors <- fractality / length(neighbors)
        
        smoothed[i, j] <- weight_original * fractal_matrix[i, j] + sum(neighbors * weight_neighbors)
      }
    }
    fractal_matrix <- smoothed
  }
  
  normalized_fractal <- range01(fractal_matrix)
  
  if (!is.null(habitat_proportion)) {
    if (habitat_proportion >= 1.0) {
      binary_landscape <- matrix(1, nrow = cells_per_row, ncol = cells_per_col)
    } else if (habitat_proportion <= 0.0) {
      binary_landscape <- matrix(0, nrow = cells_per_row, ncol = cells_per_col)
    } else {
      threshold <- quantile(normalized_fractal, 1 - habitat_proportion, na.rm = TRUE)
      binary_landscape <- matrix(0, nrow = cells_per_row, ncol = cells_per_col)
      binary_landscape[normalized_fractal >= threshold] <- 1
    }
    
    if (return_as_landscape_params) {
      return(list(habitat = binary_landscape))
    } else {
      return(binary_landscape)
    }
    
  } else {
    scaled_fractal <- min_value + (normalized_fractal * (max_value - min_value))
    
    habitat_grid <- matrix(nrow = cells_per_row, ncol = cells_per_col)
    count <- 0
    for (i in 1:cells_per_row) {
      for (j in 1:cells_per_col) {
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


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

#' Plot Landscape with World Coordinates
#' 
#' Visualizes landscapes in the world coordinate system with proper axis scaling.
#' Works with both square and rectangular landscapes. Shows the landscape as it
#' appears in the simulation, centered at the origin with correct spatial scaling.
#' 
#' @param landscape_data Either:
#'   \itemize{
#'     \item A numeric matrix of habitat values
#'     \item A list with a \code{$habitat} component (output from landscape creation functions)
#'   }
#' @param cell_size Numeric. Size of each landscape cell in world units. Should match
#'   the cell_size used in simulations. Default: 1.0
#' @param filename Character. Optional file path to save plot (e.g., "landscape.png").
#'   If NULL, displays interactively. Default: NULL
#' @param main Character. Plot title. Default: "Landscape (World Coordinates)"
#' @param colors Character. Color scheme:
#'   \itemize{
#'     \item "habitat" - Green gradient (white=matrix, dark green=high quality habitat)
#'     \item "terrain" - Built-in terrain.colors palette
#'     \item "viridis" - Viridis color scale (requires viridisLite package)
#'   }
#'   Default: "habitat"
#' @param show_legend Logical. Whether to show color legend. Default: TRUE
#' @param add_grid Logical. If TRUE, adds grid lines and crosshairs at origin. Default: TRUE
#' 
#' @return NULL, invisibly. Function is called for its side effect (plotting).
#' 
#' @details
#' The function displays landscapes in world coordinates where:
#' \itemize{
#'   \item Origin (0, 0) is at the center
#'   \item X-axis ranges from -world_width/2 to +world_width/2
#'   \item Y-axis ranges from -world_height/2 to +world_height/2
#'   \item Red crosshairs mark the origin
#'   \item Grid lines aid spatial interpretation
#' }
#' 
#' @examples
#' # Quick examples
#' fractal <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.5
#' )
#' 
#' plot_landscape(fractal, main = "Fractal Landscape")
#' 
#' # Binary landscape with custom colors
#' binary <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.4
#' )
#' 
#' plot_landscape(
#'   binary,
#'   main = "Binary Landscape",
#'   colors = "terrain"
#' )
#' 
#' \donttest{
#' # Larger examples (not run during CRAN checks)
#' 
#' # Using output from landscape function
#' landscape_list <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.6,
#'   habitat_proportion = 0.5,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' plot_landscape(landscape_list)
#' 
#' # Continuous landscape
#' continuous <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.7,
#'   min_value = 0,
#'   max_value = 1
#' )
#' 
#' plot_landscape(
#'   continuous,
#'   main = "Continuous Habitat Quality",
#'   colors = "viridis"
#' )
#' 
#' # Compare rectangular landscapes
#' par(mfrow = c(1, 2))
#' 
#' rect1 <- create_fractal_landscape(15, 20, 0.7, habitat_proportion = 0.4)
#' plot_landscape(rect1, main = "15 x 20")
#' 
#' rect2 <- create_fractal_landscape(20, 15, 0.7, habitat_proportion = 0.4)
#' plot_landscape(rect2, main = "20 x 15")
#' 
#' par(mfrow = c(1, 1))
#' }
#' 
#' @seealso \code{\link{plot_simulation_on_landscape}} to overlay simulation results,
#'   \code{\link{create_fractal_landscape}} for landscape generation
#' 
#' @export
plot_landscape <- function(landscape_data, 
                           cell_size = 1.0,
                           filename = NULL,
                           main = "Landscape (World Coordinates)",
                           colors = "habitat",
                           show_legend = TRUE,
                           add_grid = TRUE) {
  
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
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
  y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
  
  z <- t(apply(habitat_grid, 2, rev))
  
  is_binary <- all(habitat_grid %in% c(0, 1))
  z_range <- range(z, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  
  if (colors == "habitat") {
    if (is_binary) {
      # Binary landscapes: 0 = white (matrix), 1 = green (habitat)
      color_palette <- c("white", "#228B22")
      if (is_uniform) {
        color_palette <- if(z_range[1] == 0) "white" else "#228B22"
      }
    } else {
      color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
    }
  } else {
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
  
  image(x = x_coords, y = y_coords, z = z, 
        col = color_palette,
        main = main,
        xlab = "X (World Coordinates)",
        ylab = "Y (World Coordinates)",
        asp = 1,
        useRaster = TRUE)
  
  if (add_grid) {
    grid(col = "gray80", lty = 1, lwd = 0.5)
    abline(h = 0, v = 0, col = "red", lty = 2, lwd = 2)
  }
  
  invisible(NULL)
}

#' Plot Simulation Results on Landscape
#' 
#' Visualizes simulation results by overlaying survivor positions on the
#' landscape. Points can be colored by genotype, phenotype, or a single color.
#' Provides immediate visual feedback on spatial distribution and habitat selection.
#' 
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param point_size Numeric. Size of survivor points. Larger values make points more visible. Default: 2
#' @param point_color Character. Color for points when \code{color_by="none"}. 
#'   Any valid R color name or hex code. Default: "red"
#' @param point_shape Integer. Point shape following R's pch codes (e.g., 16 = filled circle, 
#'   1 = open circle, 17 = filled triangle). Default: 16
#' @param color_by Character. How to color survivor points:
#'   \itemize{
#'     \item "none" - Single color specified by \code{point_color}
#'     \item "genotype" - Color gradient based on genotype values (heat colors)
#'     \item "phenotype" - Color gradient based on phenotype values (heat colors)
#'   }
#'   Default: "none"
#' @param landscape_colors Character. Color scheme for landscape:
#'   \itemize{
#'     \item "habitat" - Green gradient for continuous, white/green for binary
#'     \item "terrain" - Built-in terrain.colors palette
#'     \item "viridis" - Viridis color scale (requires viridisLite package)
#'   }
#'   Default: "habitat"
#' @param main Character. Plot title. If NULL, automatically generated. Default: NULL
#' @param filename Character. Optional file path to save plot. If NULL, displays interactively. Default: NULL
#' @param add_stats Logical. If TRUE, adds statistics text to plot (currently not implemented). Default: TRUE
#' 
#' @return The input \code{simulation_result} object, invisibly. Allows for piping/chaining operations.
#' 
#' @details
#' Coordinate System:
#'   The function displays the landscape in world coordinates where:
#'   \itemize{
#'     \item Origin (0, 0) is at the center
#'     \item Red crosshairs mark the origin
#'     \item Grid lines aid spatial interpretation
#'   }
#' 
#' Color Mapping:
#'   \itemize{
#'     \item "none": All survivors shown in a single color (point_color)
#'     \item "genotype": Heat gradient (red = low, yellow = high) shows genetic variation
#'     \item "phenotype": Heat gradient shows expressed trait values (includes plasticity effects)
#'   }
#'   
#'   When coloring by genotype/phenotype, visual clustering of similar colors in similar
#'   habitats indicates successful habitat selection or local adaptation.
#' 
#' Interpretation:
#'   \itemize{
#'     \item Survivors clustered in habitat (green/dark areas) = effective habitat preference
#'     \item Match between point colors and background colors (when using "genotype" or "phenotype") 
#'       = phenotype-environment matching (high fitness)
#'     \item Scattered distribution = weak habitat selection or high dispersal
#'     \item Survivors in matrix (white/light areas) = tolerance to poor habitat or recent dispersal
#'   }
#' 
#' @examples
#' # Quick example
#' set.seed(400)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Run simulation
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 101
#' )
#' 
#' # Basic plot
#' plot_simulation_on_landscape(result)
#' 
#' \donttest{
#' # Larger examples (not run during CRAN checks)
#' landscape_big <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.7,
#'   habitat_proportion = 0.4,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' result_big <- twolife_simulation(
#'   landscape_params = landscape_big,
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 202
#' )
#' 
#' # Customize appearance
#' plot_simulation_on_landscape(
#'   result_big,
#'   point_size = 3,
#'   point_color = "blue",
#'   point_shape = 17,
#'   main = "Final Population Distribution"
#' )
#' 
#' # Color by genotype with variation
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape_big,
#'   individual_params = list(initial_population_size = 10),
#'   genetic_params = list(
#'     genotype_means = runif(10, 0, 1),
#'     genotype_sds = 0.1
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 303
#' )
#' 
#' plot_simulation_on_landscape(
#'   result_genetic,
#'   color_by = "genotype",
#'   point_size = 2.5
#' )
#' 
#' # Color by phenotype (requires plasticity)
#' result_plastic <- twolife_simulation(
#'   landscape_params = landscape_big,
#'   individual_params = list(initial_population_size = 40),
#'   genetic_params = list(
#'     genotype_means = runif(40, 0, 1),
#'     plasticities = 0.5,
#'     sampling_points = 5
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 404
#' )
#' 
#' plot_simulation_on_landscape(
#'   result_plastic,
#'   color_by = "phenotype",
#'   landscape_colors = "terrain"
#' )
#' }
#' 
#' @seealso \code{\link{plot.twolife_result}} for S3 method (use \code{plot(result)}),
#'   \code{\link{check_habitat_match}} for habitat-trait correlation plots
#' 
#' @export
plot_simulation_on_landscape <- function(simulation_result,
                                         point_size = 2,
                                         point_color = "red",
                                         point_shape = 16,
                                         color_by = "none",
                                         landscape_colors = "habitat",
                                         main = NULL,
                                         filename = NULL,
                                         add_stats = TRUE) {
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  if (!color_by %in% c("none", "genotype", "phenotype")) {
    stop("color_by must be one of: 'none', 'genotype', or 'phenotype'", call. = FALSE)
  }
  
  habitat_grid <- simulation_result$parameters$landscape$habitat
  cell_size <- simulation_result$parameters$landscape$cell_size
  survivors <- simulation_result$survivors
  n_survivors <- if (is.null(survivors)) 0L else nrow(survivors)
  if (is.null(main)) main <- paste("Simulation Results:", n_survivors, "Survivors")
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
  y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
  
  z <- t(apply(habitat_grid, 2, rev))
  
  is_binary <- all(habitat_grid %in% c(0, 1))
  z_range <- range(z, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  
  if (landscape_colors == "habitat") {
    if (is_binary) {
      # Binary landscapes: 0 = white (matrix), 1 = green (habitat)
      color_palette <- c("white", "#228B22")
      if (is_uniform) {
        color_palette <- if (z_range[1] == 0) "white" else "#228B22"
      }
    } else {
      color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
    }
  } else {
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
  
  image(x = x_coords, y = y_coords, z = z, 
        col = color_palette,
        main = main,
        xlab = "X (World Coordinates)",
        ylab = "Y (World Coordinates)",
        asp = 1,
        useRaster = TRUE)
  
  grid(col = "gray80", lty = 1, lwd = 0.5)
  abline(h = 0, v = 0, col = "red", lty = 2, lwd = 1)
  
  if (n_survivors > 0) {
    if (color_by == "genotype" && !is.null(survivors$genotype)) {
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
    } else if (color_by == "phenotype" && !is.null(survivors$phenotype)) {
      gpal <- heat.colors(100)
      gr <- range(survivors$phenotype, na.rm = TRUE)
      if (gr[1] != gr[2]) {
        gscaled <- (survivors$phenotype - gr[1])/(gr[2]-gr[1])
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

# ============================================================================
# S3 METHODS FOR TWOLIFE_RESULT CLASS
# ============================================================================

#' Print Method for TWoLife Results
#' 
#' Provides a clean summary when a twolife_result object is printed to console.
#' Shows key simulation outcomes without overwhelming detail.
#' 
#' @param x A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param ... Additional arguments (currently unused)
#' 
#' @return The input object \code{x}, invisibly
#' 
#' @examples
#' set.seed(500)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 123
#' )
#' 
#' # Just type the object name to see summary
#' result
#' 
#' # Or explicitly call print
#' print(result)
#' 
#' @export
print.twolife_result <- function(x, ...) {
  cat("TWoLife Simulation Result\n")
  cat("==========================\n\n")
  cat("Status:", x$summary$status, "\n")
  cat("Final population:", x$summary$final_population_size, "\n")
  cat("Total events:", x$summary$total_events, "\n")
  cat("Duration:", round(x$summary$duration, 2), "time units\n")
  
  if (x$summary$status == "surviving" && !is.null(x$survivors)) {
    cat("\nSurvivors:", nrow(x$survivors), "individuals\n")
    if ("genotype" %in% names(x$survivors)) {
      cat("Genotype range:", 
          round(range(x$survivors$genotype), 3), "\n")
    }
  }
  
  cat("\nLandscape:", 
      nrow(x$parameters$landscape$habitat), "x", 
      ncol(x$parameters$landscape$habitat), "cells\n")
  
  cat("\nUse summary() for more details\n")
  cat("Use plot() to visualize results\n")
  
  invisible(x)
}

#' Summary Method for TWoLife Results
#' 
#' Provides detailed statistical summary of simulation outcomes including
#' population dynamics, spatial patterns, and genetic characteristics.
#' 
#' @param object A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param ... Additional arguments (currently unused)
#' 
#' @return A list of class 'summary.twolife_result' containing:
#'   \describe{
#'     \item{status}{Simulation outcome ("surviving" or "extinct")}
#'     \item{n_survivors}{Final population size}
#'     \item{initial_pop}{Starting population size}
#'     \item{events}{Total number of events simulated}
#'     \item{duration}{Simulation duration in time units}
#'     \item{extinction_time}{Time of extinction (NA if surviving)}
#'     \item{landscape_size}{Dimensions of landscape (rows x columns)}
#'     \item{history_detail}{Level of event history recorded}
#'   }
#' 
#' @examples
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.7,
#'   habitat_proportion = 0.4,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 456
#' )
#' 
#' # Get detailed summary
#' summary(result)
#' 
#' # Store summary for further use
#' sim_summary <- summary(result)
#' sim_summary$n_survivors
#' 
#' @export
summary.twolife_result <- function(object, ...) {
  structure(
    list(
      status = object$summary$status,
      n_survivors = object$summary$final_population_size,
      initial_pop = object$parameters$individual$initial_population_size,
      events = object$summary$total_events,
      duration = object$summary$duration,
      extinction_time = if (object$summary$status == "extinct") 
        max(object$events$times) else NA,
      landscape_size = dim(object$parameters$landscape$habitat),
      history_detail = object$parameters$simulation$history_detail
    ),
    class = "summary.twolife_result"
  )
}

#' Print Method for TWoLife Summary
#' 
#' Prints the summary of a twolife_result object in a readable format.
#' 
#' @param x A 'summary.twolife_result' object from \code{summary.twolife_result()}
#' @param ... Additional arguments (currently unused)
#' 
#' @return The input object \code{x}, invisibly
#' 
#' @export
print.summary.twolife_result <- function(x, ...) {
  cat("TWoLife Simulation Summary\n")
  cat("===========================\n\n")
  
  cat("Population Dynamics:\n")
  cat("  Initial population:", x$initial_pop, "\n")
  cat("  Final population:", x$n_survivors, "\n")
  cat("  Status:", x$status, "\n")
  
  if (!is.na(x$extinction_time)) {
    cat("  Extinction time:", round(x$extinction_time, 2), "time units\n")
  }
  
  cat("\nSimulation Details:\n")
  cat("  Total events:", x$events, "\n")
  cat("  Duration:", round(x$duration, 2), "time units\n")
  cat("  Events per time unit:", round(x$events / x$duration, 2), "\n")
  
  cat("\nLandscape:\n")
  cat("  Dimensions:", x$landscape_size[1], "x", x$landscape_size[2], "cells\n")
  cat("  Total cells:", prod(x$landscape_size), "\n")
  
  cat("\nHistory Detail:", x$history_detail, "\n")
  
  invisible(x)
}

#' Plot Method for TWoLife Results
#' 
#' Generic plot method for twolife_result objects. Calls 
#' \code{\link{plot_simulation_on_landscape}} with convenient syntax.
#' 
#' @param x A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param ... Additional arguments passed to \code{\link{plot_simulation_on_landscape}}.
#'   Common options include: point_size, point_color, color_by, landscape_colors, main
#' 
#' @return The input object \code{x}, invisibly
#' 
#' @examples
#' set.seed(200)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 789
#' )
#' 
#' # Simple plot using S3 method
#' plot(result)
#' 
#' # With options
#' plot(result, point_size = 3, color_by = "genotype")
#' 
#' # Equivalent to:
#' plot_simulation_on_landscape(result, point_size = 3, color_by = "genotype")
#' 
#' @seealso \code{\link{plot_simulation_on_landscape}} for full documentation of plotting options
#' 
#' @export
plot.twolife_result <- function(x, ...) {
  plot_simulation_on_landscape(x, ...)
  invisible(x)
}

# ============================================================================
# VALIDATION AND ANALYSIS FUNCTIONS
# ============================================================================

#' Validate Habitat Matching for Simulation Results
#' 
#' Visualizes and quantifies whether surviving individuals are positioned
#' in habitat matching their genotype or phenotype values. Calculates
#' correlation between individual trait values and habitat environmental values to assess
#' effectiveness of habitat selection or adaptation.
#' 
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param main Character. Plot title. If NULL, automatically generated. Default: NULL
#' @param point_size Numeric. Size of points in scatter plot. Default: 2
#' @param landscape_colors Character. Color scheme for landscape background. 
#'   Options: "habitat", "terrain", "viridis". Default: "terrain"
#' @param color_by Character. Which trait to analyze:
#'   \itemize{
#'     \item "genotype" - Analyze genotype-habitat matching
#'     \item "phenotype" - Analyze phenotype-habitat matching (recommended for plasticity)
#'   }
#'   Default: "phenotype"
#' @param show_stats Logical. If TRUE, prints statistical summary to console. Default: TRUE
#' @param show_legend Logical. If TRUE, displays color legend on plot. Default: TRUE
#' 
#' @return Invisibly returns a data frame with columns:
#'   \describe{
#'     \item{id}{Integer. Individual ID}
#'     \item{x}{Numeric. X coordinate in world space}
#'     \item{y}{Numeric. Y coordinate in world space}
#'     \item{genotype}{Numeric. Genotype value}
#'     \item{phenotype}{Numeric. Phenotype value}
#'     \item{habitat_value}{Numeric. Habitat quality at individual's location}
#'     \item{row_index}{Integer. Landscape row index (for debugging)}
#'     \item{col_index}{Integer. Landscape column index (for debugging)}
#'   }
#' 
#' @details
#' What This Function Measures:
#'   This function tests whether individuals are found in habitats matching their trait
#'   values. Perfect matching would show:
#'   \itemize{
#'     \item Individuals with trait value 0.8 in habitats with quality ~0.8
#'     \item Individuals with trait value 0.2 in habitats with quality ~0.2
#'     \item Strong positive correlation (r close to 1.0)
#'   }
#' 
#' Mathematical Basis:
#'   Calculates Pearson correlation coefficient:
#'   
#'   \deqn{r = \frac{Cov(trait, habitat)}{SD_{trait} \times SD_{habitat}}}
#'   
#'   Where:
#'   \itemize{
#'     \item trait = genotype or phenotype values of survivors
#'     \item habitat = habitat quality values at survivor locations
#'     \item r ranges from -1 (perfect negative correlation) to +1 (perfect positive correlation)
#'   }
#'   
#'   Also provides Spearman rank correlation (robust to outliers and non-linear relationships).
#' 
#' Visualization Components:
#'   Creates a two-panel plot:
#'   
#'   Left panel - Spatial distribution:
#'   \itemize{
#'     \item Shows landscape with habitat quality as background color
#'     \item Survivor locations as colored points (color = trait value)
#'     \item Good matching: point colors visually match background colors
#'     \item Poor matching: points of one color scattered across all background colors
#'   }
#'   
#'   Right panel - Correlation plot:
#'   \itemize{
#'     \item X-axis: habitat value at individual location (0 to 1)
#'     \item Y-axis: trait value (genotype or phenotype)
#'     \item Each point = one surviving individual
#'     \item Best-fit line and Pearson r displayed
#'     \item Perfect matching would show points along diagonal (r = 1.0)
#'   }
#' 
#' Connection to Simulation Parameters:
#'   Correlation strength depends on:
#'   \itemize{
#'     \item sampling_points: Higher values enable better habitat selection (stronger r)
#'     \item habitat_selection_temperatures: Lower values create stronger selection (stronger r)
#'     \item mutation_rates: Higher values create more trait variation to match habitat variation
#'     \item genotype_sds: Niche width - narrower niches (lower values) favor stronger matching
#'     \item simulation duration: Longer simulations allow more adaptation (stronger r over time)
#'   }
#' 
#' Statistical Output:
#'   When show_stats = TRUE, prints to console:
#'   \itemize{
#'     \item Sample size (number of survivors)
#'     \item Pearson correlation with p-value (tests H0: r = 0)
#'     \item Spearman rank correlation (non-parametric alternative)
#'     \item Mean and range of trait values
#'     \item Mean and range of habitat values at survivor locations
#'   }
#' 
#' @examples
#' # Create landscape
#' set.seed(123)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Simulation with genetic variation
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 25,
#'     base_birth_rate = 0.5,
#'     base_mortality_rate = 0.1
#'   ),
#'   genetic_params = list(
#'     genotype_means = runif(25, 0.4, 0.6),
#'     genotype_sds = 0.2
#'   ),
#'   simulation_params = list(max_events = 250),
#'   master_seed = 789
#' )
#' 
#' # Validate genotype matching
#' validation_data <- check_habitat_match(
#'   result,
#'   color_by = "genotype"
#' )
#' 
#' # Examine data
#' head(validation_data)
#' 
#' # Check correlation
#' if (!is.null(validation_data) && nrow(validation_data) > 0) {
#'   cat("Correlation:", cor(validation_data$genotype, validation_data$habitat_value), "\n")
#' }
#' 
#' \donttest{
#' # Compare with plasticity
#' result_plastic <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(initial_population_size = 10),
#'   genetic_params = list(
#'     genotype_means = runif(10, 0, 1),
#'     plasticities = 0.5,
#'     sampling_points = 5
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 999
#' )
#' 
#' # Validate phenotype matching (should show stronger correlation)
#' check_habitat_match(
#'   result_plastic,
#'   color_by = "phenotype",
#'   main = "Phenotype Matching (with Plasticity)"
#' )
#' }
#' 
#' @seealso \code{\link{habitat_mismatch}} for fitness-based analysis
#' 
#' @export
check_habitat_match <- function(simulation_result,
                                main = NULL,
                                point_size = 2,
                                landscape_colors = "terrain",
                                color_by = "phenotype",
                                show_stats = TRUE,
                                show_legend = TRUE) {
  
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  if (!color_by %in% c("genotype", "phenotype")) {
    stop("color_by must be either 'genotype' or 'phenotype'", call. = FALSE)
  }
  
  habitat_grid <- simulation_result$parameters$landscape$habitat
  cell_size <- simulation_result$parameters$landscape$cell_size
  survivors <- simulation_result$survivors
  n_survivors <- nrow(survivors)
  
  if (is.null(main)) {
    main <- paste("Habitat Matching Validation:", n_survivors, "Survivors")
  }
  
  if (n_survivors == 0) {
    cat("No survivors to validate\n")
    return(invisible(NULL))
  }
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  
  col_indices <- round((survivors$x + world_width/2) / cell_size) + 1
  row_indices <- n_rows - round((survivors$y + world_height/2) / cell_size)
  
  col_indices <- pmax(1, pmin(n_cols, col_indices))
  row_indices <- pmax(1, pmin(n_rows, row_indices))
  
  habitat_at_survivors <- numeric(n_survivors)
  for (i in 1:n_survivors) {
    habitat_at_survivors[i] <- habitat_grid[row_indices[i], col_indices[i]]
  }
  
  validation_data <- data.frame(
    id = survivors$id,
    x = survivors$x,
    y = survivors$y,
    genotype = survivors$genotype,
    phenotype = survivors$phenotype,
    habitat_value = habitat_at_survivors,
    row_index = row_indices,
    col_index = col_indices
  )
  
  trait_values <- if (color_by == "phenotype") survivors$phenotype else survivors$genotype
  trait_name <- if (color_by == "phenotype") "Phenotype" else "Genotype"
  
  if (show_stats) {
    cat("\nHabitat Matching Statistics (colored by", trait_name, "):\n")
    cat("----------------------------\n")
    cat("Number of survivors:", n_survivors, "\n")
    cat(trait_name, "distribution:\n")
    cat("  Mean:", round(mean(trait_values), 3), "\n")
    cat("  Range:", round(min(trait_values), 3), "to", 
        round(max(trait_values), 3), "\n")
    cat("Habitat value at survivor locations:\n")
    cat("  Mean:", round(mean(habitat_at_survivors), 3), "\n")
    cat("  Range:", round(min(habitat_at_survivors), 3), "to", 
        round(max(habitat_at_survivors), 3), "\n")
    
    if (n_survivors > 1) {
      correlation <- cor(trait_values, habitat_at_survivors)
      cat("\n", trait_name, "-Habitat correlation:", round(correlation, 3), "\n", sep="")
      if (!is.na(correlation) && correlation > 0.3) {
        cat("  -> Positive match: individuals in habitat matching their", tolower(trait_name), "\n")
      }
    }
  }
  
  z_range <- range(habitat_grid, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  is_binary <- all(habitat_grid %in% c(0, 1))
  
  if (landscape_colors == "habitat") {
    if (is_binary) {
      # Binary landscapes: 0 = white (matrix), 1 = green (habitat)
      color_palette <- c("white", "#228B22")
      if (is_uniform) {
        color_palette <- if (z_range[1] == 0) "white" else "#228B22"
      }
    } else {
      color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
    }
  } else {
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
                            terrain.colors(n_colors))
  }
  
  if (!is_uniform) {
    trait_normalized <- (trait_values - z_range[1]) / (z_range[2] - z_range[1])
    trait_normalized <- pmax(0, pmin(1, trait_normalized))
    
    if (is_binary) {
      color_indices <- ifelse(trait_normalized < 0.5, 1, 2)
    } else {
      color_indices <- pmax(1, pmin(100, round(trait_normalized * 99) + 1))
    }
    survivor_colors <- color_palette[color_indices]
  } else {
    survivor_colors <- rep(color_palette, n_survivors)
  }
  
  par(mfrow = c(1, 2))
  
  x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
  y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
  z <- t(apply(habitat_grid, 2, rev))
  
  image(x = x_coords, y = y_coords, z = z, 
        col = color_palette,
        main = "Landscape (Habitat Quality)",
        xlab = "X (World Coordinates)",
        ylab = "Y (World Coordinates)",
        asp = 1,
        useRaster = TRUE)
  
  grid(col = "gray80", lty = 1, lwd = 0.5)
  abline(h = 0, v = 0, col = "red", lty = 2, lwd = 1)
  
  plot(survivors$x, survivors$y, 
       col = survivor_colors,
       pch = 16, 
       cex = point_size,
       main = paste("Survivors (Colored by", trait_name, ")"),
       xlab = "X (World Coordinates)",
       ylab = "Y (World Coordinates)",
       asp = 1,
       xlim = c(-world_width/2, world_width/2),
       ylim = c(-world_height/2, world_height/2))
  
  grid(col = "gray80", lty = 1, lwd = 0.5)
  abline(h = 0, v = 0, col = "red", lty = 2, lwd = 1)
  
  if (!is_uniform && !is_binary && show_legend) {
    legend_vals <- seq(z_range[1], z_range[2], length.out = 5)
    legend("topright", 
           legend = round(legend_vals, 2),
           col = color_palette[seq(1, 100, length.out = 5)],
           pch = 16,
           title = paste(trait_name, "Value"),
           cex = 0.8,
           bg = "white")
  }
  
  mtext("If habitat selection works: dot colors should match background colors", 
        side = 1, line = -1.5, outer = TRUE, cex = 0.9, col = "blue")
  
  mtext(main, side = 3, line = -1.5, outer = TRUE, cex = 1.2, font = 2)
  
  par(mfrow = c(1, 1))
  
  return(invisible(validation_data))
}

#' Calculate Genotype-Habitat Fitness Statistics
#' 
#' Computes fitness-based metrics quantifying how well surviving individuals match their
#' habitat. Calculates fitness using the same Gaussian function used during the simulation,
#' providing comprehensive statistics on habitat matching quality and population fitness.
#' 
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param return_individuals Logical. If TRUE, includes individual-level data in output. Default: FALSE
#' 
#' @return A list of class 'genotype_habitat_mismatch' with components:
#'   \describe{
#'     \item{n_survivors}{Integer. Number of surviving individuals analyzed}
#'     \item{mean_fitness}{Numeric. Mean fitness across all survivors (0-1 scale)}
#'     \item{median_fitness}{Numeric. Median fitness}
#'     \item{sd_fitness}{Numeric. Standard deviation of fitness}
#'     \item{min_fitness}{Numeric. Minimum fitness value}
#'     \item{max_fitness}{Numeric. Maximum fitness value}
#'     \item{high_fitness_count}{Integer. Number with fitness > 0.8}
#'     \item{medium_fitness_count}{Integer. Number with fitness 0.5-0.8}
#'     \item{low_fitness_count}{Integer. Number with fitness < 0.5}
#'     \item{percent_high_fitness}{Numeric. Percentage with fitness > 0.8}
#'     \item{percent_medium_fitness}{Numeric. Percentage with fitness 0.5-0.8}
#'     \item{percent_low_fitness}{Numeric. Percentage with fitness < 0.5}
#'     \item{correlation}{Numeric. Pearson correlation between phenotype and habitat (-1 to 1)}
#'     \item{mean_genotype}{Numeric. Mean genotype value of survivors}
#'     \item{sd_genotype}{Numeric. Standard deviation of genotype}
#'     \item{mean_phenotype}{Numeric. Mean phenotype value of survivors}
#'     \item{sd_phenotype}{Numeric. Standard deviation of phenotype}
#'     \item{mean_habitat}{Numeric. Mean habitat value at survivor locations}
#'     \item{sd_habitat}{Numeric. Standard deviation of habitat at survivor locations}
#'     \item{mean_niche_width}{Numeric. Mean niche width (genotype_sds) of survivors}
#'     \item{sd_niche_width}{Numeric. Standard deviation of niche width}
#'     \item{mean_absolute_mismatch}{Numeric. Mean |phenotype - habitat|}
#'     \item{median_absolute_mismatch}{Numeric. Median |phenotype - habitat|}
#'     \item{individuals}{Data frame with per-individual data (only if return_individuals=TRUE).
#'       Contains columns: id, x, y, genotype, phenotype, width, habitat_value, fitness,
#'       absolute_mismatch, raw_mismatch, row_index, col_index}
#'   }
#'   
#'   The object has a custom print method that displays formatted statistics.
#' 
#' @details
#' Fitness Calculation:
#'   This function uses the same fitness formula that determines demographic rates during
#'   the simulation. Fitness quantifies how well an individual's phenotype matches the
#'   habitat quality at its location.
#'   
#'   For generalists (niche_width > 0):
#'   
#'   \deqn{fitness = exp\left(-\frac{(phenotype - habitat)^2}{2 \times niche\_width^2}\right)}
#'   
#'   Where:
#'   \itemize{
#'     \item phenotype = individual's expressed trait value (genotype + plasticity)
#'     \item habitat = habitat quality value at individual's location (0 to 1)
#'     \item niche_width = genotype_sds parameter (tolerance to mismatch)
#'     \item fitness ranges from 0 (complete mismatch) to 1 (perfect match)
#'   }
#'   
#'   For perfect specialists (niche_width = 0):
#'   
#'   \deqn{fitness = 1 \text{ if } |phenotype - habitat| < 0.001, \text{ otherwise } fitness = 0}
#'   
#'   This binary fitness function means specialists only survive in exactly matching habitat.
#' 
#' Fitness Interpretation:
#'   \itemize{
#'     \item fitness = 1.0: Perfect match. Phenotype exactly equals habitat quality.
#'     \item fitness > 0.8: High fitness. Within ~1 niche width of optimum.
#'     \item fitness = 0.5-0.8: Medium fitness. 1-2 niche widths from optimum.
#'     \item fitness < 0.5: Low fitness. More than 2 niche widths from optimum.
#'     \item fitness = 0.368 (1/e): Exactly 1 niche width from optimum (Gaussian inflection point).
#'   }
#'   
#'   Example: If niche_width = 0.2 and phenotype = 0.5:
#'   \itemize{
#'     \item habitat = 0.5 -> fitness = 1.0 (perfect)
#'     \item habitat = 0.7 (1 niche width away) -> fitness = 0.368
#'     \item habitat = 0.9 (2 niche widths away) -> fitness = 0.135 (poor)
#'   }
#' 
#' Connection to Simulation Demography:
#'   During simulation, this fitness value affects demographic rates. For generalists
#'   (genotype_sds > 0), the mortality rate is interpolated based on fitness:
#'   
#'   \deqn{mortality = mortality_{max} - (fitness_{relative} \times (mortality_{max} - mortality_{min}))}
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{mortality_{max} = matrix\_mortality\_multiplier \times base\_mortality\_rate}
#'     \item \eqn{mortality_{min} = base\_mortality\_rate}
#'     \item \eqn{fitness_{relative}} = fitness at current habitat / fitness at optimal habitat
#'   }
#'   
#'   Higher fitness leads to lower mortality rate, higher birth rate (for generalists),
#'   and greater probability of leaving offspring.
#' 
#' @examples
#' # Create landscape
#' set.seed(456)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Run simulation
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 25,
#'     base_birth_rate = 0.5,
#'     base_mortality_rate = 0.1
#'   ),
#'   genetic_params = list(
#'     genotype_means = runif(25, 0.4, 0.6),
#'     genotype_sds = 0.2
#'   ),
#'   simulation_params = list(max_events = 250),
#'   master_seed = 567
#' )
#' 
#' # Calculate fitness statistics
#' fitness_stats <- habitat_mismatch(result)
#' 
#' # Print formatted output
#' print(fitness_stats)
#' 
#' # Access specific metrics
#' cat("Mean fitness:", fitness_stats$mean_fitness, "\n")
#' cat("Correlation:", fitness_stats$correlation, "\n")
#' cat("% high fitness:", fitness_stats$percent_high_fitness, "\n")
#' 
#' \donttest{
#' # Compare with neutral scenario
#' result_neutral <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(initial_population_size = 10),
#'   genetic_params = list(genotype_means = runif(10, 0, 1)),
#'   simulation_params = list(max_events = 1500, neutral_mode = TRUE),
#'   master_seed = 888
#' )
#' 
#' fitness_neutral <- habitat_mismatch(result_neutral)
#' 
#' # Compare
#' cat("With selection - Mean fitness:", fitness_stats$mean_fitness, "\n")
#' cat("Neutral - Mean fitness:", fitness_neutral$mean_fitness, "\n")
#' cat("With selection - Correlation:", fitness_stats$correlation, "\n")
#' cat("Neutral - Correlation:", fitness_neutral$correlation, "\n")
#' }
#' 
#' @seealso \code{\link{check_habitat_match}} for visual validation
#' 
#' @export
habitat_mismatch <- function(simulation_result, 
                             return_individuals = FALSE) {
  
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  habitat_grid <- simulation_result$parameters$landscape$habitat
  cell_size <- simulation_result$parameters$landscape$cell_size
  survivors <- simulation_result$survivors
  n_survivors <- nrow(survivors)
  
  if (n_survivors == 0) {
    warning("No survivors to analyze", call. = FALSE)
    return(list(
      n_survivors = 0,
      mean_fitness = NA,
      median_fitness = NA,
      sd_fitness = NA,
      min_fitness = NA,
      high_fitness_count = 0,
      medium_fitness_count = 0,
      low_fitness_count = 0,
      mean_absolute_mismatch = NA,
      correlation = NA,
      mean_genotype = NA,
      mean_phenotype = NA,
      mean_habitat = NA
    ))
  }
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  
  col_indices <- round((survivors$x + world_width/2) / cell_size) + 1
  row_indices <- n_rows - round((survivors$y + world_height/2) / cell_size)
  
  col_indices <- pmax(1, pmin(n_cols, col_indices))
  row_indices <- pmax(1, pmin(n_rows, row_indices))
  
  habitat_at_survivors <- numeric(n_survivors)
  for (i in 1:n_survivors) {
    habitat_at_survivors[i] <- habitat_grid[row_indices[i], col_indices[i]]
  }
  
  survivor_phenotypes <- survivors$phenotype
  survivor_widths <- survivors$width
  
  fitness_values <- numeric(n_survivors)
  for (i in 1:n_survivors) {
    if (survivor_widths[i] > 0) {
      deviation <- habitat_at_survivors[i] - survivor_phenotypes[i]
      variance <- survivor_widths[i]^2
      fitness_values[i] <- exp(-(deviation^2) / (2 * variance))
    } else {
      fitness_values[i] <- ifelse(
        abs(habitat_at_survivors[i] - survivor_phenotypes[i]) < 0.001,
        1.0,
        0.0
      )
    }
  }
  
  absolute_mismatch <- abs(survivor_phenotypes - habitat_at_survivors)
  raw_mismatch <- survivor_phenotypes - habitat_at_survivors
  
  correlation <- if (n_survivors > 1) {
    cor(survivor_phenotypes, habitat_at_survivors)
  } else {
    NA
  }
  
  high_fitness <- sum(fitness_values > 0.8)
  medium_fitness <- sum(fitness_values >= 0.5 & fitness_values <= 0.8)
  low_fitness <- sum(fitness_values < 0.5)
  
  results <- list(
    n_survivors = n_survivors,
    mean_fitness = mean(fitness_values),
    median_fitness = median(fitness_values),
    sd_fitness = sd(fitness_values),
    min_fitness = min(fitness_values),
    max_fitness = max(fitness_values),
    high_fitness_count = high_fitness,
    medium_fitness_count = medium_fitness,
    low_fitness_count = low_fitness,
    percent_high_fitness = (high_fitness / n_survivors) * 100,
    percent_medium_fitness = (medium_fitness / n_survivors) * 100,
    percent_low_fitness = (low_fitness / n_survivors) * 100,
    mean_absolute_mismatch = mean(absolute_mismatch),
    median_absolute_mismatch = median(absolute_mismatch),
    correlation = correlation,
    mean_genotype = mean(survivors$genotype),
    sd_genotype = sd(survivors$genotype),
    mean_phenotype = mean(survivor_phenotypes),
    sd_phenotype = sd(survivor_phenotypes),
    mean_habitat = mean(habitat_at_survivors),
    sd_habitat = sd(habitat_at_survivors),
    mean_niche_width = mean(survivor_widths),
    sd_niche_width = sd(survivor_widths)
  )
  
  if (return_individuals) {
    results$individuals <- data.frame(
      id = survivors$id,
      x = survivors$x,
      y = survivors$y,
      genotype = survivors$genotype,
      phenotype = survivor_phenotypes,
      width = survivor_widths,
      habitat_value = habitat_at_survivors,
      fitness = fitness_values,
      absolute_mismatch = absolute_mismatch,
      raw_mismatch = raw_mismatch,
      row_index = row_indices,
      col_index = col_indices
    )
  }
  
  class(results) <- c("genotype_habitat_mismatch", "list")
  return(results)
}


#' Print Method for Mismatch Statistics
#' 
#' @param x A genotype_habitat_mismatch object
#' @param ... Additional arguments (ignored)
#' 
#' @return x (invisibly)
#' 
#' @export
print.genotype_habitat_mismatch <- function(x, ...) {
  cat("Genotype-Habitat Fitness Analysis\n")
  cat("==================================\n\n")
  
  if (x$n_survivors == 0) {
    cat("No survivors to analyze\n")
    return(invisible(x))
  }
  
  cat("Sample size:", x$n_survivors, "individuals\n\n")
  
  cat("Fitness Metrics (0-1 scale, based on PHENOTYPE-environment match):\n")
  cat("------------------------------------------------------------------\n")
  cat("Mean fitness:", round(x$mean_fitness, 4), "\n")
  cat("Median fitness:", round(x$median_fitness, 4), "\n")
  cat("SD of fitness:", round(x$sd_fitness, 4), "\n")
  cat("Range:", round(x$min_fitness, 4), "to", round(x$max_fitness, 4), "\n\n")
  
  cat("Fitness Quality Distribution:\n")
  cat("-----------------------------\n")
  cat("High fitness (>0.8):", x$high_fitness_count, 
      paste0("(", round(x$percent_high_fitness, 1), "%)"), "\n")
  cat("Medium fitness (0.5-0.8):", x$medium_fitness_count, 
      paste0("(", round(x$percent_medium_fitness, 1), "%)"), "\n")
  cat("Low fitness (<0.5):", x$low_fitness_count, 
      paste0("(", round(x$percent_low_fitness, 1), "%)"), "\n\n")
  
  cat("Phenotype-Habitat Relationship:\n")
  cat("-------------------------------\n")
  cat("Correlation:", round(x$correlation, 4), "\n")
  cat("Mean genotype:", round(x$mean_genotype, 4), 
      "+/- SD:", round(x$sd_genotype, 4), "\n")
  cat("Mean phenotype:", round(x$mean_phenotype, 4), 
      "+/- SD:", round(x$sd_phenotype, 4), "\n")
  cat("Mean habitat at survivors:", round(x$mean_habitat, 4), 
      "+/- SD:", round(x$sd_habitat, 4), "\n")
  cat("Mean niche width:", round(x$mean_niche_width, 4), 
      "+/- SD:", round(x$sd_niche_width, 4), "\n")
  cat("Mean absolute mismatch:", round(x$mean_absolute_mismatch, 4), 
      "(for reference)\n\n")
  
  cat("Interpretation:\n")
  cat("---------------\n")
  
  if (x$mean_fitness > 0.8) {
    cat("Excellent fitness: Individuals well-matched to their habitat\n")
  } else if (x$mean_fitness > 0.6) {
    cat("Good fitness: Most individuals reasonably well-matched\n")
  } else if (x$mean_fitness > 0.4) {
    cat("Moderate fitness: Individuals somewhat matched to habitat\n")
  } else {
    cat("Poor fitness: Individuals poorly matched to habitat\n")
  }
  
  if (!is.na(x$correlation)) {
    if (x$correlation > 0.5) {
      cat("Strong positive correlation: Effective habitat selection\n")
    } else if (x$correlation > 0.3) {
      cat("Moderate positive correlation: Some habitat selection\n")
    } else if (x$correlation > 0) {
      cat("Weak positive correlation: Limited habitat selection\n")
    } else {
      cat("No/negative correlation: Poor or absent habitat selection\n")
    }
  } else {
    cat("Correlation: NA (insufficient data variation)\n")
  }
  
  if (x$percent_high_fitness > 50) {
    cat("Majority in high-fitness locations (>50% above 0.8 fitness)\n")
  } else if (x$percent_low_fitness > 30) {
    cat("Many individuals in low-fitness locations (>30% below 0.5 fitness)\n")
  }
  
  invisible(x)
}