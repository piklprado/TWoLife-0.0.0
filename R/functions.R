# functions.R - Complete TWoLife R interface
# UPDATED: Support for rectangular landscapes, habitat_selection_temperatures,
#          phenotypes, widths, and extended validation/analysis functions

# ============================================================================
# CORE SIMULATION FUNCTIONS
# ============================================================================

#' Run TWoLife Individual-Based Simulation
#' 
#' Runs a spatially-explicit individual-based simulation with habitat selection,
#' genetic variation, and demographic processes. Supports rectangular landscapes
#' and customizable habitat selection parameters.
#' 
#' @param landscape_params List containing landscape parameters:
#'   \itemize{
#'     \item \code{habitat}: Required matrix of habitat values
#'     \item \code{cell_size}: Size of each cell (default: 1.0)
#'     \item \code{boundary_condition}: 1=reflecting, 2=periodic (default: 1)
#'     \item \code{density_type}: 1=local, 2=global (default: 1)
#'     \item \code{matrix_mortality_multiplier}: Mortality multiplier in matrix (default: 2.0)
#'     \item \code{matrix_dispersal_multiplier}: Dispersal multiplier in matrix (default: 0.5)
#'   }
#' @param individual_params List containing individual parameters:
#'   \itemize{
#'     \item \code{neighbor_radius}: Radius for local density (default: 2.0)
#'     \item \code{vision_angle}: Angle of vision in radians (default: pi)
#'     \item \code{step_length}: Movement distance per step (default: 1.0)
#'     \item \code{base_dispersal_rate}: Base dispersal probability (default: 0.1)
#'     \item \code{base_birth_rate}: Base birth probability (default: 0.3)
#'     \item \code{base_mortality_rate}: Base mortality probability (default: 0.20)
#'     \item \code{birth_density_slope}: Birth rate density dependence (default: 0.02)
#'     \item \code{mortality_density_slope}: Mortality rate density dependence (default: 0.02)
#'     \item \code{initial_population_size}: Starting population (default: 200)
#'     \item \code{initial_placement_mode}: 1=random in habitat, 2=random anywhere, 3=custom coordinates (default: 1)
#'     \item \code{initial_x_coordinates}: X coordinates for mode 3 (optional)
#'     \item \code{initial_y_coordinates}: Y coordinates for mode 3 (optional)
#'   }
#' @param genetic_params List containing genetic parameters:
#'   \itemize{
#'     \item \code{genotype_means}: Mean genotype values (length 1 or population_size)
#'     \item \code{genotype_sds}: Standard deviations of genotypes (length 1 or population_size)
#'     \item \code{mutation_rates}: Mutation rates (length 1 or population_size)
#'     \item \code{plasticities}: Phenotypic plasticity values (length 1 or population_size)
#'     \item \code{sampling_points}: Number of sampling points for habitat selection (length 1 or population_size)
#'     \item \code{habitat_selection_temperatures}: Temperature parameter for habitat selection softmax (length 1 or population_size, must be positive)
#'   }
#' @param simulation_params List containing simulation control parameters:
#'   \itemize{
#'     \item \code{max_events}: Maximum number of events (default: 50 * population_size)
#'     \item \code{neutral_mode}: Disable selection if TRUE (default: FALSE)
#'   }
#' @param master_seed Integer seed for reproducible simulations (optional)
#' @param ... Additional named arguments (deprecated)
#' @param output_file Optional file path for detailed event output
#' 
#' @return A list of class 'twolife_result' containing:
#'   \itemize{
#'     \item \code{summary}: Summary statistics (final_population_size, total_events, duration, status)
#'     \item \code{survivors}: Data frame of surviving individuals with id, x, y, genotype, phenotype, width
#'     \item \code{spatial}: Spatial dimensions (world_width, world_height, num_patches)
#'     \item \code{events}: Event history (times, types, individual_ids, patch_ids, coordinates, genotypes)
#'     \item \code{parameters}: All input parameters organized by category
#'   }
#' 
#' @examples
#' \dontrun{
#' # Create a simple landscape
#' landscape <- create_fractal_landscape(50, 50, fractality = 0.5, 
#'                                      habitat_proportion = 0.3,
#'                                      return_as_landscape_params = TRUE)
#' 
#' # Run simulation
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(initial_population_size = 100),
#'   genetic_params = list(genotype_means = 0.5, genotype_sds = 0.1),
#'   master_seed = 42
#' )
#' 
#' # View results
#' print(result)
#' plot_simulation_on_landscape(result)
#' }
#' 
#' @export
twolife_simulation <- function(landscape_params = list(), 
                               individual_params = list(), 
                               genetic_params = list(),
                               simulation_params = list(),
                               master_seed = NULL,
                               ..., 
                               output_file = NULL) {
  
  dots <- list(...)
  if ("habitat_grid" %in% names(dots)) {
    stop("habitat_grid argument is deprecated. Please pass habitat inside landscape_params$habitat", call. = FALSE)
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
  if (is.null(simulation_params$max_events) && !("max_events" %in% names(dots))) {
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
    master_seed = master_seed,
    output_file = output_file
  )
  
  final_pop <- as.integer(length(result$survivor_x))
  total_events <- length(result$event_times)
  duration <- if(total_events > 0) max(result$event_times) else 0
  
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
    events = list(
      times = result$event_times,
      types = result$event_types,
      individual_ids = result$individual_ids,
      patch_ids = result$patch_ids,
      x_coordinates = result$x_coordinates,
      y_coordinates = result$y_coordinates,
      genotypes = result$genotypes
    ),
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
#' 
#' @param result Result object from twolife_simulation
#' 
#' @return Data frame with columns:
#'   \itemize{
#'     \item \code{time}: Event times
#'     \item \code{population_size}: Cumulative population size
#'   }
#' 
#' @examples
#' \dontrun{
#' result <- twolife_simulation(...)
#' trajectory <- compute_population_size(result)
#' plot(trajectory$time, trajectory$population_size, type = "l")
#' }
#' 
#' @export
compute_population_size <- function(result) {
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

# ============================================================================
# LANDSCAPE GENERATION FUNCTIONS
# ============================================================================

#' Create Fractal Landscape
#' 
#' Generates a fractal landscape using spatial autocorrelation. Supports both
#' square and rectangular landscapes, and can create either continuous or
#' binary (habitat/matrix) landscapes.
#' 
#' @param cells_per_row Integer. Number of cells per row
#' @param cells_per_col Integer. Number of cells per column. If NULL, creates square landscape.
#' @param fractality Numeric between 0 and 1. Higher values create more spatially
#'   autocorrelated (clumped) patterns. 0 = random, 1 = highly structured.
#' @param min_value Numeric. Minimum habitat value for continuous landscapes (default: 0.0)
#' @param max_value Numeric. Maximum habitat value for continuous landscapes (default: 1.0)
#' @param habitat_proportion Numeric between 0 and 1. If provided, creates binary
#'   landscape with this proportion of cells as habitat (value 1), rest as matrix
#'   (value 0). If NULL, creates continuous landscape.
#' @param return_as_landscape_params Logical. If TRUE, returns list(habitat=matrix)
#'   suitable for direct use in twolife_simulation. If FALSE, returns matrix only.
#' 
#' @return Either a matrix (if return_as_landscape_params=FALSE) or a list with
#'   habitat component (if return_as_landscape_params=TRUE)
#' 
#' @examples
#' \dontrun{
#' # Binary landscape (30% habitat)
#' landscape1 <- create_fractal_landscape(
#'   cells_per_row = 50,
#'   fractality = 0.6,
#'   habitat_proportion = 0.3,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Continuous landscape (habitat quality varies)
#' landscape2 <- create_fractal_landscape(
#'   cells_per_row = 50,
#'   fractality = 0.5,
#'   min_value = 0,
#'   max_value = 1
#' )
#' 
#' # Rectangular landscape
#' landscape3 <- create_fractal_landscape(
#'   cells_per_row = 50,
#'   cells_per_col = 100,
#'   fractality = 0.7,
#'   habitat_proportion = 0.4
#' )
#' }
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
  
  if (cells_per_row == cells_per_col && requireNamespace("rflsgen", quietly = TRUE)) {
    raw_fractal <- rflsgen::flsgen_terrain(
      cells_per_row,    
      cells_per_col,    
      fractality         
    )
    fractal_matrix <- as.matrix(raw_fractal)
  } else {
    fractal_matrix <- generate_rectangular_fractal(cells_per_row, cells_per_col, fractality)
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

#' Generate Rectangular Fractal Pattern
#' 
#' Internal function to generate fractal patterns for rectangular landscapes.
#' Uses iterative smoothing with neighbor averaging.
#' 
#' @param rows Number of rows
#' @param cols Number of columns
#' @param fractality Fractal parameter (0-1)
#' 
#' @return Matrix with fractal pattern
#' 
#' @keywords internal
generate_rectangular_fractal <- function(rows, cols, fractality) {
  fractal <- matrix(runif(rows * cols), nrow = rows, ncol = cols)
  n_iterations <- max(1, round(fractality * 10))
  
  for (iter in 1:n_iterations) {
    smoothed <- matrix(0, nrow = rows, ncol = cols)
    
    for (i in 1:rows) {
      for (j in 1:cols) {
        neighbors <- c()
        if (i > 1) neighbors <- c(neighbors, fractal[i-1, j])
        if (i < rows) neighbors <- c(neighbors, fractal[i+1, j])
        if (j > 1) neighbors <- c(neighbors, fractal[i, j-1])
        if (j < cols) neighbors <- c(neighbors, fractal[i, j+1])
        if (i > 1 && j > 1) neighbors <- c(neighbors, fractal[i-1, j-1])
        if (i > 1 && j < cols) neighbors <- c(neighbors, fractal[i-1, j+1])
        if (i < rows && j > 1) neighbors <- c(neighbors, fractal[i+1, j-1])
        if (i < rows && j < cols) neighbors <- c(neighbors, fractal[i+1, j+1])
        
        weight_original <- 1 - fractality
        weight_neighbors <- fractality / length(neighbors)
        
        smoothed[i, j] <- weight_original * fractal[i, j] + sum(neighbors * weight_neighbors)
      }
    }
    fractal <- smoothed
  }
  return(fractal)
}

#' Create Corner Test Landscapes
#' 
#' Creates binary landscapes with habitat confined to one corner, useful for
#' testing habitat selection and dispersal. Supports rectangular landscapes.
#' 
#' @param cells_per_row Integer >= 4. Number of rows
#' @param cells_per_col Integer >= 4. Number of columns. If NULL, creates square landscape.
#' @param corner Character. Corner location for habitat: "top-left", "top-right",
#'   "bottom-left", or "bottom-right" (default: "top-left")
#' @param corner_size Integer. Size of habitat corner in cells. If NULL,
#'   automatically set to 25% of the smallest dimension (default: NULL)
#' @param return_as_landscape_params Logical. Return as list(habitat=matrix)
#'   suitable for twolife_simulation (default: FALSE)
#' 
#' @return Either a binary matrix or list with habitat component
#' 
#' @examples
#' \dontrun{
#' # Square landscape with top-left habitat
#' landscape <- create_corner_landscape(
#'   cells_per_row = 50,
#'   corner = "top-left",
#'   corner_size = 10,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Rectangular landscape with bottom-right habitat
#' landscape2 <- create_corner_landscape(
#'   cells_per_row = 30,
#'   cells_per_col = 60,
#'   corner = "bottom-right"
#' )
#' }
#' 
#' @export
create_corner_landscape <- function(cells_per_row, 
                                    cells_per_col = NULL,
                                    corner = "top-left",
                                    corner_size = NULL,
                                    return_as_landscape_params = FALSE) {
  
  if (!is.numeric(cells_per_row) || cells_per_row < 4 || cells_per_row != round(cells_per_row)) {
    stop("cells_per_row must be an integer >= 4", call. = FALSE)
  }
  
  if (is.null(cells_per_col)) {
    cells_per_col <- cells_per_row
  }
  
  if (!is.numeric(cells_per_col) || cells_per_col < 4 || cells_per_col != round(cells_per_col)) {
    stop("cells_per_col must be an integer >= 4", call. = FALSE)
  }
  
  valid_corners <- c("top-left", "top-right", "bottom-left", "bottom-right")
  if (!corner %in% valid_corners) {
    stop("corner must be one of: ", paste(valid_corners, collapse = ", "), call. = FALSE)
  }
  
  if (is.null(corner_size)) {
    corner_size <- max(2, round(min(cells_per_row, cells_per_col) / 4))
  }
  
  if (!is.numeric(corner_size) || corner_size < 1 || 
      corner_size >= cells_per_row || corner_size >= cells_per_col) {
    stop("corner_size must be between 1 and min(cells_per_row, cells_per_col)-1", call. = FALSE)
  }
  
  landscape <- matrix(0, nrow = cells_per_row, ncol = cells_per_col)
  
  corner_coords <- switch(corner,
                          "top-left" = list(
                            rows = (cells_per_row - corner_size + 1):cells_per_row,
                            cols = 1:corner_size
                          ),
                          "top-right" = list(
                            rows = (cells_per_row - corner_size + 1):cells_per_row,
                            cols = (cells_per_col - corner_size + 1):cells_per_col
                          ),
                          "bottom-left" = list(
                            rows = 1:corner_size,
                            cols = 1:corner_size
                          ),
                          "bottom-right" = list(
                            rows = 1:corner_size,
                            cols = (cells_per_col - corner_size + 1):cells_per_col
                          )
  )
  
  landscape[corner_coords$rows, corner_coords$cols] <- 1
  
  if (return_as_landscape_params) {
    return(list(habitat = landscape))
  } else {
    return(landscape)
  }
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

#' Plot Landscape with World Coordinates
#' 
#' Visualizes landscape in world coordinate system with proper axis scaling.
#' Works with both square and rectangular landscapes.
#' 
#' @param landscape_data Matrix or list(habitat = matrix)
#' @param cell_size Numeric. Size of each cell (default: 1.0)
#' @param filename Character. File path for export (currently unused)
#' @param main Character. Plot title (default: "Landscape (World Coordinates)")
#' @param colors Character. Color scheme: "habitat" (green tones), "terrain",
#'   or "viridis" (default: "habitat")
#' @param show_legend Logical. Show color scale legend (default: TRUE, currently unused)
#' @param add_grid Logical. Add grid lines and axes (default: TRUE)
#' 
#' @return Invisibly returns NULL
#' 
#' @examples
#' \dontrun{
#' landscape <- create_fractal_landscape(50, 50, fractality = 0.5,
#'                                      habitat_proportion = 0.3)
#' plot_landscape_world_coords(landscape)
#' }
#' 
#' @export
plot_landscape_world_coords <- function(landscape_data, 
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
      color_palette <- c("#F5F5DC", "#228B22")
      if (is_uniform) {
        color_palette <- if(z_range[1] == 0) "#F5F5DC" else "#228B22"
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
    grid(col = "white", lty = 1, lwd = 0.5)
    abline(h = 0, v = 0, col = "red", lty = 2, lwd = 2)
  }
  
  invisible(NULL)
}

#' Plot Simulation Results on Landscape
#' 
#' Visualizes simulation results by overlaying survivor positions on the
#' landscape. Points can be colored by genotype, phenotype, or a single color.
#' 
#' @param simulation_result A twolife_result object
#' @param point_size Numeric. Size of points (default: 2)
#' @param point_color Character. Color of survivor points (default: "red")
#' @param point_shape Numeric. Point shape (pch value) (default: 16)
#' @param color_by Character. What to color points by: "none" (use point_color),
#'   "genotype", or "phenotype" (default: "none")
#' @param landscape_colors Character. Color scheme for background landscape:
#'   "habitat", "terrain", or "viridis" (default: "habitat")
#' @param main Character. Plot title. If NULL, auto-generates title with
#'   survivor count (default: NULL)
#' @param filename Character. Optional filename for export (currently unused)
#' @param add_stats Logical. Add population statistics (currently unused) (default: TRUE)
#' 
#' @return Invisibly returns the simulation result object
#' 
#' @examples
#' \dontrun{
#' result <- twolife_simulation(...)
#' plot_simulation_on_landscape(result, color_by = "genotype")
#' plot_simulation_on_landscape(result, color_by = "phenotype")
#' }
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
  
  # Validate color_by parameter
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
      color_palette <- c("#F5F5DC", "#228B22")
      if (is_uniform) {
        color_palette <- if(z_range[1] == 0) "#F5F5DC" else "#228B22"
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
  
  grid(col = "white", lty = 1, lwd = 0.5)
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

#' Quick Plot of Simulation Results
#' 
#' Convenience wrapper for plot_simulation_on_landscape with sensible defaults.
#' 
#' @param simulation_result A twolife_result object
#' @param ... Additional arguments passed to plot_simulation_on_landscape
#' 
#' @examples
#' \dontrun{
#' result <- twolife_simulation(...)
#' quick_plot_result(result)
#' }
#' 
#' @export
quick_plot_result <- function(simulation_result, ...) {
  plot_simulation_on_landscape(simulation_result, ...)
}

# ============================================================================
# VALIDATION AND ANALYSIS FUNCTIONS
# ============================================================================

#' Validate Habitat Matching for Simulation Results
#' 
#' Creates a two-panel visualization showing the landscape and survivor positions
#' colored by genotype or phenotype. Validates that individuals are in appropriate 
#' habitat by checking the correlation between trait and habitat quality at survivor
#' locations.
#' 
#' @param simulation_result A twolife_result object
#' @param main Character. Overall plot title. If NULL, auto-generates title
#'   with survivor count (default: NULL)
#' @param point_size Numeric. Size of survivor points (default: 2)
#' @param landscape_colors Character. Color scheme: "habitat", "terrain", or
#'   "viridis" (default: "terrain")
#' @param color_by Character. What to color survivor points by: "genotype" or
#'   "phenotype" (default: "phenotype")
#' @param show_stats Logical. Print habitat quality statistics to console
#'   (default: TRUE)
#' 
#' @return Invisibly returns data frame with columns:
#'   \itemize{
#'     \item \code{id}: Individual ID
#'     \item \code{x, y}: Spatial coordinates
#'     \item \code{genotype}: Genotype value
#'     \item \code{phenotype}: Phenotype value
#'     \item \code{habitat_value}: Habitat quality at individual's location
#'     \item \code{row_index, col_index}: Grid cell indices
#'   }
#' 
#' @details
#' The function displays two side-by-side plots:
#' \itemize{
#'   \item Left: Landscape colored by habitat quality
#'   \item Right: Survivor positions colored by selected trait (genotype or phenotype)
#' }
#' If habitat selection is working, dot colors should visually match the
#' background colors (individuals with high trait values in high-quality
#' habitat, and vice versa).
#' 
#' @examples
#' \dontrun{
#' result <- twolife_simulation(...)
#' validation_data <- validate_habitat_matching(result, color_by = "phenotype")
#' head(validation_data)
#' }
#' 
#' @export
validate_habitat_matching <- function(simulation_result,
                                      main = NULL,
                                      point_size = 2,
                                      landscape_colors = "terrain",
                                      color_by = "phenotype",
                                      show_stats = TRUE) {
  
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  # Validate color_by parameter
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
  
  # Select which trait to use for statistics
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
      if (correlation > 0.3) {
        cat("  -> Positive match: individuals in habitat matching their", tolower(trait_name), "\n")
      }
    }
  }
  
  z_range <- range(habitat_grid, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  is_binary <- all(habitat_grid %in% c(0, 1))
  
  if (landscape_colors == "habitat") {
    if (is_binary) {
      color_palette <- c("#F5F5DC", "#228B22")
      if (is_uniform) {
        color_palette <- if (z_range[1] == 0) "#F5F5DC" else "#228B22"
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
  
  grid(col = "white", lty = 1, lwd = 0.5)
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
  
  if (!is_uniform && !is_binary) {
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

#' Batch Validate Multiple Simulations
#' 
#' Validates habitat matching for multiple simulation results simultaneously,
#' creating a grid of validation plots for easy comparison.
#' 
#' @param result_list Named list of twolife_result objects
#' @param point_size Numeric. Size of points (default: 1.5)
#' @param landscape_colors Character. Color scheme (default: "terrain")
#' 
#' @return Invisibly returns named list of validation data frames
#' 
#' @examples
#' \dontrun{
#' results <- list(
#'   low_temp = result1,
#'   medium_temp = result2,
#'   high_temp = result3
#' )
#' validations <- batch_validate_habitat_matching(results)
#' }
#' 
#' @export
batch_validate_habitat_matching <- function(result_list, 
                                            point_size = 1.5,
                                            landscape_colors = "terrain") {
  
  n_results <- length(result_list)
  
  if (n_results == 0) {
    stop("result_list cannot be empty", call. = FALSE)
  }
  
  if (n_results <= 2) {
    par(mfrow = c(1, n_results * 2))
  } else if (n_results <= 4) {
    par(mfrow = c(2, 4))
  } else if (n_results <= 6) {
    par(mfrow = c(3, 4))
  } else {
    par(mfrow = c(4, 4))
  }
  
  validation_results <- list()
  
  for (i in seq_along(result_list)) {
    result_name <- if (is.null(names(result_list)[i])) paste("Result", i) else names(result_list)[i]
    result <- result_list[[i]]
    
    cat("\n", result_name, "\n", sep = "")
    validation_results[[i]] <- validate_habitat_matching(
      result, 
      main = result_name,
      point_size = point_size,
      landscape_colors = landscape_colors,
      show_stats = TRUE
    )
    names(validation_results)[i] <- result_name
  }
  
  par(mfrow = c(1, 1))
  
  return(invisible(validation_results))
}

#' Calculate Genotype-Habitat Mismatch Statistics Using Fitness
#' 
#' Computes fitness-based measures of genotype-habitat matching. Fitness is
#' calculated using PHENOTYPE (not genotype) as exp(-(habitat - phenotype)^2 / (2 * width^2)), 
#' representing how well each individual's expressed phenotype matches its environment.
#' 
#' @param simulation_result A twolife_result object
#' @param return_individuals Logical. If TRUE, includes individual-level data
#'   in the output (default: FALSE)
#' 
#' @return List of class 'genotype_habitat_mismatch' containing:
#'   \itemize{
#'     \item \code{n_survivors}: Number of surviving individuals
#'     \item \code{mean_fitness, median_fitness, sd_fitness}: Summary statistics
#'     \item \code{min_fitness, max_fitness}: Fitness range
#'     \item \code{high_fitness_count}: Count with fitness > 0.8
#'     \item \code{medium_fitness_count}: Count with fitness 0.5-0.8
#'     \item \code{low_fitness_count}: Count with fitness < 0.5
#'     \item \code{percent_high_fitness, percent_medium_fitness, percent_low_fitness}: Percentages
#'     \item \code{mean_absolute_mismatch, median_absolute_mismatch}: Absolute difference measures
#'     \item \code{correlation}: Phenotype-habitat correlation
#'     \item \code{mean_genotype, sd_genotype}: Genotype distribution
#'     \item \code{mean_phenotype, sd_phenotype}: Phenotype distribution
#'     \item \code{mean_habitat, sd_habitat}: Habitat distribution at survivors
#'     \item \code{mean_niche_width, sd_niche_width}: Niche width distribution
#'     \item \code{individuals}: Data frame (if return_individuals=TRUE)
#'   }
#' 
#' @examples
#' \dontrun{
#' result <- twolife_simulation(...)
#' mismatch_stats <- calculate_genotype_habitat_mismatch(result)
#' print(mismatch_stats)
#' 
#' # Get individual-level data
#' detailed <- calculate_genotype_habitat_mismatch(result, return_individuals = TRUE)
#' head(detailed$individuals)
#' }
#' 
#' @export
calculate_genotype_habitat_mismatch <- function(simulation_result, 
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
  
  # Use phenotypes and widths from survivors data frame
  survivor_phenotypes <- survivors$phenotype
  survivor_widths <- survivors$width
  
  # Calculate fitness based on phenotype-environment match
  fitness_values <- numeric(n_survivors)
  for (i in 1:n_survivors) {
    if (survivor_widths[i] > 0) {
      deviation <- habitat_at_survivors[i] - survivor_phenotypes[i]
      variance <- survivor_widths[i]^2
      fitness_values[i] <- exp(-(deviation^2) / (2 * variance))
    } else {
      # Perfect specialist - fitness is 1 if exact match, 0 otherwise
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

#' Compare Mismatch Statistics Across Multiple Simulations
#' 
#' Computes and compares fitness-based habitat matching statistics across
#' multiple simulation scenarios. Useful for comparing different parameter
#' settings or selection strengths.
#' 
#' @param result_list Named list of twolife_result objects
#' 
#' @return Invisibly returns list containing:
#'   \itemize{
#'     \item \code{comparison}: Data frame with summary statistics for each scenario
#'     \item \code{detailed}: List of full mismatch analysis objects
#'   }
#' 
#' @examples
#' \dontrun{
#' results <- list(
#'   no_selection = result1,
#'   weak_selection = result2,
#'   strong_selection = result3
#' )
#' comparison <- compare_mismatch_statistics(results)
#' print(comparison$comparison)
#' }
#' 
#' @export
compare_mismatch_statistics <- function(result_list) {
  
  if (!is.list(result_list) || length(result_list) == 0) {
    stop("result_list must be a non-empty list", call. = FALSE)
  }
  
  mismatch_list <- lapply(result_list, calculate_genotype_habitat_mismatch)
  
  comparison_df <- data.frame(
    Scenario = if (is.null(names(result_list))) {
      paste("Result", seq_along(result_list))
    } else {
      names(result_list)
    },
    N_Survivors = sapply(mismatch_list, function(x) x$n_survivors),
    Mean_Fitness = sapply(mismatch_list, function(x) round(x$mean_fitness, 4)),
    Median_Fitness = sapply(mismatch_list, function(x) round(x$median_fitness, 4)),
    Percent_High_Fitness = sapply(mismatch_list, function(x) round(x$percent_high_fitness, 1)),
    Percent_Low_Fitness = sapply(mismatch_list, function(x) round(x$percent_low_fitness, 1)),
    Correlation = sapply(mismatch_list, function(x) round(x$correlation, 4)),
    Mean_Mismatch = sapply(mismatch_list, function(x) round(x$mean_absolute_mismatch, 4)),
    Mean_Niche_Width = sapply(mismatch_list, function(x) round(x$mean_niche_width, 4)),
    stringsAsFactors = FALSE
  )
  
  cat("Fitness-Based Habitat Matching Comparison\n")
  cat("==========================================\n\n")
  print(comparison_df, row.names = FALSE)
  
  cat("\nKey Insights:\n")
  cat("-------------\n")
  
  if (nrow(comparison_df) > 1) {
    best_fitness <- which.max(comparison_df$Mean_Fitness)
    best_correlation <- which.max(comparison_df$Correlation)
    most_high_fitness <- which.max(comparison_df$Percent_High_Fitness)
    
    cat("Best mean fitness:", comparison_df$Scenario[best_fitness], 
        "(", comparison_df$Mean_Fitness[best_fitness], ")\n")
    cat("Best correlation:", comparison_df$Scenario[best_correlation], 
        "(r =", comparison_df$Correlation[best_correlation], ")\n")
    cat("Most high-fitness individuals:", comparison_df$Scenario[most_high_fitness],
        "(", comparison_df$Percent_High_Fitness[most_high_fitness], "%)\n")
    
    if (best_fitness == best_correlation && best_correlation == most_high_fitness) {
      cat("\n-> Same scenario excels in all fitness metrics!\n")
    }
  }
  
  return(invisible(list(
    comparison = comparison_df,
    detailed = mismatch_list
  )))
}

# ============================================================================
# PRINT METHODS
# ============================================================================

#' Print Method for TWoLife Results
#' 
#' @param x A twolife_result object
#' @param ... Additional arguments (unused)
#' 
#' @return Invisibly returns the input object
#' 
#' @export
print.twolife_result <- function(x, ...) {
  cat("TWoLife Simulation Result\n")
  cat("========================\n")
  cat("Status:", x$summary$status, "\n")
  cat("Final population:", x$summary$final_population_size, "\n") 
  cat("Duration:", round(x$summary$duration, 2), "\n")
  cat("Total events:", x$summary$total_events, "\n")
  
  if (!is.null(x$spatial$world_width) && !is.null(x$spatial$world_height)) {
    cat("World size:", round(x$spatial$world_width, 2), "x", round(x$spatial$world_height, 2), "\n")
    if (abs(x$spatial$world_width - x$spatial$world_height) < 1e-10) {
      cat("Landscape shape: Square\n")
    } else {
      cat("Landscape shape: Rectangular\n")
    }
  } else {
    cat("World size:", x$spatial$world_size, "x", x$spatial$world_size, "\n")
  }
  
  if (x$summary$final_population_size > 0) {
    cat("Survivors: Use result$survivors for details\n")
    cat("Columns: id, x, y, genotype, phenotype, width\n")
  }
  cat("\nUse compute_population_size() for trajectory analysis.\n")
  cat("Use plot_simulation_on_landscape() for visualization.\n")
  invisible(x)
}

#' Print Method for Mismatch Statistics
#' 
#' @param x A genotype_habitat_mismatch object
#' @param ... Additional arguments (unused)
#' 
#' @return Invisibly returns the input object
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
      "± SD:", round(x$sd_genotype, 4), "\n")
  cat("Mean phenotype:", round(x$mean_phenotype, 4), 
      "± SD:", round(x$sd_phenotype, 4), "\n")
  cat("Mean habitat at survivors:", round(x$mean_habitat, 4), 
      "± SD:", round(x$sd_habitat, 4), "\n")
  cat("Mean niche width:", round(x$mean_niche_width, 4), 
      "± SD:", round(x$sd_niche_width, 4), "\n")
  cat("Mean absolute mismatch:", round(x$mean_absolute_mismatch, 4), 
      "(for reference)\n\n")
  
  cat("Interpretation:\n")
  cat("---------------\n")
  
  if (x$mean_fitness > 0.8) {
    cat("✓ Excellent fitness: Individuals well-matched to their habitat\n")
  } else if (x$mean_fitness > 0.6) {
    cat("✓ Good fitness: Most individuals reasonably well-matched\n")
  } else if (x$mean_fitness > 0.4) {
    cat("⚠ Moderate fitness: Individuals somewhat matched to habitat\n")
  } else {
    cat("✗ Poor fitness: Individuals poorly matched to habitat\n")
  }
  
  if (x$correlation > 0.5) {
    cat("✓ Strong positive correlation: Effective habitat selection\n")
  } else if (x$correlation > 0.3) {
    cat("○ Moderate positive correlation: Some habitat selection\n")
  } else if (x$correlation > 0) {
    cat("○ Weak positive correlation: Limited habitat selection\n")
  } else {
    cat("✗ No/negative correlation: Poor or absent habitat selection\n")
  }
  
  if (x$percent_high_fitness > 50) {
    cat("✓ Majority in high-fitness locations (>50% above 0.8 fitness)\n")
  } else if (x$percent_low_fitness > 30) {
    cat("⚠ Many individuals in low-fitness locations (>30% below 0.5 fitness)\n")
  }
  
  invisible(x)
}