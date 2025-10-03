# functions.R - Complete TWoLife R interface
# UPDATED: Fixed NA handling in correlation comparisons

# ============================================================================
# CORE SIMULATION FUNCTIONS
# ============================================================================

#' Run TWoLife Individual-Based Simulation
#' 
#' Runs a spatially-explicit individual-based simulation with habitat selection,
#' genetic variation, and demographic processes. Supports rectangular landscapes
#' and customizable habitat selection parameters.
#' 
#' @param landscape_params List containing landscape parameters
#' @param individual_params List containing individual parameters
#' @param genetic_params List containing genetic parameters
#' @param simulation_params List containing simulation control parameters
#' @param history_detail Character. Level of event history detail to record:
#'   "minimal" for only time, event type, and individual ID (fastest, smallest memory),
#'   "standard" for adding spatial coordinates, patch ID, and genotype (default),
#'   "full" for adding phenotype and niche width (enables exact historical reconstruction)
#' @param master_seed Integer seed for reproducible simulations (optional)
#' @param ... Additional named arguments (deprecated)
#' @param output_file Optional file path for detailed event output
#' 
#' @return A list of class 'twolife_result' containing simulation results
#' 
#' @export
twolife_simulation <- function(landscape_params = list(), 
                               individual_params = list(), 
                               genetic_params = list(),
                               simulation_params = list(),
                               history_detail = "standard",
                               master_seed = NULL,
                               ..., 
                               output_file = NULL) {
  
  # Validate history_detail parameter
  valid_levels <- c("minimal", "standard", "full")
  if (!history_detail %in% valid_levels) {
    stop("history_detail must be one of: ", 
         paste(valid_levels, collapse = ", "), 
         call. = FALSE)
  }
  
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
    history_detail = history_detail,
    master_seed = master_seed,
    output_file = output_file
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
#' Works with all history_detail levels (only requires times and types).
#' 
#' @param result Result object from twolife_simulation
#' 
#' @return Data frame with columns time and population_size
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

#' Reconstruct Population State at Specific Time
#' 
#' Replays the event log from a simulation to reconstruct the exact positions
#' and states of all living individuals at a specified time point. Useful for
#' temporal analysis and creating animations.
#' 
#' Requires history_detail = "standard" or "full"
#' 
#' @param simulation_result A twolife_result object
#' @param target_time Numeric. Time point to reconstruct population state
#' @param color_by Character. What to color points by: "genotype", "phenotype", or "none"
#' @param show_plot Logical. Display visualization (default: TRUE)
#' @param point_size Numeric. Size of points in plot (default: 2)
#' 
#' @return Invisibly returns list with time, n_alive, population, events_processed, and history_level
#' 
#' @examples
#' \dontrun{
#' result <- twolife_simulation(..., history_detail = "standard")
#' state <- snapshot_at_time(result, target_time = 50)
#' head(state$population)
#' }
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

#' Check if Snapshot Capability is Available
#' 
#' Helper function to check if a result has sufficient history detail
#' for snapshot_at_time() function.
#' 
#' @param simulation_result Result object from twolife_simulation
#' 
#' @return Logical indicating if snapshot_at_time() can be used
#' 
#' @examples
#' \dontrun{
#' result <- twolife_simulation(..., history_detail = "minimal")
#' if (can_snapshot(result)) {
#'   snapshot_at_time(result, 50)
#' } else {
#'   message("Re-run with history_detail='standard' or 'full'")
#' }
#' }
#' 
#' @export
can_snapshot <- function(simulation_result) {
  history_level <- simulation_result$parameters$simulation$history_detail
  
  if (is.null(history_level)) {
    # Infer from available fields
    return("x_coordinates" %in% names(simulation_result$events))
  }
  
  return(history_level %in% c("standard", "full"))
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
#' @param cells_per_row Integer greater than or equal to 4. Number of rows
#' @param cells_per_col Integer greater than or equal to 4. Number of columns. If NULL, creates square landscape.
#' @param corner Character. Corner location for habitat: "top-left", "top-right",
#'   "bottom-left", or "bottom-right" (default: "top-left")
#' @param corner_size Integer. Size of habitat corner in cells. If NULL,
#'   automatically set to 25 percent of the smallest dimension (default: NULL)
#' @param return_as_landscape_params Logical. Return as list(habitat=matrix)
#'   suitable for twolife_simulation (default: FALSE)
#' 
#' @return Either a binary matrix or list with habitat component
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

#' Quick Plot of Simulation Results
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
      # FIXED: Check for NA before comparing
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

#' Compare Mismatch Statistics Across Multiple Simulations
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
  
  # Show history detail level
  history_level <- x$parameters$simulation$history_detail
  if (!is.null(history_level)) {
    cat("History detail:", history_level, "\n")
  }
  
  if (x$summary$final_population_size > 0) {
    cat("Survivors: Use result$survivors for details\n")
    cat("Columns: id, x, y, genotype, phenotype, width\n")
  }
  cat("\nUse compute_population_size() for trajectory analysis.\n")
  cat("Use plot_simulation_on_landscape() for visualization.\n")
  if (!is.null(history_level) && history_level %in% c("standard", "full")) {
    cat("Use snapshot_at_time() for temporal reconstruction.\n")
  }
  invisible(x)
}

#' Print Method for Mismatch Statistics
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
    cat("Excellent fitness: Individuals well-matched to their habitat\n")
  } else if (x$mean_fitness > 0.6) {
    cat("Good fitness: Most individuals reasonably well-matched\n")
  } else if (x$mean_fitness > 0.4) {
    cat("Moderate fitness: Individuals somewhat matched to habitat\n")
  } else {
    cat("Poor fitness: Individuals poorly matched to habitat\n")
  }
  
  # FIXED: Check for NA before comparing
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