#' Example Landscape: Multiple Small Habitat Patches
#'
#' A pre-computed binary landscape matrix with multiple small, disconnected
#' habitat patches. Useful for demonstrating habitat fragmentation effects,
#' dispersal limitation, and metapopulation dynamics.
#'
#' @format A 15-by-15 integer matrix (225 cells) with binary values:
#' \describe{
#'   \item{0}{Matrix (non-habitat) cells - 176 cells, 78\%}
#'   \item{1}{Habitat cells - 49 cells, 22\%}
#' }
#' 
#' The habitat is distributed across 3-4 disconnected patches of varying sizes,
#' separated by matrix habitat that creates barriers to dispersal.
#'
#' @details
#' This fragmented landscape is useful for studying:
#' \itemize{
#'   \item Population persistence in fragmented habitats
#'   \item Dispersal behavior and patch colonization
#'   \item Metapopulation dynamics and rescue effects
#'   \item Extinction-recolonization dynamics
#'   \item Source-sink dynamics across different patch sizes
#' }
#' 
#' With default parameters, populations face challenges from fragmentation
#' including increased extinction risk in small patches and isolation effects.
#'
#' @examples
#' # Load and examine the dataset
#' data(several_small)
#' dim(several_small)
#' table(several_small)
#' 
#' # Visualize the landscape
#' plot_landscape_world_coords(several_small,
#'                             main = "Fragmented Habitat",
#'                             colors = "habitat")
#' 
#' # Run a simulation
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = several_small),
#'   individual_params = list(
#'     initial_population_size = 30,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.2
#'   ),
#'   simulation_params = list(max_events = 2000)
#' )
#' 
#' # View results
#' plot_simulation_on_landscape(result)
#' 
#' # Check outcome
#' cat("Final population:", nrow(result$survivors), "\n")
#' 
#' # Population trajectory
#' pop <- compute_population_size(result)
#' plot(pop$time, pop$population_size, type = "l",
#'      xlab = "Time", ylab = "Population Size",
#'      main = "Dynamics in Fragmented Habitat")
#' 
#' \donttest{
#' # Compare with continuous habitat
#' data(single_large)
#' result2 <- twolife_simulation(
#'   landscape_params = list(habitat = single_large),
#'   individual_params = list(initial_population_size = 30),
#'   simulation_params = list(max_events = 2000)
#' )
#' 
#' cat("Fragmented final pop:", nrow(result$survivors), "\n")
#' cat("Continuous final pop:", nrow(result2$survivors), "\n")
#' }
#'
#' @seealso
#' \code{\link{single_large}} for a continuous habitat landscape
#' 
#' \code{\link{create_corner_landscape}} for creating test landscapes
#'
#' @keywords datasets
"several_small"

#' Example Landscape: Single Large Habitat Patch
#'
#' A pre-computed binary landscape matrix with a single, large, contiguous
#' habitat patch. Useful as a baseline for comparing with fragmented landscapes
#' and studying populations in continuous habitat.
#'
#' @format A 15-by-15 integer matrix (225 cells) with binary values:
#' \describe{
#'   \item{0}{Matrix (non-habitat) cells - 180 cells, 80\%}
#'   \item{1}{Habitat cells - 45 cells, 20\%}
#' }
#' 
#' The habitat forms a single unfragmented patch in the central region,
#' surrounded by matrix habitat.
#'
#' @details
#' This continuous landscape is useful for studying:
#' \itemize{
#'   \item Baseline population dynamics without fragmentation
#'   \item Density-dependent processes in high-quality habitat
#'   \item Comparing effects of fragmentation vs. continuous habitat
#' }
#' 
#' Populations typically show more stable dynamics compared to fragmented
#' landscapes, with free movement throughout the habitat patch and reduced
#' extinction risk.
#'
#' @examples
#' # Load and examine the dataset
#' data(single_large)
#' dim(single_large)
#' table(single_large)
#' 
#' # Visualize the landscape
#' plot_landscape_world_coords(single_large,
#'                             main = "Continuous Habitat",
#'                             colors = "habitat")
#' 
#' # Run a simulation
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = single_large),
#'   individual_params = list(
#'     initial_population_size = 30,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.2
#'   ),
#'   simulation_params = list(max_events = 2000)
#' )
#' 
#' # View results
#' plot_simulation_on_landscape(result)
#' 
#' # Check outcome
#' cat("Final population:", nrow(result$survivors), "\n")
#' cat("Status:", result$summary$status, "\n")
#' 
#' # Population trajectory
#' pop <- compute_population_size(result)
#' plot(pop$time, pop$population_size, type = "l",
#'      xlab = "Time", ylab = "Population Size",
#'      main = "Dynamics in Continuous Habitat")
#' 
#' \donttest{
#' # Compare with fragmented habitat
#' data(several_small)
#' result2 <- twolife_simulation(
#'   landscape_params = list(habitat = several_small),
#'   individual_params = list(initial_population_size = 30),
#'   simulation_params = list(max_events = 2000)
#' )
#' 
#' # Compare trajectories
#' pop1 <- compute_population_size(result)
#' pop2 <- compute_population_size(result2)
#' 
#' plot(pop1$time, pop1$population_size, type = "l",
#'      col = "blue", lwd = 2, ylim = c(0, max(pop1$population_size, pop2$population_size)),
#'      xlab = "Time", ylab = "Population Size",
#'      main = "Continuous vs. Fragmented")
#' lines(pop2$time, pop2$population_size, col = "red", lwd = 2)
#' legend("topright", c("Continuous", "Fragmented"),
#'        col = c("blue", "red"), lwd = 2)
#' }
#'
#' @seealso
#' \code{\link{several_small}} for a fragmented habitat landscape
#' 
#' \code{\link{create_corner_landscape}} for creating test landscapes
#'
#' @keywords datasets
"single_large"