#' TWoLife: Individual-Based Spatial Population Simulations
#'
#' @description 
#' TWoLife provides tools for running individual-based spatial population simulations
#' with genetic evolution, habitat selection, and demographic processes. The package
#' implements spatially-explicit models where individuals move, reproduce, and die
#' across heterogeneous landscapes. Simulations use efficient C++ algorithms to
#' enable large-scale studies of eco-evolutionary dynamics.
#'
#' @details
#' The main function \code{\link{twolife_simulation}} runs complete simulations.
#' Key features include:
#' \itemize{
#'   \item Spatial population dynamics in continuous or discrete landscapes
#'   \item Genetic evolution with mutation and selection
#'   \item Phenotypic plasticity and genotype-environment interactions
#'   \item Habitat selection behavior based on genetic traits
#'   \item Density-dependent demography (births and deaths)
#'   \item Multiple boundary conditions (absorbing, reflective, periodic)
#'   \item Various landscape types and spatial structures
#'   \item Fast C++ simulation engine via Rcpp for computational efficiency
#'   \item Flexible event history recording for detailed analyses
#' }
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{twolife_simulation}}}{Run complete individual-based simulation
#'     with customizable parameters for landscape, demography, and evolution}
#'   \item{\code{\link{create_fractal_landscape}}}{Generate realistic test landscapes
#'     with controlled spatial autocorrelation and fragmentation}
#'   \item{\code{\link{plot_simulation_on_landscape}}}{Visualize simulation results
#'     overlaid on habitat quality maps}
#'   \item{\code{\link{plot_landscape}}}{Plot landscape habitat quality
#'     in world coordinate system}
#'   \item{\code{\link{population_size}}}{Extract population size trajectories
#'     over time from simulation results}
#'   \item{\code{\link{snapshot_at_time}}}{Reconstruct population spatial distribution
#'     at specific time points during simulation}
#' }
#'
#' @section Typical Workflow:
#' \enumerate{
#'   \item Create or load a landscape matrix representing habitat quality
#'   \item Set simulation parameters for demography, movement, and evolution
#'   \item Run simulation with \code{twolife_simulation()}
#'   \item Analyze results using population trajectories and spatial patterns
#'   \item Visualize outcomes with plotting functions
#' }
#'
#' @section Getting Started:
#' See \code{vignette("introduction", package = "TWoLife")} for detailed examples
#' and tutorials demonstrating typical use cases.
#' 
#' Quick start example:
#' \preformatted{
#' # Create a simple landscape
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 30,
#'   fractality = 0.5,
#'   habitat_proportion = 0.4,
#'   return_as_landscape_params = TRUE
#' )
#' 
#' # Run simulation
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   simulation_params = list(max_events = 5000)
#' )
#' 
#' # Visualize results
#' plot_simulation_on_landscape(result)
#' }
#'
#' @section Parameter Categories:
#' \describe{
#'   \item{Landscape Parameters}{Habitat quality matrix, cell size, boundary
#'     conditions (absorbing/reflective/periodic), density calculation methods,
#'     and matrix habitat effects on demography and movement}
#'   \item{Individual Parameters}{Initial population size, movement behavior
#'     (neighbor radius, vision angle, step length), demographic rates
#'     (birth, mortality, dispersal), density-dependence parameters, and
#'     initial placement options}
#'   \item{Genetic Parameters}{Evolutionary mechanisms including genotype
#'     initialization, mutation rates, phenotypic plasticity, habitat selection
#'     strength, and number of genetic loci}
#'   \item{Simulation Parameters}{Runtime control including maximum events,
#'     neutral vs. adaptive mode, event history detail level, random seed
#'     for reproducibility, and optional file output}
#' }
#' 
#' @section Model Description:
#' TWoLife implements a continuous-time, spatially-explicit individual-based model
#' where demographic and movement events occur stochastically. The model combines:
#' \itemize{
#'   \item \strong{Spatial dynamics}: Individuals occupy continuous coordinates
#'     in a 2D landscape divided into discrete habitat patches
#'   \item \strong{Demography}: Birth and death rates depend on local habitat
#'     quality, individual traits, and population density
#'   \item \strong{Movement}: Individuals disperse based on habitat selection,
#'     with movement rates influenced by local conditions
#'   \item \strong{Evolution}: Genetic traits evolve through mutation and
#'     selection on habitat matching
#'   \item \strong{Plasticity}: Phenotypes can respond to local environmental
#'     conditions via reaction norms
#' }
#' 
#' The simulation uses a Gillespie algorithm for event-driven dynamics, making
#' it computationally efficient for long time scales and large populations.
#'
#' @section Performance:
#' Typical performance on modern hardware (2024):
#' \itemize{
#'   \item ~1 second per 10,000 events
#'   \item Scales linearly with population size
#'   \item Memory usage: ~10 MB per 100,000 events (standard history detail)
#'   \item Landscapes up to 200x200 cells handle efficiently
#'   \item Populations of 1,000-10,000 individuals commonly used
#' }
#'
#' @section Example Datasets:
#' The package includes two example landscapes:
#' \describe{
#'   \item{\code{\link{several_small}}}{Fragmented landscape with multiple
#'     small habitat patches for testing metapopulation dynamics}
#'   \item{\code{\link{single_large}}}{Continuous landscape with one large
#'     habitat patch for baseline comparisons}
#' }
#'
#' @references
#' Grimm, V., Berger, U., Bastiansen, F., Eliassen, S., Ginot, V., Giske, J.,
#' Goss-Custard, J., Grand, T., Heinz, S.K., Huse, G., Huth, A., Jepsen, J.U.,
#' Jørgensen, C., Mooij, W.M., Müller, B., Pe'er, G., Piou, C., Railsback, S.F.,
#' Robbins, A.M., Robbins, M.M., Rossmanith, E., Rüger, N., Strand, E.,
#' Souissi, S., Stillman, R.A., Vabø, R., Visser, U., & DeAngelis, D.L. (2006).
#' A standard protocol for describing individual-based and agent-based models.
#' \emph{Ecological Modelling}, 198(1-2), 115-126.
#' \doi{10.1016/j.ecolmodel.2006.04.023}
#' 
#' Bolnick, D.I., Svanbäck, R., Fordyce, J.A., Yang, L.H., Davis, J.M.,
#' Hulsey, C.D., & Forister, M.L. (2003). The ecology of individuals:
#' Incidence and implications of individual specialization.
#' \emph{The American Naturalist}, 161(1), 1-28. \doi{10.1086/343878}
#'
#' @author 
#' Lucas Freitas (Maintainer) \email{rodriguesdefreitas@@gmail.com}
#' (\href{https://orcid.org/0000-0002-2773-0981}{ORCID})
#' 
#' Contributors: piLaboratory
#'
#' @keywords internal
#' @aliases TWoLife-package
#' @useDynLib TWoLife, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom mathjaxr preview_rd
"_PACKAGE"
