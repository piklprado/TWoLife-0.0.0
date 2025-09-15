# TWoLife

[![R-CMD-check](https://github.com/yourusername/TWoLife/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/TWoLife/actions)

**Individual-based spatial population simulations with genetic evolution and habitat selection**

TWoLife is an R package for running spatially-explicit, individual-based population models. It simulates organisms that move, reproduce, and die across heterogeneous landscapes, with genetic evolution and behavioral habitat selection.

## Features

- **Spatial Population Dynamics**: Individuals move across realistic landscapes
- **Genetic Evolution**: Mutation, selection, and adaptation to local conditions  
- **Habitat Selection**: Active behavioral choice of optimal environments
- **Demographic Processes**: Density-dependent birth and death rates
- **Fast Simulation**: C++ backend via Rcpp for computational efficiency
- **Flexible Landscapes**: Binary, continuous, or fractal habitat patterns
- **Comprehensive Visualization**: Built-in plotting functions for results

## Installation

### From Source (Development)

```r
# Install development tools if needed
install.packages("devtools")

# Install TWoLife from local source
devtools::install()

# Or if package is on GitHub
devtools::install_github("yourusername/TWoLife")
```

### Dependencies

The package requires:
- R >= 4.0.0
- Rcpp
- rflsgen (for fractal landscape generation)

Optional for enhanced visualization:
- viridisLite

## Quick Start

```r
library(TWoLife)

# 1. Create a test landscape
habitat <- create_fractal_landscape(
  cells_per_side = 20, 
  fractality = 0.6, 
  habitat_proportion = 0.4
)

# 2. Run a basic simulation
result <- twolife_simulation(
  landscape_params = list(habitat = habitat),
  individual_params = list(initial_population_size = 50),
  simulation_params = list(max_events = 1000)
)

# 3. View results
print(result)
plot_simulation_on_landscape(result)

# 4. Analyze population trajectory
trajectory <- compute_population_size(result)
plot(trajectory$time, trajectory$population_size, 
     type = "l", main = "Population Over Time")
```

## Advanced Usage

### Genetic Evolution

```r
# Simulation with genetic diversity and evolution
result <- twolife_simulation(
  landscape_params = list(habitat = habitat),
  individual_params = list(initial_population_size = 40),
  genetic_params = list(
    genotype_means = runif(40, 0.2, 0.8),     # Diverse starting genotypes
    genotype_sds = rep(0.15, 40),             # Niche width tolerance
    mutation_rates = rep(0.03, 40),           # Evolution rate
    plasticities = rep(0.02, 40),             # Phenotypic flexibility
    sampling_points = rep(15, 40)             # Habitat selection intensity
  ),
  simulation_params = list(max_events = 2000)
)

# View genetic structure
plot_simulation_on_landscape(result, show_genotypes = TRUE)
```

### Landscape Effects

```r
# Compare different landscape types
landscapes <- list(
  continuous = create_fractal_landscape(15, fractality = 0.7),
  fragmented = create_fractal_landscape(15, fractality = 0.3, habitat_proportion = 0.3),
  corner = create_corner_landscape(15, corner = "top-left")
)

results <- lapply(landscapes, function(habitat) {
  twolife_simulation(
    landscape_params = list(habitat = habitat),
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 500)
  )
})

# Compare final population sizes
sapply(results, function(r) r$summary$final_population_size)
```

## Documentation

### Package Vignette

For detailed examples and biological interpretation:

```r
vignette("introduction", package = "TWoLife")
```

### Function Help

```r
?twolife_simulation      # Main simulation function
?create_fractal_landscape # Landscape generation
?plot_simulation_on_landscape # Result visualization
help(package = "TWoLife") # Complete function list
```

## Key Functions

| Function | Purpose |
|----------|---------|
| `twolife_simulation()` | Run complete individual-based simulation |
| `create_fractal_landscape()` | Generate realistic test landscapes |
| `create_corner_landscape()` | Create simple test landscapes |
| `plot_simulation_on_landscape()` | Visualize simulation results |
| `plot_landscape_world_coords()` | Display landscape patterns |
| `compute_population_size()` | Extract population trajectories |

## Parameter Categories

### Landscape Parameters
- **habitat**: Habitat quality matrix (required)
- **cell_size**: Spatial resolution 
- **boundary_condition**: Edge effects (absorbing/periodic/reflective)

### Individual Parameters  
- **initial_population_size**: Starting number of individuals
- **neighbor_radius**: Density dependence scale
- **step_length**: Movement distance per event
- **base_birth_rate**, **base_mortality_rate**: Demographic rates

### Genetic Parameters
- **genotype_means**: Optimal habitat values for each individual
- **genotype_sds**: Niche width (habitat tolerance)
- **mutation_rates**: Evolution rate per generation
- **sampling_points**: Habitat selection behavior intensity

## Examples Repository

See the `vignettes/` directory for:
- Basic workflow examples
- Genetic evolution demonstrations  
- Landscape effect comparisons
- Parameter sensitivity analysis

## Development

### Building from Source

```r
# Generate documentation
devtools::document()

# Run tests
devtools::test()

# Check package
devtools::check()

# Install locally  
devtools::install()
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Submit a pull request

## Citation

If you use TWoLife in your research, please cite:

```
[Add citation information when published]
```

## License

[Add license information]

## Contact

- **Issues**: [GitHub Issues](https://github.com/yourusername/TWoLife/issues)
- **Email**: your.email@domain.com

---

**Note**: This package is under active development. API may change between versions.