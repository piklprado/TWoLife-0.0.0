# TWoLife

> Individual-based spatial population simulations with genetic evolution and habitat selection

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

```r
# Install development tools if needed
install.packages("devtools")

# Install TWoLife from GitHub
devtools::install_github("DeFreitasLR/TWoLife-0.0.0")
```

## Dependencies

The package requires:
- R >= 4.0.0
- Rcpp

Optional for enhanced visualization:
- viridisLite

Note: Fractal landscape generation is implemented internally using base R functions.

## Quick Start

```r
library(TWoLife)

# 1. Create a test landscape
landscape <- create_fractal_landscape(
  cells_per_row = 20,
  fractality = 0.6,
  habitat_proportion = 0.4,
  return_as_landscape_params = TRUE
)

# 2. Run a basic simulation
result <- twolife_simulation(
  landscape_params = landscape,
  individual_params = list(initial_population_size = 50),
  simulation_params = list(max_events = 1000)
)

# 3. View results
print(result)
plot(result)

# 4. Analyze population trajectory
trajectory <- population_size(result)
plot(trajectory$time, trajectory$population_size,
     type = "l", main = "Population Over Time")
```

## Advanced Usage

### Simulation with Genetic Diversity

```r
# Create landscape
landscape <- create_fractal_landscape(
  cells_per_row = 20,
  fractality = 0.6,
  habitat_proportion = 0.5,
  return_as_landscape_params = TRUE
)

# Run simulation with genetic variation
result <- twolife_simulation(
  landscape_params = landscape,
  individual_params = list(
    initial_population_size = 40,
    base_birth_rate = 0.4,
    base_mortality_rate = 0.15
  ),
  genetic_params = list(
    genotype_means = runif(40, 0.2, 0.8),  # Diverse starting genotypes
    genotype_sds = 0.15,                    # Niche width tolerance
    mutation_rates = 0.03,                  # Evolution rate
    plasticities = 0.02,                    # Phenotypic flexibility
    sampling_points = 15                    # Habitat selection intensity
  ),
  simulation_params = list(max_events = 2000),
  master_seed = 123
)

# Visualize results with genotypes
plot(result, color_by = "genotype")

# Check habitat matching
validation <- check_habitat_match(result)
```

### Compare Different Landscapes

```r
# Generate different landscape types
continuous <- create_fractal_landscape(
  cells_per_row = 15, 
  fractality = 0.7,
  habitat_proportion = 0.6,
  return_as_landscape_params = TRUE
)

fragmented <- create_fractal_landscape(
  cells_per_row = 15, 
  fractality = 0.3, 
  habitat_proportion = 0.3,
  return_as_landscape_params = TRUE
)

# Run simulations
results <- lapply(list(continuous = continuous, fragmented = fragmented), 
  function(landscape) {
    twolife_simulation(
      landscape_params = landscape,
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 500),
      master_seed = 456
    )
  }
)

# Compare final population sizes
sapply(results, function(r) r$summary$final_population_size)
```

## Documentation

For detailed examples and biological interpretation:

```r
# View vignette
vignette("introduction", package = "TWoLife")

# Function help
?twolife_simulation        # Main simulation function
?create_fractal_landscape  # Landscape generation
?plot.twolife_result       # Result visualization (S3 method)
?population_size           # Extract population trajectories
?check_habitat_match       # Validate genotype-habitat matching

# Complete function list
help(package = "TWoLife")
```

## Main Functions

| Function | Purpose |
|----------|---------|
| `twolife_simulation()` | Run complete individual-based simulation |
| `create_fractal_landscape()` | Generate realistic fractal landscapes |
| `plot.twolife_result()` | Visualize simulation results (S3 method) |
| `plot_landscape()` | Display landscape patterns |
| `population_size()` | Extract population trajectories over time |
| `snapshot_at_time()` | Reconstruct population state at specific time |
| `check_habitat_match()` | Validate genotype-habitat matching |
| `habitat_mismatch()` | Calculate fitness statistics |
| `plot_simulation_on_landscape()` | Overlay simulation on landscape |

## Key Parameters

### Landscape Parameters
- `habitat`: Habitat quality matrix (required)
- `cell_size`: Spatial resolution (default: 1.0)
- `boundary_condition`: Edge behavior (1=reflective, 2=absorbing, 3=periodic)

### Individual Parameters
- `initial_population_size`: Starting number of individuals
- `neighbor_radius`: Density dependence spatial scale
- `step_length`: Movement distance per dispersal event
- `base_birth_rate`, `base_mortality_rate`: Demographic rates

### Genetic Parameters
- `genotype_means`: Optimal habitat values for each individual
- `genotype_sds`: Niche width (habitat tolerance)
- `mutation_rates`: Evolution rate per generation
- `plasticities`: Phenotypic plasticity
- `sampling_points`: Habitat selection behavior intensity

### Simulation Parameters
- `max_events`: Maximum number of events to simulate
- `neutral_mode`: Disable habitat selection (for null models)

## Examples and Tutorials

See the `vignettes/` directory for:
- Basic workflow examples
- Genetic evolution demonstrations
- Landscape effect comparisons
- Parameter sensitivity analysis

## Development

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

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Submit a pull request

## Citation

If you use TWoLife in your research, please cite:

```
[Citation will be added upon publication]
```

## License

[License information to be added]

## Contact

- **Issues**: [GitHub Issues](https://github.com/DeFreitasLR/TWoLife-0.0.0/issues)
- **Maintainer**: DeFreitasLR

## Status

Note: This package is under active development for CRAN submission.
