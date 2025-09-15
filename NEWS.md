# TWoLife 1.0.0

## Initial Release

### New Features

- **Core Simulation Engine**: Fast C++ implementation of individual-based spatial population models
- **Genetic Evolution**: Full genetic system with mutation, selection, and adaptation
- **Habitat Selection**: Behavioral habitat choice with sampling-based decisions
- **Landscape Generation**: Fractal and geometric landscape creation tools
- **Comprehensive Visualization**: Multiple plotting functions for landscapes and results
- **Flexible Parameter System**: Organized parameter categories for easy simulation setup

### Main Functions

- `twolife_simulation()`: Primary function for running complete simulations
- `create_fractal_landscape()`: Generate realistic fractal landscapes
- `create_corner_landscape()`: Create simple test landscapes with habitat in corners
- `plot_simulation_on_landscape()`: Visualize simulation results overlaid on landscapes
- `plot_landscape_image()`: Display landscape patterns with enhanced coloring
- `plot_landscape_world_coords()`: Plot landscapes using world coordinate system
- `compute_population_size()`: Extract and analyze population trajectories
- `quick_plot_result()`: Simplified result visualization

### Documentation

- **Comprehensive Vignette**: Step-by-step tutorial with biological interpretation
- **Complete Help System**: Detailed documentation for all functions
- **Example Workflows**: Progressive complexity examples from basic to advanced
- **Parameter Guides**: Detailed explanation of all simulation parameters

### Technical Features

- **Performance**: C++ backend for computational efficiency
- **Memory Safety**: Smart pointers and bounds checking throughout
- **Error Handling**: Comprehensive input validation and informative error messages
- **Reproducibility**: Proper random seed handling for consistent results
- **Cross-Platform**: Compatible with Windows, macOS, and Linux

### Simulation Capabilities

#### Spatial Features
- Multiple boundary conditions (absorbing, periodic, reflective)
- Continuous and discrete habitat types
- Flexible landscape sizes and resolutions
- World coordinate system matching internal simulation

#### Population Dynamics
- Individual-based birth, death, and dispersal events
- Density-dependent demographic rates (global or local)
- Stochastic event timing with exponential distributions
- Population tracking and trajectory analysis

#### Genetic Evolution
- Quantitative genetic traits with normal distributions
- Mutation with configurable rates per individual
- Phenotypic plasticity for environmental response
- Selection on habitat-matching fitness

#### Behavioral Ecology
- Active habitat selection via sampling multiple locations
- Movement with directional bias and vision angles
- Differential movement rates in habitat vs. matrix
- Adaptive dispersal strategies

### Parameter Organization

#### Landscape Parameters
- `habitat`: Habitat quality matrix (required)
- `cell_size`: Spatial resolution (default: 1.0)
- `boundary_condition`: Edge behavior (default: periodic)
- `density_type`: Density dependence scale (default: local)
- `matrix_mortality_multiplier`: Matrix death rate increase (default: 2.0)
- `matrix_dispersal_multiplier`: Matrix movement rate change (default: 0.5)

#### Individual Parameters
- `initial_population_size`: Starting population (default: 200)
- `neighbor_radius`: Density interaction scale (default: 2.0)
- `vision_angle`: Movement orientation range (default: π)
- `step_length`: Movement distance per event (default: 1.0)
- `base_dispersal_rate`: Movement rate (default: 0.1)
- `base_birth_rate`: Reproduction rate (default: 0.3)
- `base_mortality_rate`: Death rate (default: 0.2)
- `birth_density_slope`: Density effect on births (default: 0.02)
- `mortality_density_slope`: Density effect on deaths (default: 0.02)

#### Genetic Parameters
- `genotype_means`: Optimal habitat values per individual
- `genotype_sds`: Niche width (habitat tolerance)
- `mutation_rates`: Evolution rate per generation
- `plasticities`: Phenotypic flexibility
- `sampling_points`: Habitat selection behavior intensity

#### Simulation Parameters
- `max_events`: Simulation duration (auto-calculated by default)
- `neutral_mode`: Disable evolution (default: FALSE)
- `output_file`: Optional detailed output file

### Visualization Features

#### Landscape Plotting
- Multiple color schemes (terrain, viridis, plasma, habitat, quality, binary)
- Automatic legends for continuous and discrete data
- Grid lines and reference axes
- High-resolution export (PNG, PDF, JPEG)
- World coordinate system alignment

#### Result Visualization
- Survivor positions overlaid on landscapes
- Color-coding by genotype values
- Population statistics display
- Multiple export formats
- Customizable point sizes and colors

### Example Workflows

#### Basic Usage
```r
habitat <- create_fractal_landscape(10, fractality = 0.6, habitat_proportion = 0.4)
result <- twolife_simulation(
  landscape_params = list(habitat = habitat),
  individual_params = list(initial_population_size = 20),
  simulation_params = list(max_events = 200)
)
plot_simulation_on_landscape(result)
```

#### Advanced Genetics
```r
result <- twolife_simulation(
  landscape_params = list(habitat = habitat),
  individual_params = list(initial_population_size = 40),
  genetic_params = list(
    genotype_means = runif(40, 0.3, 0.7),
    genotype_sds = rep(0.15, 40),
    mutation_rates = rep(0.05, 40),
    sampling_points = rep(20, 40)
  )
)
```

### Quality Assurance

- **Input Validation**: Comprehensive parameter checking with informative errors
- **Memory Management**: Safe memory handling with automatic cleanup
- **Thread Safety**: Atomic operations for unique ID generation
- **Numerical Stability**: Robust handling of edge cases and boundary conditions
- **Documentation**: Complete roxygen2 documentation with examples

### Dependencies

- **R**: >= 4.0.0
- **Rcpp**: C++ integration
- **rflsgen**: Fractal landscape generation
- **viridisLite**: Enhanced color palettes (optional)

### Platform Support

- **Windows**: Full support with Rtools
- **macOS**: Compatible with Xcode command line tools
- **Linux**: Standard build tools

---

## Development Notes

This initial release represents a complete rewrite and modernization of the TWoLife simulation framework, with emphasis on:

1. **Performance**: C++ implementation for speed
2. **Usability**: Intuitive parameter organization and helper functions
3. **Reliability**: Comprehensive error checking and memory safety
4. **Documentation**: Detailed examples and biological interpretation
5. **Extensibility**: Modular design for future enhancements

The package is designed for researchers studying spatial population dynamics, evolutionary ecology, and conservation biology.