#include "landscape.hpp"
#include "individual.hpp"
#include "random_utils.hpp"
#include <stdexcept>
#include <algorithm>
#include <Rcpp.h>

// FIXED: Constructor with proper memory management and validation - updated for rectangular landscapes
Landscape::Landscape(double neighbor_radius,
                     int initial_population_size,
                     double vision_angle,
                     double step_length,
                     double base_dispersal_rate,
                     double base_birth_rate,
                     double base_mortality_rate,
                     double birth_density_slope,
                     double mortality_density_slope,
                     const Rcpp::NumericMatrix& habitat,
                     double cell_size,
                     int density_type,
                     double matrix_mortality_multiplier,
                     double matrix_dispersal_multiplier,
                     int initial_placement_mode,
                     int boundary_condition,
                     const vector<double>& initial_x_coordinates,
                     const vector<double>& initial_y_coordinates,
                     const vector<double>& genotype_means,
                     const vector<double>& genotype_sds,
                     const vector<double>& mutation_rates,
                     const vector<double>& plasticities,
                     const vector<int>& sampling_points,
                     const vector<double>& habitat_selection_temperatures,
                     bool neutral_mode):
  // Initialize const members - updated for rectangular landscapes
  world_width(habitat.ncol() * cell_size),
  world_height(habitat.nrow() * cell_size),
  initial_population_size(initial_population_size),
  cells_per_row(habitat.nrow()),
  cells_per_col(habitat.ncol()),
  cell_size(cell_size),
  boundary_condition(boundary_condition),
  initial_placement_mode(initial_placement_mode),
  num_habitat_patches(0),
  // FIXED: Initialize vectors with proper rectangular size
  habitat_grid(habitat.nrow(), vector<double>(habitat.ncol(), 0.0)),
  patch_grid(habitat.nrow(), vector<int>(habitat.ncol(), 0)),
  world_time(0)
{
  // Input validation - removed square matrix requirement
  if (habitat.nrow() <= 0 || habitat.ncol() <= 0) {
    throw invalid_argument("habitat matrix dimensions must be positive");
  }
  if (cell_size <= 0.0) {
    throw invalid_argument("cell_size must be positive");
  }
  if (initial_population_size <= 0) {
    throw invalid_argument("initial_population_size must be positive");
  }
  
  // FIXED: Safe copying of habitat grid from matrix input - updated for rectangular dimensions
  for (int i = 0; i < cells_per_row; i++) // Passes through the rows (y)
  {
    for (int j = 0; j < cells_per_col; j++) // Passes through the columns (x)
    {
      habitat_grid[i][j] = habitat(i, j); // Direct matrix access: i=row(y), j=col(x)
    }
  }
  
  // Attributes values to the landscape pixels, determining the patch they belong to
  int component = 0; // Counter variable used for identifying the patches
  // FIXED: Proper loop structure for rectangular landscapes
  for (int k = 0; k < cells_per_row; ++k) { // Passes through the rows
    for (int l = 0; l < cells_per_col; ++l) { // Passes through the columns
      if (!patch_grid[k][l] && habitat_grid[k][l]) {
        identify_habitat_patches(k, l, ++component); // Checks if each landscape pixel has already been assigned a patch identification value AND if it is a non-matrix patch
      }
    }
  }
  num_habitat_patches = component; // Stores the number of habitat patches in the landscape
  
  // FIXED: Initialize patch_areas vector with proper size
  patch_areas.resize(num_habitat_patches + 1, 0.0); // Points to a new vector object used for storing the area of all habitats and the matrix
  
  // Calculate patch areas - updated for rectangular landscapes
  for (int i = 0; i < cells_per_row; i++) { // Passes through the rows
    for (int j = 0; j < cells_per_col; j++) { // Passes through the columns
      int patch_id = patch_grid[i][j];
      if (patch_id >= 0 && patch_id < static_cast<int>(patch_areas.size())) {
        patch_areas[patch_id] += cell_size * cell_size; // Indexes the patch_areas vector by the patch identification numbers of each pixel in the landscape so that it adds the number of pixels of each patch multiplied by the area of a pixel
      }
    }
  }
  
  // Modification of the radius for the global density case - updated for rectangular world
  if(density_type == 0) { // Checks if the density type is set to zero (global)
    double world_area = world_width * world_height;
    neighbor_radius = sqrt(world_area / M_PI); // Changes the perception radius to the equivalent of a radius of a circle with the same area as the rectangular landscape
  }
  
  // FIXED: Exception-safe population initialization
  try {
    initialize_population(neighbor_radius, initial_population_size, vision_angle, step_length,
                          base_dispersal_rate, base_birth_rate, base_mortality_rate,
                          birth_density_slope, mortality_density_slope,
                          matrix_mortality_multiplier, matrix_dispersal_multiplier, density_type,
                          initial_x_coordinates, initial_y_coordinates, genotype_means,
                          genotype_sds, mutation_rates, plasticities, sampling_points, 
                          habitat_selection_temperatures, neutral_mode);
    
    // Update all individuals
    for (auto& individual : population) // Passes through each individual
    {
      update_neighbors(individual.get()); // Calls a function of the Individual class to update the neighbours of an individual (the individuals within a radius distance of the focal individual)
      update_habitat_value(individual.get()); // Calls a function of the Individual class to update the habitat type the individual is currently on
      update_patch_id(individual.get()); // Calls a function of the Individual class to update the patch identifier of the patch it is currently on
    }
    
    for (auto& individual : population) // Passes through each individual
    {
      double local_density = compute_local_density(individual.get()); // Calls a function of the Individual class to determine the density of individuals within a radius distance area from a focal individual
      individual->update_rates(local_density);   // Calls a function of the Individual class to update the characteristics of each individual
    }
    
  } catch (...) {
    // If initialization fails, population vector will automatically clean up
    // any successfully created individuals due to unique_ptr
    throw;  // Re-throw the exception
  }
}

// FIXED: Population initialization with smart pointers and validation - updated coordinate ranges for rectangular world
void Landscape::initialize_population(
    const double neighbor_radius,
    const int initial_population_size,
    const double vision_angle,
    const double step_length,
    const double base_dispersal_rate,
    const double base_birth_rate,
    const double base_mortality_rate,
    const double birth_density_slope,
    const double mortality_density_slope,
    const double matrix_mortality_multiplier,
    const double matrix_dispersal_multiplier,
    const int density_type,
    const vector<double>& initial_x_coordinates,
    const vector<double>& initial_y_coordinates,
    const vector<double>& genotype_means,
    const vector<double>& genotype_sds,
    const vector<double>& mutation_rates,
    const vector<double>& plasticities,
    const vector<int>& sampling_points,
    const vector<double>& habitat_selection_temperatures,
    bool neutral_mode)
{
  Individual::reset_id(); // Restarts the id counter of the individuals
  
  // FIXED: Reserve space to avoid reallocations during population growth
  population.reserve(initial_population_size * 2);  // Allow for growth
  
  // CORRECTED: Validate array sizes UNCONDITIONALLY (not dependent on neutral_mode)
  if (mutation_rates.size() != static_cast<size_t>(initial_population_size)) {
    throw invalid_argument("Mutation rates array size mismatch");
  }
  if (plasticities.size() != static_cast<size_t>(initial_population_size)) {
    throw invalid_argument("Plasticities array size mismatch");
  }
  if (sampling_points.size() != static_cast<size_t>(initial_population_size)) {
    throw invalid_argument("Sampling points array size mismatch");
  }
  if (habitat_selection_temperatures.size() != static_cast<size_t>(initial_population_size)) {
    throw invalid_argument("Habitat selection temperatures array size mismatch");
  }
  
  // Validate neutral mode specific sizes
  if (neutral_mode) {
    if (genotype_means.size() != static_cast<size_t>(initial_population_size) ||
        genotype_sds.size() != static_cast<size_t>(initial_population_size)) {
      throw invalid_argument("Genotype arrays size mismatch in neutral mode");
    }
  }
  
  vector<double> genotype, width; // Creates vectors for storing either a genotype/width value (normal case), or all of them (null case)
  
  if (!neutral_mode) { // Checks if this is not a Null simulation run
    genotype.resize(1); // Resizes the genotype vectors for containing only a value
    width.resize(1); // Resizes the width vectors for containing only a value
  } else{ // Checks if this is a Null simulation run
    genotype.reserve(initial_population_size);
    width.reserve(initial_population_size);
    for (int i = 0; i < initial_population_size; i++) { // Pass through each individual
      genotype.push_back(genotype_means[i]); // Stores each genotype value of the population in the vector
      width.push_back(genotype_sds[i]); // Stores each width value of the population in the vector
    }
  }
  
  for (int i = 0; i < initial_population_size; i++) // Goes through the amount of initial individuals selected
  {
    // CORRECTED: Validate array access including habitat_selection_temperatures
    if (i >= static_cast<int>(mutation_rates.size()) || 
        i >= static_cast<int>(plasticities.size()) ||
        i >= static_cast<int>(sampling_points.size()) ||
        i >= static_cast<int>(habitat_selection_temperatures.size())) {
      throw out_of_range("Parameter array index out of bounds");
    }
    
    if (!neutral_mode) { // Checks if this is not a Null simulation run
      if (i >= static_cast<int>(genotype_means.size()) || 
          i >= static_cast<int>(genotype_sds.size())) {
        throw out_of_range("Genotype array index out of bounds");
      }
      genotype[0] = genotype_means[i]; // Stores the genotype value of the current individual
      width[0] = genotype_sds[i]; // Stores the width value of the current individual
    }
    
    double x_coord, y_coord;
    
    // FIXED: Safe coordinate initialization based on placement mode - updated for rectangular world
    switch (initial_placement_mode) {
    case 0: // Checks if the initial_placement_mode is set to 0 (origin)
      x_coord = 0.0; // X Coordinate
      y_coord = 0.0; // Y Coordinate
      break;
      
    case 1: // Checks if the initial_placement_mode is set to 1 (Random initial positions)
      x_coord = runif(-world_width/2, world_width/2); // X Coordinate - uses world width
      y_coord = runif(-world_height/2, world_height/2); // Y Coordinate - uses world height
      break;
      
    case 2: // Checks if the initial_placement_mode is set to 2, Random initial positions with normal distribution
      x_coord = rnorm(0, sqrt(base_dispersal_rate) * step_length); // X Coordinate
      y_coord = rnorm(0, sqrt(base_dispersal_rate) * step_length); // Y Coordinate
      break;
      
    case 3: // Checks if the initial_placement_mode is set to 3, Sets given initial positions to individuals
      if (i >= static_cast<int>(initial_x_coordinates.size()) || 
          i >= static_cast<int>(initial_y_coordinates.size())) {
        throw out_of_range("Initial coordinates array index out of bounds");
      }
      x_coord = initial_x_coordinates[i]; // X Coordinate
      y_coord = initial_y_coordinates[i]; // Y Coordinate
      break;
      
    default:
      throw invalid_argument("Invalid initial placement mode");
    }
    
    // FIXED: Use make_unique for exception safety
    try {
      auto individual = make_unique<Individual>(
        x_coord, y_coord, // X and Y Coordinates
        0, // Species ID
        base_mortality_rate, // The basal death rate (The rate at which individuals die on a habitat patch without neighbours)
        runif(0, 2.0 * M_PI), // Angle that an individual is currently facing
        vision_angle, // Angle used for orientation when dispersing
        step_length, // The Length distance of a dispersal event
        base_dispersal_rate, // The rate at which individuals disperse
        neighbor_radius, // Density dependence radius
        base_birth_rate, // The basal birth rate (The rate at which individuals give birth on a habitat patch without neighbours)
        birth_density_slope, // The slope of the birth density dependence function
        mortality_density_slope, // The slope of the death density dependence function
        matrix_mortality_multiplier, // Constant that indicates how many times higher the death rate should be on non-habitat pixels
        matrix_dispersal_multiplier, // Constant that indicates how many times lower the movement rate should be on non-habitat pixels
        density_type, // Density type (0 = global, 1 = local/within an individual radius)
        mutation_rates[i], // Mutation rate
                      plasticities[i], // Plasticity
                                  genotype, // Individual's genotype trait mean
                                  width, // Individual's genotype trait width
                                  sampling_points[i], // Number of sampling points
                                                 habitat_selection_temperatures[i] // Habitat selection temperature
      );
      
      population.push_back(move(individual));
      
    } catch (const exception& e) {
      // If individual creation fails, already created individuals
      // are automatically cleaned up by unique_ptr destructors
      throw runtime_error(string("Failed to create individual ") + to_string(i) + ": " + e.what());
    }
  }
}

// Updater function
void Landscape::update(int action_code, int individual_index)
{
  if(population.size() > 0) // checks if there are any currently alive individuals
  {
    switch (action_code) // Switches between possible actions (0= death, 1= birth, 2= dispersal)
    {
    case 0: //0= death
      break;
    case 1: //1= birth
      update_habitat_value(population[population.size()-1].get()); // updates in the individual the environmental value of the pixel corresponding to its current coordinate
      update_patch_id(population[population.size()-1].get()); // updates the individual the patch identification value of the pixel corresponding to its current coordinate
      break;
    case 2: // 2= dispersal
      update_habitat_value(population[individual_index].get()); // updates the individual the environmental value of the pixel corresponding to its current coordinate
      update_patch_id(population[individual_index].get()); // updates the individual the patch identification value of the pixel corresponding to its current coordinate
      break;
    }
    
    for(auto& individual : population) //Goes through all currently alive individuals
    {
      double local_density = compute_local_density(individual.get()); // returns the density of individuals within a radius distance area from each individual
      individual->update_rates(local_density);   // Updates the rates of each individual
    }
  }
}

// Function that drafts an individual
int Landscape::select_next_individual()
{
  if (population.empty()) {
    throw runtime_error("No individuals available for selection");
  }
  
  // time for next event and simulation time update
  int selected_individual_index = 0; // Stores the first individual of the vector as the one with currently the least required time
  double min_event_time = population[0]->get_time_to_next_event(); // Drafts a time from an exponential equation used as the needed time to execute an action
  
  for(size_t i = 1; i < population.size(); i++) //Goes through all individuals
  {
    if(population[i]->get_time_to_next_event() < min_event_time) // Drafts a time for each individual and checks if it is lower than the current lower one
    {
      selected_individual_index = static_cast<int>(i); // Stores the current smaller individual
      min_event_time = population[i]->get_time_to_next_event(); // Stores the current lower time
    }
  }
  return selected_individual_index; // Returns the individual with the lowest time
}

//FIXED: Safe action performance with smart pointers and validation
bool Landscape::perform_action(int action_code, int individual_index)
{
  // Input validation
  if (individual_index < 0 || individual_index >= static_cast<int>(population.size())) {
    throw out_of_range("Individual index out of bounds");
  }
  
  bool left_boundary = false; // Initializes a boolean variable used to record if an individual has left the boundaries of the landscape as false
  
  switch(action_code) //(0= death, 1= birth, 2= dispersal)
  {
  case 0: //0= death
    update_neighbor_lists(action_code, individual_index);
    // FIXED: Smart pointer automatically handles deletion
    population.erase(population.begin() + individual_index); // Erases the position that object occupied in the vector of individuals
    break;
    
  case 1: // 1= birth
    try {
      // FIXED: Use smart pointer for new individual
      auto newborn = make_unique<Individual>(*population[individual_index]); // Calls a function to create a new individual with similar characteristics to the parent one
      population.push_back(move(newborn)); // Stores the new individual at the position of the landscape individuals vector
      update_neighbor_lists(action_code, individual_index);
    } catch (const exception& e) {
      throw runtime_error(string("Failed to create offspring: ") + e.what());
    }
    break;
    
  case 2: // 2= dispersal
    left_boundary = perform_movement(individual_index); // Performs movement action and returns if the individual left the landscape boundaries
    
    if(left_boundary) //checks if an individual left the landscape boundaries
    {
      update_neighbor_lists(0, individual_index);  // Death action - Treats an emigrating individual as if it died
      population.erase(population.begin() + individual_index);
    }
    else if(!left_boundary) //checks if an individual didn't leave the landscape boundaries
      update_neighbor_lists(action_code, individual_index); //Apply the boundary condition if the individual disperses out of the landscape boundary
    break;
    
  default:
    throw invalid_argument("Invalid action code");
  }
  return left_boundary; // Returns true if an individual left the landscape boundaries
}

// FIXED: Safe boundary condition application with rectangular logic
bool Landscape::apply_boundary_condition(Individual* individual)
{
  bool left_boundary = false; // Initializes a boolean variable used to record if an individual has left the boundaries of the landscape as false
  double half_width = world_width / 2;
  double half_height = world_height / 2;
  
  switch(boundary_condition) // Switches between conditions (0= absorbing, 1= periodic (pacman), 2= reflective)
  {
    
  case 0: //0= absorbing - updated for rectangular boundaries
    if(individual->get_x_coordinate() >= half_width || //Checks if the x coordinate is higher than the maximum value OR
       individual->get_x_coordinate() < -half_width || // if the x coordinate is lower than the minimum value OR
       individual->get_y_coordinate() >= half_height || // if the y coordinate is higher than the maximum value Or
       individual->get_y_coordinate() < -half_height) // if the y coordinate is lower than the minimum value
    {
      left_boundary = true; // Sets emigration to true
    }
    break;
    
  case 1: // 1= periodic - updated for rectangular boundaries
    if(individual->get_x_coordinate() < -half_width) //Checks if the x coordinate is lower than the minimum value
      individual->set_x_coordinate(world_width + individual->get_x_coordinate()); // Changes the x coordinate to represent the periodic condition (left edge to right edge)
    if(individual->get_x_coordinate() >= half_width) //Checks if the x coordinate is higher than the maximum value
      individual->set_x_coordinate(individual->get_x_coordinate() - world_width); // Changes the x coordinate to represent the periodic condition (right edge to left edge)
    if(individual->get_y_coordinate() < -half_height) //Checks if the y coordinate is lower than the minimum value
      individual->set_y_coordinate(world_height + individual->get_y_coordinate()); // Changes the y coordinate to represent the periodic condition (bottom edge to top edge)
    if(individual->get_y_coordinate() >= half_height) //Checks if the y coordinate is higher than the maximum value
      individual->set_y_coordinate(individual->get_y_coordinate() - world_height); // Changes the y coordinate to represent the periodic condition (top edge to bottom edge)
    break;
    
  case 2: // 2= reflective - updated for rectangular boundaries
    if(individual->get_x_coordinate() < -half_width) //Checks if the x coordinate is lower than the minimum value
      individual->set_x_coordinate(-half_width + abs(-half_width - individual->get_x_coordinate())); // Sets the x coordinate to the boundary plus its excess length (reflected)
    if(individual->get_x_coordinate() >= half_width) //Checks if the x coordinate is higher than the maximum value
      individual->set_x_coordinate(half_width - abs(half_width - individual->get_x_coordinate()));// Sets the x coordinate to the boundary minus its excess length (reflected)
    if(individual->get_y_coordinate() < -half_height) //Checks if the y coordinate is lower than the minimum value
      individual->set_y_coordinate(-half_height + abs(-half_height - individual->get_y_coordinate()));// Sets the y coordinate to the boundary plus its excess length (reflected)
    if(individual->get_y_coordinate() >= half_height) //Checks if the y coordinate is higher than the maximum value
      individual->set_y_coordinate(half_height - abs(half_height - individual->get_y_coordinate()));// Sets the y coordinate to the boundary minus its excess length (reflected)
    break;
  }
  return left_boundary; // Returns a positive boolean value if the individual emigrated
}

// Function that calculates the distance between two individuals - updated for rectangular periodic boundaries
double Landscape::compute_distance(const Individual* individual_a, const Individual* individual_b) const
{
  switch(boundary_condition) //switches between boundary condition cases (0= absorbing, 1= periodic (pacman), 2= reflective)
  {
  case 0:              // absorbing Euclidean distance
  case 2:             // reflective Euclidean distance
    return sqrt(pow(individual_a->get_x_coordinate() - individual_b->get_x_coordinate(), 2) + 
                pow(individual_a->get_y_coordinate() - individual_b->get_y_coordinate(), 2));
    
  case 1: // The Euclidean distance accounting for the periodic effect - updated for rectangular boundaries
    {
      double x1 = individual_a->get_x_coordinate(); // Gets the x coordinate of the first individual
      double x2 = individual_b->get_x_coordinate(); // Gets the x coordinate of the second individual
      double y1 = individual_a->get_y_coordinate(); // Gets the y coordinate of the first individual
      double y2 = individual_b->get_y_coordinate(); // Gets the y coordinate of the second individual
      double dx = abs(x1 - x2); // chooses the lower between x1-x2 and x2-x1
      double dy = abs(y1 - y2); //chooses the lower between y1-y2 and y2-y1
      dx = min(dx, world_width - dx); // chooses the lower between world_width-dx and dx
      dy = min(dy, world_height - dy); //chooses the lower between world_height-dy and dy
      return sqrt(dx * dx + dy * dy); // Computes and returns the distance
    }
    
  default: // defensive fallback Euclidean distance
    return sqrt(pow(individual_a->get_x_coordinate() - individual_b->get_x_coordinate(), 2) + 
                pow(individual_a->get_y_coordinate() - individual_b->get_y_coordinate(), 2));
  }
}

// Function that updates an individual's neighbourhood list
void Landscape::update_neighbors(Individual* individual) const
{
  vector<Individual*> neighbor_list; // Creates a vector of individuals to store the neighbours
  if(individual->get_density_type() == 0) // Checks if the density type is set to 0 (global)
  {
    for(auto& other_individual : population) //Goes through all the individuals
    {
      if(individual == other_individual.get()) continue; //Except if it is the focal individual itself
      neighbor_list.push_back(other_individual.get()); // Stores the copy of each other individual in a neighbours vector
    }
  }
  else //Checks if the density type is not set to 0 (1= local)
  {
    double radius = individual->get_neighbor_radius(); // Gets the radius of the focal individual
    for(auto& other_individual : population) //Goes through all currently alive individuals
    {
      if(individual == other_individual.get()) continue; //Except if it is the focal individual itself
      double distance = compute_distance(individual, other_individual.get()); // Checks if the copied individual is within a radius distance
      if(distance <= radius) {neighbor_list.push_back(other_individual.get());} // Stores the copy of each other individual within a radius distance into a neighbours vector
    }
  }
  individual->set_neighbors(neighbor_list); // Passes to an individual its neighbourhood list
}

// FIXED: Safe neighbor list updates
void Landscape::update_neighbor_lists(int action_code, int individual_index) const
{
  if (individual_index < 0 || individual_index >= static_cast<int>(population.size())) {
    throw out_of_range("Individual index out of bounds in update_neighbor_lists");
  }
  
  switch(action_code) // Switches between actions (0= death, 1= birth, 2= dispersal)
  {
  case 0: // Performs the neighbourhood update in case of a death event
  {
    uint64_t dying_id = population[individual_index]->get_id();
    for(auto& individual : population) { // Passes through each of the drafted individual's neighbours
      auto& neighbors = individual->neighbors;
      neighbors.erase( // Erases the drafted individual from the neighbourhood of its neighbours
        remove_if(neighbors.begin(), neighbors.end(),
                  [dying_id](Individual* neighbor) { 
                    return neighbor->get_id() == dying_id; 
                  }),
                  neighbors.end()
      );
    }
  }
    break; // Breaks switch
    
  case 1: // Performs the neighbourhood update in case of a birth event
    if (!population.empty()) {
      Individual* newborn = population.back().get();
      update_neighbors(newborn); // Calls a function that stores the last born individual's neighbours in its neighbourhood
      
      // Add newborn to neighbors' lists if within range
      for(auto& individual : population) { // Passes through each of the drafted individual's neighbours
        if (individual.get() == newborn) continue;
        
        if (individual->get_density_type() == 0 || 
            compute_distance(individual.get(), newborn) <= individual->get_neighbor_radius()) {
          individual->neighbors.push_back(newborn); // Adds the last born individual to its neighbours' neighbourhood
        }
      }
    }
    break; // Breaks switch
    
  case 2: // Performs the neighbourhood update in case of a migration event
    if (population[individual_index]->get_density_type() != 0) { // Checks if the density type is not set to global
      // Remove from old neighbors
      uint64_t moving_id = population[individual_index]->get_id();
      for(auto& individual : population) { // Passes through each of the drafted individual's neighbours
        auto& neighbors = individual->neighbors;
        neighbors.erase( // Erases the drafted individual from the neighbourhood of its neighbours
          remove_if(neighbors.begin(), neighbors.end(),
                    [moving_id](Individual* neighbor) { 
                      return neighbor->get_id() == moving_id; 
                    }),
                    neighbors.end()
        );
      }
      
      // Rebuild neighbor lists
      update_neighbors(population[individual_index].get()); // Calls a function that stores the neighbours in the individual's neighbourhood
      
      // Add to new neighbors
      Individual* moved_individual = population[individual_index].get();
      for(auto& individual : population) { // Passes through each living individual
        if (individual.get() == moved_individual) continue;
        
        if (compute_distance(individual.get(), moved_individual) <= individual->get_neighbor_radius()) { // Checks if the compared individual is within a radius distance
          individual->neighbors.push_back(moved_individual); // Stores the drafted individual on its neighbours' neighbourhood
        }
      }
    }
    break; // Breaks switch
  }
}

// FIXED: Safe grid access functions - updated for rectangular grid coordinates
void Landscape::update_habitat_value(Individual* individual) const {
  int grid_x = static_cast<int>(individual->get_x_coordinate() / cell_size + cells_per_col / 2); // Finds the pixel X coordinate (column)
  int grid_y = static_cast<int>((-individual->get_y_coordinate() / cell_size) + cells_per_row / 2); //Finds the pixel y coordinate (row)
  
  // FIXED: Bounds checking before array access
  if (is_valid_grid_position(grid_x, grid_y)) {
    individual->set_habitat_value(habitat_grid[grid_y][grid_x]); // FIXED: [row][col] = [grid_y][grid_x]
  } else {
    // Individual is outside grid bounds - treat as matrix
    individual->set_habitat_value(0.0);
  }
}

// Function to update the patch the individual is currently in - updated for rectangular coordinates
void Landscape::update_patch_id(Individual* individual) const {
  int grid_x = static_cast<int>(individual->get_x_coordinate() / cell_size + cells_per_col / 2); // Finds the pixel X coordinate (column)
  int grid_y = static_cast<int>((-individual->get_y_coordinate() / cell_size) + cells_per_row / 2); //Finds the pixel y coordinate (row)
  
  // FIXED: Bounds checking before array access  
  if (is_valid_grid_position(grid_x, grid_y)) {
    individual->set_patch_id(patch_grid[grid_y][grid_x]); // FIXED: [row][col] = [grid_y][grid_x]
  } else {
    // Individual is outside grid bounds
    individual->set_patch_id(0);
  }
}

// Recursive function that finds and identifies each landscape patch - updated for rectangular boundaries
void Landscape::identify_habitat_patches(int x, int y, int current_label)
{
  if (x < 0 || x >= cells_per_row) return; // Checks if the row coordinate is out of bounds
  if (y < 0 || y >= cells_per_col) return;// Checks if the column coordinate is out of bounds
  if (patch_grid[x][y] || !habitat_grid[x][y]) return; // Checks if the current pixel is already labeled AND if is a habitat pixel
  
  // marks the current cell with its respective patch ID value
  patch_grid[x][y] = current_label;
  
  // recursively mark the neighbors
  identify_habitat_patches(x + 1, y, current_label);
  identify_habitat_patches(x, y + 1, current_label);
  identify_habitat_patches(x - 1, y, current_label);
  identify_habitat_patches(x, y - 1, current_label);
  
  // 8-rule
  identify_habitat_patches(x + 1, y + 1, current_label);
  identify_habitat_patches(x - 1, y + 1, current_label);
  identify_habitat_patches(x - 1, y - 1, current_label);
  identify_habitat_patches(x + 1, y - 1, current_label);
}

// A function to calculate the density of individuals according to density type (global or local) and considering landscape boundary effects in the calculation - updated for rectangular boundaries
double Landscape::compute_local_density(const Individual* focal_individual) const
{
  double density; // Creates a temporary variable for storing density
  double radius = focal_individual->get_neighbor_radius();
  double circle_area = M_PI * radius * radius;
  double world_area = world_width * world_height; // Updated for rectangular world
  
  if(circle_area > world_area){ // check if the individual radius would be higher than the landscape
    density = focal_individual->neighborhood_size() / world_area; // Computes the standard case of density using the total area of the world
  } else {
    density = focal_individual->neighborhood_size() / circle_area; // Computes the standard case of density (when an individual radius area is completely within the landscape)
    
    // Functions for local density calculation
    
    /* 1. Circular area defining a region in which density-dependence occurs: landscape boundary effects.
     In this case, density is the number of individuals inside the circle divided by circle area.
     This is the same calculation as for global density, except by the cases in which landscape boundary affects
     the area of the circle.
     */
    
    // Condition giving the boundary effects cases - updated for rectangular boundaries
    if(focal_individual->get_density_type() == 1) // checks if the density type is set to 1 (local/ within a radius distance)
    {
      double abs_x = abs(focal_individual->get_x_coordinate()); // Stores the absolute value of x
      double abs_y = abs(focal_individual->get_y_coordinate()); // Stores the absolute value of y
      double half_width = world_width / 2; // Stores half the world width
      double half_height = world_height / 2; // Stores half the world height
      
      if(abs_x > half_width - radius || abs_y > half_height - radius) // Checks if the x coordinate is within a radius distance from the edge OR if the y coordinate is within a radius distance from the edge
      { // OBS: The absolute values retain their distance from the edge, so the computations from all quadrants can be done with a single set of equations
        
        // temporary objects
        vector<double>section_x;
        vector<double>section_y;
        
        // Functions for adjusted local density calculation, according to the specific case
        // Case 1) - Individual near vertical edge (left/right boundary)
        if(abs_x > half_width - radius && abs_y <= half_height - radius) // Checks if the x coordinate is within a radius distance from the edge AND if the y coordinate is not
        {
          section_x.push_back(half_width);
          section_x.push_back(half_width);
          section_y.push_back(abs_y + sqrt(radius * radius - ((half_width - abs_x) * (half_width - abs_x))));
          section_y.push_back(abs_y - sqrt(radius * radius - ((half_width - abs_x) * (half_width - abs_x))));
          
          double distance_sections = section_y[0] - section_y[1];
          double height = half_width - abs_x;
          double theta = acos(1 - (distance_sections * distance_sections / (2 * radius * radius))); // angle in radians
          double adjusted_area = M_PI * radius * radius - (theta * radius * radius / 2 - (distance_sections * height / 2));
          
          density = focal_individual->neighborhood_size() / adjusted_area;
        }
        
        // Case 2) - Individual near horizontal edge (top/bottom boundary)
        if(abs_x <= half_width - radius && abs_y > half_height - radius) // Checks if the y coordinate is within a radius distance from the edge AND if the x coordinate is not
        {
          section_y.push_back(half_height);
          section_y.push_back(half_height);
          section_x.push_back(abs_x + sqrt(radius * radius - ((half_height - abs_y) * (half_height - abs_y))));
          section_x.push_back(abs_x - sqrt(radius * radius - ((half_height - abs_y) * (half_height - abs_y))));
          
          double distance_sections = section_x[0] - section_x[1];
          double height = half_height - abs_y;
          double theta = acos(1 - (distance_sections * distance_sections / (2 * radius * radius))); // angle in radians
          double adjusted_area = M_PI * radius * radius - (theta * radius * radius / 2 - (distance_sections * height / 2));
          
          density = focal_individual->neighborhood_size() / adjusted_area;
        }
        
        // Case 3) & 4) - Individual near corner
        if(abs_x > half_width - radius && abs_y > half_height - radius)
        {
          
          // Case 3) - Individual far from corner
          if((abs_x - half_width) * (abs_x - half_width) + (abs_y - half_height) * (abs_y - half_height) > radius * radius)
          {
            section_x.push_back(abs_x + sqrt(radius * radius - ((half_height - abs_y) * (half_height - abs_y))));
            section_y.push_back(half_height);
            section_x.push_back(abs_x - sqrt(radius * radius - ((half_height - abs_y) * (half_height - abs_y))));
            section_y.push_back(half_height);
            section_x.push_back(half_width);
            section_y.push_back(abs_y + sqrt(radius * radius - ((half_width - abs_x) * (half_width - abs_x))));
            section_x.push_back(half_width);
            section_y.push_back(abs_y - sqrt(radius * radius - ((half_width - abs_x) * (half_width - abs_x))));
            
            double distance_sections = sqrt((section_x[3] - section_x[1]) * (section_x[3] - section_x[1]) + (section_y[1] - section_y[3]) * (section_y[1] - section_y[3]));
            double distance_sections2 = sqrt((section_x[2] - section_x[0]) * (section_x[2] - section_x[0]) + (section_y[0] - section_y[2]) * (section_y[0] - section_y[2]));
            double theta = acos(1 - (distance_sections * distance_sections / (2 * radius * radius))); // angle in radians
            double phi = acos(1 - (distance_sections2 * distance_sections2 / (2 * radius * radius))); // angle in radians
            double adjusted_area = (2 * M_PI - theta) * radius * radius / 2 + phi * radius * radius / 2 + (section_x[0] - section_x[1]) * (half_height - abs_y) / 2 + (section_y[2] - section_y[3]) * (half_width - abs_x) / 2;
            
            density = focal_individual->neighborhood_size() / adjusted_area;
          }
          // Case 4) - Individual close to corner
          else
          {
            section_x.push_back(abs_x - sqrt(radius * radius - ((half_height - abs_y) * (half_height - abs_y))));
            section_y.push_back(half_height);
            section_x.push_back(half_width);
            section_y.push_back(abs_y - sqrt(radius * radius - ((half_width - abs_x) * (half_width - abs_x))));
            
            double distance_sections = sqrt((section_x[1] - section_x[0]) * (section_x[1] - section_x[0]) + (section_y[0] - section_y[1]) * (section_y[0] - section_y[1]));
            double theta = acos(1 - (distance_sections * distance_sections / (2 * radius * radius))); // angle in radians
            double adjusted_area = theta * radius * radius / 2 + (half_width - section_x[0]) * (half_height - abs_y) + (half_height - section_y[1]) * (half_width - abs_x);
            
            density = focal_individual->neighborhood_size() / adjusted_area;
          }
        }
      }
    }
  }
  
  return density; // Returns the individual's current density
}

// FIXED: Safe movement function with proper container usage - updated for rectangular boundaries
bool Landscape::perform_movement(int individual_index){
  
  if (individual_index < 0 || individual_index >= static_cast<int>(population.size())) {
    throw out_of_range("Individual index out of bounds in perform_movement");
  }
  
  bool left_boundary = false; // Initializes boolean variable that marks if the individual emigrated
  Individual* individual = population[individual_index].get();
  double half_width = world_width / 2;
  double half_height = world_height / 2;
  
  if(individual->get_sampling_points() > 0) // Checks if the individual samples their step radius
  {
    int num_points = individual->get_sampling_points();
    
    // FIXED: Use vector instead of C-style array
    vector<vector<double>> candidate_locations(num_points, vector<double>(3)); // Creates a vector for storing the possible migration locations
    
    double sampled_angle, sampled_distance; // Temporary variables for storing angle and length
    int grid_x, grid_y; // Temporary variables for storing the possible coordinates
    
    for (int i = 0; i < num_points; i++) // Passes by each possible location
    {
      
      sampled_angle = runif(0.0, 2.0 * M_PI); // Samples an angle directly in radians
      sampled_distance = runif(0.0, individual->get_step_length()); // Samples a distance
      
      candidate_locations[i][0] = individual->get_x_coordinate() + cos(sampled_angle) * sampled_distance; // Stores the possible x
      candidate_locations[i][1] = individual->get_y_coordinate() + sin(sampled_angle) * sampled_distance; // Stores the possible y
      
      switch(boundary_condition){ // switches between boundary conditions (0= absorption, 1= periodic (pacman), 2= reflective)
      
      case 0: // 0= absorption - updated for rectangular boundaries
        
        if(candidate_locations[i][0] < -half_width || //Checks if the x coordinate is lower than the minimum value OR
           candidate_locations[i][0] >= half_width || // if the x coordinate is higher than the maximum value OR
           candidate_locations[i][1] < -half_height || // if the y coordinate is lower than the minimum value Or
           candidate_locations[i][1] >= half_height) // if the y coordinate is higher than the maximum value
        {
          
          candidate_locations[i][2] = -10000; // Sets the out of boundaries locations rank to a very small number
          
          break;
          
        }
        // Fall through to get habitat value if within bounds
        
      case 1: // 1= periodic (pacman) - updated for rectangular boundaries
        if(candidate_locations[i][0] < -half_width) //Checks if the x coordinate is lower than the minimum value
          candidate_locations[i][0] = (world_width + candidate_locations[i][0]); // Transports coordinate to the other edge
        if(candidate_locations[i][0] >= half_width) // if the x coordinate is higher than the maximum value
          candidate_locations[i][0] = (candidate_locations[i][0] - world_width); // Transports coordinate to the other edge
        if(candidate_locations[i][1] < -half_height)// if the y coordinate is lower than the minimum value
          candidate_locations[i][1] = (world_height + candidate_locations[i][1]); // Transports coordinate to the other edge
        if(candidate_locations[i][1] >= half_height ) // if the y coordinate is higher than the maximum value
          candidate_locations[i][1] = (candidate_locations[i][1] - world_height); // Transports coordinate to the other edge
        
        grid_x = static_cast<int>(candidate_locations[i][0] / cell_size + cells_per_col / 2); // Discovers the x pixel index of the coordinate (column)
        grid_y = static_cast<int>((-candidate_locations[i][1] / cell_size) + cells_per_row / 2); // Discovers the y pixel index of the coordinate (row)
        
        if (is_valid_grid_position(grid_x, grid_y)) {
          candidate_locations[i][2] = habitat_grid[grid_y][grid_x]; // FIXED: [row][col] = [grid_y][grid_x]
        } else {
          candidate_locations[i][2] = 0.0; // Default to matrix if out of bounds
        }
        
        break;
        
      case 2:// 2= reflective - updated for rectangular boundaries
        // Prevents out of boundary coordinates from being sampled
        while (candidate_locations[i][0] < -half_width || //Checks if the x coordinate is lower than the minimum value OR
               candidate_locations[i][0] >= half_width || // if the x coordinate is higher than the maximum value OR
               candidate_locations[i][1] < -half_height || // if the y coordinate is lower than the minimum value Or
               candidate_locations[i][1] >= half_height) { // if the y coordinate is higher than the maximum value
          
          sampled_angle = runif(0.0, 2.0 * M_PI); // Samples an angle directly in radians
          sampled_distance = runif(0.0, individual->get_step_length()); // re-Samples a distance
          
          candidate_locations[i][0] = individual->get_x_coordinate() + cos(sampled_angle) * sampled_distance;// Stores the possible x
          candidate_locations[i][1] = individual->get_y_coordinate() + sin(sampled_angle) * sampled_distance;// Stores the possible y
        }
        grid_x = static_cast<int>(candidate_locations[i][0] / cell_size + cells_per_col / 2); // Discovers the x pixel index of the coordinate (column)
        grid_y = static_cast<int>((-candidate_locations[i][1] / cell_size) + cells_per_row / 2); // Discovers the y pixel index of the coordinate (row)
        
        if (is_valid_grid_position(grid_x, grid_y)) {
          candidate_locations[i][2] = habitat_grid[grid_y][grid_x]; // FIXED: [row][col] = [grid_y][grid_x]
        } else {
          candidate_locations[i][2] = 0.0; // Default to matrix if out of bounds
        }
        
        break;
      }
      
    }
    
    // FIXED: Use safe habitat selection
    try {
      individual->select_habitat(candidate_locations);// Passes the possible locations for the individual to choose from
      
      if (boundary_condition == 0) {
        left_boundary = apply_boundary_condition(individual);
      }
    } catch (const exception& e) {
      throw runtime_error(string("Habitat selection failed: ") + e.what());
    }
    
  }
  
  else  // Checks if the individual doesn't sample their step radius (random)
  {
    individual->perform_random_walk(); // Calls a function that changes the X and Y coordinates of the individuals via random walk
    left_boundary = apply_boundary_condition(individual); // Applies boundary condition
  }
  return left_boundary;
}

// Placeholder for count_species - implementation depends on requirements
const int Landscape::count_species() const {
  // TBI - To Be Implemented
  return 1; // Placeholder
}