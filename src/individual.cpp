#include "individual.hpp"
#include "random_utils.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

// Initialize atomic counter
atomic<uint64_t> Individual::MAX_ID{0};

/** \brief Individual exception
 Displays warning messages for impossible values
 (\ref TBI). */
class individual_exception: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Exception occurred in Individual class!";
  }
} individual_ex;

// Constructor with comprehensive input validation - Class Constructor
Individual::Individual(double x_coordinate,
                       double y_coordinate,
                       int species_id,
                       double base_mortality_rate,
                       double current_bearing,
                       double vision_angle,
                       double step_length,
                       double base_dispersal_rate,
                       double neighbor_radius,
                       double base_birth_rate,
                       double birth_density_slope,
                       double mortality_density_slope,
                       double matrix_mortality_multiplier,
                       double matrix_dispersal_multiplier,
                       int density_type,
                       double mutation_rate,
                       double plasticity,
                       const vector<double>& genotype_means,
                       const vector<double>& genotype_sds,
                       int sampling_points,
                       double habitat_selection_temperature):
  // Thread-safe ID assignment with overflow check - This section is similar to regular attribution and is run before brackets
  id(++MAX_ID),
  x_coordinate(x_coordinate),
  y_coordinate(y_coordinate),
  species_id(species_id),
  base_mortality_rate(base_mortality_rate),
  current_bearing(current_bearing),
  vision_angle(vision_angle),
  step_length(step_length),
  base_dispersal_rate(base_dispersal_rate),
  current_dispersal_rate(base_dispersal_rate),
  neighbor_radius(neighbor_radius),
  base_birth_rate(base_birth_rate),
  birth_density_slope(birth_density_slope),
  mortality_density_slope(mortality_density_slope),
  matrix_mortality_multiplier(matrix_mortality_multiplier),
  matrix_dispersal_multiplier(matrix_dispersal_multiplier),
  density_type(density_type),
  mutation_rate(mutation_rate),
  plasticity(plasticity),
  habitat_selection_temperature(habitat_selection_temperature),
  sampling_points(sampling_points),
  // Initialize other members - Consistent order
  current_birth_rate(0.0),
  current_death_rate(0.0),
  time_to_next_event(0.0),
  max_density(0.0),
  equilibrium_rate(0.0),
  habitat_value(0.0),
  patch_id(0)
{
  // Input validation with detailed error messages - Checks for impossible values
  if (base_mortality_rate < 0) {
    throw invalid_argument("Base mortality rate cannot be negative");
  }
  if (step_length < 0) {
    throw invalid_argument("Step length cannot be negative");
  }
  if (base_dispersal_rate < 0) {
    throw invalid_argument("Base dispersal rate cannot be negative");
  }
  if (neighbor_radius < 0) {
    throw invalid_argument("Neighbor radius cannot be negative");
  }
  if (base_birth_rate < 0) {
    throw invalid_argument("Base birth rate cannot be negative");
  }
  if (sampling_points < 0) {
    throw invalid_argument("Sampling points cannot be negative");
  }
  if (habitat_selection_temperature <= 0) {
    throw invalid_argument("Habitat selection temperature must be positive");
  }
  if (genotype_means.empty()) {
    throw invalid_argument("Genotype means vector cannot be empty");
  }
  if (genotype_sds.empty()) {
    throw invalid_argument("Genotype standard deviations vector cannot be empty");
  }
  if (genotype_means.size() != genotype_sds.size()) {
    throw invalid_argument("Genotype means and sds vectors must have same size");
  }
  
  // Check for ID overflow (extremely unlikely but good practice)
  if (id == 0) {
    throw runtime_error("Individual ID counter overflow occurred");
  }
  
  // Safe copying of genotype data
  this->genotype_means = genotype_means;
  this->genotype_sds = genotype_sds;
  
  // Safe environmental optimum calculation
  environmental_optimum.reserve(genotype_means.size());
  
  try {
    for (size_t i = 0; i < genotype_means.size(); i++) { // passes by the genotype means of the individual
      double phenotype = rnorm(genotype_means[i], plasticity); // Generates phenotype by adding a random value to the genotype value
      
      // Validate the generated phenotype
      if (!isfinite(phenotype)) {
        throw runtime_error("Generated phenotype is not finite");
      }
      
      environmental_optimum.push_back(phenotype);
    }
  } catch (const std::exception& e) {
    throw runtime_error(string("Failed to initialize environmental optimum: ") + e.what());
  }
  
  // Safe density calculations with validation
  if(birth_density_slope != 0 && mortality_density_slope != 0) // Checks if the birth and death rates are density dependent
  {
    double denominator = birth_density_slope + mortality_density_slope;
    if (denominator != 0) {
      max_density = (base_birth_rate - base_mortality_rate) / denominator; // Sets the maximum density
      equilibrium_rate = base_mortality_rate + mortality_density_slope * max_density; // Sets the point where rates are equal
    } else {
      max_density = 0.0;
      equilibrium_rate = base_mortality_rate;
    }
  } else {
    max_density = 0.0;
    equilibrium_rate = base_mortality_rate;
  }
}

/** Safe copy constructor with proper error handling and CONSISTENT INITIALIZATION ORDER
 All the characteristics of the parent individual will be copied, except:
 - id (will be set to new value)
 - list of Neighbours (will be updated)
 - time until next event (will be drafted)
 @param rhs Parent individual */
Individual::Individual(const Individual& rhs)
  :
  id(++MAX_ID),
  x_coordinate(rhs.x_coordinate),
  y_coordinate(rhs.y_coordinate),
  species_id(rhs.species_id),
  base_mortality_rate(rhs.base_mortality_rate),
  vision_angle(rhs.vision_angle),
  step_length(rhs.step_length),
  base_dispersal_rate(rhs.base_dispersal_rate),
  current_dispersal_rate(rhs.base_dispersal_rate),
  neighbor_radius(rhs.neighbor_radius),
  base_birth_rate(rhs.base_birth_rate),
  birth_density_slope(rhs.birth_density_slope),
  mortality_density_slope(rhs.mortality_density_slope),
  matrix_mortality_multiplier(rhs.matrix_mortality_multiplier),
  matrix_dispersal_multiplier(rhs.matrix_dispersal_multiplier),
  density_type(rhs.density_type),
  mutation_rate(rhs.mutation_rate),
  plasticity(rhs.plasticity),
  habitat_selection_temperature(rhs.habitat_selection_temperature),
  sampling_points(rhs.sampling_points),
  // Initialize in SAME ORDER as main constructor
  current_birth_rate(0.0),
  current_death_rate(0.0),
  time_to_next_event(0.0),
  max_density(rhs.max_density),
  equilibrium_rate(rhs.equilibrium_rate),
  habitat_value(rhs.habitat_value),
  patch_id(0)
{
  // Check for ID overflow
  if (id == 0) {
    throw runtime_error("Individual ID counter overflow occurred during reproduction");
  }
  
  // Sets a random orientation in radians
  current_bearing = runif(0, 2.0 * M_PI);
  
  try {
    if (rhs.genotype_means.size() == 1) { // Checks if this is not a null simulation
      // Normal case: single genotype with mutation
      double mutated_genotype = rnorm(rhs.genotype_means[0], mutation_rate); // Generates the genotype of the offspring by drafting a random value from a normal distribution with mean= father genotype and sd= mutation rate
      
      if (!isfinite(mutated_genotype)) {
        throw runtime_error("Mutated genotype is not finite");
      }
      
      genotype_means.push_back(mutated_genotype);
      genotype_sds = rhs.genotype_sds; // Sets the individual's width value (or values)
      
      // Generate phenotype
      double phenotype = rnorm(mutated_genotype, plasticity); // Generates the phenotype of the offspring by adding a random value to its genotype value
      if (!isfinite(phenotype)) {
        throw runtime_error("Generated offspring phenotype is not finite");
      }
      environmental_optimum.push_back(phenotype);
      
    } else { // Checks if this is a null simulation
      // Neutral mode: multiple genotypes
      genotype_means.reserve(rhs.genotype_means.size());
      environmental_optimum.reserve(rhs.genotype_means.size());
      
      for (size_t i = 0; i < rhs.genotype_means.size(); i++) { // passes by the genotype means of the individual
        double mutated_genotype = rnorm(rhs.genotype_means[i], mutation_rate); // Generates the genotype of the offspring by drafting a random value from a normal distribution with mean= father genotype and sd= mutation rate
        
        if (!isfinite(mutated_genotype)) {
          throw runtime_error("Mutated genotype is not finite at position " + to_string(i));
        }
        
        genotype_means.push_back(mutated_genotype);
        
        double phenotype = rnorm(mutated_genotype, plasticity); // Generates phenotype by adding a random value to the genotype values
        if (!isfinite(phenotype)) {
          throw runtime_error("Generated phenotype is not finite at position " + to_string(i));
        }
        environmental_optimum.push_back(phenotype);
      }
      
      genotype_sds = rhs.genotype_sds; // Sets the individual's width values
    }
    
  } catch (const std::exception& e) {
    throw runtime_error(string("Failed to create offspring: ") + e.what());
  }
}

/** Function that updates both death and birth rates of individuals based on the number of individuals within the area of their density dependence radius
 It then calls a function to draft the time needed for execution of that individual based on its birth, death and dispersal rates*/
void Individual::update_rates(double local_density)
{
  double density = local_density; // Density of neighbour individuals (includes focal one)
  
  if (genotype_sds[0] == 0) { // Checks if individuals are completely specialists
    double optimal_habitat = environmental_optimum[0]; // Individual's optimal habitat value
    
    if (habitat_value == optimal_habitat) {
      // Individual is in exactly optimal habitat
      current_birth_rate = base_birth_rate - birth_density_slope * density; // Computes the actual birth rate on habitat patch (that is influenced by the density of neighbours)
      current_death_rate = base_mortality_rate + mortality_density_slope * density; // Computes the actual death rate on habitat patch (that is influenced by the density of neighbours)
      current_dispersal_rate = base_dispersal_rate; // Sets the standard movement rates
    } else {
      // Individual is in suboptimal habitat - treat as matrix
      current_birth_rate = 0; // Sets birth rate to 0
      current_death_rate = matrix_mortality_multiplier * base_mortality_rate + mortality_density_slope * density; // Sets the higher mortality rate attributed to non-optimal habitat
      current_dispersal_rate = matrix_dispersal_multiplier * base_dispersal_rate; // Sets the higher/lower movement rates attributed to non-optimal habitat
    }
  }
  else{
    
    current_birth_rate = base_birth_rate - birth_density_slope * density; // Computes the actual birth rate on habitat patch (that is influenced by the density of neighbours)
    current_death_rate = ((matrix_mortality_multiplier * base_mortality_rate) - ((sum_normal_densities(habitat_value, environmental_optimum, genotype_sds) / normal_density(environmental_optimum[0], environmental_optimum[0], genotype_sds[0])) * ((matrix_mortality_multiplier * base_mortality_rate) - base_mortality_rate))); // Computes the actual death rate on habitat patch (that is influenced by the suitability of its current habitat)
    
  }
  
  if(current_birth_rate < 0) // Checks if the birth rate is lower than possible
  {
    current_birth_rate = 0; // Sets to the lowest possible value for Birth
  }
  
  draw_next_event_time(); // Calls the function to draft the time needed to execute the next action
}

// Function that drafts the time for the next event based on all rates
void Individual::draw_next_event_time()
{
  time_to_next_event = rexp(1.0 / (current_dispersal_rate + current_birth_rate + current_death_rate)); // Sets the time to the individual
}

// Function that drafts one of the three possible actions (accounting for their respective rates) and returns the drafted action to the landscape
int Individual::select_action()
{
  vector<double> cumulative_probabilities; // creates a vector for storing the probabilities
  double total_rate = current_death_rate + current_birth_rate + current_dispersal_rate; // Sums all rates
  cumulative_probabilities.push_back(current_death_rate / total_rate); // stores the death event probability at the first position
  cumulative_probabilities.push_back(cumulative_probabilities[0] + (current_birth_rate / total_rate)); // stores the birth event probability at the second position
  cumulative_probabilities.push_back(cumulative_probabilities[1] + (current_dispersal_rate / total_rate)); // stores the dispersal event probability at the third position
  double random_value; // temporary variable for storing a random value
  random_value = runif(0.0, 1.0); // Samples between 0 and 1
  
  int selected_action; // temporary variable for storing drafted event
  for(unsigned int i = 0; i < cumulative_probabilities.size() - 1; i++) // Goes through all events
  {
    if(cumulative_probabilities[i] > random_value) // Finds the first action probability higher than the drafted number
    {
      selected_action = i; // Stores the drafted action
      return selected_action; // Returns the drafted action to the landscape (0 = death, 1 = birth, 2 = dispersal)
      break;
    }
  }
  return cumulative_probabilities.size() - 1; // returns sorted action
}

// Function that dislocates the XY coordinates of an individual by a "step_length" distance.
void Individual::perform_random_walk()
{
  // Drafts a random turning angle within the vision range (centered on 0) and updates the bearing
  current_bearing += runif(-vision_angle / 2.0, vision_angle / 2.0);
  
  // Wraps negative/overflow values to keep the bearing within [0, 2π)
  if (current_bearing < 0)
    current_bearing += 2.0 * M_PI;
  else if (current_bearing >= 2.0 * M_PI)
    current_bearing -= 2.0 * M_PI;
  
  // Computes displacement using the updated bearing
  double dx = cos(current_bearing) * step_length;
  double dy = sin(current_bearing) * step_length;
  
  // Updates the position
  x_coordinate += dx;
  y_coordinate += dy;
}

// Function that returns the value of the probability density function for the normal distribution
double Individual::normal_density(double x, double mean, double sd) const {
  
  return (1 / (sd * sqrt(2 * M_PI)) * (exp(-1 * (pow(x - mean, 2) / (2 * pow(sd, 2))))));
  
}

// Function that returns the value of the sum of several probability density functions for the normal distribution
double Individual::sum_normal_densities(double x, const vector<double>& mean_values, const vector<double>& sd_values) const {
  
  double probability_cumsum = 0; // Temporary variable for storing the sum
  
  for (size_t i = 0; i < mean_values.size(); i++) { // Passes by each normal distribution
    probability_cumsum = probability_cumsum + normal_density(x, mean_values[i], sd_values[i]); // Adds it to the counter
  }
  return probability_cumsum / mean_values.size(); // Returns the mean probability density
  
}

// Safe habitat selection with proper bounds checking - Function that selects from within a given set of coordinates based on the dispersing individual's preference ranks via a softmax function
void Individual::select_habitat(const vector<vector<double>>& candidate_locations) {
  if (candidate_locations.empty()) {
    throw invalid_argument("No candidate locations provided");
  }
  
  const size_t num_candidates = candidate_locations.size();
  
  // Validate that all candidate locations have correct format [x, y, environmental_value]
  for (size_t i = 0; i < num_candidates; ++i) {
    if (candidate_locations[i].size() != 3) {
      throw invalid_argument("Each candidate location must have exactly 3 values [x, y, env_value]");
    }
  }
  
  // Use vector instead of C-style array
  vector<double> fitness_scores(num_candidates); // Temporary vector for storing the scores
  double cumsum = 0.0; // Temporary variables for storing the drafting components
  
  // Calculate fitness scores for each candidate location
  for (size_t i = 0; i < num_candidates; ++i) { // Passes by the points
    double env_value = candidate_locations[i][2];
    
    try {
      // Calculate raw fitness based on environmental suitability
      double raw_fitness = sum_normal_densities(env_value, environmental_optimum, genotype_sds);
      
      // Apply temperature to the fitness calculation
      double fitness = exp(raw_fitness / habitat_selection_temperature);
      
      // Handle potential numerical issues
      if (!isfinite(fitness)) {
        fitness = 0.0;  // Set to minimum fitness if not finite
      }
      
      fitness_scores[i] = fitness;
      cumsum += fitness; // Sums that score to the total
      
    } catch (const std::exception& e) {
      // If fitness calculation fails, set to minimum
      fitness_scores[i] = 0.0;
    }
  }
  
  // Handle edge case where all fitness scores are zero
  if (cumsum <= 0.0) {
    // Fallback to uniform selection
    for (size_t i = 0; i < num_candidates; ++i) {
      fitness_scores[i] = 1.0;
    }
    cumsum = static_cast<double>(num_candidates);
  }
  
  // Select location using weighted random selection
  double choice = runif(0.0, 1.0); // Drafts between 0 and 1
  double accumulated_probability = 0.0; // Temporary variables for storing the drafting components
  size_t selected_index = 0; // Temporary variable for storing the chosen index
  
  for (size_t i = 0; i < num_candidates; ++i) { // Passes by the points
    accumulated_probability += fitness_scores[i] / cumsum; // Computes the comparative value of the ranks
    
    if (accumulated_probability > choice) { // Checks the comparative value of the ranks has exceeded the drafted value
      selected_index = i; // Chooses the index of the selected coordinate
      break;
    }
  }
  
  // Ensure we have a valid selection (fallback to last candidate)
  if (selected_index >= num_candidates) {
    selected_index = num_candidates - 1;
  }
  
  // Safe coordinate assignment with validation
  try {
    double new_x = candidate_locations[selected_index][0]; // Sets new x coordinate
    double new_y = candidate_locations[selected_index][1]; // Sets new y coordinate
    
    // Validate coordinates are finite
    if (!isfinite(new_x) || !isfinite(new_y)) {
      throw runtime_error("Selected coordinates are not finite");
    }
    
    x_coordinate = new_x;
    y_coordinate = new_y;
    
  } catch (const std::exception& e) {
    throw runtime_error(string("Failed to set new coordinates: ") + e.what());
  }
}