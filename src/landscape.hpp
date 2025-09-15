#ifndef LANDSCAPE_H
#define LANDSCAPE_H
#include "individual.hpp"
#include <vector>
#include <memory>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <Rcpp.h>

using namespace std;

/** \brief The Landscape class implements the world where the agents exist.
 *
 * The Landscape class is responsible for generating the domain where the agents (individuals) reside.
 * This class is responsible for creating individuals, updating the world clock, and establishing communication between individuals
 * \sa \ref Landscape::update */

class Landscape
{
private:
  
  // Private properties - REORDERED to match constructor initialization exactly
  /** Length of a square landscape side */
  const double world_side_length;
  /** Number of individuals at the start of the simulation */
  const unsigned long initial_population_size;
  /** Number of pixels in a square landscape side */
  const int cells_per_side;
  /** Resolution Length of a square pixel side */
  const double cell_size;
  /** The boundary condition type affects how individuals interact with the edges of the landscape (0= absorbing, 1= periodic (pacman), 2= reflective) */
  const int boundary_condition;
  /** The initial position of individuals (0 = origin, 1 = random, 2 = normally distributed with mean on origin, 3 = input initial coordinates)*/
  const int initial_placement_mode;
  /** number of non-matrix patches */
  int num_habitat_patches;
  /** FIXED: Dynamic 2D vector for habitat grid (safe memory management) - Matrix containing the environmental values of the landscape pixels (0= matrix, 1= habitat)*/
  vector<vector<double>> habitat_grid;
  /** FIXED: Dynamic 2D vector for patch grid (safe memory management) - Matrix containing the fragment ID that each pixel resides (0 for matrix; 1, 2, 3, ... for fragments) */
  vector<vector<int>> patch_grid;
  /** FIXED: Vector instead of raw pointer for patch areas - Pointer for a vector storing the area of each patch */
  vector<double> patch_areas;
  /** FIXED: Smart pointer vector for safe population management - Vector for storing the individuals of the population */
  vector<unique_ptr<Individual>> population;
  
  //Private methods
  
  /** Function responsible for calling the constructor of the Individual class to create N individuals at the start of the simulation */
  void initialize_population(
      /** Radius of density dependent influence */
      const double neighbor_radius,
      /** Number of individuals at the start of the simulation */
      const int initial_population_size,
      /** Angle used for orientation when dispersing */
      const double vision_angle,
      /** The Length distance of a dispersal event (random walk) or maximum dispersal radius (Selection step) */
      const double step_length,
      /** The rate at which individuals disperse */
      const double base_dispersal_rate,
      /** The basal birth rate (The rate at which individuals give birth on a habitat patch without neighbours) */
      const double base_birth_rate,
      /** The basal death rate (The rate at which individuals die on a habitat patch without neighbours) */
      const double base_mortality_rate,
      /** The slope of the birth density dependence function */
      const double birth_density_slope,
      /** The slope of the death density dependence function*/
      const double mortality_density_slope,
      /** Constant that indicates how many times higher the death rate should be on non-habitat pixels */
      const double matrix_mortality_multiplier,
      /** Constant that indicates how many times lower the movement rate should be on non-habitat pixels */
      const double matrix_dispersal_multiplier,
      /** Density type (0 = global, 1 = local/within an individual radius) */
      const int density_type,
      /** FIXED: Safe vector parameters - Vector containing the x coordinates of initial individuals */
      const vector<double>& initial_x_coordinates,
      /** FIXED: Safe vector parameters - Vector containing the y coordinates of initial individuals */
      const vector<double>& initial_y_coordinates,
      /** FIXED: Safe vector parameters - Vector containing the genotypical trait means of the initial individuals */
      const vector<double>& genotype_means,
      /** FIXED: Safe vector parameters - Vector containing the standard deviations of the initial individuals */
      const vector<double>& genotype_sds,
      /** FIXED: Safe vector parameters - Vector containing the mutation rates of the initial individuals */
      const vector<double>& mutation_rates,
      /** FIXED: Safe vector parameters - Vector containing the plasticities of the initial individuals */
      const vector<double>& plasticities,
      /** FIXED: Safe vector parameters - Vector containing the number of points initial individuals sample when selecting habitats */
      const vector<int>& sampling_points,
      /** Boolean value switching simulations to neutral state (all individuals acting as an average individual) */
      bool neutral_mode);
  
  /** Creates an individual's neighbourhood list (each of the individuals within a radius distance of the focal individual), or all of them (global)
   @param individual an object of the Individual class */
  void update_neighbors(Individual* individual) const;
  
  /** Updates an individual's neighbourhood list (each of the individuals within a radius distance of the focal individual), or all of them (global)
   @param action_code a variable containing the action performed by the individual
   @param individual_index index of the individual that performed the action */
  void update_neighbor_lists(int action_code, int individual_index) const;
  
  /** Updates in the individual the environmental value of the pixel corresponding to its current coordinate (if binary: 0= matrix, 1= habitat)
   @param individual an object of the Individual class */
  void update_habitat_value(Individual* individual) const;
  
  /** Updates in the individual the patch identification value of the pixel corresponding to its current coordinate
   @param individual an object of the Individual class */
  void update_patch_id(Individual* individual) const;
  
  /** Applies the boundary condition after a dispersal event
   @param individual an object of the Individual class */
  bool apply_boundary_condition(Individual* individual);
  
  /** FIXED: Helper function for safe grid access */
  bool is_valid_grid_position(int x, int y) const {
    return x >= 0 && x < cells_per_side && y >= 0 && y < cells_per_side;
  }
  
public:
  
  /** The world counter used for storing how much time has already passed */
  double world_time;
  
  //Public methods
  /** FIXED: Constructor with safe parameter types - Constructor of the Landscape class */
  Landscape(
    /** Radius of density dependent influence */
    const double neighbor_radius,
    /** Number of individuals at the start of the simulation */
    const int initial_population_size,
    /** Angle used for orientation when dispersing */
    const double vision_angle,
    /** The Length distance of a dispersal event (random walk) or maximum dispersal radius (Selection step) */
    const double step_length,
    /** The rate at which individuals disperse */
    const double base_dispersal_rate,
    /** The basal birth rate (The rate at which individuals give birth on a habitat patch without neighbours) */
    const double base_birth_rate,
    /** The basal death rate (The rate at which individuals die on a habitat patch without neighbours) */
    const double base_mortality_rate,
    /** The slope of the birth density dependence function */
    const double birth_density_slope,
    /** The slope of the death density dependence function */
    const double mortality_density_slope,
    /** Square habitat matrix containing environmental values (0= matrix, 1= habitat) */
    const Rcpp::NumericMatrix& habitat,
    /** Resolution Length of a square pixel side*/
    const double cell_size,
    /** Density type (0 = global, 1 = local/within an individual radius) */
    const int density_type,
    /** Constant that indicates how many times higher the death rate should be on non-habitat pixels*/
    const double matrix_mortality_multiplier,
    /** Constant that indicates how many times lower the movement rate should be on non-habitat pixels */
    const double matrix_dispersal_multiplier,
    /** The initial position of individuals (0 = origin, 1 = random, 2 = normally distributed with mean on origin)*/
    const int initial_placement_mode,
    /** The boundary condition type affects how individuals interact with the edges of the landscape (0= absorbing, 1= periodic (pacman), 2= reflective)*/
    const int boundary_condition,
    /** FIXED: Safe vector parameters - Vector containing the x coordinates of initial individuals */
    const vector<double>& initial_x_coordinates,
    /** FIXED: Safe vector parameters - Vector containing the y coordinates of initial individuals */
    const vector<double>& initial_y_coordinates,
    /** FIXED: Safe vector parameters - Vector containing the genotypical trait means of the initial individuals */
    const vector<double>& genotype_means,
    /** FIXED: Safe vector parameters - Vector containing the standard deviations of the initial individuals */
    const vector<double>& genotype_sds,
    /** FIXED: Safe vector parameters - Vector containing the mutation rates of the initial individuals */
    const vector<double>& mutation_rates,
    /** FIXED: Safe vector parameters - Vector containing the plasticities of the initial individuals */
    const vector<double>& plasticities,
    /** FIXED: Safe vector parameters - Vector containing the number of points initial individuals sample when selecting habitats */
    const vector<int>& sampling_points,
    /** Boolean value switching simulations to neutral state */
    bool neutral_mode
  );
  
  /** FIXED: Automatic destructor with smart pointers */
  ~Landscape() = default;
  
  /** Function that calls other functions of the Individual class to update the vector of individuals of the Landscape object (update_neighbor_lists, update_habitat_value, update_patch_id, update_rates()) 
   @param action_code A value representing the selected action (0= death, 1= birth, 2= dispersal)
   @param individual_index Index of the Individual performing the action
   */
  void update(int action_code, int individual_index);
  
  /** Function that uses the get_time_to_next_event function of the Individual class to draft times for each individual within the landscape and returns the individual with the least amount of time required for the execution of the next action */
  int select_next_individual();
  
  /** Function that receives the selected individual and calls a function member of the Individual class to randomly select one of the three possible actions for that individual to perform (0= death, 1= birth, 2= dispersal)
   @param individual_index The position of the individual with the lowest drafted time */
  int select_action(const int individual_index){return population[individual_index]->select_action();}
  
  /** Function that executes the selected action. It also returns a positive boolean value if the individual disperses out of the landscape boundary
   @param action_code A value representing the selected action (0= death, 1= birth, 2= dispersal)
   @param individual_index The position of the individual with the lowest drafted time */
  bool perform_action(int action_code, int individual_index);
  
  /** Function that updates the world clock, adding time required for the action execution by the selected individual
   @param individual_index The position of the individual with the lowest drafted time */
  void advance_world_time(const int individual_index){world_time += population[individual_index]->get_time_to_next_event();}
  
  /** FIXED: Safe population count - Function that calls a function of the Individual class to return the total number of individuals currently in the landscape */
  size_t count_individuals() const{return population.size();}
  
  /** FIXED: Safe individual access - Function that returns an individual from within the landscape vector of individuals
   @param index The position (in vector) of the individual to be returned */
  Individual* get_individual(int index) const {
    if (index < 0 || index >= static_cast<int>(population.size())) {
      throw out_of_range("Individual index out of bounds");
    }
    return population[index].get();
  }
  
  /** Function that returns the current number of species in the landscape
   \ref TBI */
  const int count_species() const;
  
  /** Function that returns the length of square landscape side */
  const double get_world_size() const {return world_side_length;}
  
  /** Function that returns the current world time */
  const double get_world_time() const {return world_time;}
  
  /** Function that computes and returns the distance between a pair of individuals
   @param individual_a First individual
   @param individual_b Second individual */
  double compute_distance(const Individual* individual_a, const Individual* individual_b) const;
  
  /** Function that computes and returns the density of individuals within a radius distance area from a focal individual
   @param focal_individual Focal individual */
  double compute_local_density(const Individual* focal_individual) const;
  
  /** Function that returns a negative boolean value if the inserted individual was present at the start of the simulation and a positive value if the individual was born during the simulation run time (used for coloring original individuals differently)
   @param individual Individual to be tested */
  const bool is_born_during_simulation(const Individual* individual) const {
    return individual->get_id() > initial_population_size;
  }
  
  /** Recursive function used to find and identify which pixels belong to each of the landscape patches
   @param x The initial value of the x coordinate
   @param y The initial value of the y coordinate
   @param current_label The initial patch identification number */
  void identify_habitat_patches(int x, int y, int current_label);
  
  /** Function that returns the number of fragments on the landscape */
  int get_num_patches() const {return num_habitat_patches;}
  
  /** FIXED: Safe patch area access - Function that computes and returns the area of a fragment in the landscape
   @param patch_index The patch number */
  double get_patch_area(int patch_index) const {
    if (patch_index < 0 || patch_index >= static_cast<int>(patch_areas.size())) {
      throw out_of_range("Invalid patch index");
    }
    return patch_areas[patch_index];
  }
  
  /** Function that decides upon, and calls other functions to perform, the desired method of dispersion (Random walk or habitat selection).
   In the habitat selection case this function also samples x Points within the individual's step_length radius and the individual's relative fitness on that location.
   @param individual_index The index of the individual with the lowest drafted time */
  bool perform_movement(int individual_index);
  
};

#endif // LANDSCAPE_H