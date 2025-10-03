#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <vector>
#include <atomic>
#include <stdexcept>
#include <cmath>

using namespace std;

/** \brief The Individual class represents an agent of the simulation
 *
 *  This class contains information pertinent to an individual, including its location, state, 
 *  a vector containing pointers to close neighbours, length distance of a dispersal event, 
 *  rate of dispersal, etc. This class DOES NOT contain methods of landscape responsibility, 
 *  such as the update of the neighbour vector of each individual. */
class Individual
{
private:
  
  // Private properties - REORDERED to match constructor initialization order
  
  /** FIXED: Thread-safe unique individual identifier */
  static atomic<uint64_t> MAX_ID;
  const uint64_t id;
  /** The X coordinate of the individual */
  double x_coordinate;
  /** The Y coordinate of the individual */
  double y_coordinate;
  /** The species identifier number of the individual (TBI)*/
  const int species_id;
  /** The basal death rate (The rate at which individuals die on a habitat patch without neighbours) */
  const double base_mortality_rate;
  /** Angle that an individual is currently facing */
  double current_bearing;
  /** Angle used for orientation when dispersing */
  const double vision_angle;
  /** The Length distance of a dispersal event (random walk) or maximum dispersal radius (Selection step) */
  const double step_length;
  /** The basal Dispersal rate of the individual  */
  const double base_dispersal_rate;
  /** Final/current Dispersal rate of the individual (all effects considered) */
  double current_dispersal_rate;
  /** Radius of density dependent influence */
  double neighbor_radius;
  /** The basal birth rate (The rate at which individuals give birth on a habitat patch without neighbours)*/
  const double base_birth_rate;
  /** The slope of the birth density dependence function*/
  const double birth_density_slope;
  /** The slope of the death density dependence function */
  const double mortality_density_slope;
  /** Constant that indicates how many times higher the death rate should be on non-habitat pixels */
  const double matrix_mortality_multiplier;
  /** Constant that indicates how many times lower the movement rate should be on non-habitat pixels*/
  const double matrix_dispersal_multiplier;
  /** Density type (0 = global, 1 = local/within an individual radius)*/
  const int density_type;
  /**The mutation rate for genetic variation between generations */
  const double mutation_rate;
  /** The plasticity for phenotypic variation from genotype */
  const double plasticity;
  /** Temperature parameter for softmax habitat selection (controls selectivity) */
  const double habitat_selection_temperature;
  /** The number of sampling coordinates drafted at each dispersal event */
  const int sampling_points;
  /** Final/current birth rate of the individual (all effects considered) */
  double current_birth_rate;
  /** Final/current death rate of the individual (all effects considered) */
  double current_death_rate;
  /** Drafted time required for the individual to execute an action */
  double time_to_next_event;
  /** Maximum density the individual is capable of bearing*/
  double max_density;
  /** Birth and death rates when they are at equilibrium (populational K) */
  double equilibrium_rate;
  /** Identifier of the Type of habitat the individual is currently on (0 = matrix; 1 = habitat) */
  double habitat_value;
  /** Identifier of the patch of habitat the individual is currently on  */
  int patch_id;
  /** FIXED: Safe genetic properties using vectors - The expressed phenotypic optimum environmental value for an individual (where its mortality rate is the lowest)*/
  vector<double> environmental_optimum;
  /** FIXED: Safe genetic properties using vectors - The genetic optimum environmental value for an individual (underlying genetic potential) */
  vector<double> genotype_means;
  /** FIXED: Safe genetic properties using vectors - The standard deviation of environmental usage by an individual, how generalist it is*/
  vector<double> genotype_sds;
  
  // Private Methods
  /** Function that generates random numbers, following an exponential distribution, corresponding to the time needed to execute the next action, taking into account the birth, death and dispersal rates. */
  void draw_next_event_time();
  
public:
  /** Vector for storing the neighbours of the individual */
  vector<Individual*> neighbors;
  
  /** Constructor of the Individual class
   Must be called by the Landscape class to position individuals at the start of the simulation. */
  Individual(
    /** The X coordinate of the individual */
    double x_coordinate,
    /** The Y coordinate of the individual */
    double y_coordinate,
    /** The species identifier number of the individual  \ref TBI */
    const int species_id,
    /** The basal death rate (The rate at which individuals die on a habitat patch without neighbours) */
    const double base_mortality_rate,
    /** Angle that an individual is currently facing */
    double current_bearing,
    /** Angle used for orientation when dispersing */
    const double vision_angle,
    /** The Length distance of a dispersal event (random walk) or maximum dispersal radius (Selection step) */
    const double step_length,
    /** The basal Dispersal rate of the individual */
    const double base_dispersal_rate,
    /** Radius of density dependent influence */
    const double neighbor_radius,
    /** The basal birth rate (The rate at which individuals give birth on a habitat patch without neighbours) */
    const double base_birth_rate,
    /** The slope of the birth density dependence function */
    const double birth_density_slope,
    /** The slope of the death density dependence function */
    const double mortality_density_slope,
    /** Constant that indicates how many times higher the death rate should be on non-habitat pixels */
    const double matrix_mortality_multiplier,
    /** Constant that indicates how many times lower the movement rate should be on non-habitat pixels */
    const double matrix_dispersal_multiplier,
    /** Density type (0 = global, 1 = local/within an individual radius)*/
    const int density_type,
    /** The mutation rate for genetic variation between generations */
    double mutation_rate,
    /** The plasticity for phenotypic variation from genotype */
    double plasticity,
    /** FIXED: Safe genetic parameters using const references - The genetic optimum environmental value for an individual (underlying genetic potential) */
    const vector<double>& genotype_means,
    /** FIXED: Safe genetic parameters using const references - The standard deviation of environmental usage by an individual, how generalist it is*/
    const vector<double>& genotype_sds,
    /** The number of sampling coordinates drafted at each dispersal event */
    int sampling_points,
    /** Temperature parameter for softmax habitat selection; lower = more selective, higher = more random */
    double habitat_selection_temperature = 1.0
  );
  
  /** Copy constructor, used for generating new individuals by asexual reproduction
   All the characteristics of the parent individual will be copied, except:
   - id (will be set to new value)
   - list of Neighbours (will be updated)
   - time until next event (will be drafted)
   @param rhs Parent individual */
  Individual(const Individual& rhs);
  
  /** FIXED: Proper destructor (automatic with vectors) */
  ~Individual() = default;
  
  /** FIXED: Delete problematic operations for safety */
  Individual& operator=(const Individual&) = delete;
  Individual(Individual&&) = delete;
  Individual& operator=(Individual&&) = delete;
  
  /** Function that sets the maximum id to 0 */
  static void reset_id() {MAX_ID = 0;}
  
  /** FIXED: Thread-safe ID getter - Function that returns the ID of an individual */
  uint64_t get_id() const {return id;}
  
  /** Function that passes the input vector of neighbours of an individual (individuals within a radius distance) to the individual.
   Must be called at each time step by the landscape
   @param neighbors_list The list of individuals within a radius distance of the focal individual*/
  void set_neighbors(const vector<Individual*>& neighbors_list){neighbors = neighbors_list;}
  
  /** Function that updates the habitat type of the pixel the individual is currently on
   Must be called at each time step by the landscape (0= matrix 1= habitat)
   @param habitat_type Pixel address on the landscape matrix */
  void set_habitat_value(const double habitat_type){habitat_value = habitat_type;}
  
  /** Function that updates the fragment identifier of the pixel the individual is currently on
   @param patch_label Pixel address on the patch id matrix */
  void set_patch_id(const int patch_label){patch_id = patch_label;}
  
  /** Function that updates the X coordinate of the individual
   @param new_x The new x Coordinate*/
  void set_x_coordinate(double new_x){x_coordinate = new_x;}
  
  /** Function that updates the Y coordinate of the individual
   @param new_y The new y Coordinate */
  void set_y_coordinate(double new_y){y_coordinate = new_y;}
  
  /** Function that returns the x coordinate of the individual */
  inline const double get_x_coordinate() const {return x_coordinate;}
  
  /** Function that returns the y coordinate of the individual */
  inline const double get_y_coordinate() const {return y_coordinate;}
  
  /** Function that returns the density dependence radius of the individual */
  const double get_neighbor_radius() const {return neighbor_radius;}
  
  /** Function that returns the density type affecting the individual (0= global, 1=local) */
  const int get_density_type() const {return density_type;}
  
  /** FIXED: Safe neighborhood size calculation - Function that returns the number of individuals within the radius of the focal individual for density calculations (it also includes the focal individual) */
  size_t neighborhood_size() const {return neighbors.size() + 1;}
  
  /** Function that returns the environmental value of the pixel the individual is currently on */
  const int get_patch_id() const {return patch_id;}
  
  /** Function that returns the step length of the individual*/
  const double get_step_length() const {return step_length;}
  
  /** Function that returns the number of sampling coordinates an individual has */
  const int get_sampling_points() const {return sampling_points;}
  
  /** FIXED: Safe genotype access with validation - Function that returns the genotype mean of the individual */
  const double get_genotype_mean() const {
    if (genotype_means.empty()) {
      throw runtime_error("No genotype means available");
    }
    return genotype_means[0];
  }
  
  /** FIXED: Safe phenotype access with validation - Function that returns the environmental optimum (phenotype) of the individual */
  const double get_environmental_optimum() const {
    if (environmental_optimum.empty()) {
      throw runtime_error("No environmental optimum available");
    }
    return environmental_optimum[0];
  }
  
  /** FIXED: Safe genotype SD access with validation - Function that returns the genotype standard deviation (niche width) of the individual */
  const double get_genotype_sd() const {
    if (genotype_sds.empty()) {
      throw runtime_error("No genotype sds available");
    }
    return genotype_sds[0];
  }
  
  // Other Public Methods
  /** Function that Returns the drafted time for executing an action by the individual based on its birth, death and dispersal rates. */
  const double get_time_to_next_event() const {return time_to_next_event;}
  
  /** Function that drafts one of the three possible actions (0 = death, 1 = birth, 2 = dispersal) (accounting for their respective rates) and returns the drafted action to the landscape */
  int select_action();
  
  /** Function that updates both death and birth rates of individuals based on the number of individuals within the area of their density dependence radius
   It then calls a function to draft the time needed for execution of that individual based on its birth, death and dispersal rates
   @param local_density The number of individuals pondered by the area of the density dependence radius */
  void update_rates(double local_density);
  
  /** Function that dislocates the XY coordinates of an individual by a "step_length" distance and oriented by the angle the individual is currently facing added of a random component (vision_angle) (if vision_angle == 360 -> Random Walk) */
  void perform_random_walk();
  
  /** Function that returns the value of the probability density function for the normal distribution
   @param x Quantile of interest
   @param mean Mean value of the distribution
   @param sd Standard deviation of the distribution*/
  double normal_density(double x, double mean, double sd) const;
  
  /** Function that returns the value of the sum of several probability density functions for the normal distribution
   @param x Quantile of interest
   @param mean_values vector containing the Mean values of the distribution
   @param sd_values vector containing the Standard deviations of the distribution*/
  double sum_normal_densities(double x, const vector<double>& mean_values, const vector<double>& sd_values) const;
  
  /** FIXED: Safe habitat selection with proper container - Function that selects from within a given set of coordinates based on the dispersing individual's preference ranks via a softmax function with configurable temperature
   @param candidate_locations Vector containing X and Y coordinates (first two cols), and the environmental values of that coordinate's pixel (third col) */
  void select_habitat(const vector<vector<double>>& candidate_locations);
  
};

#endif // INDIVIDUAL_H