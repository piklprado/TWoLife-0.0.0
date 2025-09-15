#include <Rcpp.h>
#include <fstream>
#include "landscape.hpp"
#include "individual.hpp"
#include "random_utils.hpp"
using namespace Rcpp;

//' Run TWoLife Simulation
 //' 
 //' @param neighbor_radius Radius of density dependence influence
 //' @param initial_population_size Number of individuals at start
 //' @param vision_angle Angle used for orientation when dispersing
 //' @param step_length Distance of dispersal event
 //' @param base_dispersal_rate Rate at which individuals disperse
 //' @param base_birth_rate Basal birth rate
 //' @param base_mortality_rate Basal death rate
 //' @param birth_density_slope Slope of birth density dependence
 //' @param mortality_density_slope Slope of death density dependence
 //' @param habitat Square matrix of environmental values (0=matrix, 1=habitat)
 //' @param cell_size Resolution length of pixel side
 //' @param density_type Density type (0=global, 1=local)
 //' @param matrix_mortality_multiplier Mortality multiplier for matrix
 //' @param matrix_dispersal_multiplier Dispersal multiplier for matrix
 //' @param initial_placement_mode Initial position mode
 //' @param boundary_condition Boundary condition type
 //' @param max_events Maximum number of events
 //' @param initial_x_coordinates X coordinates of initial individuals
 //' @param initial_y_coordinates Y coordinates of initial individuals
 //' @param genotype_means Genotypical trait means
 //' @param genotype_sds Standard deviations
 //' @param mutation_rates Mutation rates
 //' @param plasticities Plasticities
 //' @param sampling_points Number of sampling points
 //' @param neutral_mode Neutral simulation mode
 //' @param output_file Optional output file path
 //' @return List containing simulation results
 //' @export
 // [[Rcpp::export]]
 List run_twolife_simulation(
     double neighbor_radius,
     int initial_population_size,
     double vision_angle,
     double step_length,
     double base_dispersal_rate,
     double base_birth_rate,
     double base_mortality_rate,
     double birth_density_slope,
     double mortality_density_slope,
     NumericMatrix habitat,
     double cell_size,
     int density_type,
     double matrix_mortality_multiplier,
     double matrix_dispersal_multiplier,
     int initial_placement_mode,
     int boundary_condition,
     double max_events,
     NumericVector initial_x_coordinates,
     NumericVector initial_y_coordinates,
     NumericVector genotype_means,
     NumericVector genotype_sds,
     NumericVector mutation_rates,
     NumericVector plasticities,
     IntegerVector sampling_points,
     bool neutral_mode,
     Nullable<String> output_file = R_NilValue
 ) {
   
   try {
     // Convert Rcpp vectors to std::vectors
     std::vector<double> initial_x_vec = as<std::vector<double>>(initial_x_coordinates);
     std::vector<double> initial_y_vec = as<std::vector<double>>(initial_y_coordinates);
     std::vector<double> genotype_means_vec = as<std::vector<double>>(genotype_means);
     std::vector<double> genotype_sds_vec = as<std::vector<double>>(genotype_sds);
     std::vector<double> mutation_rates_vec = as<std::vector<double>>(mutation_rates);
     std::vector<double> plasticities_vec = as<std::vector<double>>(plasticities);
     std::vector<int> sampling_points_vec = as<std::vector<int>>(sampling_points);
     
     // Create landscape
     auto world = std::make_unique<Landscape>(
       neighbor_radius, initial_population_size, vision_angle, step_length,
       base_dispersal_rate, base_birth_rate, base_mortality_rate,
       birth_density_slope, mortality_density_slope, habitat,
       cell_size, density_type, matrix_mortality_multiplier,
       matrix_dispersal_multiplier, initial_placement_mode,
       boundary_condition, initial_x_vec, initial_y_vec,
       genotype_means_vec, genotype_sds_vec, mutation_rates_vec,
       plasticities_vec, sampling_points_vec, neutral_mode
     );
     
     // Storage for simulation events
     std::vector<double> event_times;
     std::vector<int> event_types;
     std::vector<uint64_t> individual_ids;
     std::vector<int> patch_ids;
     std::vector<double> x_coordinates;
     std::vector<double> y_coordinates;
     std::vector<double> genotypes;
     
     // Optional file output
     std::unique_ptr<std::ofstream> file_output;
     if (output_file.isNotNull()) {
       String filename = output_file.get();
       file_output = std::make_unique<std::ofstream>(filename);
       if (!file_output->is_open()) {
         stop("Could not open output file: " + std::string(filename));
       }
     }
     
     // Record initial state
     for (size_t i = 0; i < world->count_individuals(); i++) {
       Individual* ind = world->get_individual(static_cast<int>(i));
       event_times.push_back(world->world_time);
       event_types.push_back(-1); // Initial state marker
       individual_ids.push_back(ind->get_id());
       patch_ids.push_back(ind->get_patch_id());
       x_coordinates.push_back(ind->get_x_coordinate());
       y_coordinates.push_back(ind->get_y_coordinate());
       genotypes.push_back(ind->get_genotype_mean());
       
       if (file_output) {
         *file_output << world->world_time << " " << ind->get_id() << " " 
                      << ind->get_patch_id() << " " << ind->get_x_coordinate() 
                      << " " << ind->get_y_coordinate() << " " 
                      << ind->get_genotype_mean() << std::endl;
       }
     }
     
     // Run simulation
     while (world->world_time < max_events && 
            world->count_individuals() > 0 && 
            world->count_individuals() < 1000000) {
       
       // Check for user interrupt periodically
       if (event_times.size() % 1000 == 0) {
         Rcpp::checkUserInterrupt();
       }
       
       int selected_individual_index = world->select_next_individual();
       int action_code = world->select_action(selected_individual_index);
       world->advance_world_time(selected_individual_index);
       
       Individual* selected_ind = world->get_individual(selected_individual_index);
       uint64_t individual_id = selected_ind->get_id();
       double individual_x = selected_ind->get_x_coordinate();
       double individual_y = selected_ind->get_y_coordinate();
       int individual_patch = selected_ind->get_patch_id();
       double individual_genotype = selected_ind->get_genotype_mean();
       
       bool left_boundary = world->perform_action(action_code, selected_individual_index);
       
       if (action_code == 2 && left_boundary) {
         world->update(0, selected_individual_index);
         action_code = 3; // Emigration marker
       } else {
         world->update(action_code, selected_individual_index);
         
         // Update coordinates if movement occurred and individual stayed
         if (action_code == 2 && !left_boundary && 
             selected_individual_index < static_cast<int>(world->count_individuals())) {
           Individual* moved_ind = world->get_individual(selected_individual_index);
           individual_x = moved_ind->get_x_coordinate();
           individual_y = moved_ind->get_y_coordinate();
           individual_patch = moved_ind->get_patch_id();
         }
       }
       
       // Record event
       event_times.push_back(world->world_time);
       event_types.push_back(action_code);
       individual_ids.push_back(individual_id);
       patch_ids.push_back(individual_patch);
       x_coordinates.push_back(individual_x);
       y_coordinates.push_back(individual_y);
       genotypes.push_back(individual_genotype);
       
       if (file_output) {
         *file_output << world->world_time << " " << action_code << " " 
                      << individual_id << " " << individual_patch << " " 
                      << individual_x << " " << individual_y << " " 
                      << individual_genotype << std::endl;
       }
     }
     
     if (file_output) {
       *file_output << "EOF" << std::endl;
       file_output->close();
     }
     
     // Collect final survivors
     NumericVector survivor_x(world->count_individuals());
     NumericVector survivor_y(world->count_individuals());
     NumericVector survivor_genotypes(world->count_individuals());
     IntegerVector survivor_ids(world->count_individuals());
     
     for (size_t i = 0; i < world->count_individuals(); i++) {
       Individual* survivor = world->get_individual(static_cast<int>(i));
       survivor_x[i] = survivor->get_x_coordinate();
       survivor_y[i] = survivor->get_y_coordinate();
       survivor_genotypes[i] = survivor->get_genotype_mean();
       survivor_ids[i] = static_cast<int>(survivor->get_id());
     }
     
     return List::create(
       Named("final_population_size") = static_cast<int>(world->count_individuals()),
       Named("survivor_x") = survivor_x,
       Named("survivor_y") = survivor_y,
       Named("survivor_genotypes") = survivor_genotypes,
       Named("survivor_ids") = survivor_ids,
       Named("event_times") = NumericVector(event_times.begin(), event_times.end()),
       Named("event_types") = IntegerVector(event_types.begin(), event_types.end()),
       Named("individual_ids") = NumericVector(individual_ids.begin(), individual_ids.end()),
       Named("patch_ids") = IntegerVector(patch_ids.begin(), patch_ids.end()),
       Named("x_coordinates") = NumericVector(x_coordinates.begin(), x_coordinates.end()),
       Named("y_coordinates") = NumericVector(y_coordinates.begin(), y_coordinates.end()),
       Named("genotypes") = NumericVector(genotypes.begin(), genotypes.end()),
       Named("world_size") = world->get_world_size(),
       Named("num_patches") = world->get_num_patches()
     );
     
   } catch (const std::exception& e) {
     stop("Simulation failed: " + std::string(e.what()));
   }
 }