#include <Rcpp.h>
#include <fstream>
#include "landscape.hpp"
#include "individual.hpp"
#include "random_utils.hpp"
using namespace Rcpp;

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
    NumericVector habitat_selection_temperatures,
    bool neutral_mode,
    String history_detail = "standard",
    Nullable<int> master_seed = R_NilValue,
    Nullable<String> output_file = R_NilValue
) {
  
  try {
    // Validate history_detail parameter - convert Rcpp::String to std::string
    std::string detail_level = std::string(history_detail);
    if (detail_level != "minimal" && detail_level != "standard" && detail_level != "full") {
      stop("history_detail must be 'minimal', 'standard', or 'full'");
    }
    
    // Set the random seed if provided
    if (master_seed.isNotNull()) {
      int seed_value = as<int>(master_seed.get());
      Rcpp::Environment base_env("package:base");
      Rcpp::Function set_seed = base_env["set.seed"];
      set_seed(seed_value);
    }
    
    // Convert Rcpp vectors to std::vectors
    std::vector<double> initial_x_vec = as<std::vector<double>>(initial_x_coordinates);
    std::vector<double> initial_y_vec = as<std::vector<double>>(initial_y_coordinates);
    std::vector<double> genotype_means_vec = as<std::vector<double>>(genotype_means);
    std::vector<double> genotype_sds_vec = as<std::vector<double>>(genotype_sds);
    std::vector<double> mutation_rates_vec = as<std::vector<double>>(mutation_rates);
    std::vector<double> plasticities_vec = as<std::vector<double>>(plasticities);
    std::vector<int> sampling_points_vec = as<std::vector<int>>(sampling_points);
    std::vector<double> habitat_selection_temperatures_vec = as<std::vector<double>>(habitat_selection_temperatures);
    
    // Create landscape
    auto world = std::make_unique<Landscape>(
      neighbor_radius, initial_population_size, vision_angle, step_length,
      base_dispersal_rate, base_birth_rate, base_mortality_rate,
      birth_density_slope, mortality_density_slope, habitat,
      cell_size, density_type, matrix_mortality_multiplier,
      matrix_dispersal_multiplier, initial_placement_mode,
      boundary_condition, initial_x_vec, initial_y_vec,
      genotype_means_vec, genotype_sds_vec, mutation_rates_vec,
      plasticities_vec, sampling_points_vec, 
      habitat_selection_temperatures_vec, neutral_mode
    );
    
    // Storage for simulation events - ALWAYS tracked
    std::vector<double> event_times;
    std::vector<int> event_types;
    std::vector<uint64_t> individual_ids;
    
    // Storage for standard and full detail levels
    std::vector<int> patch_ids;
    std::vector<double> x_coordinates;
    std::vector<double> y_coordinates;
    std::vector<double> genotypes;
    
    // Storage for full detail level only
    std::vector<double> phenotypes;
    std::vector<double> widths;
    
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
      
      // Always record these
      event_times.push_back(world->world_time);
      event_types.push_back(-1); // Initial state marker
      individual_ids.push_back(ind->get_id());
      
      // Standard and full detail levels
      if (detail_level != "minimal") {
        patch_ids.push_back(ind->get_patch_id());
        x_coordinates.push_back(ind->get_x_coordinate());
        y_coordinates.push_back(ind->get_y_coordinate());
        genotypes.push_back(ind->get_genotype_mean());
      }
      
      // Full detail level only
      if (detail_level == "full") {
        phenotypes.push_back(ind->get_environmental_optimum());
        widths.push_back(ind->get_genotype_sd());
      }
      
      if (file_output) {
        *file_output << world->world_time << " -1 " << ind->get_id();
        if (detail_level != "minimal") {
          *file_output << " " << ind->get_patch_id() << " " 
                       << ind->get_x_coordinate() << " " 
                       << ind->get_y_coordinate() << " " 
                       << ind->get_genotype_mean();
        }
        if (detail_level == "full") {
          *file_output << " " << ind->get_environmental_optimum() 
                       << " " << ind->get_genotype_sd();
        }
        *file_output << std::endl;
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
      
      // Always capture these
      uint64_t individual_id = selected_ind->get_id();
      
      // Standard/full detail captures
      double individual_x = 0.0;
      double individual_y = 0.0;
      int individual_patch = 0;
      double individual_genotype = 0.0;
      double individual_phenotype = 0.0;
      double individual_width = 0.0;
      
      if (detail_level != "minimal") {
        individual_x = selected_ind->get_x_coordinate();
        individual_y = selected_ind->get_y_coordinate();
        individual_patch = selected_ind->get_patch_id();
        individual_genotype = selected_ind->get_genotype_mean();
      }
      
      if (detail_level == "full") {
        individual_phenotype = selected_ind->get_environmental_optimum();
        individual_width = selected_ind->get_genotype_sd();
      }
      
      bool left_boundary = world->perform_action(action_code, selected_individual_index);
      
      if (action_code == 2 && left_boundary) {
        world->update(0, selected_individual_index);
        action_code = 3; // Emigration marker
      } else {
        world->update(action_code, selected_individual_index);
        
        // Update coordinates if movement occurred and individual stayed
        if (action_code == 2 && !left_boundary && 
            selected_individual_index < static_cast<int>(world->count_individuals())) {
          if (detail_level != "minimal") {
            Individual* moved_ind = world->get_individual(selected_individual_index);
            individual_x = moved_ind->get_x_coordinate();
            individual_y = moved_ind->get_y_coordinate();
            individual_patch = moved_ind->get_patch_id();
          }
        }
      }
      
      // Record event - always record minimal data
      event_times.push_back(world->world_time);
      event_types.push_back(action_code);
      individual_ids.push_back(individual_id);
      
      // Standard and full detail levels
      if (detail_level != "minimal") {
        patch_ids.push_back(individual_patch);
        x_coordinates.push_back(individual_x);
        y_coordinates.push_back(individual_y);
        genotypes.push_back(individual_genotype);
      }
      
      // Full detail level only
      if (detail_level == "full") {
        phenotypes.push_back(individual_phenotype);
        widths.push_back(individual_width);
      }
      
      if (file_output) {
        *file_output << world->world_time << " " << action_code << " " << individual_id;
        if (detail_level != "minimal") {
          *file_output << " " << individual_patch << " " 
                       << individual_x << " " << individual_y << " " 
                       << individual_genotype;
        }
        if (detail_level == "full") {
          *file_output << " " << individual_phenotype << " " << individual_width;
        }
        *file_output << std::endl;
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
    NumericVector survivor_phenotypes(world->count_individuals());
    NumericVector survivor_widths(world->count_individuals());
    IntegerVector survivor_ids(world->count_individuals());
    
    for (size_t i = 0; i < world->count_individuals(); i++) {
      Individual* survivor = world->get_individual(static_cast<int>(i));
      survivor_x[i] = survivor->get_x_coordinate();
      survivor_y[i] = survivor->get_y_coordinate();
      survivor_genotypes[i] = survivor->get_genotype_mean();
      survivor_phenotypes[i] = survivor->get_environmental_optimum();
      survivor_widths[i] = survivor->get_genotype_sd();
      survivor_ids[i] = static_cast<int>(survivor->get_id());
    }
    
    // Build result list conditionally based on history_detail
    List result = List::create(
      Named("final_population_size") = static_cast<int>(world->count_individuals()),
      Named("survivor_x") = survivor_x,
      Named("survivor_y") = survivor_y,
      Named("survivor_genotypes") = survivor_genotypes,
      Named("survivor_phenotypes") = survivor_phenotypes,
      Named("survivor_widths") = survivor_widths,
      Named("survivor_ids") = survivor_ids,
      Named("event_times") = NumericVector(event_times.begin(), event_times.end()),
      Named("event_types") = IntegerVector(event_types.begin(), event_types.end()),
      Named("individual_ids") = NumericVector(individual_ids.begin(), individual_ids.end())
    );
    
    // Add standard detail fields
    if (detail_level != "minimal") {
      result["patch_ids"] = IntegerVector(patch_ids.begin(), patch_ids.end());
      result["x_coordinates"] = NumericVector(x_coordinates.begin(), x_coordinates.end());
      result["y_coordinates"] = NumericVector(y_coordinates.begin(), y_coordinates.end());
      result["genotypes"] = NumericVector(genotypes.begin(), genotypes.end());
    }
    
    // Add full detail fields
    if (detail_level == "full") {
      result["phenotypes"] = NumericVector(phenotypes.begin(), phenotypes.end());
      result["widths"] = NumericVector(widths.begin(), widths.end());
    }
    
    // Add spatial and metadata
    result["world_width"] = world->get_world_width();
    result["world_height"] = world->get_world_height();
    result["world_size"] = world->get_world_size();
    result["num_patches"] = world->get_num_patches();
    result["history_detail"] = history_detail;
    
    return result;
    
  } catch (const std::exception& e) {
    stop("Simulation failed: " + std::string(e.what()));
  }
}