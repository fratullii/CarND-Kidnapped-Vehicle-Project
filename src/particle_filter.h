/**
 * particle_filter.h
 * 2D particle filter class.
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <string>
#include <vector>
#include "helper_functions.h"

struct Particle {
  int id;
  double x;
  double y;
  double theta;
  double weight;
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;
};


class ParticleFilter {  
 public:
  // Constructor
  // @param num_particles Number of particles
  ParticleFilter() : num_particles(0), is_initialized(false) {}

  // Destructor
  ~ParticleFilter() {}

  /**
   * init Initializes particle filter by initializing particles to Gaussian
   *   distribution around first position and all the weights to 1.
   * @param x Initial x position [m] (simulated estimate from GPS)
   * @param y Initial y position [m]
   * @param theta Initial orientation [rad]
   * @param std[] Array of dimension 3 [standard deviation of x [m], 
   *   standard deviation of y [m], standard deviation of yaw [rad]]
   */
  void init(double x, double y, double theta, double std[]);

  /**
   * prediction Predicts the state for the next time step
   *   using the process model.
   * @param delta_t Time between time step t and t+1 in measurements [s]
   * @param std_pos[] Array of dimension 3 [standard deviation of x [m], 
   *   standard deviation of y [m], standard deviation of yaw [rad]]
   * @param velocity Velocity of car from t to t+1 [m/s]
   * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
   */
  void prediction(double delta_t, double std_pos[], double velocity, 
                  double yaw_rate);
  
  /**
   * dataAssociation Finds which observations correspond to which landmarks 
   *   (likely by using a nearest-neighbors data association).
   * @param particle_obs LandmarkObst to be associated to a landmark
   * @param map_landmarks Landmark in a given maps
   */
  void dataAssociation(LandmarkObs observation, const vector<Map::single_landmark_s> &landmarks);
  
  /**
   * updateWeights Updates the weights for each particle based on the likelihood
   *   of the observed measurements. 
   * @param sensor_range Range [m] of sensor
   * @param std_landmark[] Array of dimension 2
   *   [Landmark measurement uncertainty [x [m], y [m]]]
   * @param observations Vector of landmark observations
   * @param map Map class containing map landmarks
   */
  void updateWeights(double sensor_range, double std_landmark[], 
                     const std::vector<LandmarkObs> &observations,
                     const Map &map_landmarks);
  
  /**
   * resample Resamples from the updated set of particles to form
   *   the new set of particles.
   */
  void resample();

  /**
   * Set a particles list of associations, along with the associations'
   *   calculated world x,y coordinates
   * This can be a very useful debugging tool to make sure transformations 
   *   are correct and assocations correctly connected
   */
  void SetAssociations(Particle& particle, const std::vector<int>& associations,
                       const std::vector<double>& sense_x, 
                       const std::vector<double>& sense_y);
  /**
   * initialized Returns whether particle filter is initialized yet or not.
   */
  const bool initialized() const {
    return is_initialized;
  }

  /**
   * Switch from the vehicle's coordinate system to the map's cooordinate system
   * 
   * 
   * @param obs_map Reference to the new LandmarkObs object
   * @param obs_veh Reference to the LandmarkObs object with coordinates in the vehicle's system
   * @param particle Reference to the particle objects for which the map coordinates must be computed
   */
  void changeCoordinates(LandmarkObs &obs_map, const LandmarkObs &obs_veh, const Particle &particle);

  /**
   * Assign to a vector only the landmarks within the sensor range of the vehicle, assuming it is
   * in the particle position
   * 
   * @param range_landmarks Vector of LandmarksObs objects to be initialize within this method
   * @param map_landmarks Map class containing map landmarks
   * @param range Sensor range [m]
   */
  void selectLandmarks(vector<Map::single_landmark_s> &range_landmarks, const Particle &particle, 
                        const Map &map_landmarks, double range);

  /**
   * Used for obtaining debugging information related to particles.
   */
  std::string getAssociations(Particle best);
  std::string getSenseCoord(Particle best, std::string coord);

  // Set of current particles
  std::vector<Particle> particles;

 private:
  // Number of particles to draw
  int num_particles; 
  
  // Flag, if filter is initialized
  bool is_initialized;
  
  // Vector of weights of all particles
  std::vector<double> weights; 
};

#endif  // PARTICLE_FILTER_H_