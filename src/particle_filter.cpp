/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <unordered_set>

#include "helper_functions.h"

//debug using
using std::cout;
using std::endl;

using std::string;
using std::vector;
using std::default_random_engine;
using std::normal_distribution;
using std::discrete_distribution;
using std::unordered_set;

void ParticleFilter::init(double x, double y, double theta, double std[]){
  /**
   * TODO: Set the number of particles. Initialize all particles to 
  *   first position (based on estimates of x, y, theta and their uncertainties./c
  *   from GPS) and all weights to 1. 
  * TODO: Add random Gaussian noise to each particle.
  * NOTE: Consult particle_filter.h for more information about this method 
  *   (and others in this file).
  */

  // Define normal distributions around gps initial measurements
  default_random_engine gen;
  normal_distribution<double> dist_x (x, std[0]);
  normal_distribution<double> dist_y (y, std[1]);
  normal_distribution<double> dist_theta (theta, std[2]);

  num_particles = 400;  // TODO: Set the number of particles
  for (int i = 0; i < num_particles; ++i) {
    
    // Sample from the normal distribution
    Particle temp_particle;
    temp_particle.x = dist_x(gen);
    temp_particle.y = dist_y(gen);
    temp_particle.theta = dist_theta(gen);
    
    particles.push_back(temp_particle);

  }
  
  is_initialized = true;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   **/

double vel_yaw_ratio;
yaw_rate != 0 ? vel_yaw_ratio = velocity / yaw_rate : vel_yaw_ratio = velocity;

default_random_engine gen;

  for (int i = 0; i < num_particles; ++i){

    // Apply bicycle motion model to predict particle state
    particles[i].x += vel_yaw_ratio * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
    particles[i].y += vel_yaw_ratio * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
    particles[i].theta += yaw_rate * delta_t;

    // Add random gaussian noise by sampling from a normal distribution 
    
    normal_distribution<double> dist_x (0., std_pos[0]);
    normal_distribution<double> dist_y (0., std_pos[1]);
    normal_distribution<double> dist_theta (0., std_pos[2]);

    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);

  }

}

void ParticleFilter::dataAssociation(LandmarkObs &observation, const Map &landmarks_map) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   **/

   // Implementing naive nearest neighbor search
   observation.id = 0;
   double min_dist = dist2(observation.x, observation.y, 
                          landmarks_map.landmark_list[0].x_f, landmarks_map.landmark_list[0].y_f);
   for(int i = 1; i < landmarks_map.landmark_list.size(); ++i){
     double land_obs_dist = dist2(observation.x, observation.y, 
                                landmarks_map.landmark_list[0].x_f, landmarks_map.landmark_list[0].y_f);
     if (land_obs_dist < min_dist){
       observation.id = i;
       min_dist = land_obs_dist;
     }
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  weights.clear();

  for(int i = 0; i < particles.size(); ++i){

    particles[i].weight = 1.;

    // Filter landmarks out of particle range
    Map particle_landmarks;
    selectLandmarks(particle_landmarks, particles[i], map_landmarks, sensor_range);

    for(int j=0; j < observations.size(); ++j){

      LandmarkObs particle_obs;
      changeCoordinates(particle_obs, observations[j], particles[i]);
      dataAssociation(particle_obs, particle_landmarks);
      particles[i].weight *= multiv_prob(particle_obs.x, particle_obs.y,
                                         particle_landmarks.landmark_list[particle_obs.id].x_f,
                                         particle_landmarks.landmark_list[particle_obs.id].y_f,
                                         std_landmark[0], std_landmark[1]);

    }

    weights.push_back(particles[i].weight);
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */


  default_random_engine gen;
  discrete_distribution<int> dist(weights.begin(), weights.end());
  vector<int> sampled_particles_index;

  for(int i = 0; i < num_particles; ++i){
    sampled_particles_index.push_back(dist(gen));
  }

  vector<int>::iterator it;
  it = std::unique(sampled_particles_index.begin(), sampled_particles_index.end());
  sampled_particles_index.resize(std::distance(sampled_particles_index.begin(), it));

  vector<Particle> new_particles;
  for(int ind : sampled_particles_index){
    new_particles.push_back(particles[ind]);
  }
  particles = new_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

void ParticleFilter::changeCoordinates(LandmarkObs &obs_map, const LandmarkObs &obs_veh, const Particle &particle){

  double cos_theta = cos(particle.theta);
  double sin_theta = sin(particle.theta);
  obs_map.x = cos_theta*obs_veh.x - sin_theta*obs_veh.y;
  obs_map.y = sin_theta*obs_veh.x + cos_theta*obs_veh.y;

}

void ParticleFilter::selectLandmarks(Map &range_landmarks, const Particle &particle,
                    const Map &map_landmarks, double range){

  for(Map::single_landmark_s landmark : map_landmarks.landmark_list){

    if(dist2(particle.x, particle.y, landmark.x_f, landmark.y_f) <= pow(range, 2)){
      range_landmarks.landmark_list.push_back(landmark);
    }

  }
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}