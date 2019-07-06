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

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  //std::default_random_engine gen;
  normal_distribution<double> noise_x(0.0, std[0]);
  normal_distribution<double> noise_y(0.0, std[1]);
  normal_distribution<double> noise_theta(0.0, std[2]);
  num_particles = 100;  // TODO: Set the number of particles
 
  for(int i=0;i<num_particles;++i){
    Particle p = Particle();
    p.id = i;
    p.x = x;
    p.y = y;
    p.theta = theta;
    p.weight = 1;
    
    p.x += noise_x(gen);
    p.y += noise_y(gen);
    p.theta += noise_theta(gen);
    
    particles.push_back(p);
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
   */
  //std::default_random_engine gen;
  normal_distribution<double> noise_x(0, std_pos[0]);
  normal_distribution<double> noise_y(0, std_pos[1]);
  normal_distribution<double> noise_theta(0, std_pos[2]);
  
  for (int i=0; i<num_particles; ++i) {
    if (fabs(yaw_rate) < 0.00001) {
      particles[i].x += velocity * cos(particles[i].theta) * delta_t;
      particles[i].y += velocity * sin(particles[i].theta) * delta_t;
    } else {
      particles[i].x += velocity / yaw_rate * ( sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta) );
      particles[i].y += velocity / yaw_rate * ( cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t) );
      particles[i].theta += yaw_rate * delta_t;
    }
    particles[i].x += noise_x(gen);
    particles[i].y += noise_y(gen);
    particles[i].theta += noise_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  /*std::cout << "predictions:" << std::endl;
  for(unsigned int i=0;i<predicted.size(); ++i){
    std::cout << predicted[i].id << " " <<  predicted[i].x << " "<<  predicted[i].y << std::endl;
  }*/
  
  for(unsigned int i=0; i<observations.size(); ++i) {
    double min_d = 1000000;
    int map_id = -1;
    for(unsigned int j=0; j<predicted.size(); ++j) {
      double d = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
      if(d < min_d) {
        min_d = d;
        map_id = predicted[j].id;
      }
    }
    observations[i].id = map_id;
  }
  
  /*std::cout << "observations:" << std::endl;
  for(unsigned int i=0;i<observations.size(); ++i){
    std::cout << observations[i].id << " " <<  observations[i].x << " "<<  observations[i].y << std::endl;
  }*/
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a multi-variate Gaussian 
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

  for (int i=0; i<num_particles; ++i) {
    // Create observation predictions based on particle position.
    // Use landmarks whose distance are within sensor range.
    std::vector<LandmarkObs> predictions;
    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); ++j) {
      float landmark_x = map_landmarks.landmark_list[j].x_f;
      float landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;
      if (dist(particles[i].x, particles[i].y, landmark_x, landmark_y) < sensor_range) {
        predictions.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
      }
    }
    
    // Transform from vehicle coordinates to map coordinates
    std::vector<LandmarkObs> map_observations;
    for (unsigned int j=0; j<observations.size(); j++) {
      double transformed_x = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta);
      double transformed_y = particles[i].y + observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta);

      map_observations.push_back(LandmarkObs{observations[j].id, transformed_x, transformed_y});
    }

    //Nearest-neighbors
    dataAssociation(predictions, map_observations);

    // Update weights
    particles[i].weight = 1.0;
    for (unsigned int j=0; j<map_observations.size(); j++) {
      LandmarkObs m = map_observations[j];
      LandmarkObs p;
      // Find the predicted landmark whose id is same as the observed landmark
      for (unsigned int k = 0; k < predictions.size(); k++) {
        if (predictions[k].id == m.id) {
          p = predictions[k];
        }
      }
      
      //std::cout <<"m:"<<m.id <<","<<m.x<<","<<m.y<< " p:"<<p.id<<","<<p.x<<","<<p.y<<std::endl;
      
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      
      double f = pow(m.x-p.x,2)/(2.0*pow(s_x,2))+pow(m.y-p.y,2)/(2.0*pow(s_y,2));
      double normalizer = 1.0/(2.0*M_PI*s_x*s_y);
      //std::cout <<   f << " " << std::endl;
      particles[i].weight *= normalizer*exp(-f);                                                          
    }
    //std::cout << std::endl;

  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  // Create discrete distribution based on particles weights
  std::vector<double> w;
  for (int i=0; i<num_particles; ++i) {
    w.push_back(particles[i].weight);
  }
  std::discrete_distribution<int> dist{w.begin(), w.end()};

  // Resample
  std::vector<Particle> new_particles;
  for (int i=0; i<num_particles; ++i) {
    new_particles.push_back(particles[dist(gen)]);
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