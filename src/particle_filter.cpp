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
#include <limits>
#include <cassert>

#include "constants.h"
#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    /**
    std = std_x; std_y, std_theta
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */

    num_particles = my_constants::num_particles;
    std::default_random_engine generator;
    std::normal_distribution<float> dist_x(x, std[0]);
    std::normal_distribution<float> dist_y(y, std[1]);
    std::normal_distribution<float> dist_theta(theta, std[2]);

    for (int n=0; n<num_particles; ++n) {
        Particle particle;
        particle.id = n;
        particle.x = dist_x(generator);
        particle.y = dist_y(generator);
        particle.theta = dist_theta(generator);
        particle.weight = 1.0;
        particles.push_back(particle);
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

    // generate gaussian distributions
    std::default_random_engine gen;
    std::normal_distribution<float> dist_x(0, std_pos[0]);
    std::normal_distribution<float> dist_y(0, std_pos[1]);
    std::normal_distribution<float> dist_theta(0, std_pos[2]);

    // for each particle...
    for (auto& p : particles) {
        // move
        if (std::abs(yaw_rate) > my_constants::yaw_rate_min) {
            double theta = yaw_rate * delta_t;
            p.x += velocity / yaw_rate * (sin(p.theta + theta) - sin(p.theta));
            p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + theta));
            p.theta += theta;
        }
        else {
            p.x += velocity * delta_t * cos(p.theta);
            p.y += velocity * delta_t * sin(p.theta);
        }

        // add measurement noise
        p.x += dist_x(gen);
        p.y += dist_y(gen);
        p.theta += dist_theta(gen);
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
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations_in_VCS,
                                   const Map &map_landmarks) {
    /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   *
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

    // transform Map to vector of LandmarkObs
    vector<LandmarkObs> landmarks;
    for(auto map_landmark : map_landmarks.landmark_list) {
        LandmarkObs landmark;
        landmark.id = map_landmark.id_i;
        landmark.x = map_landmark.x_f;
        landmark.y = map_landmark.y_f;
        landmarks.push_back(landmark);
    }

    // for each particle
    for (auto& p : particles) {

        // transform observations from VCS into MCS
        vector<LandmarkObs> observations_in_MCS;
        for (auto& obs_in_VCS : observations_in_VCS) {
            LandmarkObs obs_in_MCS;
            obs_in_MCS.id = obs_in_VCS.id;
            obs_in_MCS.x = obs_in_VCS.x * cos(p.theta) - obs_in_VCS.y * sin(p.theta) + p.x;
            obs_in_MCS.y = obs_in_VCS.x * sin(p.theta) + obs_in_VCS.y * cos(p.theta) + p.y;
            observations_in_MCS.push_back(obs_in_MCS);
        }

        // update particle weight with likelihood of observations
        for (auto& obs_in_MCS : observations_in_MCS) {
            LandmarkObs landmark_nearest = get_nearest_landmark(obs_in_MCS, landmarks);
            double likelihood = calc_likelihood_of_observation(obs_in_MCS, landmark_nearest, std_landmark);
            p.weight *= likelihood;
        }
    }

    get_debug_info();
    normalize_weights();
}

void ParticleFilter::resample() {
    /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    /* particles_resampled[N];
        cumsum[num_particles] = ...
        sum = sum(cumsum)
        while(n < num_particles):
            index_random = random()
            index_particle = ...
            particles_resampled[n] = particles[index_particle]
        particles = particles_resampled;
        normalize_weights() # so that sum_weights = 0
        */

    // define weight distribution
    std::vector<float> weights;
    for (auto& p:particles) {
        weights.push_back(p.weight);
    }
    std::discrete_distribution<> distribution(weights.begin(), weights.end());
    std::default_random_engine gen;

    // create new particle list based on existing one
    std::vector<Particle> particles_new;
    while (particles_new.size() < num_particles) {
        int id_picked = distribution(gen);
        particles_new.push_back(particles[id_picked]);
    }
    particles = particles_new;
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



LandmarkObs ParticleFilter::get_nearest_landmark(LandmarkObs observation, std::vector<LandmarkObs> landmarks)
{
    double dist_final = std::numeric_limits<float>::max();
    LandmarkObs landmark_nearest;
    for (auto& lm:landmarks) {
        double dist = std::sqrt( std::pow(lm.x - observation.x,2) + std::pow(lm.y - observation.y,2));
        if (dist < dist_final) {
            dist_final = dist;
            landmark_nearest = lm;
        }
    }
    return landmark_nearest;
}

double ParticleFilter::calc_likelihood_of_observation(LandmarkObs observation,
                                                      LandmarkObs landmark,
                                                      double sigma_landmark[2])
{
    double likelihood = ( 1 / (2 * M_PI * sigma_landmark[0]* sigma_landmark[1])
            * std::exp( -1 *
                        ( 0.5 * std::pow(observation.x - landmark.x,2) / std::pow(sigma_landmark[0],2)
                        + 0.5 * std::pow(observation.y - landmark.y,2) / std::pow(sigma_landmark[1],2)
            )
            )
            );
    return likelihood;
}


void ParticleFilter::get_debug_info()
{
    // get highest weight and sum of weight
    double weights_sum = 0;
    double weight_max = 0;
    for (auto& p : particles) {
        weights_sum += p.weight;
        if (p.weight > weight_max) {
            weight_max = p.weight;
        }
    }

    // only for breakpoint setting
    int dummy_breakpoint = 0;
    if (weight_max < 0.0001) {
        dummy_breakpoint += 1;
    }
}

void ParticleFilter::normalize_weights()
{
    double sum_weights = 0;
    for (auto& particle : particles) {
        sum_weights += particle.weight;
    }
    assert(sum_weights>0);
    for (auto& particle : particles) {
        particle.weight /= sum_weights;
    }
}
