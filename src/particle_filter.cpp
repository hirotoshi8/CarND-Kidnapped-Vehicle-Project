/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

//#include <Eigen>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double init_x, double init_y, double init_theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// debug
	cout << "****************Initialized in particle filter !********************" << endl;

	// Creates a normal (Gaussian) distribution for x, y, theta.
	default_random_engine gen;
	normal_distribution<double> dist_x(init_x, std[0]);
	normal_distribution<double> dist_y(init_y, std[1]);
	normal_distribution<double> dist_theta(init_theta, std[2]);

	// Particles for initialization
	Particle initial_particle;

	for (int i = 0; i < num_particles; ++i) {

		initial_particle.x = dist_x(gen);
		initial_particle.y = dist_y(gen);
		initial_particle.theta = dist_theta(gen);
		initial_particle.weight = 1;
		// Store the initial particles
		this->particles.push_back(initial_particle);

		//debug
		//cout << "debug initial position" << endl;
		//cout << this->particles[i].x <<" "<< this->particles[i].y <<" "<< this->particles[i].theta << endl;
	}
	this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//debug
	cout << "**************** Predict in particle filter !********************" << endl;

	// Creates a normal (Gaussian) distribution for x, y, theta.
	default_random_engine gen;
	normal_distribution<double> noise_x(0, std_pos[0]);
	normal_distribution<double> noise_y(0, std_pos[1]);
	normal_distribution<double> noise_theta(0, std_pos[2]);

	// When yaw rate is zero
	if (yaw_rate != 0) {
		// Predict the vehicle position and heading
		double ratio = velocity / yaw_rate;

		for (int i = 0; i < this->num_particles; i++) {
			double pre_theta = this->particles[i].theta;
			this->particles[i].x += ratio*(sin(pre_theta + yaw_rate*delta_t) - sin(pre_theta)) + noise_x(gen);
			this->particles[i].y += ratio*(-cos(pre_theta + yaw_rate*delta_t) + cos(pre_theta)) + noise_y(gen);
			this->particles[i].theta += yaw_rate*delta_t + noise_theta(gen);

			//debug
			cout << "debug prediction" << endl;
			cout << this->particles[i].x << this->particles[i].y << this->particles[i].theta << endl;
		}
	}
	else
	{
		// Predict the vehicle position and heading
		double ratio = velocity*delta_t;

		for (int i = 0; i < this->num_particles; i++) {
			double pre_theta = this->particles[i].theta;
			this->particles[i].x += ratio*(sin(pre_theta)) + noise_x(gen);
			this->particles[i].y += ratio*(-cos(pre_theta)) + noise_y(gen);
			this->particles[i].theta += noise_theta(gen);

			//debug
			cout << "debug prediction in zero divide" << endl;
			cout << this->particles[i].x << this->particles[i].y << this->particles[i].theta << endl;
		}

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void transform_local_to_map(const Particle p, const LandmarkObs &obs, LandmarkObs *out) {
	// Transform local map to global map
	double x_particle = p.x;
	double y_particle = p.y;
	double heading_particle = p.theta;

	double obs_x = obs.x;
	double obs_y = obs.y;

	// Transform
	out->x = x_particle + cos(heading_particle)*obs_x - sin(heading_particle)*obs_y;
	out->y = y_particle + sin(heading_particle)*obs_x + cos(heading_particle)*obs_y;
}

double calculate_prob(double k, double std_x2, double std_y2, double map_x, double map_y, double x_f, double y_f) {
	// Calcilate Guassian Probability
	double x_term = (x_f - map_x)*(x_f - map_x) / (2 * std_x2);
	double y_term = (y_f - map_y)*(y_f - map_y) / (2 * std_y2);
	double prob = k * exp( -( (x_term + y_term) ));

	return prob;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	//debug
	cout << "**************** Update in particle filter !********************" << endl;

	// Update the particles weight for each particles 
		// Transform the particles as vehicle candidates from local_map to global_map
		// Data Association with Nearest Neighbors
		// Update the weight
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];
	double std_x2 = std_x*std_x;
	double std_y2 = std_y*std_y;
	double k = 1 / (2 * M_PI* std_x* std_y);

	std::vector<double> landmark_matching_probability;


	// Loop each particles
	for (int i = 0; i < this->num_particles; i++) {
		
		// Reset
		this->particles[i].sense_x.clear();
		this->particles[i].sense_y.clear();
		landmark_matching_probability.clear();

		// Transform the measuremsnts(obserbations) from local map to global map
		for (auto iter_obs = observations.begin(); iter_obs != observations.end(); ++iter_obs) {
			// Return value
			LandmarkObs landmark_found_by_particle_in_map;
			// transfere the particles position from local-map(particle-centered-map) to global-map
			transform_local_to_map(this->particles[i], *iter_obs, &landmark_found_by_particle_in_map);
			this->particles[i].sense_x.push_back(landmark_found_by_particle_in_map.x);
			this->particles[i].sense_y.push_back(landmark_found_by_particle_in_map.y);
		}

		// Data association with Nearest Neighbors
		for (int index_sensed = 0; index_sensed < particles[i].sense_x.size(); index_sensed++) {
			// Reset
			double min_dist = sensor_range;
			double prob = 0;

			// Search neighbors
			for (int index_map = 0; index_map < map_landmarks.landmark_list.size(); index_map++) {
				// Obsersation(Measurements)
				double x_map = this->particles[i].sense_x.at(index_sensed);
				double y_map = this->particles[i].sense_y.at(index_sensed);
				// Prediction
				double x_f = map_landmarks.landmark_list.at(index_map).x_f;
				double y_f = map_landmarks.landmark_list.at(index_map).y_f;
				// Calculate the distance
				double distance = dist(x_map, y_map, x_f, y_f);

				// Compare the distance
				if (distance < min_dist) {
					min_dist = distance;
					//update probability
					prob = calculate_prob(k, std_x2, std_y2, x_map, y_map, x_f, y_f);
					//Assign the "landmark ID" to measurements
					//this->particles[i].associations.push_back(map_landmarks.landmark_list.at(index_map).id_i);
				}
				else { /* Nothing to do */ }
			}
			// the probability is determined, owned by this particle
			landmark_matching_probability.push_back(prob);
			//debug
			//cout << "debug probability against the NN" << endl;
			//cout << prob << endl;
		}

		// Update the weight of this particles
		double final_prob = 1;

		for (int i = 0; i <  landmark_matching_probability.size(); i++) {
			// Caluculate final prob of the particle
			final_prob  *= landmark_matching_probability.at(i);
		}
		this->particles[i].weight = final_prob;
		//debug
		cout << "debug final prob of the particle" << endl;
		cout << this->particles[i].weight << endl;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	//debug
	cout << "**************** Resample in particle filter !********************" << endl;

	double weight_sum = 0;
	double max_weight = 0;

	// calculate the sum
	for (int i = 0; i < this->num_particles; i++) {
		weight_sum += this->particles[i].weight;
	}
	//debug
	cout << "debug weight sum" << endl;
	cout << weight_sum << endl;

	// Normalize the weights and Search the max weight
	for (int i = 0; i < this->num_particles; i++) {

		this->particles[i].weight /= weight_sum;

		//debug
		cout << "debug weight normalize" << endl;
		cout << this->particles[i].weight << endl;

		if (this->particles[i].weight > max_weight) {
			max_weight = this->particles[i].weight;
		}
	}
	//debug
	cout << "debug max weight" << endl;
	cout << max_weight << endl;

	// Wheel resampling
	// Creates a normal (Gaussian) distribution for x, y, theta.
	default_random_engine gen;
	std::uniform_int_distribution<int> U_dist_index(0, this->num_particles - 1);
	std::uniform_real_distribution<double> U_dist_weight(0, 2 * max_weight);
	//discrete_distribution<double> dist_x(init_x, std[0]);
	
	int index = U_dist_index(gen);
	//debug
	cout << "debug sampling index based on ""uniform_int_distribution""" << endl;
	cout << index << endl;

	double beta = 0;
	std::vector<Particle> re_particles;

	for (int count = 0; count < this->num_particles; count++){
		double random = U_dist_weight(gen);
		beta += random;
		//debug
		//cout << "debug sampling beta" << endl;
		//cout << random << endl;

		while (beta > this->particles[index].weight) {
			beta -= this->particles[index].weight;
			index = (index + 1) % this->num_particles;
			//debug
			//cout << "debug wheel resampling" << endl;
			//cout << index << endl;
		}
		// Resample the particle
		re_particles.push_back(this->particles[index]);
		//debug
		cout << "debug resampling index" << endl;
		cout << index << endl;
	}
	// Update the particles
	this->particles = re_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
