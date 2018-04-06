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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	default_random_engine gen;

	//Normal distribution for noise
	normal_distribution<double> x_dist(x, std[0]);
	normal_distribution<double> y_dist(y, std[1]);
	normal_distribution<double> theta_dist(theta, std[2]);

	//init particles
	for(int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = x_dist(gen);
		particle.y = y_dist(gen);
		particle.theta = theta_dist(gen);
		particle.weight = 1.0;

		particles.push_back(particle);
		weights.push_back(particle.weight);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	for (int i = 0; i < num_particles; i++) 
	{
		//calculate particles
		// when yaw_rate is zero
		if (fabs(yaw_rate) < 0.000001)
		{
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}
		else 
		{
			particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}

		//Normal distribution for noise
		normal_distribution<double> x_dist(particles[i].x, std_pos[0]);
		normal_distribution<double> y_dist(particles[i].y, std_pos[1]);
		normal_distribution<double> theta_dist(particles[i].theta, std_pos[2]);
		//add random noise
		particles[i].x = x_dist(gen);
		particles[i].y = y_dist(gen);
		particles[i].theta = theta_dist(gen);
	}	

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) 
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i < observations.size(); i++) 
	{
		LandmarkObs current_obervation = observations[i];
		int closest_id = -1;
		double min_distance = numeric_limits<double>::max();

		for (int j = 0; j < predicted.size(); j++)
		{
			LandmarkObs current_predicted = predicted[j];
			double current_distance = dist(current_obervation.x, current_obervation.y,
														current_predicted.x, current_predicted.y);

			if (current_distance < min_distance)
			{
				min_distance = current_distance;
				closest_id = current_predicted.id; 
			}
		}
		observations[i].id = closest_id;
	}

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

	double sigma_x = std_landmark[0];
	double sigma_y = std_landmark[1];
	for(int i = 0; i < num_particles; i++)
	{
		Particle particle = particles[i];
		
		// convert map_landmarks to predicted landmarks
		vector<LandmarkObs> predicted_landmarks;

		for(int k = 0; k < map_landmarks.landmark_list.size(); k++)
		{
			LandmarkObs current_predicted_landmark;
			current_predicted_landmark.id = map_landmarks.landmark_list[k].id_i;
			current_predicted_landmark.x = map_landmarks.landmark_list[k].x_f;
			current_predicted_landmark.y = map_landmarks.landmark_list[k].y_f;
			if (fabs(current_predicted_landmark.x - particle.x) <= sensor_range &&
				fabs(current_predicted_landmark.y - particle.y) <= sensor_range)
				predicted_landmarks.push_back(current_predicted_landmark);
		}

		// convert observed landmark into map coodinates
		vector<LandmarkObs> observations_in_map_coodinates;

		for (int j = 0 ; j < observations.size(); j++)
		{
			LandmarkObs current_observation = observations[j];
			LandmarkObs current_observation_in_map_coordinates;

			current_observation_in_map_coordinates.id = current_observation.id;
			current_observation_in_map_coordinates.x = particle.x + current_observation.x * cos(particle.theta)
																								- current_observation.y * sin(particle.theta);
			current_observation_in_map_coordinates.y = particle.y + current_observation.x * sin(particle.theta)
																								+ current_observation.y * cos(particle.theta);
			observations_in_map_coodinates.push_back(current_observation_in_map_coordinates);
		}


		dataAssociation(predicted_landmarks, observations_in_map_coodinates);
		particles[i].weight = 1.0;

    for (int j = 0; j < observations_in_map_coodinates.size(); j++) 
    {
      
      double predicted_x = 0.0;
      double predicted_y = 0.0;
      double observed_x = observations_in_map_coodinates[j].x;
      double observed_y = observations_in_map_coodinates[j].y;

      int associated_prediction = observations_in_map_coodinates[j].id;

      //get predited x,y of the prediction-associated with observation
      for (unsigned int k = 0; k < predicted_landmarks.size(); k++) 
      {
        if (predicted_landmarks[k].id == associated_prediction) 
        {
          predicted_x = predicted_landmarks[k].x;
          predicted_y = predicted_landmarks[k].y;
        }
      }

      //calculate weight by using multivariate Gaussian
      double sigma_x = std_landmark[0];
      double sigma_y = std_landmark[1];
      double gauss_norm = ( 1.0 / (2 * M_PI * sigma_x * sigma_y));
      double exponet = exp(-(pow(predicted_x - observed_x, 2) / (2.0 * pow(sigma_x, 2)) 
		 						+ (pow(predicted_y - observed_y, 2) / (2.0 * pow(sigma_y, 2)) )));
      double update_weight = gauss_norm * exponet;

      //update weight
      particles[i].weight *= update_weight;
    }
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> resample_particles;
	vector<double> weights;
	double max_weight = -1.0;
	for (int i = 0; i < num_particles; i++) 
	{
		weights.push_back(particles[i].weight);
		if (particles[i].weight > max_weight)
			max_weight = particles[i].weight;
	}

	//generate a random number based on weights
	default_random_engine gen;
	
	double beta = 0.0;
	uniform_real_distribution<double> uniform_real_distribution(0, max_weight);
	uniform_int_distribution<int> uniform_int_distribution(0, num_particles - 1);
  int index = uniform_int_distribution(gen);
	//resample wheel
	for (int i = 0; i < particles.size(); i++)
	{
		beta += uniform_real_distribution(gen) * 2.0;
		
		while (weights[index] < beta)
		{
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resample_particles.push_back(particles[index]);
	}
  particles = resample_particles;

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
