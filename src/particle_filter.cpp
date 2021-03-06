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

#define EPS 0.00001

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
    
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    if (is_initialized) {
        return;
    }
    
    // Initializing the number of particles
    num_particles = 100;
    
    // Set the standard deviations for x, y, and theta
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];
    
    // Creating normal (Gaussian) distributions for x, y, and theta.
    normal_distribution<double> dist_x(x, std_x);
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);
    
    // Generate particles with normal distribution with mean on GPS values.
    for (int i = 0; i < num_particles; i++) {
        
        Particle particle;
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0;
        
        particles.push_back(particle);
    }
    
    // The filter is now initialized.
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    //
    double std_x = std_pos[0];
    double std_y = std_pos[1];
    double std_theta = std_pos[2];
    
    // Creating normal distributions
    normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_y);
    normal_distribution<double> dist_theta(0, std_theta);
    
   //Add measurements to each particle and add random Gaussian noise.
    for (int i = 0; i < num_particles; i++) {
        
        double theta = particles[i].theta;
        
        if ( fabs(yaw_rate) < EPS ) {
            particles[i].x += velocity * delta_t * cos( theta );
            particles[i].y += velocity * delta_t * sin( theta );
            
        } else {
            particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
            particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
            particles[i].theta += yaw_rate * delta_t;
        }
        
        // Adding Gaussian sensor noise.
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    int nObservations = observations.size();
    int nPredictions = predicted.size();
    
    for (int i = 0; i < nObservations; i++) {
        
        double minDist = numeric_limits<double>::max();
        
        // Initialize the mapId.
        int mapId = -1;
        
        for (int j = 0; j < nPredictions; j++ ) {
            
            double xDiff = observations[i].x - predicted[j].x;
            double yDiff = observations[i].y - predicted[j].y;
            
            double distance = sqrt(xDiff * xDiff + yDiff * yDiff);
            
           
            if ( distance < minDist ) {
                minDist = distance;
                mapId = predicted[j].id;
            }
        }
        
        // Update the observation identifier.
        observations[i].id = mapId;
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
    
    double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    
    for (int i = 0; i < num_particles; i++) {
        
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;
        
        vector<LandmarkObs> inRangeLandmarks;
        for(int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            float landmarkX = map_landmarks.landmark_list[j].x_f;
            float landmarkY = map_landmarks.landmark_list[j].y_f;
            int id = map_landmarks.landmark_list[j].id_i;
            double dX = x - landmarkX;
            double dY = y - landmarkY;
            double distance = sqrt(dX*dX + dY*dY);
            if ( distance <= sensor_range ) {
                inRangeLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY });
            }
        }
        
        // Transform observation coordinates.
        vector<LandmarkObs> transformedObservations;
        for( int k = 0;  k< observations.size(); k++) {
            double xx = x + cos(theta)*observations[k].x - sin(theta)*observations[k].y ;
            double yy = y + sin(theta)*observations[k].x + cos(theta)*observations[k].y ;
            transformedObservations.push_back(LandmarkObs{ observations[k].id, xx, yy });
        }
        
        //  associate observations to landmark.
        dataAssociation(inRangeLandmarks, transformedObservations);
        
        
        particles[i].weight = 1.0;
        // Calculate weights.
        for(int j = 0; j < transformedObservations.size(); j++) {
            double observationX = transformedObservations[j].x;
            double observationY = transformedObservations[j].y;
            
            int landmarkId = transformedObservations[j].id;
            
            double landmarkX, landmarkY;
            int k = 0;
             int nLandmarks = inRangeLandmarks.size();
            bool found = false;
            while( !found && k < nLandmarks ) {
                if ( inRangeLandmarks[k].id == landmarkId) {
                    found = true;
                    landmarkX = inRangeLandmarks[k].x;
                    landmarkY = inRangeLandmarks[k].y;
                }
                k++;
            }
            
            // Calculating weight.
            double dX = observationX - landmarkX;
            double dY = observationY - landmarkY;
            
            double weight = ( 1/(2*M_PI*std_x*std_y)) * exp( -( dX*dX/(2*std_x*std_y) + (dY*dY/(2*std_x*std_y)) ) );
            if (weight == 0) {
                particles[i].weight *= EPS;
            } else {
                particles[i].weight *= weight;
            }
        }
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    
    double beta = 0.0;

    
    vector<double> weights;
    double maxWeight = numeric_limits<double>::min();
    for(int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
        if ( particles[i].weight > maxWeight ) {
            maxWeight = particles[i].weight;
        }
    }
    
    // Creating uniform distributions.
    uniform_real_distribution<double> distDouble(0.0, maxWeight);
    uniform_int_distribution<int> distInt(0, num_particles - 1);
    
    // Generating index.
    int index = distInt(gen);

    
        // Resampling wheel
    vector<Particle> resampledParticles;
    for(int i = 0; i < num_particles; i++) {
        beta += distDouble(gen) * 2.0;
        while( beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        resampledParticles.push_back(particles[index]);
    }
    
    particles = resampledParticles;
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
