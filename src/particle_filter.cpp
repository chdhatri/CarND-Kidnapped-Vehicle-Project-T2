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
    if(is_initialized) {
        return;
    }
    num_particles = 100;
    
    //Set standard deviations for x, y, and theta.
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];
    
    // creates a normal (Gaussian) distribution for x  y and theta.
    normal_distribution<double> dist_x(x, std_x);
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);
    
    for (int i = 0; i < num_particles; ++i) {
        Particle particle;
        
        particle.id = 1;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1;
        
        particles.push_back(particle);
    }
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
    
    // creates a normal (Gaussian) distribution for x  y and theta.
    normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_y);
    normal_distribution<double> dist_theta(0, std_theta);
    
    //Add measurements to each particle and add random Gaussian noise.
    for(int i=0; i<num_particles;i++) {
        
        double theta = particles[i].theta;
        
        if(fabs(yaw_rate) < EPS){
            particles[i].x += velocity * delta_t * cos( theta );
            particles[i].y += velocity * delta_t * sin( theta );
            
        } else {
            
            particles[i].x += velocity/yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
            particles[i].y += velocity/yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
            
        }
        
        //Adding Gaussian sensor noise.
        
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
    
    for(int i=0; i < observations.size(); i++) {
        
        // Initialize min distance
        double minDist = numeric_limits<double>::max();
        
        //initialize map id
        int mapId = -1;
        
        for(int j=0; j < predicted.size(); j++) {
            
            double xDiff = observations[i].x - predicted[j].x;
            double yDiff = observations[i].y - predicted[j].y;
            
            double dist = xDiff * xDiff + yDiff * yDiff;
            
            if ( dist < minDist) {
                minDist = dist;
                mapId = predicted[j].id;
            }

    }
        // Update the observation map id to the predicted id.
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
    
    /* find landmarks in each particle range
    ** Transform observation coordinates 
    */
    for(int i =0; i < num_particles;i++) {
        double x_part = particles[i].x;
        double y_part = particles[i].y;
        double theta_part = particles[i].theta;
        
        //find landmarks in each particle range
        vector<LandmarkObs> inRangeLandmarks;
        for(int j = 0 ; j < map_landmarks.landmark_list.size(); j++) {
             float x_landmark = map_landmarks.landmark_list[j].x_f;
             float y_landmark = map_landmarks.landmark_list[j].y_f;
             int id_landmark = map_landmarks.landmark_list[j].id_i;
             
             //distance between land mark and particle
             double dist_x = x_part - x_landmark;
             double dist_y = y_part - y_landmark;
             double dist = sqrt(dist_x * dist_x + dist_y * dist_y);
             
             if(dist <= sensor_range) {
                 inRangeLandmarks.push_back(LandmarkObs{ id_landmark, dist_x, dist_y });
             }
        }
        
        // Transform observation coordinates (TOB's).
        vector<LandmarkObs> transformedObservations;
        for(int k = 0 ; k < observations.size(); k++) {
            double map_x = x_part + cos(theta_part)*observations[k].x - sin(theta_part)*observations[k].y;
            double map_y = y_part + sin(theta_part)*observations[k].x + cos(theta_part)*observations[k].y;
            transformedObservations.push_back(LandmarkObs{ observations[k].id, map_x, map_y });
            
        }
        // Observation association to landmark.
        dataAssociation(inRangeLandmarks, transformedObservations);
        
        particles[i].weight = 1.0;
        
        for(int m=0; m < transformedObservations.size();m++) {
            double tobs_x = transformedObservations[m].x;
            double tobs_y = transformedObservations[m].y;
            
            //associate each transformed observation with a land mark identifier.
            int landmarkId = transformedObservations[i].id;
            double landmarkx, landmarky;
            int i = 0;
            int nLandmarks = inRangeLandmarks.size();
            bool found = false;
            while(!found && i < nLandmarks) {
                if(inRangeLandmarks[i].id == landmarkId) {
                    found = true;
                    landmarkx = inRangeLandmarks[i].x;
                    landmarky = inRangeLandmarks[i].y;
                }
                
                i++;
            }
            
            //calculating weight

            double dX = tobs_x - landmarkx;
            double dY = tobs_y - landmarky;
            
            double weight = ( 1/(2*M_PI*std_x*std_y)) * exp( -( dX*dX/(2*std_x*std_x) + (dY*dY/(2*std_y*std_y)) ) );
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
