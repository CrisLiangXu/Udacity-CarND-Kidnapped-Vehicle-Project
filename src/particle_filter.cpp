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
#include <sstream>

#include "helper_functions.h"

//using std::string;
//using std::vector;
using namespace std;

// declare a random engine
static default_random_engine engine;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * Set the number of particles. Initialize all particles to first
   *   position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */

  if(is_initialized){return;}
  // Initializing the number of particles
  num_particles = 100;

  //Extracting standard deviations
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  // normal distributions for given GPS
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);


  // Creating particals with normal distributions for given GPS
  for(int i = 0; i<num_particles; i++)
  {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(engine);
    particle.y = dist_y(engine);
    particle.theta = dist_theta(engine);
    particle.weight = 1.0;

    particles.push_back(particle);
  }

  /* After finished initialization, is_initialized wil return ture*/
  is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  //normal distribution for sensor noise
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	for(int i = 0; i < num_particles; i++) {
		if(fabs(yaw_rate) < 0.00001) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
      /*yaw rate = 0*/
		}
		else{
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate * delta_t;
      /*yaw rate !=0 */
		}

		//Adding Noise
		particles[i].x += dist_x(engine);
		particles[i].y += dist_y(engine);
		particles[i].theta += dist_theta(engine);
	}

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
  for(unsigned int i = 0; i < observations.size(); i++) {
    unsigned int nObs = observations.size();
    unsigned int nPred = predicted.size();
    for(unsigned int i = 0; i < nObs; i++)
    {
      double minDist = numeric_limits<double>::max();
      int mapId = -1;
      for(unsigned j = 0; j< nPred; j++)
      {
        double xDist = observations[i].x - predicted[j].x;
        double yDist = observations[i].y - predicted[j].y;
        double distance = xDist * xDist + yDist * yDist;
        if(distance < minDist)
        {
          minDist = distance;
          mapId = predicted[j].id;
        }
        observations[i].id = mapId;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * Update the weights of each particle using a mult-variate Gaussian
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

  for(int i = 0; i<num_particles; i++)
  {
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;

    //vector for landmark prediction
    vector<LandmarkObs> predictions;

    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      float landmark_x = map_landmarks.landmark_list[j].x_f;
      float landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;
      if(fabs(landmark_x - particle_x) <= sensor_range && fabs(landmark_y - particle_y) <= sensor_range)
      {
        predictions.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
      }
    }

    //vector for landmark observation
    vector<LandmarkObs> mapped_obs;
    for(unsigned int k = 0; k < observations.size(); k++)
    {
      double mapped_x = cos(particle_theta)*observations[k].x-sin(particle_theta)*observations[k].y + particle_x;
      double mapped_y = sin(particle_theta)*observations[k].x+cos(particle_theta)*observations[k].y + particle_y;
      mapped_obs.push_back(LandmarkObs{observations[k].id, mapped_x, mapped_y});
    }

    /*Data association for the predictions and transformed observations on current particle */
    dataAssociation(predictions, mapped_obs);
    particles[i].weight = 1.0;
    for(unsigned int m = 0; m < mapped_obs.size(); m++)
    {
      double obs_x, obs_y, predict_x, predict_y;
      obs_x = mapped_obs[m].x;
      obs_y = mapped_obs[m].y;

      //the prediction coordinate base on the current observation
      for(unsigned int k = 0; k < predictions.size(); k++)
      {
        if(predictions[k].id == mapped_obs[m].id)
        {
          predict_x = predictions[k].x;
          predict_y = predictions[k].y;
        }
     }

      //Weight for this observation with multivariate Gaussian
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(predict_x-obs_x,2)/(2*pow(s_x, 2)) + (pow(predict_y-obs_y,2)/(2*pow(s_y, 2))) ) );

      //Product of this obersvation weight with total observations weight
      particles[i].weight *= obs_w;
    }
  }
}

void ParticleFilter::resample()
{
  /**
   * Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<double> weights;
  double weight_max = numeric_limits<double>::min();   //max weight
  for(int i = 0; i < num_particles; i++)
  {
    weights.push_back(particles[i].weight);
    if(particles[i].weight  > weight_max)
    {
      weight_max = particles[i].weight;
    }
  }
  uniform_real_distribution<double> distDouble(0.0, weight_max);
  uniform_int_distribution<int> distInt(0, num_particles - 1);
    int index = distInt(engine);
    double beta = 0.0;
    vector<Particle> resampledParticles;
    for(int j = 0; j < num_particles; j++) {
        beta += distDouble(engine) * 2.0;
        while(beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        resampledParticles.push_back(particles[index]);
    }
  particles = resampledParticles;

}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y)
{
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
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
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord)
{
  vector<double> v;

  if (coord == "X")
  {
    v = best.sense_x;
  }
  else
  {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
