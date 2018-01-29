 /**
 * Copyright (C) 2018 Greenweaves Software Limited
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>
 */

#ifndef _GALAXY_H
#define _GALAXY_H

#include "particle.h"

 /**
  * Create all bodies needed at start of run
  */
std::vector<Particle*>  createParticles(int numbodies,double inivel,double ini_radius,double mass);
 
 /**
  * Execute simulation
  */
// void simulate(int start_iter,int max_iter,std::vector<Body*> bodies, double theta, double G, double dt, int img_iter,std::string path,std::string config_file_name,int check_energy);
 
 /**
  * Generate help text
  */
void help(int numbodies,double inivel,double ini_radius,double mass,int max_iter,double theta, double G, double dt, int img_iter,std::string path,std::string config_file_name);
 
 /**
  * Restore configuration from saved file
  */
// bool restore_config(std::string path,std::string name,std::vector<Body*>& bodies, int &iter, double &theta, double &G, double& dt);

 /**
  * Restore value stored by encode
  */
// void save_config( std::vector<Body*>& bodies, int iter, double theta, double G, double dt, std::string path,std::string name="config.txt");

/**
 * Retrieve position, mass, and velocities stored for one Body
 */
// Body * extract_body(std::string line);

/**
 * Calculate kinetic energy for bodies
 */
// double get_kinetic_energy(std::vector<Body*> bodies);

/**
 * Calculate gravitational potential energy for bodies
 */
// double get_potential_energy(std::vector<Body*> bodies,double G);

/**
 * Print total energy
 */
// void print_energy(double energy, double total_energy,double initial_energy);

#endif

 
 
 
 
 
 
 
