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
 
#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <unistd.h>
#include <getopt.h>
#include <random>
#include <algorithm>
#include "tree.h"
#include "barnes_hut.h"
#include "galaxy.h"
#include "utils.h"

/**
 *  Long version of command line options.
 */

struct option long_options[] = {
	{"dt",  		required_argument,	0, 'd'},
	{"config",  	required_argument,	0, 'c'},
	{"G",  			required_argument, 	0, 'G'},
    {"help",  		no_argument, 		0, 'h'},
	{"img_iter",	required_argument, 	0, 'i'},
	{"max_iter",  	required_argument, 	0, 'm'},
	{"numbodies",  	required_argument, 	0, 'n'},
	{"path",  		required_argument, 	0, 'p'},
	{"ini_radius",  required_argument, 	0, 'r'},
	{"mass",  		required_argument, 	0, 's'},
	{"theta",  		required_argument, 	0, 't'},
	{"inivel",  	required_argument, 	0, 'v'},
	{0, 			0, 					0, 0}
};	
	
/**
 * Main program. Parse command line options, create bodies, then run simulation.
 */
int main(int argc, char **argv) {

    // Theta-criterion of the Barnes-Hut algorithm.
    double theta = 0.5;
    // Mass of a body.
    double mass = 1.0;
    // Initially, the bodies are distributed inside a circle of radius ini_radius.
    double ini_radius = 0.1;
    // Initial maximum velocity of the bodies.
    double inivel = 0.1;
    // The "gravitational constant" is chosen so as to get a pleasant output.
    double G = 4.e-6;
    // Discrete time step.
    double dt = 1.e-3;
    // outside the initial radius are removed).
    int numbodies = 1000;
    // Number of time-iterations executed by the program.
    int max_iter = 10000;
    // Frequency at which PNG images are written.
    int img_iter = 20;
	
	std::string config_file_name="config.txt";
	
	std::string path = "./configs/";
	
	int option_index = 0;
	int c;
	
	while ((c = getopt_long (argc, argv, "c:d:G:hi:m:n:p:r:Ss:t:v:",long_options, &option_index)) != -1)
    switch (c){
		case 'c':{
			std::stringstream param(optarg);
			param>>config_file_name;
			std::cout<<"Configuration File:="<<config_file_name<<std::endl;
			break;
		}
		
		case 'd':{
			std::stringstream param(optarg);
			param>>dt;
			std::cout<<"dt="<<dt<<std::endl;
			break;
		}
		
		case 'G':{
			std::stringstream param(optarg);
			param>>G;
			std::cout<<"G="<<G<<std::endl;
			break;
		}
		
		case 'h':{
			help( numbodies, inivel, ini_radius, mass, max_iter, theta,  G,  dt,  img_iter, path,config_file_name);
			return 0;
		}
		
		case 'i':{
			std::stringstream param(optarg);
			param>>img_iter;
			std::cout<<"Frequency at which PNG images are written="<<img_iter<<std::endl;
			break;
		}
		
		case 'm':{
			std::stringstream param(optarg);
			param>>max_iter;
			std::cout<<"Number of iterations="<<max_iter<<std::endl;
			break;
		}
		
		case 'n':{
			std::stringstream param(optarg);
			param>>numbodies;
			std::cout<<"Number of bodies="<<numbodies<<std::endl;
			break;
		}
		
		case 'p':{
			std::stringstream param(optarg);
			param>>path;
			std::cout<<"Path="<<path<<std::endl;
			break;
		}
		
		case 'r':{
			std::stringstream param(optarg);
			param>>ini_radius;
			std::cout<<"Initial radius="<<ini_radius<<std::endl;
			break;
		}
		
		case 'S':{
			std::cout<<"Seed random number generator"<<std::endl;
			std::srand(1);
			break;
		}
		
		case 's':{
			std::stringstream param(optarg);
			param>>mass;
			std::cout<<"mass="<<mass<<std::endl;
			break;
		}
			
		case 't':{
			std::stringstream param(optarg);
			param>>theta;
			std::cout<<"Theta="<<theta<<std::endl;
			break;
		}
		
		case 'v':{
			std::stringstream param(optarg);
			param>>inivel;
			std::cout<<"Velocity="<<inivel<<std::endl;
			break;
		}
	}
	std::vector<Body*> bodies0;
	int iter=0;
    if (restore_config(path,config_file_name, bodies0,  iter,  theta,  G,  dt)) {
		std::cout <<"Resume at "<<iter <<  ", theta="<<theta<<", G="<< G <<", dt="<<  dt << ", size="<< bodies0.size() << std::endl;
		simulate(iter, max_iter, bodies0,  theta,  G,  dt,  img_iter, path,config_file_name);
	} else {
		std::cout << "Configuration file not found: starting from a new configuration" << std::endl;
		std::vector<Body*> bodies=createBodies(numbodies, inivel, ini_radius, mass );
		simulate(0, max_iter, bodies,  theta,  G,  dt,  img_iter, path,config_file_name);
	}
}

 /**
  * Create all bodies needed at start of run
  */
 std::vector<Body*>  createBodies(int numbodies,double inivel,double ini_radius,double mass ){
	std::cout << "Initializing " << numbodies << " bodies" << std::endl;
	std::vector<std::vector<double>> positions=direct_sphere(3,numbodies);
	std::vector<Body*> product;
	for (std::vector<std::vector<double>>::iterator it = positions.begin() ; it != positions.end(); ++it) {
        const double px = (*it)[0]* 2.*ini_radius + 0.5-ini_radius;
        const double py = (*it)[1]* 2.*ini_radius + 0.5-ini_radius;
		const double pz = (*it)[2]* 2.*ini_radius + 0.5-ini_radius;
        const double rpx = px-0.5;
        const double rpy = py-0.5;
		const double rpz = pz-0.5;
        const double rnorm = std::sqrt(sqr(rpx)+sqr(rpy)+sqr(rpz));
        const double vx = -rpy * inivel * rnorm / ini_radius;
        const double vy =  rpx * inivel * rnorm / ini_radius;
		const double vz = std::rand()%2==0 ? rpx : -rpx;
        product.push_back( new Body(mass, px, py, pz, vx, vy,vz) );
    }
	return product;
 }
 
  /**
  * Execute simulation
  */
 void simulate(int start_iter,int max_iter,std::vector<Body*> bodies, double theta, double G, double dt, int img_iter,std::string path,std::string config_file_name) {
	bool exiting=false;
    for (int iter=start_iter; iter<max_iter+start_iter && !exiting; ++iter) {
        Node* root = NULL;    // The oct-tree is recomputed at each iteration.
        for (unsigned i=0; i<bodies.size(); ++i) {
            bodies[i] -> resetToZerothQuadrant();
            root = add(bodies[i], root);
        }
 
        verlet(bodies, root, theta, G, dt); // Compute forces and advance bodies.
 
        delete root;  // De-allocate the oct-tree.
		
		exiting=killed();
        if (iter%img_iter==0||exiting) {
            std::cout << "Writing images at iteration " << iter << std::endl;
            save_bodies(bodies, iter/img_iter,path);
			save_config(bodies, iter, theta, G, dt,path,config_file_name);
        }
		
    }
}

double config_version=0.0;


 /**
  * Restore configuration from saved file
  */
bool restore_config(std::string path,std::string name,std::vector<Body*>& bodies, int& iter, double &theta, double &G, double &dt) {
	std::stringstream file_name;
    file_name << path<< name;
	std::ifstream config_file(file_name.str().c_str());

    if(! config_file.is_open())
        return false;
	enum State{expect_version, expect_iteration, expect_theta, expect_g, expect_dt, expect_body, expect_eof};
	State state=State::expect_version;
	while(! config_file.eof())   {
		std::string line;
        getline(config_file,line);
        std::stringstream ss(line);
		std::string token;
		switch(state) {
			case State::expect_version:
				token = line.substr(1+line.find("="));
				std::cout << token << std::endl;
				state=State::expect_iteration;
				break;
			case State::expect_iteration:
				token = line.substr(1+line.find("="));
				iter=atoi(token.c_str());
				std::cout << iter << std::endl;
				state=State::expect_theta;
				break;
			case State::expect_theta:
				token = line.substr(1+line.find("="));
				std::cout << decode(token) << std::endl;
				theta=decode(token);
				state=State::expect_g;
				break;
			case State::expect_g:
				token = line.substr(1+line.find("="));
				std::cout << decode(token) << std::endl;
				G=decode(token);
				state=State::expect_dt;
				break;
			case State::expect_dt:
				token = line.substr(1+line.find("="));
				std::cout << decode(token) << std::endl;
				dt=decode(token);
				state=State::expect_body;
				break;
			case State::expect_body:
				if (line.find("End")==0)
					state=State::expect_eof;
				else
					bodies.push_back(extract_body(line));
				break;
			case State::expect_eof:
				if (line.length()>0){
					std::cout<<"Unexpected text following end"<<std::endl;
					std::cout<<line<<std::endl;
					return false;
				}
			default:
				if (line.length()>0){
					std::cout<<"Unexpected state: "<<state<<std::endl;
					return false;
				}
		}
    }
	if (state!=State::expect_eof) {
		std::cout<<"Unexpected state: "<<state<<"-" <<State::expect_eof <<std::endl;
		return false;
	}
	return true;
}

/**
 * Retrieve position, mass, and velocities stored for one Body
 */
Body * extract_body(std::string line){
	enum State {expect_i,expect_x,expect_y,expect_z,expect_m,expect_vx,expect_vy,expect_vz,end_line};
	State state=expect_i;
	double px, py, pz, m, vx,vy,vz;
	while (state!=end_line) {
		int pos=line.find(",");
		std::string token=pos>=0 ? line.substr(0,pos) : line;
		line=line.substr(pos+1);
	
		switch (state){
			case expect_i:
				state=expect_x;
				break;
			case expect_x:
				state=expect_y;
				px=decode(token);
				break;
			case expect_y:
				state=expect_z;
				py=decode(token);
				break;
			case expect_z:
				state=expect_m;
				pz=decode(token);
				break;
			case expect_m:
				state=expect_vx;
				m=decode(token);
				break;
			case expect_vx:
				vx=decode(token);
				state=expect_vy;
				break;
			case expect_vy:
				vy=decode(token);
				state=expect_vz;
				break;
			case expect_vz:
				vz=decode(token);
				state=end_line;
		}
	}
	return new Body(m,px,py,pz,vx,vy,vz);
}


void save_config( std::vector<Body*>& bodies, int iter, double theta, double G, double dt, std::string path,std::string name) {
	std::stringstream file_name;
    file_name << path<< name;
/**
 * Save configuration so it can be restarted later
 */	backup(file_name.str().c_str());
    std::ofstream ofile(file_name.str().c_str());
	ofile << "Version="<<config_version<<"\n";
	ofile << "iteration=" << iter  << "\n";
	ofile << "theta=" << encode(theta)  << "\n";
	ofile << "G=" << encode(G)  << "\n";
	ofile << "dt=" << encode(dt)  << "\n";
    for (unsigned i=0; i<bodies.size(); ++i) {
		double px, py, pz;
		bodies[i] -> getPos(px, py,pz);
		double vx, vy,vz;
		bodies[i] -> getVel(vx, vy,vz);
		double m=bodies[i]->getMass();
		ofile <<i<<","<< encode(px)<<","<< encode(py)<<","<< encode(pz)<<","<< encode(m) <<","<< encode(vx)<<","<< encode(vy)<<","<< encode(vz)<<"\n";
	}
	ofile << "End\n";
}

/**
  * Generate help text
  */
void help(int numbodies,double inivel,double ini_radius,double mass,int max_iter,double theta, double G, double dt, int img_iter,std::string path,std::string config_file_name) {
	std::cout << "Galaxy Simulator based on Barnes Hut code from University of Geneva." << std::endl<<std::endl;
	std::cout << "Parameters, showing default values" <<std::endl;
	std::cout << "\t-c,--config\t\tConfugyration file [" << config_file_name<<"]"<< std::endl;
	std::cout << "\t-d,--dt\t\tTime Step for Integration [" << dt<<"]"<< std::endl;
	std::cout << "\t-G,--G\t\tGravitational Constant [" << G << "]"<<std::endl;
	std::cout << "\t-h,--help\tShow help text" << std::endl;
	std::cout << "\t-i,--img_iter\tFrequency for writing positions [" << img_iter << "]"<< std::endl;
	std::cout << "\t-m,--max_iter\tMaximum number of iterations [" << max_iter << "]"<< std::endl;
	std::cout << "\t-n,--numbodies\tNumber of bodies [" << numbodies<< "]"<<std::endl;
	std::cout << "\t-p,--path\tPath for writing configurations [" << path << "]"<< std::endl;
	std::cout << "\t-r,--ini_radius\tInitial Radius [" << ini_radius << "]"<<std::endl;
	std::cout << "\t-s,--mass\tMass of bodies [" << mass << "]"<<std::endl;
	std::cout << "\t-t,--theta\tTheta-criterion of the Barnes-Hut algorithm [" << theta << "]"<< std::endl;
	std::cout << "\t-v,--inivel\tInitial velocities [" << inivel << "]"<<std::endl;
}


 
 