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
#include "tree.h"
#include "barnes_hut.h"
#include "galaxy.h"

struct option long_options[] = {
	{"dt",  required_argument, 0, 'd'},
	{"G",  required_argument, 0, 'G'},
    {"help",  no_argument, 0, 'h'},
	{"img_iter",  required_argument, 0, 'i'},
	{"max_iter",  required_argument, 0, 'm'},
	{"numbodies",  required_argument, 0, 'n'},
	{"path",  required_argument, 0, 'p'},
	{"ini_radius",  required_argument, 0, 'r'},
	{"mass",  required_argument, 0, 's'},
	{"theta",  required_argument, 0, 't'},
	{"inivel",  required_argument, 0, 'v'},
	{0, 0, 0, 0}
};	
	
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
	
	std::string path = "./configs/";
	
	int option_index = 0;
	int c;
	

	while ((c = getopt_long (argc, argv, "d:G:hi:m:n:p:r:s:t:v:",long_options, &option_index)) != -1)
    switch (c){
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
			help();
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
    
    std::vector<Body*> bodies;
	createBodies(numbodies, inivel, ini_radius, mass,bodies );
	simulate( max_iter, bodies,  theta,  G,  dt,  img_iter, path);
 
}


 void createBodies(int numbodies,double inivel,double ini_radius,double mass,std::vector<Body*>& bodies ){
	    std::srand(1);
    // x- and y-pos are initialized to a square with side-length 2*ini_radius.
    std::vector<double> posx(numbodies), posy(numbodies), posz(numbodies);
	
    for (int i=0; i<numbodies; ++i) {
        posx[i] = ((double) std::rand() / (double)RAND_MAX) * 2.*ini_radius + 0.5-ini_radius;
        posy[i] = ((double) std::rand() / (double)RAND_MAX) * 2.*ini_radius + 0.5-ini_radius;
		posz[i] = ((double) std::rand() / (double)RAND_MAX) * 2.*ini_radius + 0.5-ini_radius;
    }
	
	for (int i=0; i<numbodies; ++i) {
        const double px = posx[i];
        const double py = posy[i];
		const double pz = posz[i];
        const double rpx = px-0.5;
        const double rpy = py-0.5;
        const double rnorm = std::sqrt(sqr(rpx)+sqr(rpy));
        if ( rnorm < ini_radius ) {
            const double vx = -rpy * inivel * rnorm / ini_radius;
            const double vy =  rpx * inivel * rnorm / ini_radius;
			const double vz = 0;
            bodies.push_back( new Body(mass, px, py, pz, vx, vy,vz) );
        }
    }
 }
 
 void simulate(int max_iter,std::vector<Body*> bodies, double theta, double G, double dt, int img_iter,std::string path) {

    for (int iter=0; iter<max_iter; ++iter) {
        // The quad-tree is recomputed at each iteration.
        Node* root = 0;
        for (unsigned i=0; i<bodies.size(); ++i) {
            bodies[i] -> resetToZerothQuadrant();
            root = add(bodies[i], root);
        }
        // Computation of forces, and advancement of bodies.
        verlet(bodies, root, theta, G, dt);
        // De-allocate the quad-tree.
        delete root;
        // Output.
        if (iter%img_iter==0) {
            std::cout << "Writing images at iteration " << iter << std::endl;
            save_bodies(bodies, iter/img_iter,path);
        }
    }
}

void help() {
	std::cout << "Galaxy Simulator" << std::endl;
	std::cout << "\t-d,--dt\t\tTime Step for Integration" << std::endl;
	std::cout << "\t-G,--G\t\tGravitational Constant" << std::endl;
	std::cout << "\t-h,--help\tHelp text" << std::endl;
	std::cout << "\t-i,--img_iter\tFrequency for writing positions" << std::endl;
	std::cout << "\t-m,--max_iter\tMaximum number of iterations" << std::endl;
	std::cout << "\t-n,--numbodies\tNumber of bodies" << std::endl;
	std::cout << "\t-p,--path\tPath for writing configurations" << std::endl;
	std::cout << "\t-r,--ini_radius\tInitial Radius" << std::endl;
	std::cout << "\t-s,--mass\tMass of bodies" << std::endl;
	std::cout << "\t-t,--theta\tTheta-criterion of the Barnes-Hut algorithm" << std::endl;
	std::cout << "\t-v,--inivel\tInitial velocities" << std::endl;
}
