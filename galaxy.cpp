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


struct option long_options[] =
{
  {"numbodies",  required_argument, 0, 'n'},
  {0, 0, 0, 0}
};	
	
int main(int argc, char **argv) {


    // Theta-criterion of the Barnes-Hut algorithm.
    double theta = 0.5;
    // Mass of a body.
    const double mass = 1.0;
    // Initially, the bodies are distributed inside a circle of radius ini_radius.
    const double ini_radius = 0.1;
    // Initial maximum velocity of the bodies.
    const double inivel = 0.1;
    // The "gravitational constant" is chosen so as to get a pleasant output.
    const double G = 4.e-6;
    // Discrete time step.
    const double dt = 1.e-3;
    // outside the initial radius are removed).
    int numbodies = 1000;
    // Number of time-iterations executed by the program.
    const int max_iter = 10000;
    // Frequency at which PNG images are written.
    const int img_iter = 20;
	
	std::string path = ".\\configs\\";
	
	int option_index = 0;
	int c;
	while ((c = getopt_long (argc, argv, "n:t:p:",long_options, &option_index)) != -1)
    switch (c){
		case 'n':{
			std::stringstream param(optarg);
			param>>numbodies;
			std::cout<<"Number of bodies="<<numbodies<<std::endl;
			break;
		}
			
		case 't':{
			std::stringstream param(optarg);
			param>>theta;
			std::cout<<"Theta="<<theta<<std::endl;
			break;
		}
		
		case 'p':{
			std::stringstream param(optarg);
			param>>path;
			std::cout<<"Path="<<path<<std::endl;
			break;
		}
	}
    
    
    // The pseudo-random number generator is initialized at a deterministic
    // value, for proper validation of the output for the exercise series.
    std::srand(1);
    // x- and y-pos are initialized to a square with side-length 2*ini_radius.
    std::vector<double> posx(numbodies), posy(numbodies), posz(numbodies);
    for (int i=0; i<numbodies; ++i) {
        posx[i] = ((double) std::rand() / (double)RAND_MAX) * 2.*ini_radius + 0.5-ini_radius;
        posy[i] = ((double) std::rand() / (double)RAND_MAX) * 2.*ini_radius + 0.5-ini_radius;
		posz[i] = ((double) std::rand() / (double)RAND_MAX) * 2.*ini_radius + 0.5-ini_radius;
    }
    // Initially, the bodies have a radial velocity of an amplitude proportional to
    // the distance from the center. This induces a rotational motion creating a
    // "galaxy-like" impression.
    std::vector<Body*> bodies;
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

    // Principal loop over time iterations.
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
