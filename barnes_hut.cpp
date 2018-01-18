/* Simple C++ implementation of a Barnes-Hut galaxy simulator.
 * This file is part of the exercise series of the University of Geneva's
 * MOOC "Simulation and Modeling of Natural Processes".
 *
 * Author: Jonas Latt
 * E-mail contact: jonas.latt@unige.ch
 * Important: don't send questions about this code to the above e-mail address.
 * They will remain unanswered. Instead, use the resources of the MOOC.
 * 
 * Copyright (C) 2016 Université de Genève
 * 24 rue du Général-Dufour
 * CH - 1211 Genève 4
 * Switzerland
 *
 * This code is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

#include "tree.h"


// Barnes-Hut algorithm: Creation of the quad-tree. This function adds
// a new body into a quad-tree node. Returns an updated version of the node.
Node* add(Body* body, Node* node) {
    // To limit the recursion depth, set a lower limit for the size of quadrant.
  static const double smallest_quadrant = 1.e-4;
  Node* new_node = 0;
  // 1. If node n does not contain a body, put the new body b here.
  if (node==0)
    new_node = body;
  else {
    if (node->getS() < smallest_quadrant)
      return node;
    
    // 3. If node n is an external node, then the new body b is in conflict
    //    with a body already present in this region. ...
    if (node->isEndnode()) {
      //    ... Subdivide the region further by creating an internal node.
      new_node = new InternalNode(node);
      //    ... And to start with, insert the already present body recursively
      //        into the appropriate quadrant.
      const int quadrant = node->intoNextQuadrant();
      new_node->setChild(quadrant, node);
    }
    else // 2. If node n is an internal node, we don't to modify its child. 
      new_node = node;

    // 2. and 3. If node n is or has become an internal node ...
    //           ... update its mass and "center-of-mass times mass"
    new_node->addMassCom(body);
    // ... and recursively add the new body into the appropriate quadrant.
    const int quadrant = body->intoNextQuadrant();
    new_node->setChild (
			quadrant,
			add(body, new_node->extractChild(quadrant)) );
  }
  return new_node;
}

// Compute the force of all other nodes onto a given node, divided by
// the node's mass, and divided by the gravitational constant G.
// This amounts to a recursive evaluation of the quad-tree created by
// the Barnes-Hut algorithm.
void accelerationOn( Body const* body, Node const* node, double theta,
                     double& ax, double& ay, double& az){
  // 1. If the current node is an external node, 
  //    calculate the force exerted by the current node on b.
  const double dsqr = node->distSqr(body);
  if (node->isEndnode())
    node->accelerationOn(body, ax, ay, az,dsqr);

  // 2. Otherwise, calculate the ratio s/d. If s/d < Î¸, treat this internal
  //    node as a single body, and calculate the force it exerts on body b.
  else if (sqr(node->getS()) < dsqr*sqr(theta))
    node->accelerationOn(body, ax, ay, az,dsqr);
  
  else { // 3. Otherwise, run the procedure recursively on each child.
    ax = 0.;
    ay = 0.;
    az = 0.;
    for (int i=0; i<n_children; ++i) {
      Node const* c = node->getChild(i);
      if (c!=0) {
	double ax_, ay_,az_;
	accelerationOn(body, c, theta, ax_, ay_,az_);
	ax += ax_;
	ay += ay_;
	az += az_;
      }
    }
  } 
}

// Execute a time iteration according to the Verlet algorithm.
void verlet( std::vector<Body*>& bodies, Node* root,
             double theta, double G, double dt ) {
    for(size_t i=0; i<bodies.size(); ++i) {
      double ax, ay, az;
      accelerationOn(bodies[i], root, theta, ax, ay, az);
        ax *= G;
        ay *= G;
	az *= G;
        bodies[i]->advance(ax, ay, az,dt);
    }
}


// Write the position of all bodies into a text file.
// The text file can be converted into an image with the
// Python script make_img.py
// Batch-processing of all text files is achieved with the
// shell script dat2img.
void save_bodies( std::vector<Body*>& bodies, int i, std::string path){
    std::stringstream fNameStream;
    fNameStream << path<< "body_" << std::setfill('0') << std::setw(6) << i << ".dat";
    std::ofstream ofile(fNameStream.str().c_str());
    for (unsigned i=0; i<bodies.size(); ++i) {
      double px, py, pz;
        bodies[i] -> getPos(px, py,pz);
        ofile << std::setprecision(12)
              << std::setw(20) << px
	      << std::setw(20) << py
              << std::setw(20) << pz << "\n";
    }
}

int main() {
    // Theta-criterion of the Barnes-Hut algorithm.
    const double theta = 0.5;
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
    const int numbodies = 1000;
    // Number of time-iterations executed by the program.
    const int max_iter = 10000;
    // Frequency at which PNG images are written.
    const int img_iter = 20;

    std::string path = ".\\configs\\";
    
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

