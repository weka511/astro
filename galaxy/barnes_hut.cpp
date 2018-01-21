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
#include "barnes_hut.h"

// Barnes-Hut algorithm: Creation of the quad-tree. This function adds
// a new body into a quad-tree node. Returns an updated version of the node.
Node* add(Body* body, Node* node) {
    // To limit the recursion depth, set a lower limit for the size of quadrant.
  static const double smallest_quadrant = 1.e-4;
  Node* new_node = NULL;
  // 1. If node n does not contain a body, put the new body b here.
  if (node==NULL)
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


