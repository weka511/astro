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
#include "tree.h"
#include "catch.hpp"
#include "integrators.h"
#include "barnes_hut.h"

 TEST_CASE( "Tests for Verlet integrator", "[verlet]" ) {
    
	/**
	 * Test bsed on assignment
	 * https://www.coursera.org/learn/modeling-simulation-natural-processes/exam/HFO5I/particles-and-point-like-objects
	 */
	SECTION("Test Euler starting step"){
		std::vector<Body*> bodies;
		bodies.push_back(new Body (1e20,0,0,0,0,0,0));
		bodies.push_back(new Body (1e12,1e6,0,0,0,0,0));
		double theta=0.5;
		double G = 6.6740e-11;
		double dt=1000;
		Node* root = NULL;
		for (unsigned i=0; i<bodies.size(); ++i) {
            bodies[i] -> resetToZerothOctant();
            root = add(bodies[i], root);
        }
		double ax1, ay1, az1;
		accelerationOn(bodies[0], root, theta, ax1, ay1, az1);
		ax1*=G; ay1*=G; az1*=G;
		REQUIRE(ax1==6.67400E-11);
		double ax2, ay2, az2;
		accelerationOn(bodies[1], root, theta, ax2, ay2, az2);
		ax2*=G; ay2*=G; az2*=G;
		REQUIRE(ax2==Approx(-6.67400E-3));
		Euler euler(bodies,G,root,theta);
		euler.step(1000);
		double vx1, vy1, vz1;
		bodies[0]->getVel(vx1, vy1, vz1);
		REQUIRE(vx1==Approx(6.6740800E-08).epsilon(0.0001));
		double vx2, vy2, vz2;
		bodies[1]->getVel(vx2, vy2, vz2);
		REQUIRE(vx2==Approx(-6.6740800E+00).epsilon(0.0001));
		double x1,y1,z1;
		bodies[0]->getPos(x1, y1, z1);
		REQUIRE(x1==Approx(0.00006674));
		double x2,y2,z2;
		bodies[1]->getPos(x2, y2, z2);
		REQUIRE(x2==Approx(993325.92000000));
	}
 }