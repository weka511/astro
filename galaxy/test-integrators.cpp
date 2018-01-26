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
    
	SECTION("Simple manipulators for Body"){
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
		REQUIRE(ay1==0);
		REQUIRE(az1==0);
		double ax2, ay2, az2;
		accelerationOn(bodies[1], root, theta, ax2, ay2, az2);
		ax2*=G; ay2*=G; az2*=G;
		REQUIRE(ax2==Approx(-6.67400E-3));
		REQUIRE(ay2==0);
		REQUIRE(az2==0);
	}
 }