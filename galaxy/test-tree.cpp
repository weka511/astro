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

TEST_CASE( "Tests for Node and its descendants", "[node]" ) {
    
	SECTION("Simple manipulators for Body"){
		Body b(6.5,1,2,3,0,0,0);
		double x,y,z;
		b.getPos(x,y,z);
		REQUIRE(x==1);
		REQUIRE(y==2);
		REQUIRE(z==3);
		double mx,my,mz;
		b.get_m_pos(mx,my,mz);
		REQUIRE(mx==6.5);
		REQUIRE(my==13);
		REQUIRE(mz==19.5);
	}

	SECTION("Verlet Integration"){
		Body b(6.5,1,2,3,4,5,6);
		b.advance(0.1,0.01,0.001, 3);
		double vx,vy,vz;
		b.getVel(vx,vy,vz);
		REQUIRE(vx==4+3*0.1);
		REQUIRE(vy==5+3*0.01);
		REQUIRE(vz==6+3*0.001);
		double x,y,z;
		b.getPos(x,y,z);
		REQUIRE(x==1+3*(4+3*0.1));
		REQUIRE(y==2+3*(5+3*0.01));
		REQUIRE(z==3+3*(6+3*0.001));
	
	}
}