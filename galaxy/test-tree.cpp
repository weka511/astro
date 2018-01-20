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