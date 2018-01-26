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
#include "integrators.h"
#include "barnes_hut.h"

void Euler::step(double dt){
	for (std::vector<Body*>::iterator it = _bodies.begin() ; it != _bodies.end(); ++it) {
		double ax, ay, az;
		accelerationOn(*it, _root, _theta, ax, ay, az);
		ax *= _G;
		ay *= _G;
		az *= _G;
		double vx, vy, vz;
		(*it)->getVel(vx, vy, vz);
		vx+=(ax*dt); vy=(ay*dt); vz+=(az*dt);
		// std::cout << ax << ", " << vx << std::endl;
		(*it)->setVel(vx, vy, vz);
		double x, y, z;
		(*it)->getPos(x, y, z);
		x+=vx*dt; y+=vy*dt; z+=vz*dt;
		(*it)->setPos(x, y, z);
    }
}