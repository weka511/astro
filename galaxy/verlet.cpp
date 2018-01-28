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
 
 #include <cstdlib>
 #include <iostream>
 #include <algorithm>
 
 #include "verlet.h"



void  euler(Particle* particle,double dt){
	std::cout << "Euler: " <<dt<< std::endl;
	double vx,vy,vz;
	particle->getVel( vx, vy,  vz);
	double ax,ay,az;
	particle->getAcc( ax, ay,  az);
	particle->setVel( vx+dt*ax, vy+dt*ay,  vz+dt*az);
}

void  verlet(Particle* particle,double dt) {
	std::cout << "Verlet: " <<dt<< std::endl;
	double x,y,z;
	particle->getPos( x, y,  z);
	double vx,vy,vz;
	particle->getVel( vx, vy,  vz);
	double ax,ay,az;
	particle->getAcc( ax, ay,  az);
	particle->setPos( x+dt*vx, y+dt*vy,  z+ dt*vz);
	particle->setVel( vx+dt*ax, vy+dt*ay,  vz+dt*az);
}

