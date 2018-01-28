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
 
 #include <iostream>
 
 
 #include "verlet.h"

void get_acceleration_shm(std::vector<Particle*> particles){
	double x,y,z;
	particles[0]->getPos( x, y,  z);
	particles[0]->setAcc( -x, 0,  0);
}


bool print_values(std::vector<Particle*> particles) {
	double x,y,z;
	particles[0]->getPos( x, y,  z);
	double vx,vy,vz;
	particles[0]->getVel( vx, vy,  vz);
	std::cout<<x<<","<<vx<<"," << x*x + vx*vx<<std::endl;
	return true;
}

int main(int argc, char** argv){
    int max_iter = 100000; // Number of time-iterations executed by the program.
	double dt=0.01;
	std::vector<Particle*> particles;
	particles.push_back(new Particle(1,0));
	run_verlet(&get_acceleration_shm, max_iter, dt,	particles,	&print_values);
	for (std::vector<Particle*>::iterator it = particles.begin() ; it != particles.end(); ++it) 
		delete (*it);

	 return EXIT_SUCCESS;
}