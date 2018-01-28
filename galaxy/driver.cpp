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

void get_acceleration(std::vector<Particle*> particles){

}


int main(int argc, char** argv){
    int max_iter = 10000; // Number of time-iterations executed by the program.
	double dt=0.0001;
	std::vector<Particle*> particles;
	particles.push_back(new Particle(0,0,1));
	particles.push_back(new Particle(0,0,2));
	get_acceleration(particles);
	std::for_each(particles.begin(),particles.end(),[dt](Particle* p){euler(p,0.5*dt);});
	for (int i=1;i<max_iter;i++) {
		get_acceleration(particles);
		std::for_each(particles.begin(),particles.end(),[dt](Particle* p){verlet(p,dt);});
	}
	for (std::vector<Particle*>::iterator it = particles.begin() ; it != particles.end(); ++it) 
		delete (*it);

	 return EXIT_SUCCESS;
}