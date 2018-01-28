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
 
#include "treecode.h"
#include <algorithm>

Node * Node ::create(std::vector<Particle*> particles){
	Node * product=new Node();
	for (int index=0;index<particles.size();index++)
		product->insert(particles[index],index,particles);
	return product;
}

void Node::insert(Particle * particle,int particle_index,std::vector<Particle*> particles) {
	switch(_particle_index){
		case Unused:
			_particle_index=particle_index;
			return;
		case Internal:
			_child[_get_child_index(particle)]->insert(particle,particle_index,particles);
			return;
		default:
			int incumbent=_particle_index;
			_particle_index=Internal;
			for (int i=0;i<N_Children;i++)
				_child[i]=new Node();
			int child_index_new=_get_child_index(particle);
			int child_index_incumbent=_get_child_index(particles[incumbent]);
			if (child_index_new==child_index_incumbent)
				_child[child_index_incumbent]->_pass_down(particle_index,incumbent,particles);
			else {
				_child[child_index_new]->insert(particle,particle_index,particles);
				_child[child_index_incumbent]->insert(particles[incumbent],incumbent,particles);
			}
			

	}
}
