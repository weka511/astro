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
	Node * product=new Node(-1);
	std::for_each(	particles.begin(),
					particles.end(),
					[product](Particle* particle){product->insert(particle);});
	return product;
}

void Node::insert(Particle * particle) {
	
}