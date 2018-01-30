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
 
#include <algorithm>
#include "barnes-hut.h"
#include "center-of-mass.h"

void get_acceleration_bh(std::vector<Particle*> particles,double theta) {
	Node * root=Node::create(particles);
	CentreOfMassCalculator calculator(particles);
	root->visit(calculator);
	calculator.display();
	std::for_each(particles.begin(),
				particles.end(),
				[root,theta](Particle*me){
					BarnesHutVisitor visitor(me,theta);
					root->visit(visitor);
					visitor.store_accelerations();
				});
	delete root;
}

bool BarnesHutVisitor::visit(Node * node) {
}

void BarnesHutVisitor::propagate(Node * node,Node * child){
}

bool BarnesHutVisitor::depart(Node * node) {
}

void BarnesHutVisitor::store_accelerations() {
	
}