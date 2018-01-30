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
#include <iostream>
#include <limits>

Node::Node(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
: _particle_index(Unused),
	_xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax),
	_xmean(0.5*(xmin+ xmax)), _ymean(0.5*(ymin+ ymax)), _zmean(0.5*(zmin+ zmax)),
	_m(0.0d),_x(0.0d),_y(0.0d),_z(0.0d) {
	for (int i=0;i<N_Children;i++)
		_child[i]=NULL;
}

Node::get_limits(std::vector<Particle*> particles,double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax){
	xmin=std::numeric_limits<double>::max();
	xmax=-xmin;
	ymin=std::numeric_limits<double>::max();
	ymax=-ymin;
	zmin=std::numeric_limits<double>::max();
	zmax=-zmin;
	std::for_each(particles.begin(),
					particles.end(),
					[&xmin,&xmax,&ymin,&ymax,&zmin,&zmax](Particle* particle){
						double x,y,z;
						particle->getPos(x,y,z);
						if (x<xmin) xmin=x;
						if (x>xmax) xmax=x;
						if (y<ymin) ymin=y;
						if (y>ymax) ymax=y;
						if (z<zmin) zmin=z;
						if (z>zmax) zmax=z;
					});
	xmin=ymin=zmin=std::min(xmin,std::min(ymin,zmin));  // because Barnes Hut requires cubes
	xmax=ymax=zmax=std::max(xmax,std::max(ymax,zmax));
}
Node * Node::create(std::vector<Particle*> particles){
	double xmin, xmax, ymin, ymax, zmin, zmax;
	Node::get_limits(particles,xmin, xmax, ymin, ymax, zmin, zmax);
	Node * product=new Node(xmin,xmax,ymin,ymax,zmin,zmax);
	for (int index=0;index<particles.size();index++)
		product->insert(index,particles);
	return product;
}

void Node::insert(int particle_index,std::vector<Particle*> particles) {
	switch(_particle_index){
		case Unused:
			_particle_index=particle_index;
			return;
		case Internal:
			_child[_get_child_index(particles[particle_index])]->insert(particle_index,particles);
			return;
		default:     //oops - we already have a particle here, so have to move it
			int incumbent=_particle_index;
			_split_node();
			_insert_or_propagate(particle_index,incumbent,particles);
	}
}

void Node::_pass_down(int particle_index,int incumbent,std::vector<Particle*> particles) {
	_split_node();
	_insert_or_propagate(particle_index,incumbent,particles);
} 

void Node::_insert_or_propagate(int particle_index,int incumbent,std::vector<Particle*> particles) {
	int child_index_new=_get_child_index(particles[particle_index]);
	int child_index_incumbent=_get_child_index(particles[incumbent]);
	if (child_index_new==child_index_incumbent)
		_child[child_index_incumbent]->_pass_down(particle_index,incumbent,particles);
	else {
		_child[child_index_new]->insert(particle_index,particles);
		_child[child_index_incumbent]->insert(incumbent,particles);
	}
}

int Node::_get_child_index(Particle * particle) {
	double x,y,z;
	particle->getPos(x,y,z);
	const int i=x>_xmean;
	const int j=y>_ymean;
	const int k=z>_zmean;
	return _get_child_index(i,j,k);
}

void Node::_split_node() {
	_particle_index=Internal;ode:
	double xmin, xmax, ymin, ymax, zmin, zmax;
	for (int i=0;i<2;i++) {
		if (i==0) {
			xmin=_xmin;
			xmax=_xmean;
		} else {
			xmin=_xmean;
			xmax=_xmax;
		}
		for (int j=0;j<2;j++) {
			if (j==0) {
				ymin=_ymin;
				ymax=_ymean;
			} else {
				ymin=_ymean;
				ymax=_ymax;
			}
			for (int k=0;k<2;k++) {
				if (k==0) {
					zmin=_zmin;
					zmax=_zmean;
				} else {
					zmin=_zmean;
					zmax=_zmax;
				}
				_child[_get_child_index(i,j,k)]=new Node(xmin, xmax, ymin, ymax, zmin, zmax);
			}	// k
		}		// j
	}			// i
} 

bool Node::visit(Visitor & visitor) {
	Node::Visitor::Status status=visitor.visit(this);
	bool should_continue=status != Node::Visitor::Status::Stop;
	if (_particle_index==Internal && status==Node::Visitor::Status::Continue)
		for (int i=0;i<N_Children&&should_continue;i++) {
		should_continue=_child[i]->visit(visitor);
		visitor.propagate(this,_child[i]);
	}
	return should_continue ? visitor.depart(this) : false;
}


