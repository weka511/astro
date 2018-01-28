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
 
#ifndef _TREECODE_H
#define _TREECODE_H

#include <cstdlib>
#include <vector>
#include "verlet.h"
#include "particle.h"

class Node {
  public:
  	enum Status {Internal=-2, Unused=-1};
	enum {N_Children=8};
    Node(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
	: _particle_index(Unused),
		_xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax),
		_xmean(0.5*(xmin+ xmax)), _ymean(0.5*(ymin+ ymax)), _zmean(0.5*(zmin+ zmax))	{
		for (int i=0;i<N_Children;i++)
			_child[i]=NULL;
	}
	
	void insert(int particle_index,std::vector<Particle*> particles);
	
	virtual ~Node() {
		for (int i=0;i<N_Children;i++)
			if (_child[i]!=NULL)
				delete _child[i];
	}
	
	static Node * create(std::vector<Particle*> particles);
  private:
	int _get_child_index(int i, int j, int k) {return 4*i+2*j+k;}
	
	int _get_child_index(Particle * particle);
	
	void _pass_down(int particle_index,int incumbent,std::vector<Particle*> particles);
	
	void _insert_or_propagate(int particle_index,int incumbent,std::vector<Particle*> particles);
	
	void _split_node();
	
	int _particle_index;
	
	const double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax, _xmean, _ymean, _zmean;
	
	Node * _child[N_Children];
};

#endif
