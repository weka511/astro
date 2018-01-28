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
    Node()
	: _particle_index(Unused) {
		for (int i=0;i<N_Children;i++)
			_child[i]=NULL;
	}
	
	void insert(Particle * particle,int particle_index,std::vector<Particle*> particles);
	
	virtual ~Node() {
		for (int i=0;i<N_Children;i++)
			if (_child[i]!=NULL)
				delete _child[i];
	}
	
	static Node * create(std::vector<Particle*> particles);
  private:
	int _get_child_index(Particle * particle) {return -1;}     //TODO
	void _pass_down(int particle_index,int incumbent,std::vector<Particle*> particles) {;}    //TODO
	int _particle_index;
	Node * _child[N_Children];
};

#endif
