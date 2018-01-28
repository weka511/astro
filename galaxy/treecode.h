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
    Node(int body_index)
	: _body_index(body_index) {
		for (int i=0;i<8;i++)
			_child[i]=NULL;
	}
	
	void insert(Particle * particle);
	
	virtual ~Node() {
		for (int i=0;i<8;i++)
			if (_child[i]!=NULL)
				delete _child[i];
	}
	
	static Node * create(std::vector<Particle*> particles);
  private:
	int _body_index;
	Node * _child[8];
};

#endif
