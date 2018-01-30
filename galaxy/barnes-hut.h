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
 

#ifndef _BARNES_HUT_H
#define _BARNES_HUT_H
#include <vector>
#include "particle.h"
#include "treecode.h"

void get_acceleration_bh(std::vector<Particle*>,double theta);


class BarnesHutVisitor :  public Node::Visitor{
  public:
	BarnesHutVisitor(Particle* me,const double theta) : _me(me),_theta(theta) {}
	virtual bool visit(Node * node);
	virtual void propagate(Node * node,Node * child);
	virtual bool depart(Node * node);
	void store_accelerations();
  private:
	Particle * _me;
	const double _theta;
};

#endif