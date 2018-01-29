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
	  class Visitor {
		public:
		  virtual bool visit(Node * node)=0;
		  virtual void propagate(Node * node,Node * child){;}
	  };
  	enum Status {Internal=-2, Unused=-1};
	enum {N_Children=8};
    Node(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);
	
	void insert(int particle_index,std::vector<Particle*> particles);
	
	virtual ~Node() {
		for (int i=0;i<N_Children;i++)
			if (_child[i]!=NULL)
				delete _child[i];
	}
	
	bool visit(Visitor& visitor);
	
	int getStatus() { return _particle_index;}
	
	void getPhysics(double& m, double& x, double& y, double &z) {m=_m;x=_x;y=_y;z=_z;}
	
	void setPhysics(double m, double x, double y, double z) {_m=m,_x=x;_y=y;_z=z;}
	
	void accumulatePhysics(Node* other) {_m+=other->_m,_x+=other->_x;_y+=other->_y;_z+=other->_z;}
	
	static Node * create(std::vector<Particle*> particles);
	
	static get_limits(std::vector<Particle*> particles,double& xmin,double& xmax,double& ymin,double& ymax,double& zmin,double& zmax);
	
  private:
	int _get_child_index(int i, int j, int k) {return 4*i+2*j+k;}
	
	int _get_child_index(Particle * particle);
	
	void _pass_down(int particle_index,int incumbent,std::vector<Particle*> particles);
	
	void _insert_or_propagate(int particle_index,int incumbent,std::vector<Particle*> particles);
	
	void _split_node();
	
	int _particle_index;
	
	const double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax, _xmean, _ymean, _zmean;
	
	Node * _child[N_Children];
	
	double _m, _x, _y, _z;
};

class CentreOfMassCalculator : public Node::Visitor {
  public:
	CentreOfMassCalculator(std::vector<Particle*> particles);
	bool visit(Node * node);
	virtual void propagate(Node * node,Node * child);
	void display();
	
  private:
	std::vector<Particle*> _particles;
	std::vector<bool> indices;
};

#endif
