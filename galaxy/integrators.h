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
 
#ifndef _INTEGRATORS_H
#define _INTEGRATORS_H

 class Euler {
	public:
		Euler(std::vector<Body*>& bodies,double G,Node* root,double theta) : _bodies(bodies),_G(G),_root(root),_theta(theta) {;}
		void step(double dt);
	private:
		std::vector<Body*>& _bodies;
		double 				_G;
		Node* 				_root;
		double 				_theta;
 };
 
 class Verlet {
	 
 };
 
 #endif
 