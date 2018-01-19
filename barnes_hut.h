/* Simple C++ implementation of a Barnes-Hut galaxy simulator.
 * This file is part of the exercise series of the University of Geneva's
 * MOOC "Simulation and Modeling of Natural Processes".
 *
 * Author: Jonas Latt
 * E-mail contact: jonas.latt@unige.ch
 * Important: don't send questions about this code to the above e-mail address.
 * They will remain unanswered. Instead, use the resources of the MOOC.
 * 
 * Copyright (C) 2016 Université de Genève
 * 24 rue du Général-Dufour
 * CH - 1211 Genève 4
 * Switzerland
 *
 * This code is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

Node* add(Body* body, Node* node);

void accelerationOn( Body const* body, Node const* node, double theta,
                     double& ax, double& ay, double& az);

void verlet( std::vector<Body*>& bodies, Node* root,
             double theta, double G, double dt );

void save_bodies( std::vector<Body*>& bodies, int i, std::string path);
			 