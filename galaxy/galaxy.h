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

 /**
  * Create all bodies needed at start of run
  */
 void createBodies(int numbodies,double inivel,double ini_radius,double mass,std::vector<Body*>&);
 
 /**
  * Execute simulation
  */
 void simulate(int max_iter,std::vector<Body*> bodies, double theta, double G, double dt, int img_iter,std::string path);
 
 /**
  * Generate help text
  */
 void help();
 
 /**
  * Check for presence of killfile
  */
 bool killed(std::string killfile="kill");
 
 /**
  *  Sample points from hypersphere
  */
 std::vector<std::vector<double>> direct_sphere(int d=3,int n=1,double mean=0);
 
 
 
 
