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
 
 std::string encode(const double small);
 
 double decode(std::string str);
 
 /**
  * Check for presence of killfile
  */
 bool killed(std::string killfile="kill");
 
 /**
 * If file exists, copy to backup
 */
 void backup(std::string file_name, std::string backup="~");
 
 /**
  *  Sample points from hypersphere
  */
 std::vector<std::vector<double>> direct_sphere(int d=3,int n=1,double mean=0);
 
 
 
 
 