# Copyright (C) 2019 Greenweaves Software Limited

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this software.  If not, see <http://www.gnu.org/licenses/>

def get_average_interval():
     pass

if __name__=='__main__':
     with open('data/earth_mars.dat') as data_file:
          data = []
          for line in data_file:
               parts = line.strip().split(',')
               name  = parts[0]
               a     = abs(float(parts[1]))
               data.append((name,a))
               
          get_average_interval()   
            