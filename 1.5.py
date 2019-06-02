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

def identify_commensurables(periods,tolerance=0.001):
     ratios = sorted([(i,j,periods[i]/periods[j]) for j in range(len(periods)) for i in range(j)],
                     key=lambda x:x[2])
     print (ratios)
     for p in range(1,11):
          for i,j,ratio in ratios:
               if abs(ratio-p/(p+1))<tolerance:
                    print (i,j,ratio,p/(p+1))
    
if __name__=='__main__':
     with open('data/1.5.txt') as data:
          periods = []
          for line in data:
               periods.append(float(line.strip()))
          identify_commensurables(periods)
            