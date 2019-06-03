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

import random, matplotlib.pyplot as plt, math


def get_mean_motion_ratios(data,max_p=100):
     cs = []
     for j in range(len(data)):
          _,T1 = data[j]
          n1   = 360/T1
          for i in range(j):
               _,T2    = data[i]
               n2      = 360/T2
               n_ratio = n1/n2
               p0      = 0
               p1      = None
               for p in range(max_p):
                    r = p/(p+1)
                    if r<n_ratio:
                         p0 = r
                    if r>n_ratio:
                         p1=r
                         break

               a = (n_ratio-p0)/(p1-p0)
               b = 0 if a<=0.5 else 1
               c = 2*math.pi * (a-b)
               print (p0,n_ratio,p1,a,b,c)
               cs.append((p,c))
     plt.hist([c for _,c in cs])
     plt.hist([abs(c) for _,c in cs if abs(c)>0.5])
     plt.title('Distribution of c')
     print (len(set([p for p,_ in cs])))
     
if __name__=='__main__':
     with open('data/saturn.dat') as data_file:
          data = []
          for line in data_file:
               parts = line.strip().split(',')
               name  = parts[0]
               T     = abs(float(parts[1]))
               data.append((name,T))
          get_mean_motion_ratios(data)
     plt.show()     
            