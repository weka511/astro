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
     
     # get_bounds
     #
     # calculate p and p' from Murray and Dermott, Section 1.7.
     # We use r0 for the lower bound (MD p'/(p'+1)) and r1 for the upper 
     def get_bounds(n_ratio):
          if 1/3<n_ratio and n_ratio < 1/2: # See MD Section 1.7
               return (1/3,1/2)
          r0      = 0
          r1      = None
          for p in range(max_p):
               r = p/(p+1)
               if r<n_ratio:
                    r0 = r
               if r>n_ratio:
                    r1 = r
                    return (r0,r1)
               
     rcs = []
     for j in range(len(data)):
          _,T1 = data[j]
          n1   = 360/T1
          for i in range(j):
               _,T2    = data[i]
               n2      = 360/T2
               (r0,r1) =get_bounds(n1/n2)
               a = (n1/n2-r0)/(r1-r0)   # Murray & Dermott (1.19)
               b = 0 if a<=0.5 else 1   # Murray & Dermott (1.20)
               c = 2 * math.pi * (a-b)    # Murray & Dermott (1.21)
               rcs.append(((r0,r1),c))
     cs = sorted([abs(c) for _,c in rcs])
     for i in range(len(cs)):
          print (i,cs[i])
     plt.hist(cs,bins=100)
     #n,bins,_=plt.hist([abs(c) for _,c in cs],bins=100,cumulative=True)     
     plt.title('Distribution of c')
     print (len(set([r for r,_ in rcs])))
     
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
            