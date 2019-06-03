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

import random, matplotlib.pyplot as plt


def identify_commensurabilities(periods,tolerance=0.001,maxp=10):
     ratios = sorted([(i,j,periods[i]/periods[j]) for j in range(len(periods)) for i in range(j)],
                     key=lambda x:x[2])
     return [(i,j,ratio, p) for p in range(1,maxp+1) for i,j,ratio in ratios if abs(ratio-p/(p+1))<tolerance]

def monte_carlo(n,N,target_commensurabilities,T=20,seed=None,maxp=10):
     plt.figure()
     random.seed(seed)
     print (n,N,target_commensurabilities[0])
     _,_,_,target = target_commensurabilities[0]
     counts = []
     ps = []
     for i in range(N):
          periods            = [random.random()*T for i in range(n)] 
          commensurabilities = identify_commensurabilities(periods)
          counts.append(len(commensurabilities))
          if len(commensurabilities)==1:
               _,_,_,p=commensurabilities[0]
               ps.append(p)
     
     plt.subplot(211)         
     plt.hist(counts,bins=[x for x in range(0,max(counts)+3)])
     plt.title('Number of commensurabilities')

     plt.subplot(212)
     plt.title('Index of commensurabilities assuming number there is only one.')
     
     n,bins, _=plt.hist(ps,bins=[x + 0.1 for x in range(maxp+2)])
     i=0
     while bins[i]<target:
          i+=1
     return n[i-1]/sum(n)
     
if __name__=='__main__':
     N = 100
     from numpy import mean,std
     with open('data/1.5.txt') as data:
          periods            = [float(line.strip()) for line in data]
          commensurabilities = identify_commensurabilities(periods)
          for i,j,ratio,p in commensurabilities:
               print ( i,j,ratio,p,p/(p+1) )
          
          probs = [monte_carlo(len(periods),1000,commensurabilities,seed=None) for i in range(N)]
          plt.figure()
          plt.hist(probs)
          plt.title('N={0}, mean= {1:.3f}, std= {2:.3f}'.format(N,mean(probs),std(probs)))
          
     plt.show()
            