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

import random, matplotlib.pyplot as plt, numpy as np


def identify_commensurabilities(periods,tolerance=0.001,maxp=10):
     ratios = sorted([(i,j,periods[i]/periods[j]) for j in range(len(periods)) for i in range(j)],
                     key=lambda x:x[2])
     return [(i,j,ratio, p) for p in range(1,maxp+1) for i,j,ratio in ratios if abs(ratio-p/(p+1))<tolerance]

def monte_carlo(n,N,target_commensurabilities,T=20,seed=None,maxp=10):
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
     
     n,bins = np.histogram(ps,bins=[x + 0.1 for x in range(maxp+2)])
     i=0
     while bins[i]<target:
          i+=1
     return n[i-1]/sum(n)
     
if __name__=='__main__':
     import argparse
     parser = argparse.ArgumentParser(description='Calculate probability for commensurabilities for problem 1.5.')    
     parser.add_argument('-N',type=int, help='Number of Monte Carlo calculations',default=100)
     parser.add_argument('--seed','-s',help='Seed for random number generator',default=None)
     args = parser.parse_args()
     random.seed(args.seed)
     
     with open('data/commensurability.dat') as data:
          periods            = [float(line.strip()) for line in data]
          commensurabilities = identify_commensurabilities(periods)
 
          probs = [monte_carlo(len(periods),1000,commensurabilities,seed=None) for i in range(args.N)]
          plt.figure()
          plt.hist(probs)
          plt.title('N={0}, mean= {1:.3f}, std= {2:.3f}'.format(args.N,np.mean(probs),np.std(probs)))
          
     plt.show()
            