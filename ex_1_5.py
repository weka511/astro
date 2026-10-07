#!/usr/bin/env python

# Copyright (C) 2019-2026 Greenweaves Software Limited

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

'''
Murray & Dermott, Exercise 1.5. We are given the periods of six planets orbiting a star. 
By considering the fifteen possible ratios and the ten first order commensurabilities, 
identify the pairs of planets such that ratios ae within 0.0001 of the commensurabilities. 
Estimate probability of this occurring by chance if the periods were randomly distributed. 
'''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import strftime
from matplotlib.pyplot import figure, show
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'

class Logger: 
     instance = None

     @staticmethod
     def create(path_name:str):
          '''
          Set up console logger and file logger

          Parameters:
              path_name    Path name for log file
          '''
          def add_handler(handler,logger,formatter=Formatter('%(message)s'),level:int=INFO):
               handler.setLevel(level)
               handler.setFormatter(formatter)
               logger.addHandler(handler)    
          Logger.instance = getLogger(__name__)
          Logger.instance.setLevel(DEBUG)
          add_handler(FileHandler(path_name),Logger.instance,level=DEBUG)
          add_handler(StreamHandler(),Logger.instance)   
          np.set_printoptions(linewidth=np.nan) # Prevent lines being split when we log numpy arrays
          
def identify_commensurabilities(periods, atol=0.001, p_max=10):
     '''
     Find the ratios that are close to integer ratios p/(p+1)
     
     Parameters:
         periods  Orbital periods, used to calculate ratios
         atol     Tolerance: two numbers are deemed equal if they are within atol of each other
         p_max    Maximum value of p in p(p+1)
         
     Returns:
        List of tuples (i, j, r, p), where i and j are the indices of the periods that define the ratio,
                                     r is the ratio, and p the corresponding integer such that p(p+1)
                                     is within atol of r.
     '''   
     ratios = sorted([(i, j, periods[i] / periods[j]) 
                         for j in range(len(periods)) 
                         for i in range(j)],
                     key=lambda x: x[2])
     for i,j,r in ratios:
          Logger.instance.debug(f'{i}, {j}, {r}')
     commensurabilities = [p/(p+1) for p in range(1,p_max+1)]
     for c in commensurabilities:
          Logger.instance.debug(f'{c}')
                               
     for i, j, r in ratios:
          distances = [abs(r-c) for c in commensurabilities]
          index = np.argmin(distances)     
 
     return [(i, j, r, p) for p in range(1, p_max + 1) 
               for i, j, r in ratios
               if abs(r - p / (p + 1)) < atol]


def sample(m, N, target, 
                T=20, p_max=10, rng=np.random.default_rng()):
     '''
     Sample possible orbits
     
     Parameters:
         m         Number of periods
         N         Number of samples
         target    The commensurabilities identified in dataset
         T         Upper bound on orbital periods
         p_max     Largest integer to be considered in ratios p/(p+1)
         rng       Random number generator
     '''
     def sample_single():
          ps = []
          for i in range(N):
               commensurabilities = identify_commensurabilities(rng.uniform(low=0,high=T,size=m))
               if len(commensurabilities) == 1:
                    _, _, _, p = commensurabilities[0]
                    ps.append(p)
          return ps     
     
     _, _, _, target_p = target[0]
 
     ps = sample_single()
     counts = np.zeros((p_max+1))
     for p in ps:
          counts[p] += 1
     return counts[target_p] / counts.sum()    


def parse_args():
     parser = ArgumentParser(description='Calculate probability for commensurabilities for problem 1.5.')
     parser.add_argument('-M', type=int, help='Number of samples per Monte Carlo run', default=100)
     parser.add_argument('-N', type=int, help='Number of Monte Carlo runs', default=100)
     parser.add_argument('--seed', '-s', help='Seed for random number generator', default=None)
     parser.add_argument('--figs', default='./figs', help=f'Path to plots')
     parser.add_argument('--data', default='./data', help=f'Path to data files')
     parser.add_argument('--commensurability',default='commensurability',help='File name for planetary data')
     parser.add_argument('--show', default=False, action='store_true', help='Controls whether plot will be displayed')
     parser.add_argument('--logs',default='./logs', help=f'Path to log files')
     return parser.parse_args()

def get_data(path):
     '''
     Load periods from file
     
     Parameters:
         path     Identifies data file
     '''
     with open(path) as data:
          data_reader = reader(data)
          return np.array([float(row[0]) for row in data_reader])
 
def get_title(commensurabilities):
     '''
     Display the commensurabilities
     '''
     (i, j, r, p) = commensurabilities[0]
     return f'Planets {i+1} & {j+1} have ratio {r:.5f}, within {abs(r-p/(p+1)):.5f} of {p}:{p+1}'

def main():
     '''
     Identify commensurability and estimate probability of value occurring by chance.
     '''
     args = parse_args()
     Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
     rng = np.random.default_rng(args.seed)
     periods = get_data((Path(args.data)/args.commensurability).with_suffix('.csv'))
     commensurabilities = identify_commensurabilities(periods)

     probs = [sample(len(periods), args.M, commensurabilities, rng=rng) for i in range(args.N)]
     fig = figure()
     fig.suptitle(f'Exercise 1.5: There is {len(commensurabilities)} commensurability')
 
     ax = fig.add_subplot(1,1,1)
     ax.hist(probs,
             color='xkcd:blue',
             bins='fd',
             density=True,
             label=f'N={args.N},mean={np.mean(probs):.4f},std={np.std(probs):.4f}')
             
     ax.set_title(get_title(commensurabilities))
     ax.set_xlabel(r'$\frac{p}{p+1}$')
     ax.legend()
     fig.savefig(Path(args.figs)/Path(__file__).stem)

     if args.show:
          show()

if __name__ == '__main__':
     main()

