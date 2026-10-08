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
Murray & Dermott, Exercise 1.4. Taking the orbital periods listed in Table A.9 but excluding Epimetheus,
Telesto, Calypso and Helen, use the criteria given in Section 1.7 to show there are 28 ratios of mean 
motions to consider in the Saturn System.
'''

from argparse import ArgumentParser
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import strftime,time
from matplotlib.pyplot import figure, show
import numpy as np
from md_data import create_data
from ratios import get_bounds,get_abc
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
          
def generate_pairs(data):
     '''
     Used to iterate through data, returning each pair once (unordered).
     '''
     def generate_periods():
          for i in range(len(data)):
               name, T = data[i]
               yield name,360/T
     
     for name1,n1 in generate_periods():
          for name2,n2 in generate_periods():
               if n1 < n2:
                    yield name1,name2,n1,n2
                    
def create_mean_motion_ratios(names_and_periods):
     '''
     Calculate ratios between periods for each pair of orbits
     
     Parameters:
         names_and_periods  A list of pairs, name and period.
         
     Returns:
        A list of tuples for each pair:
          First satellite in pair
          Second satellite in pair
          Integer p s.t. p/(p+1) is upper bound of ratio of mean motions
          Integer p' s.t. p'/(p'+1) is upper bound of ratio of mean motions
          Number calculate by M&D (1.22)
     '''
     product = []
     counts = dict()
     
     for name_i,name_j,n1,n2 in generate_pairs(names_and_periods):
          (p,p_prime) = get_bounds(n1,n2)
          try:
               counts[f'{p}-{p_prime}'] += 1
          except KeyError:
               counts[f'{p}-{p_prime}'] = 1          
          if p_prime == 0: continue       
          a,b,c = get_abc(n1,n2,p,p_prime)
          if abs(c) < 0.15:   
               product.append((name_i,name_j, p,p_prime, c))

     return product,len(counts)


def parse_args():
     parser = ArgumentParser(description='Calculate probability of delta_n<0.1 for problem 1.4.')
     parser.add_argument('--data', default='./data', help=f'Path to data files')
     parser.add_argument('--figs', default='./figs', help=f'Path to plots')
     parser.add_argument('--logs',default='./logs', help=f'Path to log files')
     parser.add_argument('--tolerance',default=0.15,type=float)
     parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
     parser.add_argument('--exclude',
                         default=['Epimetheus','Telesto','Calypso','Helene'],
                         nargs='*',
                         help='A list of satellites whose data are to be ignored') 
     parser.add_argument('--planet',default='saturn',help='File name for planetary data')
     parser.add_argument('--bins',default=100,type=int,help='Number of bins for histogram')
     return parser.parse_args()

def get_bar(mean_motion_ratios,tolerance=0.15):
     '''
     Used to display coloured bars in lowest panel
     '''
     c_indices = np.argsort([abs(c) for _,_,_,_,c in mean_motion_ratios])
     
     labels = []
     cc = []
     for i in c_indices:
          name_i,name_j,r0, r1, c = mean_motion_ratios[i]
          if r0 > 0 and abs(c) < tolerance:
               labels.append(f'{name_i}-{name_j}')
               cc.append(abs(c))
               #print (mean_motion_ratios[i])
     return labels,cc,['xkcd:red','xkcd:green','xkcd:blue','xkcd:yellow']

def main():
     args = parse_args()
     Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')

     mean_motion_ratios,count = create_mean_motion_ratios(
                                   create_data(
                                        (Path(args.data)/args.planet).with_suffix('.csv'),
                                        exclude=args.exclude ))

     labels,cc,bar_colours = get_bar(mean_motion_ratios,tolerance=args.tolerance)
     cs = sorted([abs(c) for _,_,_,_,c in mean_motion_ratios])

     fig = figure(figsize=(12, 12))
     ax1 = fig.add_subplot(2, 1, 1)
     ax1.hist(cs,bins=args.bins,color='skyblue', edgecolor='white' )
     ax1.set_title(f'Distribution of c. There are {count} distinct ratios')
     ax1.set_xlabel('c')
     ax2 = fig.add_subplot(2, 1, 2)

     ax2.bar(labels,cc,color=bar_colours)
     
     fig.tight_layout(h_pad=2)
     fig.savefig(Path(args.figs)/Path(__file__).stem)

     if args.show:
          show()

if __name__ == '__main__':
     main()
