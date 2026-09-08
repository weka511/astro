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
Murray & Dermott, Exercise 1.4. Taking the orbital periods lisyted in Table A.9 but excluding Epimetheus,
Telesto, Calypso and Helen, use the criteria geiven in Section 1.7 to show there are 28 ratios of mean 
motions to consider in the Saturn System.
'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from matplotlib.pyplot import figure, show
import numpy as np


def get_bounds(n_ratio):
     '''
     calculate p and p' from Murray and Dermott, Section 1.7.
     We use r0 for the lower bound (MD p'/(p'+1)) and r1 for the upper
     
     Parameters:
         n_ratio
         
     Returns:
         r0,r1, integer ratios, r0 < n_ration < r1
     '''  
     r0 = 0
     r1 = 0
     p = 1

     while r1 < n_ratio:
          r0 = r1
          r1 = max(r1,p/(p + 1))
          p += 1
          
     if 1/3 < n_ratio and n_ratio < 1/2: 
          r0 = 1/3
     return (r0, r1)

def generate_pairs(data):
     for i in range(len(data)):
          name_i, T1 = data[i]
          n1 = 360/T1
          for j in range(len(data)):
               name_j, T2 = data[j]
               n2 = 360/T2
               if n1 < n2:
                    yield name_i,name_j,n1,n2
                    
def create_mean_motion_ratios(data):
     product = []
     for name_i,name_j,n1,n2 in generate_pairs(data):
          (r0, r1) = get_bounds(n1/n2)
          a = (n1/n2 - r0)/(r1 - r0)   # Murray & Dermott (1.19)
          b = 0 if a <= 0.5 else 1   # Murray & Dermott (1.20)
          c = 2*np.pi*(a - b)    # Murray & Dermott (1.21)
          product.append((name_i,name_j,r0, r1, c))
    
     return product


def parse_args():
     parser = ArgumentParser(description='Calculate probability of delta_n<0.1 for problem 1.4.')
     parser.add_argument('--data', default='./data', help=f'Path to data files')
     parser.add_argument('--figs', default='./figs', help=f'Path to plots')
     parser.add_argument('--tolerance',default=0.15,type=float)
     return parser.parse_args()

def create_data(data_path,
                exclude=[
                     'Epimetheus',
                     'Telesto',
                     'Calypso',
                     'Helene']):
     '''
     Read data file
     
     Parameters:
         data_path
         exclude
     '''
     with open(data_path) as data_file:
          product = []
          data_reader = reader(data_file)
          for row in data_reader:
               if row[0] not in exclude:
                    product.append([row[0],abs(float(row[1]))])     
     return product

def get_bar(mean_motion_ratios,tolerance=0.15):
     c_indices = np.argsort([abs(c) for _,_,_,_,c in mean_motion_ratios])
     labels = []
     cc = []
     for i in c_indices:
          name_i,name_j,r0, r1, c = mean_motion_ratios[i]
          if r0 > 0 and abs(c) < tolerance:
               labels.append(f'{name_i}-{name_j}')
               cc.append(abs(c))
               print (mean_motion_ratios[i])
     return labels,cc,['xkcd:red','xkcd:green','xkcd:blue','xkcd:yellow']

def main():
     args = parse_args()

     mean_motion_ratios = create_mean_motion_ratios(
               create_data(
                    (Path(args.data)/Path(__file__).stem).with_suffix('.csv')))

     labels,cc,bar_colours = get_bar(mean_motion_ratios,tolerance=args.tolerance)
 
     
     cs = sorted([abs(c) for _,_,_,_,c in mean_motion_ratios])
     fig = figure(figsize=(12, 12))
     ax1 = fig.add_subplot(2, 1, 1)
     ax1.hist(cs,bins=100)
     ax1.set_title('Distribution of c')
     ax1.set_xlabel('c')
     ax2 = fig.add_subplot(2, 1, 2)

     ax2.bar(labels,cc,color=bar_colours)
     
     fig.tight_layout(h_pad=2)
     fig.savefig(Path(args.figs)/Path(__file__).stem)

     show()

if __name__ == '__main__':
     main()
