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
Murray & Dermott, Exercise 1.3. Setellites of Saturn
'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from matplotlib.pyplot import figure, show
import numpy as np


def get_bounds(n_ratio, max_p=100):
     '''
     calculate p and p' from Murray and Dermott, Section 1.7.
     We use r0 for the lower bound (MD p'/(p'+1)) and r1 for the upper
     '''
     if 1/3 < n_ratio and n_ratio < 1/2: # See MD Section 1.7
          return (1/3, 1/2)
     r0 = 0
     r1 = None
     for p in range(max_p):
          r = p/(p + 1)
          if r < n_ratio:
               r0 = r
          if r > n_ratio:
               r1 = r
               return (r0, r1)


def get_mean_motion_ratios(data, max_p=100):
     rcs = []
     for j in range(len(data)):
          _, T1 = data[j]
          n1 = 360/T1
          for i in range(j):
               _, T2 = data[i]
               n2 = 360/T2
               (r0, r1) = get_bounds(n1/n2, max_p=max_p)
               a = (n1/n2 - r0)/(r1 - r0)   # Murray & Dermott (1.19)
               b = 0 if a <= 0.5 else 1   # Murray & Dermott (1.20)
               c = 2*np.pi*(a - b)    # Murray & Dermott (1.21)
               rcs.append(((r0, r1), c))
     return rcs


def parse_args():
     parser = ArgumentParser(description='Calculate probability of delta_n<0.1 for problem 1.4.')
     parser.add_argument('--data', default='./data', help=f'Path to data files')
     parser.add_argument('--figs', default='./figs', help=f'Path to plots')
     return parser.parse_args()

def create_data(data_path):
     '''
     Read data file
     
     Parameters:
         data_path
     '''
     with open(data_path) as data_file:
          product = []
          data_reader = reader(data_file)
          for row in data_reader:
               product.append([row[0],abs(float(row[1]))])     
     return product

def main():
     args = parse_args()

     rcs = get_mean_motion_ratios(
               create_data(
                    (Path(args.data)/'saturn').with_suffix('.csv')))
     cs = sorted([abs(c) for _, c in rcs])
     for i in range(len(cs)):
          print(i, cs[i])
     fig = figure(figsize=(6, 6))
     ax = fig.add_subplot(1, 1, 1)
     ax.hist(cs, bins=100,density=True)
     ax.set_title('Distribution of c')
     ax.set_xlabel('c')
     print(len(set([r for r, _ in rcs])))
     fig.savefig(Path(args.figs)/Path(__file__).stem)

     show()

if __name__ == '__main__':
     main()

