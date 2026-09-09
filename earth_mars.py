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
Murray and Dermott, Exercise 2.2

 1. Determine average times or orbital conjunction between earth & Mars
 2. Show that the minimum distance varies by factor of almost 2.
 3. Determine prbital motions over the period 1982-2002. Neglecting the 
    relative orbital inclinations, show that the closest opposition occurred
    in September 1988, and the furthest in February 1995, and determine
    the minimum distances at these times.
'''

from argparse import ArgumentParser
from pathlib import Path
from math import floor
import numpy as np
from matplotlib.pyplot import figure, show
from mpl_toolkits.mplot3d import Axes3D
from utilities import get_date, get_planetary_data
from orbital import get_mean_longitude,Calendar,create_orbit,is_minimum,get_distance

_version__ = '1.0'
__author__ = 'Simon Crase'

class Conjunctions:
     '''
     Find conjunctions in the orbits of two planets
     '''     
     def __init__(self,earth,mars,
                  From=Calendar.get_julian_date(1985,1,1),
                  To=Calendar.get_julian_date(2002,12,31),
                  Incr=10,
                  is2D=False):
  
          self.Xs,self.Ys,self.Zs,self.ts = create_orbit(earth,
                                                         lambda_dot=1293740.63,Nr=99,From=From,To=To,Incr=Incr,is2D=is2D)
          self.Xm,self.Ym,self.Zm,_ = create_orbit(mars,
                                                   lambda_dot=217103.78,Nr=53,From=From,To=To,Incr=Incr,is2D=is2D)
          self.distances = [get_distance(self.Xs[i],self.Ys[i],self.Zs[i],self.Xm[i],self.Ym[i],self.Zm[i]) 
                            for i in range(len(self.Xs))]
          self.conjunctions = [(i,self.ts[i],self.distances[i]) for i in range(1,len(self.distances)-1)
                          if is_minimum(self.distances[i-1],self.distances[i],self.distances[i+1])]
          
     def get_min_max(self):
          distances = [d for _,_,d in self.conjunctions]
          i0 = np.argmin(distances)
          _,t0,_ = self.conjunctions[i0]
          i1 = np.argmax(distances)
          _,t1,_ = self.conjunctions[i1]      
          return i0,t0,distances[i0],i1,t1,distances[i1]
          

def parse_args():
     parser = ArgumentParser('Find Conjunctions between Earth and Mars')
     parser.add_argument('--from',dest='from_date',default='1985-1-1',help='First date in range')
     parser.add_argument('--to',dest='to_date',default='2002-12-31',help='Last date in range')
     parser.add_argument('--incr',type=int,default=1,help='step size')
     parser.add_argument('--2D',dest='is2D',action='store_true',help='Ignore inclindations')
     parser.add_argument('--figs', default='./figs', help=f'Path to plots')
     parser.add_argument('--data', default='./data', help=f'Path to data files')
     return parser.parse_args()
     
def main():
     args = parse_args()
     data = get_planetary_data((Path(args.data)/Path(__file__).stem).with_suffix('.csv'))
     f1,f2,f3 = get_date(args.from_date)
     t1,t2,t3 = get_date(args.to_date)
     fig = figure(figsize=(20, 20), dpi=80)
     conjunctions = Conjunctions(data['Earth'],data['Mars'],
                         From=Calendar.get_julian_date(f1,f2,f3),
                         To=Calendar.get_julian_date(t1,t2,t3),
                         Incr=args.incr,
                         is2D=args.is2D)
         
     i0,t0,distances0,i1,t1,distances1  = conjunctions.get_min_max()
     Y0,M0,D0 = Calendar.get_calendar_date(t0)
     Y1,M1,D1 = Calendar.get_calendar_date(t1)

     # Plot conjunctions
     ax5 = fig.add_subplot(111) 
     ax5.plot(conjunctions.ts,conjunctions.distances,'g',label='Distance between planets')
     ax5.scatter([t for _,t,_ in conjunctions.conjunctions],[d for _,_,d in conjunctions.conjunctions],c='m',label='Conjunction')
     ax5.scatter(t0,distances0,marker=11,s=144,c='b',label=f'Closest at {Y0}/{M0:02}/{int(D0):02}, {distances0:.4} AU')
     ax5.axvline(x=t0, color='m', linestyle='-.')
     ax5.axhline(y=distances0, color='g', linestyle='-.')
     ax5.scatter(t1,distances1,marker=10,s=144,c='r',label=f'Furthest at {Y1}/{M1:02}/{int(D1):02}, {distances1:.4} AU')
     ax5.axvline(x=t1, color='m', linestyle=':')
     ax5.axhline(y=distances1, color='g', linestyle=':')
     ax5.set_xticklabels([])
     ax5.set_ylabel('Distance(AU)')
     ax5.set_xlabel('t')
     ax5.legend(loc='upper center')
     
     fig.suptitle('Murray and Dermott, Exercise 2.2')
     fig.tight_layout(pad=2,h_pad=5,w_pad=3)
     fig.savefig((Path(args.figs)/Path(__file__).stem).with_suffix('.png'))
     
     show()
     
if __name__=='__main__':
     main()
     