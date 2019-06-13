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

# Murray and Dermott, Exercise 2.2
#
# 1. Determine average times or orbital conjunction between earth & Mars
# 2. Show that the minimum distance varies by factor of almost 2.
# 3. Determine prbital motions over the period 1982-2002. Neglecting the 
#    relative orbital inclinations, show that the closest opposition occurred
#    in September 1988, and the furthest in February 1995, and determine
#    the minimum distances at these times.

from orbital import get_xy,get_mean_longitude,compose,get_julian_date,get_calendar_date,create_orbit,is_minimum,get_distance
from math import pi,radians,sqrt,floor
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import matmul, mean, std
from utilities import get_date
    
# find_conjunctions
#
# Find conjunctions in the orbits of two planets

def find_conjunctions(earth,
                      mars,
                      From=get_julian_date(1985,1,1),
                      To=get_julian_date(2002,12,31),
                      Incr=10,
                      is2D=False):
     Xs,Ys,Zs,ts  = create_orbit (earth,lambda_dot=1293740.63,Nr=99,From=From,To=To,Incr=Incr,is2D=is2D)
     Xm,Ym,Zm,_   = create_orbit (mars,lambda_dot=217103.78,Nr=53,From=From,To=To,Incr=Incr,is2D=is2D)
     distances    = [get_distance(Xs[i],Ys[i],Zs[i],Xm[i],Ym[i],Zm[i]) for i in range(len(Xs))]
     conjunctions = [(i,ts[i],distances[i]) for i in range(1,len(distances)-1)
                     if is_minimum(distances[i-1],distances[i],distances[i+1])]
     times        = [t for _,t,_ in conjunctions]
     intervals =    [12*(times[i]-times[i-1])/365 for i in range(1,len(times))]
     print ('Average time between conjunctions = {0:.2f} months, sigma = {1:.2f} months.'.format(mean(intervals),std(intervals)))
     print ('Ratio(Largest/Smallest)={0:.2f}'.format(max([d for _,_,d in conjunctions])/min([d for _,_,d in conjunctions])))
     for _,t,d in conjunctions:
          Y,M,D = get_calendar_date(t)
          print ('{0:02d}-{1:02d}-{2}: {3:.4f} AU'.format(int(floor(D)),M,Y,d))
               
     fig = plt.figure(figsize=(20, 20), dpi=80)
     
     # 3D Plot
     
     ax1 = fig.add_subplot(231, projection='3d',aspect='equal')     
     ax1.scatter(Xs, Ys, Zs, c='b', edgecolor='face', s=1)
     ax1.scatter(Xm, Ym, Zm, c='r', edgecolor='face', s=1) 
     ax1.set_xlabel('X')
     ax1.set_ylabel('Y')
     ax1.set_zlabel('Z') 
     
     # 2D plot X & Y
     
     ax2 = fig.add_subplot(232,aspect='equal') 
     p21 = ax2.scatter(Xs, Ys, c='b', edgecolor='face',s=1,label='Earth')
     p22 = ax2.scatter(Xm, Ym,c='r',edgecolor='face',s=1,label='Mars') 
     ax2.set_xlabel('X')
     ax2.set_ylabel('Y')
     
     # 2D plot Y & Z
     
     ax3 = fig.add_subplot(233) 
     ax3.scatter(Ys, Zs, c='b', edgecolor='face',s=1,label='Earth')
     ax3.scatter(Ym,Zm,c='r',edgecolor='face',s=1,label='Mars')
     ax3.set_xlabel('Y')
     ax3.set_ylabel('Z')
     
     # 2D plots - X & Z
     
     ax4 = fig.add_subplot(234) 
     ax4.scatter(Xs, Zs, c='b', edgecolor='face',s=1,label='Earth')
     ax4.scatter(Xm,Zm,c='r',edgecolor='face',s=1,label='Mars')
     ax4.set_xlabel('X')
     ax4.set_ylabel('Z')
     
     # Plot conjunctions
     
     ax5 = fig.add_subplot(235) 
     p51 = ax5.plot(ts,distances,'g',label='Distance')
     p52 = ax5.scatter([t for _,t,_ in conjunctions],[d for _,_,d in conjunctions],c='m',label='Conjunction')
     ax5.set_xticklabels([])
     ax5.set_ylabel('Distance(AU)')
     ax5.set_xlabel('t')
     ax5.legend()
     
     fig.legend([p21,p22],['Earth','Mars'],'upper center')


             
if __name__=='__main__':
     from utilities import get_data_file_name,get_planetary_data
     from argparse import ArgumentParser
     parser = ArgumentParser('Find Conjunctions between Earth and Mars')
     parser.add_argument('--from',dest='from_date',default='1985-1-1',help='First date in range')
     parser.add_argument('--to',dest='to_date',default='2002-12-31',help='Last date in range')
     parser.add_argument('--incr',type=int,default=10,help='step size')
     parser.add_argument('--2D',dest='is2D',action='store_true',help='Ignore inclindations')
     args     = parser.parse_args()
     data     = get_planetary_data(get_data_file_name()) 
     f1,f2,f3 = get_date(args.from_date)
     t1,t2,t3 = get_date(args.to_date)
     find_conjunctions(data['Earth'],data['Mars'],
                          From = get_julian_date(f1,f2,f3),
                          To   = get_julian_date(t1,t2,t3),
                          Incr = args.incr,
                          is2D = args.is2D) 
     plt.savefig(get_data_file_name(path='images',ext='png'))
     
     plt.show()       