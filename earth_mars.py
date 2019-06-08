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
# ...

from orbital import get_xy,get_lambda,compose
from math import pi,radians,sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import matmul

def orbit(planet,N=10000,lambda_dot=1293740.63,Nr=99):
     a,e,I,varpi,Omega,lambda0 = planet
     Xs = []
     Ys = []
     Zs = []
     
     for i in range(N):
          T = 0.0001 * i
          x,y = get_xy(T=T, 
                       lambdaT = get_lambda(T,lambda0=lambda0,lambda_dot=lambda_dot,Nr=Nr), 
                       eccentricity = e,a=a,varpi=radians(varpi))
          P = compose(omega=radians(varpi - Omega), I=radians(I), Omega=radians(Omega))
          W = matmul(P,[[x],[y],[0]])
          Xs.append(W[0])
          Ys.append(W[1])
          Zs.append(W[2])
     return (Xs,Ys,Zs)
 
 
def is_minimum(a,b,c):
     return a>b and b < c

 
def get_average_interval(earth,mars):
     Xs,Ys,Zs = orbit (earth)
     Xm,Ym,Zm = orbit (mars,lambda_dot=217103.78,Nr=53)
     distances = [sqrt((Xs[i]-Xm[i])**2 + (Ys[i]-Ym[i])**2 + (Zs[i]-Zm[i])**2) for i in range(len(Xs))]
     conjunctions = [(i,distances[i]) for i in range(1,len(distances)-1) if is_minimum(distances[i-1],distances[i],distances[i+1])]
     print ('Ratio={0:2f}'.format(max([d for _,d in conjunctions])/min([d for _,d in conjunctions])))
     
     fig = plt.figure(figsize=(20, 20), dpi=80)
     
     ax1 = fig.add_subplot(231, projection='3d',aspect='equal')     
     ax1.scatter(Xs, Ys, Zs, c='b', edgecolor='face',s=1,label='Earth')
     ax1.scatter(Xm,Ym,Zm,c='r',edgecolor='face',s=1,label='Mars') 
     ax1.set_xlabel('X')
     ax1.set_ylabel('Y')
     ax1.set_zlabel('Z') 
     
     ax2 = fig.add_subplot(232,aspect='equal') 
     p21 = ax2.scatter(Xs, Ys, c='b', edgecolor='face',s=1,label='Earth')
     p22 = ax2.scatter(Xm, Ym,c='r',edgecolor='face',s=1,label='Mars') 
     ax2.set_xlabel('X')
     ax2.set_ylabel('Y')
     
     ax3 = fig.add_subplot(233) 
     ax3.scatter(Ys, Zs, c='b', edgecolor='face',s=1,label='Earth')
     ax3.scatter(Ym,Zm,c='r',edgecolor='face',s=1,label='Mars')
     ax3.set_xlabel('Y')
     ax3.set_ylabel('Z')
     
     ax4 = fig.add_subplot(234) 
     ax4.scatter(Xs, Zs, c='b', edgecolor='face',s=1,label='Earth')
     ax4.scatter(Xm,Zm,c='r',edgecolor='face',s=1,label='Mars')
     ax4.set_xlabel('X')
     ax4.set_ylabel('Z')
     
     ax5 = fig.add_subplot(235) 
     p51 = ax5.plot(distances,'g',label='Distance')
     p52 = ax5.scatter([i for i,_ in conjunctions],[d for _,d in conjunctions],c='m',label='Conjunction')
     ax5.legend()
     
     fig.legend([p21,p22],['Earth','Mars'],'upper center')
     
if __name__=='__main__':
     from utilities import get_data_file_name,get_planetary_data
    
     data = get_planetary_data(get_data_file_name())    
     get_average_interval(data['Earth'],data['Mars']) 
     plt.savefig(get_data_file_name(path='images',ext='png'))
     
     plt.show()       