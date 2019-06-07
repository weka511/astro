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
     
def get_average_interval(earth,mars):
     fig = plt.figure()
     ax = fig.add_subplot(111, projection='3d')     
     Xs,Ys,Zs = orbit (earth)
     ax.scatter(Xs, Ys, Zs, c='b', edgecolor='face',s=1,label='Earth')
     Xm,Ym,Zm = orbit (mars,lambda_dot=217103.78,Nr=53)
     ax.scatter(Xm,Ym,Zm,c='r',edgecolor='face',s=1,label='Mars') 
     ax.set_xlabel('X')
     ax.set_ylabel('Y')
     ax.set_zlabel('Z') 
     ax.legend()
     
     plt.figure()
     ds = [sqrt((Xs[i]-Xm[i])**2 + (Ys[i]-Ym[i])**2 + (Zs[i]-Zm[i])**2) for i in range(len(Xs))]
     plt.plot(ds,'g',label='Distance')
  
     mins = [(i,ds[i]) for i in range(1,len(ds)-1) if ds[i-1] > ds[i] and ds[i]<ds[i+1]]

     plt.scatter([i for i,_ in mins],[d for _,d in mins],label='Distance at Opposition')
     plt.legend()
     print ('Ratio={0:2f}'.format(max([d for _,d in mins])/min([d for _,d in mins])))

     
if __name__=='__main__':
     from utilities import get_data_file_name,get_planetary_data
    
     data = get_planetary_data(get_data_file_name())    
     get_average_interval(data['Earth'],data['Mars'])   
     plt.show()       