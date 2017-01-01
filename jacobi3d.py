# Copyright (C) 2016 Greenweaves Software Pty Ltd

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

from mpl_toolkits.mplot3d import Axes3D
import numpy as np, matplotlib.pyplot as plt, math as m, jacobi,matplotlib.colors as clrs

def plot_3d(limit=2,eps=0.001,minZ=-3.9,maxZ=-2.84,steps=1000):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')    
    xlist = np.linspace(-limit, limit+eps, steps)   
    ylist = np.linspace(-limit, limit+eps, steps)
    X,Y = np.meshgrid(xlist, ylist)   
    Z = -jacobi.jacobi(X,Y)
    Z[Z<minZ]= np.nan
    Z[Z>maxZ]= np.nan
    norm = clrs.Normalize(vmin = minZ, vmax = maxZ, clip = False)
    surf=ax.plot_surface(X, Y, Z,cmap=plt.cm.jet,norm=norm)
    ax.set_xlim(-limit,limit)
    ax.set_ylim(-limit,limit)
    fig.colorbar(surf, shrink=0.5, aspect=5)

if __name__=='__main__':    
    plot_3d(limit=1)
