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

import numpy as np, matplotlib.pyplot as plt, math as m 
from matplotlib.backends.backend_pdf import PdfPages

@np.vectorize
def jacobi(x,y,n=1,mu2=0.2,Cj=3.9):
    n2  = n*n
    mu1 = 1-mu2
    x2  = x*x
    y2  = y*y
    r1  = m.sqrt((x+mu2)*(x+mu2) + y2)
    r2  = m.sqrt((x-mu1)*(x-mu1) + y2)
    return  n2*(x2 + y2) + 2*(mu1/r1 + mu2/r2) - Cj

def bounds(Z):
    minZ = float('inf')
    maxZ = - minZ
    for zz in Z:
        for z in zz:
            if z < minZ:
                minZ = z
            if z > maxZ:
                maxZ = z
    return (m.floor(minZ),m.ceil(maxZ))     

def plot_jacobi(fig=1,n=1,mu2=0.2,Cj=3.9,limit=5,origin='lower',pdf_pages=None,close=True): 
    plt.figure(fig)
    xlist = np.linspace(-limit, limit+0.001, 100)   
    ylist = np.linspace(-limit, limit+0.001, 100)
    X,Y = np.meshgrid(xlist, ylist)   
    Z = jacobi(X,Y,n,mu2,Cj)
    (z0,z1) = bounds(Z)
    
    levels=list(range(z0,0,10))+list(range(0,z1+1,10))

    ticks = [z0,0,z1]
    c = plt.pcolormesh(X,Y,Z)
    cbar = plt.colorbar(c,orientation='vertical',ticks=ticks)
    cbar.ax.set_yticklabels(['min', '0', 'max']) 
    CS3=plt.contourf(X, Y, Z, levels, cmap=plt.cm.jet,origin=origin)
    CS2 = plt.contour(X,Y,Z,levels=[0],colors='w',origin=origin,hold='on',linewidths=(1,))
    cbar.add_lines(CS2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'Zero velocity surfaces for $n={0},\mu_2 = {1},C_j={2}$'.format(n,mu2,Cj))
    if pdf_pages==None:
        plt.savefig('Jacobi{0}.png'.format(fig))
    else:
        pdf_pages.savefig()
    if close:
        plt.close(fig)
    
if __name__=='__main__':
    with PdfPages('jacobi.pdf') as pdf_pages:
        fig=1
        for i in range(1,11):
            for j in range(25):
                plot_jacobi(fig=fig,n=1,mu2=i/20.0,Cj=3.0+j/10,limit=2.5,pdf_pages=pdf_pages)
                fig+=1
