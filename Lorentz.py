#!/usr/bin/env python

# Copyright (C) 2015-2026 Greenweaves Software Limited

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

# This tests the ImplicitRungeKutts Inntegrator by calculating the
# evolution of the Lorentz Attractor

from argparse import ArgumentParser
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import rki

__version__ = '1.0'
__author__ = 'Simon Crase'

class Lorentz:
    def __init__(self,sigma,rho,beta):
        self.sigma = sigma
        self.rho   = rho
        self.beta  = beta
        
    def dx(self,x):
        return np.array([
            self.sigma*(x[1]-x[0]),      \
            x[0]*(self.rho-x[2])-x[1],   \
            x[0]*x[1]-self.beta*x[2]
        ])
    
    def __str__(self):
        return 'sigma={0},rho={1}.beta={2}'.format(self.sigma,self.rho,self.beta)
 
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('-n', type=int,default=10000)
    parser.add_argument('--step',type=float,default=0.01)
    return parser.parse_args()
      
def main():
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
     
    lorentz = Lorentz(10,28,8.0/3.0)
    rk = rki.ImplicitRungeKutta4(lambda x: lorentz.dx(x),10,0.0000001)
    XYZ = np.zeros((args.n+1,3))
    XYZ[0,:] = np.array([1,0,1])
    for i in range(args.n):
        XYZ[i+1,:]  = rk.step(args.step,XYZ[i,:])
        
    fig = figure(figsize=(12,12))   
    ax  = fig.add_subplot(111, projection='3d')    
    ax.scatter(XYZ[:,0], XYZ[:,1], XYZ[:,2],s=1,edgecolors ='face')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(lorentz)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()  
    
if __name__=='__main__':
    main()