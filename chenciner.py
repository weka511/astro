#!/usr/bin/env python

#   Copyright (C) 2026 Simon Crase

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''
Chenciner choreography using symplectic integrator
'''

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from rki import ImplicitRungeKutta4,Driver
from threebody import Hamiltonian

__version__ = '1.0'
__author__ = 'Simon Crase'


      
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('file_name',help='File name for initial conditions')
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    parser.add_argument('-N','--Iterations',default=100,type=int,help='Number of steps')
    return parser.parse_args()

def read_data(file_name,dim=2):
    '''
    Read initial conditions from a file
    '''
    state = -1
    i = 0
    with open(file_name) as f:
        data_reader = reader(f)
        for row in data_reader:
            values = np.array([float(m) for m in row])
            match(state):
                case -1:
                    masses = values
                    Q = np.zeros((len(masses),dim))
                    P = np.zeros((len(masses),dim))
                    state = 0
                case 0:
                    Q[i,:] = values
                    i += 1
                    if i == len(masses):
                        i = 0
                        state = 1
                case 1:
                    P[i,:] = values
                    i += 1
 
        return Q,P,masses
    


def main():
    '''
    Read initial values and integrate equations of motion
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    R,R_dot,m = read_data(Path(args.data)/args.file_name)
    hamiltonian = Hamiltonian(m) 
    y = hamiltonian.create_initial_values(R,R_dot,m)
    integrator = ImplicitRungeKutta4(lambda y: hamiltonian.dH(y), 10, 1e-9)
    driver = Driver(integrator,h=0.1)
    XY = np.zeros((args.Iterations+1,6))

    XY[0,:] = hamiltonian.get_coordinates(y)
    for i in range(args.Iterations):
        y = driver.step(y)
        XY[i,:] = hamiltonian.get_coordinates(y)
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    ax1.plot(XY[0,:],XY[1,:],'r')
    ax1.scatter(R[0,0],R[0,1],c='r')
    ax1.scatter(XY[0,0],XY[0,1],c='r',marker='X')
    ax1.scatter(R[1,0],R[1,1],c='g')
    ax1.scatter(XY[0,2],XY[0,3],c='g',marker='X')
    ax1.scatter(R[2,0],R[2,1],c='b')
    ax1.scatter(XY[0,4],XY[0,5],c='b',marker='X')
    #ax2 = fig.add_subplot(2,2,2,adjustable='box',aspect=1.0)
    #ax2.plot(XY[2,:],XY[3,:],'g')
    #ax2.scatter(R[1,0],R[1,1],c='g')
    #ax3 = fig.add_subplot(2,2,3,adjustable='box',aspect=1.0)
    #ax3.plot(XY[4,:],XY[5,:],'b')
    #ax3.scatter(R[2,0],R[2,1],c='b')
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__=='__main__':
    main()
