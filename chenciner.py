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
Calculate and plot Chenciner & Montgomery choreography using symplectic integrator
'''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter
from pathlib import Path
from time import time,strftime
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from rki import ImplicitRungeKutta4, ImplicitRungeKutta2
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
    parser.add_argument('-n','--n',default=12,type=float,help='Number of orbits')
    parser.add_argument('--freq',default=100,type=int,help='Print progress every freq steps')
    parser.add_argument('--step',type=float,default=0.001)
    parser.add_argument('--logs', default='./logs', help=f'Path to log files')
    parser.add_argument('--integrator',default='ImplicitRungeKutta4')
    groupImplicitRungeKutta = parser.add_argument_group('ImplicitRungeKutta','Used with ImplicitRungeKutta')
    groupImplicitRungeKutta.add_argument('--atol',type=float,default=1e-16)
    groupImplicitRungeKutta.add_argument('--max_iterations',type=int,default=1000)
    return parser.parse_args()

def read_data(file_name,dim=2):
    '''
    Read initial conditions from a file
    
    Parameters:
        file_name
        dim
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

def plot_orbits(R,XY,ax=None):
    '''
    Plot orbits
    
    Parameters:
        R     Initial positions
        XY    Computed positions
        ax    Axis for plotting
    '''
    ax.scatter(XY[:,0],XY[:,1],c='r',marker='*',s=1)
    ax.scatter(XY[:,2],XY[:,3],c='g',marker='*',s=1)
    ax.scatter(XY[:,4],XY[:,5],c='b',marker='*',s=1)
    # Plot initial condition for all three masses
    ax.scatter(R[0,0],R[0,1],c='r',marker='+',s=200)
    ax.scatter(R[1,0],R[1,1],c='g',marker='+',s=200)
    ax.scatter(R[2,0],R[2,1],c='b',marker='+',s=200)
    
def create_logger(path_name):
    '''
    Set up console logger and file logger
    '''
    def add_handler(handler,logger,formatter=Formatter('%(message)s')):
        handler.setLevel(INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)    
    product = getLogger(__name__)
    product.setLevel(INFO)
    add_handler(FileHandler(path_name),product)
    add_handler(StreamHandler(),product)   
    np.set_printoptions(linewidth=np.nan) # Prevent lines being split when we log numpy arrays
    return product

def create_integrator(name,hamiltonian):
    '''
    Factory method for setting up integrator
    
    Parameters:
        name           Name of integrator to be used
        himiltonian    The aomiltonian that evolves system
    '''
    match name:
        case 'ImplicitRungeKutta4':
            return ImplicitRungeKutta4(lambda y: hamiltonian.dH(y), 
                                       max_iterations=1000, 
                                       atol=1e-16)  
        case 'ImplicitRungeKutta2':
            return ImplicitRungeKutta2(lambda y: hamiltonian.dH(y), 
                                       max_iterations=1000, 
                                       atol=1e-16)            
def main():
    '''
    Read initial values, integrate equations of motion and plot results
    '''
    args = parse_args()
    logger = create_logger(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    rcParams['text.usetex'] = True
    start  = time()
    
    R,R_dot,m = read_data(Path(args.data)/args.file_name)
    hamiltonian = Hamiltonian(m) 
    y = hamiltonian.create_initial_values(R,R_dot,m)
    integrator = create_integrator(args.integrator,hamiltonian)

    T = 6.32591398
    N = int (args.n*T/args.step)  
    XY = np.zeros((N+1,6))

    XY[0,:] = hamiltonian.get_coordinates(y)

    for i in range(N):
        y = integrator.step(args.step,y)
        XY[i+1,:] = hamiltonian.get_coordinates(y,atol=1.0e-4)
        if i%args.freq == 0:
            logger.info ((f'Step {i}, E={hamiltonian.get_total_energy(y)}'))
            
    fig = figure(figsize=(8,8))
    fig.suptitle(f'{Path(__file__).stem}, Orbits={args.n:,}, Integrator={integrator}')
    plot_orbits(R,XY,ax = fig.add_subplot(1,1,1))
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
