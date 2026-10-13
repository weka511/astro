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

'''Integrate equations (3.16) & (3.17) for a particle near a Lagrange Point '''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from rki import ImplicitRungeKutta4

__version__ = '1.0'
__author__ = 'Simon Crase'

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--logs',default='./logs', help=f'Path to log files')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('--mu2', default=0.02, type=float,help='Smaller of the two main masses')
    parser.add_argument('-N','--NIterations',type=int,default=1000)
    parser.add_argument('--step',type=float,default=0.01)
    return parser.parse_args()


class Logger: 
    '''
    Singleton Logger
    '''
    instance = None
    
    @staticmethod
    def create(path_name:str):
        '''
        Set up console logger and file logger
        
        Parameters:
            path_name    Path name for log file
        '''
        def add_handler(handler,logger,formatter=Formatter('%(message)s'),level:int=INFO):
            handler.setLevel(level)
            handler.setFormatter(formatter)
            logger.addHandler(handler)    
        Logger.instance = getLogger(__name__)
        Logger.instance.setLevel(DEBUG)
        add_handler(FileHandler(path_name),Logger.instance,level=DEBUG)
        add_handler(StreamHandler(),Logger.instance)   
        np.set_printoptions(linewidth=np.nan) # Prevent lines being split when we log numpy arrays
 
class F:
    def __init__(self,mu2):
        self.mu1 = 1 - mu2
        self.mu2 = mu2
        self.n = 1
        
    def __call__(self,x):
        dx = np.zeros_like(x)
        r1 = np.sqrt((x[0]+self.mu2)**2 + x[1]**2)         # (3.8)
        r2 = np.sqrt((x[0]-self.mu1)**2 + x[1]**2)         # (3.9)
        dx[0] = x[2]
        dx[1] = x[3]
 
        dx[2] = (2*self.n*x[3] + self.n**2*x[0] 
                 - (self.mu1*(x[0]+self.mu2)/r1**3 + self.mu2*(x[0]-self.mu1)/r2**3))   # (3.16)
        
        dx[3] = (- 2*self.n*x[2] + self.n**2 * x[1]                                  # (3.17)
                 -(self.mu1/r1**3 + self.mu2/r2**3))
        return dx

def main():
    '''
    Integrate equations (3.16) & (3.17) for a particle near a Lagrange Point
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    
    f = F(args.mu2)
    integrator = ImplicitRungeKutta4(f,atol=1e-8)
    X = np.zeros((args.NIterations+1,4))
    X[0,0] = f.mu1 - 1/2
    X[0,1] = np.sqrt(3)/2
    for i in range(args.NIterations):
        X[i+1,:] = integrator.step(args.step,X[i,:])
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    
    ax1 = fig.add_subplot(1,2,1,adjustable='box',aspect='equal')
    ax1.scatter(X[:,0],X[:,1],s=1,c='xkcd:blue')
    ax1.scatter(X[0,0],X[0,1],marker='x',c='xkcd:blue')
    
    ax2 = fig.add_subplot(1,2,2,adjustable='box',aspect='equal')
    ax2.scatter(X[:,2],X[:,3],s=1,c='xkcd:blue')
  
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__ == '__main__':
    main()
