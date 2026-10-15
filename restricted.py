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
from matplotlib.pyplot import figure, show, colorbar
from matplotlib import rcParams
import numpy as np
from rki import ImplicitRungeKutta4,ImplicitRungeKutta2

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
    parser.add_argument('--mu2', default=0.01, type=float,help='Smaller of the two main masses')
    parser.add_argument('-N','--NIterations',type=int,default=1000,help='Number of iterations for integrator')
    parser.add_argument('--step',type=float,default=0.01,help='Stepsize for integrator')
    parser.add_argument('--integrator',
                        choices=IntegratorFactory.get_choices(),
                        default=IntegratorFactory.get_default(),
                        help='Choice of integrator')
    parser.add_argument('-n','--n',type=int,default=1,help='Number of iterations for Implicit Euler')
    parser.add_argument('--dx',type=float,default=[0,0],nargs=2,help='Initial perturbation (position)')
    parser.add_argument('--dv',type=float,default=[0,0],nargs=2,help='Initial perturbation (velocity)')
    parser.add_argument('--atol',type=float,default=1.0e-12,help='Tolerance for solving equation in implicit Runge Kutta')
    parser.add_argument('-m','--max_iterations',type=int,default=250,
                        help='Number of iterations for solving equation in implicit Runge Kutta')
    parser.add_argument('--freq',type=int,default=1000,help='Used to record progress')
    parser.add_argument('--stride',type=int,default=1,help='Used to skip over some points while plotting')
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
    '''
    Used in equations of motion: dx/dt = F(x)
    '''
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
                 -(self.mu1/r1**3 + self.mu2/r2**3)*x[1])
        return dx

class Euler:
    '''
    Euler method of integration
    
    A simple work-horse integrator
    '''
    def __init__(self,dy,n=1):
        self.dy = dy
        self.n = n
        
    def __str__(self):
        if self.n == 1:
            return 'Euler'
        else:
            return f'Implicit Euler, n={self.n}'
        
    def step(self, h, y):
        y1 = y.copy()
        for i in range(self.n):
            y1 = y + h*self.dy(y1)
        return y1

class IntegratorFactory: 
    '''
    This class sets up an Integrator.
    '''
    @staticmethod
    def get_choices():
        '''
        This method sets the options for the command line
        '''
        return [
            'euler',
            'rki2',
            'rki4'
        ]
    
    @staticmethod
    def create(f,args):
        '''
        Create an inegrator
        '''
        match args.integrator:
            case 'euler':
                return Euler(f,args.n)
            case 'rki2':
                return ImplicitRungeKutta2(f,
                                           atol=args.atol,
                                           max_iterations=args.max_iterations)
            case 'rki4':
                return ImplicitRungeKutta4(f,
                                           atol=args.atol,
                                           max_iterations=args.max_iterations)
            
    @staticmethod
    def get_default():
        return 'rki4'
     
def main():
    '''
    Integrate equations (3.16) & (3.17) for a particle near a Lagrange Point
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    
    f = F(args.mu2)
    integrator = IntegratorFactory.create(f,args)
    
    X = np.zeros((args.NIterations+1,4))
    X[0,0] = f.mu1 - 1/2
    X[0,1] = np.sqrt(3)/2
    X[0,0:2] += args.dx
    X[0,2:] += args.dv    
    for i in range(args.NIterations):
        X[i+1,:] = integrator.step(args.step,X[i,:])
        if i > 0 and i%args.freq == 0:
            Logger.instance.info(f'Step {i:,} of {args.NIterations:,}')
   
    fig = figure(figsize=(12,4))
    fig.suptitle(rf'{Path(__file__).stem}: ({integrator}), $\mu_2=${args.mu2}, N={args.NIterations:,}, step={args.step}')
    
    ax1 = fig.add_subplot(1,2,1,adjustable='box',aspect='equal')
    colorbar(
        ax1.scatter(X[::args.stride,0],X[::args.stride,1],
                    s=1,
                    c=np.sqrt(np.sqrt(X[::args.stride,2]**2 + X[::args.stride,3]**2)),
                    cmap='viridis'),
        label='Speed',
        location='bottom',
        shrink=0.5)
    
    ax1.scatter(X[0,0],X[0,1],marker='x',c='xkcd:blue')
    ax1.set_title('Position')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    
    ax2 = fig.add_subplot(1,2,2,adjustable='box',aspect='equal')
    ax2.scatter(X[::args.stride,2],X[::args.stride,3],s=1,c='xkcd:blue')
    ax2.set_title('Velocity')
    ax2.set_xlabel('u')
    ax2.set_ylabel('v')    
  
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        Logger.instance.exception(e,stack_info=True)
