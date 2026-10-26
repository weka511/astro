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
Figure 3-10 Plot the radial distances of the five Lagrangian points 
from the centre of mass as a function of log mu2
'''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from lagrange import L1,L2,L3

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
    
def main():
    '''
    Calculate and plot the Lagrangian points
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')

    xs = np.linspace(-1,-12,100)
    Ls = np.full((6,len(xs)),float('nan'))
    for i in range(len(xs)):
        mu2 = 10**xs[i]#np.exp(xs[i])
        mu1 = 1 - mu2
        solver1 = L1(mu2)
        solver2 = L2(mu2)
        solver3 = L3(mu2)
        alpha = solver1.get_alpha()
        try:
            r2,_ = solver1.solve(alpha)
            Ls[0,i] = mu1 - r2
        except RuntimeError:
            pass
        
        try:
            r2,_ = solver2.solve(alpha)
            Ls[1,i] = r2 + mu1
        except RuntimeError as e:
            pass
        
        beta = (-(7/12)*solver3.mu2/solver3.mu1 
                + (7/12)*(solver3.mu2/solver3.mu1)**2 
                - (13223/20736)*(solver3.mu2/solver3.mu1)**3)              
        try:      
            r1,_ = solver3.solve(beta + 1,bound=0.01)
            Ls[2,i] = r1 + mu2
        except RuntimeError:
            pass   
        
        Ls[4,i] = alpha + (1/3)*alpha**2 -(1/9) *alpha**3 -(31/81)*alpha**4 + mu1   # (3.88)
        Ls[3,i] = np.sqrt(1+mu1**2-mu1)
        Ls[5,i] = mu1 -(alpha - (1/3)*alpha**2-(1/9)*alpha**3-(23/81)*alpha**4)    #  (3.83)
        
    fig = figure(figsize=(8,8))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1) 
    ax1.plot(xs,Ls[0,:],label='L1',c='xkcd:purple')
    ax1.plot(xs,Ls[1,:],label='L2',c='xkcd:green')
    ax1.plot(xs,Ls[2,:],label='L3',c='xkcd:blue')
    ax1.plot(xs,Ls[3,:],label='L4,L5',c='xkcd:pink')
    ax1.plot(xs,Ls[4,:],label='L2 (3.88)',c='xkcd:brown',ls=':')
    ax1.plot(xs,Ls[5,:],label='L1 (3.83)',c='xkcd:red',ls=':')
    ax1.xaxis.set_inverted(True)
    ax1.set_xticks(np.arange(-2,-12,-2))
    ax1.set_yticks(np.arange(0.9, 1.1, 0.05))
    ax1.set_xlim((-1,-11))
    ax1.set_ylim((0.9,1.1))
    ax1.set_xlabel(r'$\log{\mu_2}$')
    ax1.set_ylabel('Radial Distance')

    ax1.legend(title='Lagrangian Points',loc='upper right')
    
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
