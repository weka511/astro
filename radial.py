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
    Do whatever...
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')

    xs = np.linspace(-3,-10,10)
    mu2s = np.exp(xs)
    Ls = np.full((3,len(xs)),float('nan'))
    for i in range(len(xs)):
        solver1 = L1(mu2s[i])
        solver2 = L2(mu2s[i])
        solver3 = L3(mu2s[i])
        alpha = solver1.get_alpha()
        try:
            Ls[0,i],_ = solver1.solve(alpha)
        except RuntimeError:
            pass
        try:
            Ls[1,i],_ = solver2.solve(alpha)
        except RuntimeError:
            pass
        try:
            beta = (-(7/12)*solver3.mu2/solver3.mu1 
                    + (7/12)*(solver3.mu2/solver3.mu1)**2 
                    - (13223/20736)*(solver3.mu2/solver3.mu1)**3)            
            Ls[2,i],_ = solver3.solve(beta + 1,bound=0.01)
        except RuntimeError:
            pass                
        
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1) 
    ax1.plot(xs,Ls[0,:],label='L1')
    ax1.plot(xs,Ls[1,:],label='L2')
    ax1.plot(xs,Ls[2,:],label='L3')
    ax1.xaxis.set_inverted(True)
    ax1.legend()
    
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
