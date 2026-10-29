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
Exercise 1.1 from Murray and Dermott. Calculate major exes using Kepler's model (Platonic solids).
Fit best possible permuation of Platonic Solids, neglecting eccentricity.
'''

from argparse import ArgumentParser
from csv import reader
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
import pandas as pd
from sympy.utilities.iterables import multiset_permutations

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
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    parser.add_argument('source', choices=['kepler','planets'])
    return parser.parse_args()

def get_data(path):
    data = pd.read_csv(path)
    return data['a'].tolist()[:6]  

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
        
def create_platonic_ratios():
    product = {
        12:np.sqrt(15-6*np.sqrt(5)),
        6: np.sqrt(3),
        4: 3
    }
    product[8] = product[6]
    product[20] = product[12]
    return product

def fit(a,platonic_solids,platonic_ratios):
    '''
    Calculate Major axis for each planet assuming a Platonic sold model
    '''
    a_kepler = np.zeros((6))
    a_kepler[2] = a[2]                              # Earth
    a_kepler[1] = a_kepler[2] / platonic_ratios[platonic_solids[0]] # Earth/Venus icosahedron
    a_kepler[0] = a_kepler[1] / platonic_ratios[platonic_solids[1]]  # Venus/Mercury octahedron
    a_kepler[3] = a_kepler[2] * platonic_ratios[platonic_solids[2]]  # Earth/Mars dodecadedron
    a_kepler[4] = a_kepler[3] * platonic_ratios[platonic_solids[3]]  # Mars/Jupiter tetrahedron
    a_kepler[5] = a_kepler[4] * platonic_ratios[platonic_solids[4]]  # Jupiter/Saturn cube
    rms = np.sqrt(np.mean((a-a_kepler)**2))
    return a_kepler,rms

def main():
    '''
    Calculate major exes using Kepler's model (Platonic solids).
    Fit best possible permuation of Platonic Solids, neglecting eccentricity.
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
 
    platonic_ratios = create_platonic_ratios()
    a = np.array(get_data((Path(args.data)/args.source).with_suffix('.csv')))
    platonic_solids = [20,8,12,4,6]
    a_kepler,rms = fit(a,platonic_solids,platonic_ratios)
    rms_best = rms
    a_best = a_kepler.copy()
    model = platonic_solids.copy()
    for solids in multiset_permutations(platonic_solids):
        aa,rms0 = fit(a,solids,platonic_ratios)
        if rms0 < rms_best:
            rms_best = rms0
            a_best = aa
            model = solids.copy()

    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    
    ax1 = fig.add_subplot(2,2,1)
    ax1.scatter(a,a_kepler,c='xkcd:purple',label='Predicted')
    ax1.plot(a,a,c='xkcd:green',label='Actual')
    ax1.legend()
    ax1.set_xlabel('Semi-major axis')
    ax1.set_ylabel('Semi-major axis (predicted)')
    ax1.set_title(f'Data source={args.source}, RMS error = {rms:.3f}')
    
    ax2 = fig.add_subplot(2,2,2)
    ax2.scatter(a,a_best,c='xkcd:purple',label='Predicted')
    ax2.plot(a,a,c='xkcd:green',label='Actual')
    ax2.legend()
    ax2.set_xlabel('Semi-major axis')
    ax2.set_ylabel('Semi-major axis (predicted)')
    ax2.set_title(f'Model={model}, RMS error = {rms_best:.3f}')
    
    fig.tight_layout(h_pad=2,pad=3)
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
