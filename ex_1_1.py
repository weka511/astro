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

class Platonic:
    TETRAHEDRON = 4
    CUBE = 6
    OCTAHEDRON = 8
    DODECAHEDRON = 12
    ICOSAHEDRON = 20
    
    @staticmethod
    def create_ratios():
        '''
        Used to create ratios of radii of circumsphere to insphere
        '''
        product = {
            Platonic.DODECAHEDRON:np.sqrt(15-6*np.sqrt(5)),
            Platonic.CUBE: np.sqrt(3),
            Platonic.TETRAHEDRON: 3
        }
        
        product[Platonic.OCTAHEDRON] = product[Platonic.CUBE]
        product[Platonic.ICOSAHEDRON] = product[Platonic.DODECAHEDRON]
        return product

class PlatonicSolidModel:
    
    MERCURY = 0
    VENUS = 1
    EARTH = 2
    MARS = 3
    JUPITER = 4
    SATURN = 5
    N = 6
    
    ratios = Platonic.create_ratios()
    
    def __init__(self,solids):
        self.solids = solids.copy()
    
    def fit(self,a):
        '''
        Calculate Major axis for each planet assuming a Platonic solid model
        '''
        self.a_kepler = np.zeros((PlatonicSolidModel.N))
        self.a_kepler[PlatonicSolidModel.EARTH] = a[PlatonicSolidModel.EARTH]
        for i in [PlatonicSolidModel.VENUS,PlatonicSolidModel.MERCURY]:
            self.a_kepler[i] = self.a_kepler[i+1]/PlatonicSolidModel.ratios[self.solids[1-i]] 
   
        for i in [PlatonicSolidModel.MARS,PlatonicSolidModel.JUPITER,PlatonicSolidModel.SATURN]:
            self.a_kepler[i] = self.a_kepler[i-1]*PlatonicSolidModel.ratios[self.solids[i-1]]         
   
        self.rms = np.sqrt(np.mean((a-self.a_kepler)**2))
        return self.rms

def main():
    '''
    Calculate major exes using Kepler's model (Platonic solids).
    Fit best possible permuation of Platonic Solids, neglecting eccentricity.
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
 
    a = np.array(get_data((Path(args.data)/args.source).with_suffix('.csv')))
    platonic_solids = [
        Platonic.ICOSAHEDRON,
        Platonic.OCTAHEDRON,
        Platonic.DODECAHEDRON,
        Platonic.TETRAHEDRON,
        Platonic.CUBE
    ]
    kepler_model = PlatonicSolidModel(platonic_solids)
    kepler_rms = kepler_model.fit(a)

    models = [PlatonicSolidModel(solids) for solids in multiset_permutations(platonic_solids)]
    rmss = [m.fit(a) for m in models]
    i_min = np.argmin(rmss)
    i_max = np.argmax(rmss)
    
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    
    ax1 = fig.add_subplot(2,2,1)
    ax1.scatter(a,kepler_model.a_kepler,c='xkcd:purple',label='Predicted')
    ax1.plot(a,a,c='xkcd:green',label='Actual')
    ax1.legend()
    ax1.set_xlabel('Semi-major axis')
    ax1.set_ylabel('Semi-major axis (predicted)')
    ax1.set_title(f'Data source={args.source}, RMS error = {kepler_rms:.3f}')
    
    ax2 = fig.add_subplot(2,2,2)
    ax2.scatter(a,models[i_min].a_kepler,c='xkcd:purple',label='Predicted')
    ax2.plot(a,a,c='xkcd:green',label='Actual')
    ax2.legend()
    ax2.set_xlabel('Semi-major axis')
    ax2.set_ylabel('Semi-major axis (predicted)')
    ax2.set_title(f'Model={models[i_min].solids}, RMS error = {models[i_min].rms:.3f}')
    
    ax3 = fig.add_subplot(2,2,3)
    ax3.scatter(a,models[i_max].a_kepler,c='xkcd:purple',label='Predicted')
    ax3.plot(a,a,c='xkcd:green',label='Actual')
    ax3.legend()
    ax3.set_xlabel('Semi-major axis')
    ax3.set_ylabel('Semi-major axis (predicted)')
    ax3.set_title(f'Model={models[i_max].solids}, RMS error = {models[i_max].rms:.3f}')   
    
    ax4 = fig.add_subplot(2,2,4)
    ax4.hist(rmss,color='xkcd:blue',density=True)
    ax4.set_xlabel('RMS error')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Distribution of Error')
    ax4.axvline(x=kepler_rms, color='xkcd:red', linestyle='--', linewidth=2, label='Kepler')
    ax4.legend()
    
    fig.tight_layout(h_pad=5,pad=3)
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
