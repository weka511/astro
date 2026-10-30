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
Exercise 1.2 from Murray and Dermott. Calculate major axes using Kepler's model (Platonic solids),
allowing for  eccentricity.
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
from platonic import Platonic,PlatonicSolidModel

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
    return data['Planet'].tolist()[:6],data['a'].tolist()[:6],data['e'].tolist()[:6]   

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

def get_aphelion(a,e):
    return a*(1+e)

def get_perihelion(a,e):
    return a*(1-e)

def get_major_axis(x,e,from_perihelion=True):
    return x/(1-e) if from_perihelion else x/(1+e)

def build_polyhedra(a,e,platonic_solids = PlatonicSolidModel.create_kepler()):        
    perihelion_kepler = np.zeros((PlatonicSolidModel.N))
    aphelion_kepler = np.zeros((PlatonicSolidModel.N))    
 
    ratios = Platonic.create_ratios()
    aphelion_kepler[PlatonicSolidModel.EARTH] = get_aphelion(a[PlatonicSolidModel.EARTH],e[PlatonicSolidModel.EARTH])
    for i in range(PlatonicSolidModel.MARS,PlatonicSolidModel.SATURN+1):
        perihelion_kepler[i] = aphelion_kepler[i-1]*ratios[platonic_solids[i-1]]
        major_axis = get_major_axis(perihelion_kepler[i],e[i])
        aphelion_kepler[i] = get_aphelion(major_axis,e[i])
        
    perihelion_kepler[PlatonicSolidModel.EARTH] = get_perihelion(a[PlatonicSolidModel.EARTH],e[PlatonicSolidModel.EARTH])
    for i in range(PlatonicSolidModel.VENUS,PlatonicSolidModel.MERCURY-1,-1): 
        aphelion_kepler[i] = perihelion_kepler[i+1]/ratios[platonic_solids[i] ]
        major_axis = get_major_axis(aphelion_kepler[i],e[i],from_perihelion=False)
        perihelion_kepler[i] = get_perihelion(major_axis,e[i])
        
    return aphelion_kepler,perihelion_kepler

def get_rms_error(perihelion_data,perihelion_kepler,aphelion_data,aphelion_kepler):
    error = np.concatenate([
        perihelion_data[:PlatonicSolidModel.EARTH]-perihelion_kepler[:PlatonicSolidModel.EARTH],
        aphelion_data[PlatonicSolidModel.MARS:]-aphelion_kepler[PlatonicSolidModel.MARS:],
                            ])
    return np.sqrt(np.mean(error**2))

def search(a,e,perihelion_data,aphelion_data):
    rms_s = []
    platonic_solids = PlatonicSolidModel.create_kepler()
    for solids in multiset_permutations(platonic_solids):
        aphelion_kepler_x,perihelion_kepler_x = build_polyhedra(a,e,platonic_solids=solids) 
        rms_s.append(get_rms_error(perihelion_data,perihelion_kepler_x,aphelion_data,aphelion_kepler_x))
    return rms_s

def main():
    '''
    Do whatever...
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    Logger.create(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    names,a,e = get_data((Path(args.data)/args.source).with_suffix('.csv'))
    perihelion_data = np.zeros((PlatonicSolidModel.N))
    aphelion_data = np.zeros((PlatonicSolidModel.N))

    for i in range(len(names)):
        perihelion_data[i] = get_perihelion(a[i],e[i])
        aphelion_data[i] = get_aphelion(a[i],e[i])
        if i < 3:
            Logger.instance.info(f'{names[i]} Aphelion: {aphelion_data[i]}')
        else:
            Logger.instance.info(f'{names[i]} Perihelion: {perihelion_data[i]}') 
    
    aphelion_kepler,perihelion_kepler = build_polyhedra(a,e)
    rms_s = search(a,e,perihelion_data,aphelion_data)
    rms_k = get_rms_error(perihelion_data,perihelion_kepler,aphelion_data,aphelion_kepler)
    fig = figure(figsize=(12,12))
    
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(2,2,1)
    ax1.scatter(aphelion_data[PlatonicSolidModel.MARS:],aphelion_kepler[PlatonicSolidModel.MARS:],
                label='Aphelion (outer planets)')
    ax1.plot(aphelion_data[PlatonicSolidModel.MARS:],aphelion_data[PlatonicSolidModel.MARS:],
                label='Aphelion (outer planets)')
    ax1.scatter(perihelion_data[:PlatonicSolidModel.EARTH],perihelion_kepler[:PlatonicSolidModel.EARTH],
                label='Perihelion (inner planets)')
    ax1.plot(perihelion_data[:PlatonicSolidModel.EARTH],perihelion_data[:PlatonicSolidModel.EARTH],
                label='Perihelion (inner planets)')
    ax1.legend()
    ax1.set_xlabel('Perihelion/Aphelion')
    ax1.set_ylabel('Perihelion/Aphelion')
    ax1.set_title('Kepler Platonic Solids Model, with eccentricities.'
                  f' RMS={rms_k:.3f}')
    
    ax2 = fig.add_subplot(2,2,2)
    ax2.hist(rms_s,color='xkcd:blue',label='RMS error',density=True)
    ax2.axvline(x=rms_k, color='xkcd:red', linestyle='--', linewidth=2, label='RMS Kepler')
    ax2.legend()
    ax2.set_title('Analysis of RMS error')
    ax2.set_xlabel('RMS')
    ax2.set_ylabel('Frequency')
    
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
