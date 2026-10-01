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

'''Access ansd validate Leighton and Murray data'''

from argparse import ArgumentParser
from csv import DictReader
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    parser.add_argument('--name', default='saturn', help=f'Path to data files')
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    return parser.parse_args()

def create_data(data_path,
                exclude=[
                    'Epimetheus',
                     'Telesto',
                     'Calypso',
                     'Helene'
                     ],
                field='T'):
    '''
    Read data from one of the tables from Murray and Dermott, Appendex A

    Parameters:
        data_path     Path to scv file containing data
        exclude       A list of satellites whose data are to be ignored
        field         The name of the field whose data is to be read
    '''
    with open(data_path) as data_file:
        product = []
        data_reader = DictReader(data_file)
        for row in data_reader:
            if row['Satellite'] not in exclude:
                product.append([row['Satellite'],abs(float(row[field]))])     
    return product

def create_ratios(path):

    names_and_periods = create_data(path)
    names_and_semimajor_axes = create_data(path,field='a')
    n = len(names_and_periods)
    names = []
    ratios = np.zeros((n))
    for i in range(n):
        name,T = names_and_periods[i]
        _,a = names_and_semimajor_axes[i]
        ratios[i] = T**2/a**3
        names.append(name)
    return names,ratios/ratios.mean()

def main():
    '''
    Do whatever...
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1,adjustable='box',aspect=1.0)
    names,ratios = create_ratios((Path(args.data)/args.name).with_suffix('.csv'))
    ax1.bar(names,ratios)
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
