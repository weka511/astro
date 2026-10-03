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
Exercise 1.6 (e) Use the data in Appendix A to find the periods of all possible pairs of
periods among the planates and the prograde salellites of Mars, Jupiter, Saturn, Uranus, 
and Neptune, with mean radii > 100 km and orbital eccentricites < 0.15. Taking i_max =7, show
that thirty pairs of objects have ratios of orbital periods witin epsilon-max of a permitted
commensurability.
'''

from argparse import ArgumentParser
from pathlib import Path
from sys import float_info
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

__version__ = '1.0'
__author__ = 'Simon Crase'

class SimpleRatios:
    '''
    This class keeps track of the allowable commensurabilities.
    
    Attributes:
        imax
        pairs
        ratios
        eps_max
    '''
    def __init__(self,imax:int=7):
        self.imax = imax
        self.pairs = sorted([(i1,i2) for i2 in range(1,imax+1) for i1 in range(1,i2)],key=lambda x:x[0]/x[1])
        self.ratios = [i1/i2 for i1,i2 in self.pairs]
        self.eps_max = 0.5*min([a-b for a in self.ratios for b in self.ratios if a > b])
        
    def get_match(self,target:float):
        '''
        Find the pair whose ratio is the best match to one that has been supplied
        
        Parameters:
            target  
        '''
        best = np.argmin(abs(target - self.ratios))
        return self.pairs[best],abs(target-self.ratios[best])
    
    def get_sequence(self,i1:int,i2:int):
        '''
        A function used to assign a sequence number to each pair,
        e.g. for selecting colours
        
        Parameters:
            i1
            i2
        '''
        return self.pairs.index((i1,i2))
    
class Resonance:
    '''
    This class keeps track of the actual ratios between a pair of orbits
    
    Attributes:
        name1
        name2
        T1
        T2
        i1
        i2
    '''
    
    @staticmethod
    def create(names:[str],periods:[float]):
        '''
        Populate a list with all possible pairs of orbits 
        from a group of planets or satellites
        
        Parameters:
           names     List of names of planets or satellites
           periods   List of periods for orbits, one for each name
        '''
        return [Resonance(names[i],names[j],periods[i],periods[j]) 
                          for i in range(len(names)) 
                          for j in range(i+1,len(names))
                          if periods[i] > 0 and periods[j] > 0]
    
    def __init__(self,name1:str,name2:str,T1:float,T2:float):
        '''
        Create a Resonance
        
        Parameters:
            name1
            name2
            T1
            T2
        '''
        self.name1 = name1
        self.name2 = name2
        self.T1 = T1
        self.T2 = T2
        self.i1 = None
        self.i2 = None
        
    def can_match(self,simple_ratios):
        '''
        Used to find the simple integer ratio that best matches this resonance
        
        Parameters:
            simple_ratios   Object keepting track of commensurabilities   
        '''
        best_pair,distance = simple_ratios.get_match(self.T1/self.T2)
        if distance < simple_ratios.eps_max:
            self.i1,self.i2 = best_pair
            return True
        return False
    
def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('names',nargs='*',help='File name for satellite data')
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    parser.add_argument('--planets', default='planets',help='File name for planetary data')
    parser.add_argument('--R_min', default=100, type=float,help='We require data to have R greater than this value')
    parser.add_argument('--e_max', default=0.15, type=float,help='We reject data unless eccentricity less than this value')
    parser.add_argument('--imax',default=7,type=int,help='Limit on denominator for acceptable ratios')
    return parser.parse_args()

def get_periods(path:Path,key:str='Planet',e_max:float=0.15,R_min:float=0):
    '''
    Retrieve periods for acceptable planets and satellites. If the column for the period has no data,
    use Kepler's Third Law to calculate periods
    
    Parameters:
        path      Path to planetary data
        key       Indicates whther we are dealing with planets or satellites
        e_max     We reject data unless eccentricity less than this value
        R_min     We require data to have R greater than this value
        
    Returns:
       Names of planets or satellites, and their orbital periods, 
       either read of calculated
    '''
    df = pd.read_csv(path)
    df_acceptable = df[(df['e'] < e_max) & (df['R'] > R_min)]
    names = df_acceptable[key].to_list()
    try:
        return names,df_acceptable['T'].to_numpy()
    except KeyError:
        return names,df_acceptable['a'].to_numpy()**(3/2)

def main():
    '''
    Use the data in Appendix A to find the periods of all possible pairs of
    periods among the planates and the prograde salellites of Mars, Jupiter, Saturn, Uranus, 
    and Neptune, with mean radii > 100 km and orbital eccentricites < 0.1. Taking i_max =7, show
    that thirty pairs of objects have ratios of orbital periods witin epsilon-max of a permitted
    commensurability.
    '''
    rcParams['text.usetex'] = True
    rcParams['figure.constrained_layout.use'] = True
    start  = time()
    args = parse_args()
    simple_ratios = SimpleRatios()
    names,periods = get_periods((Path(args.data)/args.planets).with_suffix('.csv'),
                                e_max=args.e_max)
       
    resonances = Resonance.create(names,periods)
    
    for primary in args.names:
        names,periods = get_periods((Path(args.data)/primary).with_suffix('.csv'),
                                    key='Satellite',
                                    e_max=args.e_max,
                                    R_min=args.R_min)
        
        resonances += Resonance.create(names,periods)
    
    acceptable_resonances = [resonance for resonance in resonances if resonance.can_match(simple_ratios)]
    
    unique_ratios = {f'{resonance.i1}:{resonance.i2}' for resonance in acceptable_resonances}
    
    colours = [
        'xkcd:purple','xkcd:green','xkcd:blue','xkcd:pink','xkcd:brown','xkcd:red',
        'xkcd:light blue','xkcd:teal','xkcd:orange','xkcd:light green','xkcd:magenta','xkcd:yellow',
        'xkcd:sky blue','xkcd:grey','xkcd:lime green','xkcd:light purple','xkcd:violet','xkcd:dark green',
        'xkcd:turquoise','xkcd:lavender','xkcd:dark blue','xkcd:tan','xkcd:cyan','xkcd:aqua',
        'xkcd:forest green','xkcd:mauve','xkcd:dark purple','xkcd:bright green','xkcd:maroon','xkcd:olive'
    ]  
        
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(1,1,1)
    for r in acceptable_resonances:
        ax1.scatter(r.i1,r.i2,
                    label=f'{r.name1}-{r.name2} {r.i2}:{r.i1}',
                    c=colours[simple_ratios.get_sequence(r.i1,r.i2)])
    ax1.legend(ncols=max(1,len(acceptable_resonances)//12))
    ax1.set_xlabel('$i_1$')
    ax1.set_ylabel('$i_2$')
    ax1.set_ylim((0,args.imax+1))
    ax1.set_xlim((0,args.imax))
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_title(f'$R>${args.R_min}, $e<${args.e_max}, $imax=${args.imax}. ' 
                  f'There are {len(acceptable_resonances)} acceptable resonances with {len(unique_ratios)} unique ratios')
 
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
    
if __name__=='__main__':
    main()

#ax1.bar(observed_ratios,5,
        #width=0.1,
        #facecolor=colours,
        #edgecolor='xkcd:white'
#)
#ax2 = fig.add_subplot(2,1,2)
#ax2.bar(matched_ratios,5,
        #width=0.1,
        #facecolor=colours,
        #edgecolor='xkcd:white'
#)