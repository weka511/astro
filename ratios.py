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
from logging import basicConfig,getLogger,INFO,FileHandler,StreamHandler,Formatter,DEBUG
from pathlib import Path
from sys import float_info
from time import strftime,time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

__version__ = '1.0'
__author__ = 'Simon Crase'

class SimpleRatios:
    '''
    This class keeps track of the allowable commensurabilities.
    
    Attributes:
        imax     Ratios are of from i1/i2 where i1<i2, i2 <= imax
        pairs    Allowable ratios expressed as pairs
        ratios   Allowable ratios as floating poinbt values
    '''
    def __init__(self,imax:int=7):
        self.imax = imax
        self.pairs = sorted([(i1,i2) for i2 in range(1,imax+1) for i1 in range(1,i2)],key=lambda x:x[0]/x[1])
        self.ratios = [i1/i2 for i1,i2 in self.pairs]
        
    def get_eps_max(self):
        '''
        Half the separation of the two closest ratios
        '''
        return 0.5/(self.imax*(self.imax-1))
        
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
        name1    Name of first planet or satellite
        name2    Name of second planet or satellite
        T1       Period of first planet or satellite
        T2       Period of second planet or satellite
        i1       First member of pair representing T1/T2
        i2       Second member of pair representing T1/T2
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
            name1    Name of first planet or satellite
            name2    Name of second planet or satellite
            T1       Period of first planet or satellite
            T2       Period of second planet or satellite
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
            simple_ratios   Object keeping track of commensurabilities   
        '''
        best_pair,distance = simple_ratios.get_match(self.T1/self.T2)
        if distance < simple_ratios.get_eps_max():
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
    parser.add_argument('--planets',default='planets',help='File name for planetary data')
    parser.add_argument('--R_min',default=100, type=float,help='We require data for satellites to have R greater than this value')
    parser.add_argument('--e_max',default=0.15, type=float,help='We require data for satellites to have eccentricity less than this value')
    parser.add_argument('--imax',default=7,type=int,help='Limit on denominator for acceptable ratios')
    parser.add_argument('--logs',default='./logs', help=f'Path to log files')
    return parser.parse_args()

def get_periods(path:Path,key:str='Planet',e_max:float=1.0,R_min:float=0,logger=None):
    '''
    Retrieve periods for acceptable planets and satellites. If the column for the period has no data,
    use Kepler's Third Law to calculate periods. (I'd like to remain consistent with Murray and Dermott, 
    so I am  not using data from third parties except as a check).
    
    Parameters:
        path      Path to planetary data
        key       Indicates whther we are dealing with planets or satellites
        e_max     We reject data unless eccentricity less than this value
        R_min     We require data to have R greater than this value
        
    Returns:
       Names of planets or satellites, and their orbital periods, 
       either read of calculated
     
    I have compared my calcukated values of T with values 
    from https://en.wikipedia.org/wiki/Orbital_period  
    
    Mercury   0.24084236579905632     0.240846
    Venus     0.6151859910713233      0.615
    Earth     1.0                     1
    Mars      1.8807584230745054      1.881
    Jupiter   11.869327586344829     11.86
    Saturn    29.452516340193995     29.46
    Uranus    84.0727581815471       84.01
    Neptune  164.88365834412897     164.8
    Pluto    248.08098430718422     248.1
    '''
    
    df = pd.read_csv(path)
    df_acceptable = df[(df['e'] < e_max) & (df['R'] > R_min)]
    names = df_acceptable[key].to_list()
    try:
        return names,df_acceptable['T'].to_numpy()
    except KeyError:
        periods = df_acceptable['a'].to_numpy()**(3/2)
        m = len(names)
        
        try:
            for i in range(len(names)):
                logger.debug(f'{names[i]},{periods[i]/periods[2]}')
        except AttributeError:
            pass
        return names,periods

def create_xkcd_colours():
    '''
    48 colours from https://blog.xkcd.com/2010/05/03/color-survey-results/
    '''
    return [
        'xkcd:purple','xkcd:green','xkcd:blue','xkcd:pink','xkcd:brown','xkcd:red',
        'xkcd:light blue','xkcd:teal','xkcd:orange','xkcd:light green','xkcd:magenta','xkcd:yellow',
        'xkcd:sky blue','xkcd:grey','xkcd:lime green','xkcd:light purple','xkcd:violet','xkcd:dark green',
        'xkcd:turquoise','xkcd:lavender','xkcd:dark blue','xkcd:tan','xkcd:cyan','xkcd:aqua',
        'xkcd:forest green','xkcd:mauve','xkcd:dark purple','xkcd:bright green','xkcd:maroon','xkcd:olive',
        'xkcd:salmon','xkcd:beige','xkcd:royal blue','xkcd:navy','xkcd:lilac','xkcd:black',
        'xkcd:hot pink','xkcd:light brown','xkcd:pale green','xkcd:peach','xkcd:olive green','xkcd:dark pink',
        'xkcd:periwinkle','xkcd:sea green','xkcd:lime','xkcd:indigo','xkcd:mustard','xkcd:light pink'
    ]  

def create_legend(simple_ratios,unique_ratios,colours,ax,loc='lower left',max_rows=4):
    '''
    Used to create a second legend showing the unique ratios
    
    Parameters:
        simple_ratios
        unique_ratios
        colours
        ax
        loc
    ''' 
    return ax.legend(handles = [Patch(color=colours[simple_ratios.get_sequence(ratio[0],ratio[1])],
                                      label=f'{ratio[1]}:{ratio[0]}') 
                                for ratio in unique_ratios], 
                     ncols=max(1,len(unique_ratios)//max_rows),
                     loc=loc,
                     title=f'The {len(unique_ratios)} unique ratios')

def create_logger(path_name:str):
    '''
    Set up console logger and file logger
    
    Parameters:
        path_name    Path name for log file
    '''
    def add_handler(handler,logger,formatter=Formatter('%(message)s'),level:int=INFO):
        handler.setLevel(level)
        handler.setFormatter(formatter)
        logger.addHandler(handler)    
    product = getLogger(__name__)
    product.setLevel(DEBUG)
    add_handler(FileHandler(path_name),product,level=DEBUG)
    add_handler(StreamHandler(),product)   
    np.set_printoptions(linewidth=np.nan) # Prevent lines being split when we log numpy arrays
    return product

def main():
    '''
    Use the data in Appendix A to find the periods of all possible pairs of
    periods among the planates and the prograde salellites of Mars, Jupiter, Saturn, Uranus, 
    and Neptune, with mean radii > 100 km and orbital eccentricites < 0.1. Taking i_max =7, show
    that thirty pairs of objects have ratios of orbital periods within epsilon-max of a permitted
    commensurability.
    '''
    rcParams['text.usetex'] = True
    rcParams['figure.constrained_layout.use'] = True
    start  = time()
    args = parse_args()
    logger = create_logger(f'{Path(args.logs)/Path(__file__).stem}{strftime('%Y%m%d%H%M%S')}.log')
    simple_ratios = SimpleRatios()
    names,periods = get_periods((Path(args.data)/args.planets).with_suffix('.csv'),logger=logger) 
    resonances = Resonance.create(names,periods)
    
    for primary in args.names:
        names,periods = get_periods((Path(args.data)/primary).with_suffix('.csv'),
                                    key='Satellite',
                                    e_max=args.e_max,
                                    R_min=args.R_min)
        
        resonances += Resonance.create(names,periods)
    
    acceptable_resonances = [resonance for resonance in resonances if resonance.can_match(simple_ratios)]
    logger.info(f'There were {len(acceptable_resonances)} out of {len(resonances)} possible')
    unique_ratios = {(resonance.i1,resonance.i2) for resonance in acceptable_resonances}
    logger.info('Unique resonances')
    for ratio in unique_ratios:
        logger.info(f'{ratio[1]}:{ratio[0]}')
    colours = create_xkcd_colours()  

    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax = fig.add_subplot(1,1,1)
    logger.info('All resonances')
    for r in acceptable_resonances:
        ax.scatter(r.i1,r.i2,
                    label=f'{r.name1}-{r.name2} {r.i2}:{r.i1}',
                    c=colours[simple_ratios.get_sequence(r.i1,r.i2)])
        logger.info(f'{r.name1}-{r.name2} {r.i2}:{r.i1}')
    
    legend1 = create_legend(simple_ratios,unique_ratios,colours,ax)
    ax.legend(loc='lower right',ncols=max(1,len(acceptable_resonances)//8),
              title=f'The {len(acceptable_resonances)} acceptable resonances out of out of {len(resonances)} ')
    fig.gca().add_artist(legend1)
    ax.set_xlabel('$i_1$')
    ax.set_ylabel('$i_2$')
    ax.set_ylim((0,args.imax+1))
    ax.set_xlim((0,args.imax))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    names = ', '.join(name.title() for name in args.names)
    ax.set_title(f'Planets, plus satellites of {names}: $R>${args.R_min}, $e<${args.e_max}, $imax=${args.imax}.')
 
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()
    
if __name__=='__main__':
    main()
