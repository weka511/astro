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


from pathlib import Path
from unittest import main,TestCase
import numpy as np
from scipy.special import binom
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
    


def get_periods(path:Path,key:str='Planet',e_max:float=1.0,R_min:float=-1):
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
     
    I have compared my calculated values of T with values 
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
    df_acceptable = df[(df['e'] < e_max) & (df['R'] > R_min)] if R_min>0 else df[(df['e'] < e_max)]
    names = df_acceptable[key].to_list()
    try:
        return names,df_acceptable['T'].to_numpy()
    except KeyError:
        periods = df_acceptable['a'].to_numpy()**(3/2)
        m = len(names)
        
        #try:
            #for i in range(len(names)):
                #Logger.instance.debug(f'{names[i]},{periods[i]/periods[2]}')
        #except AttributeError:
            #pass
        return names,periods







class Ratios:
    def __init__(self):
        self.ratios = [(1,2)]
        self.Nrs = []
        self.i_maxen = []
        self.Nr_naive = []
        self.primes = []
        self.eps_max = []
        
    def build(self,Max_imax):
        primes = create_primes(Max_imax)
        for imax in range(3,Max_imax+1):
            self.i_maxen.append(imax)
            self.Nr_naive.append(imax*(imax-1)//2)
            self._add_ratios(imax,primes)
            self.eps_max.append(0.5/(imax*(imax-1)))
            
    def get_p(self):
        '''
        Probability that a given ratio is commensurable
        '''
        Nr = max(self.Nrs)
        imax = max(self.i_maxen)
        eps_max = min(self.eps_max)
        return (2*eps_max*Nr) / ((imax-2)/imax+ eps_max)
    
    def get_P(self,Np,Nobs,p):
        '''
        Probability that ensemble of ratios is commensurable.
        This is just the binomial distribution.
        '''        
        return binom(Np,Nobs) * p**Np * (1-p)**(Np-Nobs)
                 
    def _add_ratios(self,imax,primes):
        for i in range(1,imax):
            self.ratios.append((i,imax))
        factors = factorize(imax,primes)
        if len(factors) > 1:
            self.ratios = sorted(self._purge_duplicates(self.ratios),
                                 key=lambda ratio:ratio[1]*imax+ratio[0])
        self.Nrs.append((len(self.ratios)))
            
    def _purge_duplicates(self,ratios):
        result = set()
        for a,b in ratios:
            factor = 2
            while factor <= a:
                if a%factor == 0 and b%factor == 0:
                    a //= factor
                    b //= factor
                factor += 1
            result.add((a,b))
                       
        return list(result)    

def create_primes(N):
    '''
    Create a list of prime numbers using the sieve of Eratosthenes
    
    Parameters:
        N          Largest condidate to be considered
        
    Returns:
        List of primes up to (and, if relevant), N
    '''
    canditates = list(range(2,N+1))
    i = 0
    while i < len(canditates):
        sieved = canditates[0:i+1]
        for j in range(i+1,len(canditates)):
            if canditates[j] % canditates[i] != 0:
                sieved.append(canditates[j])
        canditates = sieved
        i += 1
    return canditates
            
def factorize(n,primes):
    '''
    Factorize a number into a list of primes
    
    Parameters:
        n       Number to be factorized
        primes  List of primes to be considered
        
    Returns:
      List of factors (repeated as many times as necessary)
    '''
    factors = []
    for p in primes:
        while n%p == 0:
            factors.append(p)
            n //= p
    return factors   

def get_bounds(n1,n2):
    '''
    Establish bounds for a ratio (M&D Section 1.7)
    
    Parameters:
        n1     Mean motion of one orbit
        n2     Mean motion of second orbit: n1 < n2
        
    Returns:
        Closest set of inegers, p,p_prime, such that 
        p_prime/(p_prime+1) < n1/n2 < p(p+1)
    '''
    
    assert n1 < n2
    r = n1/n2
    #if 1/3 < r < 1/2:
        #return 3,2
    p = int(np.ceil(r/(1-r)))
    p_prime = p-1
    assert p_prime/(p_prime+1) < r < p/(p+1)
    return p,p_prime
        
def get_abc(n1,n2,p,p_prime):
    '''
    Calculate metrics from Murray & Dermott Section 1.7
    
    Parameters:
        n1
        n2
        p
        p_prime
    '''
    r_prime = 1/3 if 1/3 < n1/n2 and n1/n2 < 1/2 else p_prime/(p_prime+1)
    a = (n1/n2 - r_prime) / ( p/(p+1) - r_prime)    # M & D (1.19)     
    b = 0 if a <= 0.5 else 1                 # M & D (1.20)
    c = 2*np.pi*(a - b)                      # M & D (1.21)
    return a,b,c
    
class TestBounds(TestCase):
    def test32(self):
        p,p_prime = get_bounds(598, 626)
        self.assertEqual(22,p)
        self.assertEqual(21,p_prime)
        
    def test_atlas_pan(self):
        n1 = 360/0.6019
        n2 = 360/0.575
        p,p_prime = get_bounds(n1, n2)
        self.assertLess(p_prime/(p_prime+1), n1/n2)
        self.assertLess(n1/n2, p/(p+1))
        
        
    def test_enceladus_pan(self):
        n1 = 360/1.370218
        n2 = 360/0.575
        p,p_prime = get_bounds(n1, n2)
        self.assertLess(p_prime/(p_prime+1), n1/n2)
        self.assertLess(n1/n2, p/(p+1))   

        
if __name__ == '__main__':
    main()

