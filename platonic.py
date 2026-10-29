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

'''Common code needed by Murray & Dermott, Exercises 1 and 2'''

import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'

class Platonic:
    '''
    This class represents the 5 Platonic Solids
    '''
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
    
    @staticmethod
    def to_str(i):
        match i:
            case Platonic.TETRAHEDRON:
                return 'Tetrahedron'
            case Platonic.CUBE:
                return 'Cube'
            case Platonic.OCTAHEDRON:
                return 'Octahedron'
            case Platonic.DODECAHEDRON:
                return 'Dodecahedron'
            case Platonic.ICOSAHEDRON:
                return 'Icosahedron'           
            

class PlatonicSolidModel:
    '''
    This class models the Solar system using Platonic solids,
    in the manner of Kepler.
    '''
    MERCURY = 0
    VENUS = 1
    EARTH = 2
    MARS = 3
    JUPITER = 4
    SATURN = 5
    N = 6
    
    ratios = Platonic.create_ratios()

    @staticmethod   
    def create_kepler():
        return  [
            Platonic.ICOSAHEDRON,
            Platonic.OCTAHEDRON,
            Platonic.DODECAHEDRON,
            Platonic.TETRAHEDRON,
            Platonic.CUBE
        ]
    
    def __init__(self,solids):
        self.solids = solids.copy()
            
    def generate_annotations(self,a):
        j = 0
        inner = True
        for i in range(PlatonicSolidModel.N):
            if i == PlatonicSolidModel.EARTH: 
                inner = False
            else:
                yield a[i], self.a_kepler[i], Platonic.to_str(self.solids[j]),inner
                j += 1
    
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
