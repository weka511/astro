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
This module is responsible for instantiating an integrator
'''

from rki import ImplicitRungeKutta,ImplicitRungeKutta4,ImplicitRungeKutta2,NaiveSolver
from rk import Rk4, Nystrom
from euler import Euler
from verlet import VelocityVerlet

__version__ = '1.0'
__author__ = 'Simon Crase'

class IntegratorFactory: 
    '''
    This class is primaraily responsible for instantiating an Integrator. It also advertises
    the available choices, so they can appear in help text, and establishes a default.
    '''
    @staticmethod
    def get_choices():
        '''
        This method sets the options for the command line
        
        Returns:
           List of integrators that are available
        '''
        return [
            'euler',
            'rki2',
            'rki4',
            'verlet',
            'rk4',
            'nystrom'
        ]
    
    @staticmethod
    def create(dy,integrator,atol,max_iterations):
        '''
        Instantiate an integrator
        
        Parameters:
           dy              Function f for ODE: dy/dt=f(y)
           integrator      A key used to specify the inegrator that will be created
           max_iterations  Maximum number of iterations for solving implicit equations
           atol            Maximum error we can tolerate when solving implicit equations
           
        Returns:
           Newly created integrator
        '''
  
        match integrator:
            case 'euler':
                return  Euler(dy,n=max_iterations,atol=atol)
            
            case 'rki2':
                product = ImplicitRungeKutta2(dy)
 
            case 'rki4':
                product = ImplicitRungeKutta4(dy) 
                
            case 'verlet':
                return VelocityVerlet(dy)
            
            case 'rk4':
                return Rk4(dy)
            
            case 'nystrom':
                return Nystrom(dy)            
             
            case _:
                raise ValueError(f'Integrator {integrator} not found.')
   
        product.attach(NaiveSolver(atol=atol, max_iterations=max_iterations))           
  
        return product        
            
    @staticmethod
    def get_default():
        '''
        Default if no integrator has been specified
        '''
        return 'rki4'