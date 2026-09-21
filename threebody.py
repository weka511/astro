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
Chenciner choreography using symplectic integrator
'''

from unittest import TestCase,main
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'

class Hamiltonian:
    '''
    Attributes:
        M 
        mu 
        m
        g1 
        g2
    '''
    index_r = 0         # Index into y: [r,theta,rho,Theta,p,l,P,L]
    index_theta = 1     # Index into y:
    index_rho = 2       # Index into y:
    index_Theta = 3     # Index into y:
    index_p = 4         # Index into y:
    index_l = 5         # Index into y:
    index_P = 6         # Index into y:
    index_L = 7         # Index into y:
    
    def __init__(self,m,G=1):
        '''
        Parameters:
            m
            G
        '''
        self.G = G
        self.M = m.sum()        # Section 5
        self.mu = m[0] + m[1]   # Section 5
        self.m = m
        self.g1 = m[0]*m[1]/self.mu    # Reduced mass - Section 5 - just after (28) 
        self.g2 = m[2]*self.mu/self.M  # Reduced mass - Section 5 - just after (28) 
        
    def create_initial_values(self,R,R_dot,m):
        '''
        Determine the canonical coordinates for a configuration
        
        Parameters:
            R        Positions of the 3 points, shape =(3,2)
            R_dot    Velocities of the 3 points
            m        The three masses
        '''
        r_polar,theta = Geometry.get_polar_coordinates( R[1,:] - R[0,:] ) # Equation (27a)
        r_dot =  R_dot[1,:] - R_dot[0,:]
        p = self.g1 * Geometry.get_r_velocity(r_dot, theta)                            # Linear momentum
        l = self.g1 * r_polar**2 * Geometry.get_theta_dot(r_dot, theta, r_polar)       # angular momentum
        centre_of_mass01 = np.dot(m[:2],R[:2,:])/self.mu
        rho_dot = R_dot[2,:] - centre_of_mass01  # Parallel to formula for rho
        rho_polar,Theta = Geometry.get_polar_coordinates( R[2,:] - centre_of_mass01)
        P = self.g2 * Geometry.get_r_velocity(rho_dot, Theta)                           # Linear momentum
        L = self.g2 * rho_polar**2 * Geometry.get_theta_dot(rho_dot, Theta, rho_polar)  # angular momentum
        
        return np.array([r_polar,theta,rho_polar,Theta,p,l,P,L])    
    
    @staticmethod
    def create_vectors(y):
        r = y[Hamiltonian.index_r]
        theta = y[Hamiltonian.index_theta]
        rho = y[Hamiltonian.index_rho]
        Theta = y[Hamiltonian.index_Theta]
        return r*np.array([np.cos(theta),np.sin(theta)]),rho*np.array([np.cos(Theta),np.sin(Theta)])
    
    def get_coordinates(self,y,atol=1e-12):
        '''
        Convert canonical coordinates back to a position in the original space
        
        Parameters:
            y         A configuration in canonical coordinates
        '''
        r,rho = Hamiltonian.create_vectors(y) 
        result= np.hstack((
                (-rho - (2*self.m[1] + self.m[0])*r/self.mu),
                (-rho + (self.m[1] + 2*self.m[0])*r/self.mu),
                (2*rho + (self.m[1] - self.m[0])*r/self.mu)
            ))/3
        assert np.abs(result.sum()) < atol
        return result
      
    def dH(self,y):  # r theta p l R Theta P L
        r = y[Hamiltonian.index_r]
        theta = y[Hamiltonian.index_theta]
        rho = y[Hamiltonian.index_rho]
        Theta = y[Hamiltonian.index_Theta]  
        p = y[Hamiltonian.index_p]
        l = y[Hamiltonian.index_l]
        P = y[Hamiltonian.index_P]
        L = y[Hamiltonian.index_L]          
        dV =self.dV(r,theta,rho,Theta)
        dH = np.zeros((8))
        dH[Hamiltonian.index_r] = p/self.g1                                             # (30a)
        dH[Hamiltonian.index_theta] = l/(self.g1*r**2)                                  # (30a)
        dH[Hamiltonian.index_rho] = P/self.g2                                           # (30c)
        dH[Hamiltonian.index_Theta] = L/(self.g2*rho**2)                                # (30c)
        dH[Hamiltonian.index_p] = l**2 / (self.g1*r**3) - dV[Hamiltonian.index_r]       # (30b)
        dH[Hamiltonian.index_l] = - dV[Hamiltonian.index_theta]                         # (30b)  
        dH[Hamiltonian.index_P] = L**2 / (self.g2*rho**3) - dV[Hamiltonian.index_rho]   # (30d)
        dH[Hamiltonian.index_L] = - dV[Hamiltonian.index_Theta]                         # (30d)
        return dH
 
    def get_r(self,r,theta,rho,Theta):
        r12 = r
        r23 = np.sqrt(rho**2 
                          - 2*(self.m[0]/self.mu)*r*rho*np.cos(Theta-theta)
                          + (self.m[0]/self.mu)**2 * r**2)
        r13 = np.sqrt(rho**2 + 
                          2*(self.m[1]/self.mu)*r*rho*np.cos(Theta-theta)
                          + (self.m[1]/self.mu)**2 * r**2)
        return r12,r23,r13
    
    def dV(self,r,theta,rho,Theta):  
        '''
        Calculate derivatives of potential energy
        '''
        r12,r23,r13 = self.get_r(r,theta,rho,Theta)

        T = np.array([
            self.m[0]*self.m[1]/r12**2,
            self.m[1]*self.m[2]/r23**2,
            self.m[0]*self.m[2]/r13**2
        ])
 
        cos_difference = np.cos(Theta-theta)
        sin_difference = np.sin(Theta-theta)
        S = np.c_[np.array([1,0,0,0]),
                  np.array([
                      -(self.m[0]/self.mu) * (cos_difference*rho - (self.m[0]/self.mu)*r),
                      -(self.m[0]/self.mu) * sin_difference*r*rho,
                      rho - (self.m[0]/self.mu)*r*cos_difference,
                      (self.m[0]/self.mu) *r*rho * sin_difference
                  ])  / r23                  ,
                  np.array([
                      (self.m[1]/self.mu) * (cos_difference*rho + (self.m[1]/self.mu)*r),
                      (self.m[1]/self.mu) * sin_difference*r*rho,
                      rho + (self.m[1]/self.mu)*r*cos_difference,
                      -(self.m[1]/self.mu) *r*rho * sin_difference
                ]) / r13                  ]    
        return self.G * np.dot(S,T)
    
    def get_total_energy(self,y):
        r = y[Hamiltonian.index_r]
        theta = y[Hamiltonian.index_theta]
        rho = y[Hamiltonian.index_rho]
        Theta = y[Hamiltonian.index_Theta]  
        p = y[Hamiltonian.index_p]
        l = y[Hamiltonian.index_l]
        P = y[Hamiltonian.index_P]
        L = y[Hamiltonian.index_L]  
        r12,r23,r13 = self.get_r(r,theta,rho,Theta)
        T = (p**2/(2*self.g1) + P**2/(2*self.g2) 
             + l**2/(2*self.g1*r**2) + L**2/(2*self.g2*rho**2))  #eq (29)
        V = -self.G * (self.m[0]*self.m[1]/r12* 
                       + self.m[1]*self.m[2]/r23 
                       + self.m[0]*self.m[2]/r13)
        return T + V

class Geometry: 
       
    @staticmethod
    def get_angle(vector):
        '''
        Determine the angle for a vector
        
        Parameters:
            velocity   Velocity vector
        '''
        return np.arctan2(vector[1], vector[0])

    @staticmethod
    def get_r_velocity(velocity, theta):
        '''
        Get the radial component of velocity
        
        Parameters:
            velocity   Velocity vector
            theta
        '''
        return np.dot(np.array([np.cos(theta),  np.sin(theta)]), velocity)
    
    @staticmethod
    def get_theta_dot(velocity, theta, r):
        '''
        Get the angular component of velocity
        
        Parameters:
            velocity   Velocity vector
            theta
            r
        '''    
        return np.dot(np.array([-np.sin(theta),  np.cos(theta)]), velocity)/r
    
    @staticmethod
    def get_polar_coordinates(r):                      
        return np.linalg.norm(r),Geometry.get_angle(r)
    
class Test1(TestCase):
    def test1(self):
        m = np.array([1.0, 1.0, 1.0])
        R = np.array([
            [0.97000436, -0.24308753],
            [0, 0],
            [-0.97000436, 0.24308753]
            ])
        R_dot = np.array([
            [0.46620369, 0.43236573],
            [-0.93240737, -0.86473146],
            [0.46620369, 0.43236573]            
            ])
        hamiltonian = Hamiltonian(m) 
        canonical = hamiltonian.create_initial_values(R,R_dot,m) 
        converted_back = hamiltonian.get_coordinates(canonical)
        np.testing.assert_almost_equal(R[0,:],converted_back[0:2],decimal=16)
        np.testing.assert_almost_equal(R[1,:],converted_back[2:4],decimal=16)
        np.testing.assert_almost_equal(R[2,:],converted_back[4:6],decimal=16)

if __name__ == '__main__':
    main()
    