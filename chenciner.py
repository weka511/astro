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

from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from time import time
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
import numpy as np
from rki import ImplicitRungeKutta4,Driver

__version__ = '1.0'
__author__ = 'Simon Crase'


class Hamiltonian:
    index_r = 0
    index_theta = 1
    index_rho = 2
    index_Theta = 3
    index_p = 4
    index_l = 5    
    index_P = 6
    index_L = 7
    def __init__(self,m,G=1,clone=False,atol=1e-16):
        self.G = G
        self.M = m.sum()        # Section 5
        self.mu = m[0] + m[1]   # Section 5
        self.m = m
        self.g1 = m[0]*m[1]/self.mu    # Reduced mass - Section 5 - just after (28) 
        self.g2 = m[2]*self.mu/self.M  # Reduced mass - Section 5 - just after (28) 
        
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
        return np.array([
            p/self.g1,
            l/(self.g1*p**2),
            l**2 / (self.g1*r**3) - dV[Hamiltonian.index_r],
            - dV[Hamiltonian.index_theta],
            P/self.g1,
            L/(self.g1*rho**2),
            L**2 / (self.g1*rho**3) - dV[Hamiltonian.index_rho],
            - dV[Hamiltonian.index_Theta],            
        ])
    
    def dV(self,r,theta,rho,Theta):  
 
        r12 = r
        r23 = np.sqrt(rho**2 - 2*(self.m[0]/self.mu)*r*rho*np.cos(Theta-theta) + (self.m[0]/self.mu)**2*r**2)
        r13 = np.sqrt(rho**2 + 2*(self.m[1]/self.mu)*r*rho*np.cos(Theta-theta) + (self.m[1]/self.mu)**2*r**2)
        T = np.array([
            self.m[0]*self.m[1]/r12**2,
            self.m[1]*self.m[2]/r23**2,
            self.m[0]*self.m[2]/r13**2
        ])
 
        S = np.c_[np.array([1,0,0,0]),
                  np.array([
                      -(self.m[0]/self.mu)*(np.cos(Theta-theta)*rho-(self.m[0]/self.mu)*r),
                      -(self.m[0]/self.mu)*np.sin(Theta-theta)*r*rho,
                      rho - (self.m[0]/self.mu)*r*np.sin(Theta-theta),
                      (self.m[0]/self.mu)*np.sin(Theta-theta)*r*rho])/r23                  ,
                  np.array([
                      (self.m[1]/self.mu)*(np.cos(Theta-theta)*rho+(self.m[1]/self.mu)*r),
                      (self.m[1]/self.mu)*np.sin(Theta-theta)*r*rho,
                      rho + (self.m[1]/self.mu)*r*np.sin(Theta-theta),
                      -(self.m[1]/self.mu)*np.sin(Theta-theta)*r*rho])/r13                  ]    
        return self.G * np.dot(S,T)

class Geometry: 
    @staticmethod
    def adjust_quadrant(vector):
        '''
        Determine the correct quadrant for a vector.
        
        Parameters:
            vector      Velocity vectpr
        Returns:
            Angle for start of quadrant
        '''
        if   vector[0] >= 0 and vector[1] >= 0:  return 0
        elif vector[0] < 0 and vector[1] >= 0: return np.pi / 2
        elif vector[0] < 0 and vector[1] < 0:  return np.pi
        else:
            return 3 * np.pi / 2 
    
    @staticmethod
    def get_angle(vector):
        '''
        Determine the angle for a vector
        
        Parameters:
            velocity   Velocity vector
        '''
        return Geometry.adjust_quadrant(vector) + (np.atan(vector[1] / vector[0]) if vector[0] != 0 else 0)
    
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
        
def create_initial_values(hamiltonian,R,R_dot,m):
    r_polar,theta = Geometry.get_polar_coordinates( R[1,:] - R[0,:] ) # Equation (27a)
    r_dot =  R_dot[1,:] - R_dot[0,:]
    p = hamiltonian.g1 * Geometry.get_r_velocity(r_dot, theta)                            # Linear momentum
    l = hamiltonian.g1 * r_polar**2 * Geometry.get_theta_dot(r_dot, theta, r_polar)       # angular momentum
    centre_of_mass01 = np.dot(m[:2],R[:2,:])/hamiltonian.mu
    rho_dot = R_dot[2,:] - np.dot(m[:2],R_dot[:2,:])/hamiltonian.mu  # Parallel to formula for rho
    rho_polar,Theta = Geometry.get_polar_coordinates( R[2,:] - centre_of_mass01)
    P = hamiltonian.g2 * Geometry.get_r_velocity(rho_dot, Theta)                           # Linear momentum
    L = hamiltonian.g2 * rho_polar**2 * Geometry.get_theta_dot(rho_dot, Theta, rho_polar)  # angular momentum
    
    return np.array([r_polar,theta,rho_polar,Theta,p,l,P,L])

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('file_name',help='File name for initial conditions')
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('--data', default='./data', help=f'Path to data files')
    parser.add_argument('-N','--Iterations',default=100,type=int,help='Number of steps')
    return parser.parse_args()

def read_data(file_name,dim=2):
    '''
    Read initial conditions from a file
    '''
    state = -1
    i = 0
    with open(file_name) as f:
        data_reader = reader(f)
        for row in data_reader:
            values = np.array([float(m) for m in row])
            match(state):
                case -1:
                    masses = values
                    Q = np.zeros((len(masses),dim))
                    P = np.zeros((len(masses),dim))
                    state = 0
                case 0:
                    Q[i,:] = values
                    i += 1
                    if i == len(masses):
                        i = 0
                        state = 1
                case 1:
                    P[i,:] = values
                    i += 1
 
        return Q,P,masses
    
def create_vectors(y):
    r = y[Hamiltonian.index_r]
    theta = y[Hamiltonian.index_theta]
    rho = y[Hamiltonian.index_rho]
    Theta = y[Hamiltonian.index_Theta]
    return r*np.array([np.cos(theta),np.sin(theta)]),rho*np.array([np.cos(Theta),np.sin(Theta)])

def main():
    '''
    Read initial values and integrate equations of motion
    '''
    rcParams['text.usetex'] = True
    start  = time()
    args = parse_args()
    R,R_dot,m = read_data(Path(args.data)/args.file_name)
    hamiltonian = Hamiltonian(m) 
    y = create_initial_values(hamiltonian,R,R_dot,m)
    integrator = ImplicitRungeKutta4(lambda y: hamiltonian.dH(y), 10, 1e-9)
    driver = Driver(integrator, h=0.1)
    XY = np.zeros((args.Iterations+1,6))
    r,rho = create_vectors(y) 
    XY[0,0:2] = (-rho - (2*hamiltonian.m[1]+hamiltonian.m[0])*r/hamiltonian.mu)/3
    XY[0,2:4] = (-rho + (hamiltonian.m[1]+2*hamiltonian.m[0])*r/hamiltonian.mu)/3
    XY[0,4:] = (2*rho + (hamiltonian.m[1]-hamiltonian.m[0])*r/hamiltonian.mu)/3    
    for i in range(args.Iterations):
        y = driver.step(y)
        r,rho = create_vectors(y) 
        XY[i+1,0:2] = (-rho - (2*hamiltonian.m[1]+hamiltonian.m[0])*r/hamiltonian.mu)/3
        XY[i+1,2:4] = (-rho + (hamiltonian.m[1]+2*hamiltonian.m[0])*r/hamiltonian.mu)/3
        XY[i+1,4:] = (2*rho + (hamiltonian.m[1]-hamiltonian.m[0])*r/hamiltonian.mu)/3
    fig = figure(figsize=(12,12))
    fig.suptitle(Path(__file__).stem)
    ax1 = fig.add_subplot(2,2,1,adjustable='box',aspect=1.0)
    ax1.plot(XY[0,:],XY[1,:],'r')
    ax1.scatter(R[0,0],R[0,1],c='r')
    ax1.scatter(XY[0,0],XY[0,1],c='r',marker='X')
    ax2 = fig.add_subplot(2,2,2,adjustable='box',aspect=1.0)
    ax2.plot(XY[2,:],XY[3,:],'g')
    ax2.scatter(R[1,0],R[1,1],c='g')
    ax3 = fig.add_subplot(2,2,3,adjustable='box',aspect=1.0)
    ax3.plot(XY[4,:],XY[5,:],'b')
    ax3.scatter(R[2,0],R[2,1],c='b')
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
