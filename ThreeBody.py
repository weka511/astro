#!/usr/bin/env python

# Copyright (C) 2015,2026 Greenweaves Software Pty Ltd

# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this software.  If not, see <http://www.gnu.org/licenses/>

'''
Integrate 2 dimensional, but otherwise general, 3 body problem
'''

from argparse import ArgumentParser
from pathlib import Path
from time import time
import numpy as np
from matplotlib.pyplot import figure, show
from matplotlib import rcParams
from integrators import Hamiltonian, KotovychBowman
from utilities import get_r,newton_raphson,guarded_sqrt

__version__ = '1.0'
__author__ = 'Simon Crase'

class ThreeBody(Hamiltonian):
    '''
    Hamiltonian for 2 Dimensional, but otherwise general, 3 body problem
    '''
    def __init__(self,R,R_dot,m,G=1,clone=False):
        self.G = G
        self.M = m.sum()        # Section 5
        self.mu = m[0] + m[1]   # Section 5
        self.m = m
        self.g1 = m[0] * m[1] / self.mu    # Reduced masse - Section 5
        self.g2 = m[2] * self.mu / self.M  # Reduced masse - Section 5
        if clone:
            return

        r = R[1,:] - R[0,:]                 # Equation (27a)
        r_dot =  R_dot[1,:] - R_dot[0,:]
        r_polar = np.linalg.norm(r)
        theta = get_angle(r)
        p = self.g1 * get_r_velocity(r_dot, theta)
        l = self.g1 * r_polar**2 * get_theta_dot(r_dot, theta, r_polar)
        rho = self.M*R[2,:]/self.mu                       # Section 5
        assert np.linalg.norm(R[2,:] - R[0,:] - rho - (m[1]/self.mu)*r) < 1e-16,'Equation (27b)'    
        assert np.linalg.norm(R[2,:] - R[1,:] - rho + (m[0]/self.mu)*r) < 1e-16,'Equation (27c)'  
        rho_dot = self.M*R_dot[2,:]/self.mu         # Parallel to formula for rho
        rho_polar = np.linalg.norm(rho)
        Theta = get_angle(rho)
        P = self.g2 * get_r_velocity(rho_dot, Theta)                           # Linear momentum
        L = self.g2 * rho_polar**2 * get_theta_dot(rho_dot, Theta, rho_polar)  # angular mementum
        self.x = [r_polar, theta, rho_polar, Theta, p, l, P, L]
        self.transform()

    def dx(self):
        '''
        Calculate derivatives of r, theta, rho, Theta, p, l, P, L
        '''
        [r, theta, rho, Theta, p, l, P, L] = self.x
        [Vr, Vtheta, Vrho, VTheta] = self.dV(r, theta, rho, Theta)
        return np.array([
            p/self.g1,                           # Equation (30a)
            l/(self.g1*r**2),                    # Equation (30a)
            P/self.g2,                           # Equation (30c)
            L/(self.g2*rho**2),                  # Equation (30c)
            l**2/(self.g1*r**3) - Vr,            # Equation (30b)
            -Vtheta,                             # Equation (30b)
            L**2/(self.g2*rho**3) - Vrho,        # Equation (30d)
            -VTheta                              # Equation (302)
        ])

    def d_xi(self):
        [r_dot, theta_dot, rho_dot, Theta_dot, p_dot, l_dot, P_dot, L_dot] = self.dx()
        dH = self.dH()
        return np.array([dH[0], dH[1], dH[2], rho_dot, l_dot, L_dot, theta_dot, Theta_dot])

    def create(self, x):
        product = ThreeBody(np.zeros((3)), np.zeros((3)), self.m, self.G, True)
        product.g1 = self.g1
        product.g2 = self.g2
        product.x = self.x.copy()
        product.transform()
        return product

    def transform(self):
        '''
        Transform to the coordinates we will use for integration
        '''
        [r, theta, rho, Theta, p, l, P, L] = self.x
        self.xi = [
            (p**2 + (l/r)**2)/(2*self.g1),        # Equation (31a)
            (P**2 + (L/rho)**2)/(2*self.g2),      # Equation (31a)
            self.V(r, theta, rho, Theta),         # Equation (31b)
            rho,                                  # Equation (31b)
            l,                                    # Equation (31b)
            L,                                    # Equation (31b)
            theta,                                # Equation (31b)
            Theta,                                # Equation (31b)
         ]

    def _invert(self, hamiltonian):
        [r, _, _, _, p, _, P, _] = self.x
        [_, _, _, rho, l, L, theta, Theta] = self.xi
        r = self.g(r, self.xi[2], rho, theta, Theta)
        p = np.sign(p) * guarded_sqrt(2 * self.g1 * (self.xi[0] - l * l / (2.0 * self.g1 * r * r)))
        P = np.sign(P) * guarded_sqrt(2 * self.g2 * (self.xi[1] - L * L / (2.0 * self.g2 * rho * rho)))

        self.x = np.array([r, theta, rho, Theta, p, l, P, L])

    def get_energy(self):
        '''
        Calculate total energy    equation (29)
        '''
        [r, theta, rho, Theta, p, l, P, L] = self.x
        T = (p**2/self.g1 + P**2/self.g2 + l**2/(self.g1 * r**2) + L**2/(self.g2 * rho**2)) / 2 
        return (T + self.V(r, theta, rho, Theta))

    def dV(self, r, theta, rho, Theta):
        r23_sq = rho**2 - 2 * (self.m[0] / self.mu) * rho * r * np.cos(Theta - theta) + (self.m[0] / self.mu)**2 * r**2
        r23 = np.sqrt(r23_sq)
        dr23_drho = (rho - (self.m[0] / self.mu) * r * np.cos(Theta - theta)) / r23
        dr23_dr = (-(self.m[0] / self.mu) * rho * np.cos(Theta - theta) + (self.m[0] / self.mu)**2 * r) / r23
        dr23_dTheta = 2 * (self.m[0] / self.mu) * r * rho * np.sin(Theta - theta) / r23
        dr23_dtheta = - dr23_dTheta

        r31_sq = (rho**2 
                  + 2 * (self.m[2] / self.mu) * rho * r * np.cos(Theta - theta)  
                  +(self.m[2] / self.mu) **2 * r**2
                  )
        r31 = np.sqrt(r31_sq)
        dr31_drho = (rho + (self.m[1] / self.mu) * r * np.cos(Theta - theta)) / r31
        dr31_dr = ((self.m[1] / self.mu) * rho * np.cos(Theta - theta) 
                   + (self.m[1] / self.mu) * (self.m[1] / self.mu) * r) / r31

        dr31_dTheta = -2 * (self.m[1] / self.mu) * r * rho * np.sin(Theta - theta) / r31
        dr31_dtheta = - dr31_dTheta

        V23 = self.G * self.m[1] * self.m[2] / r23_sq
        V31 = self.G * self.m[2] * self.m[0] / r31_sq
        V12 = self.G * self.m[0] * self.m[1] / r**2

        return [
            V12 + V31 * dr31_dr + V23 * dr23_dr,
            V31 * dr31_dtheta + V23 * dr23_dtheta,
            V31 * dr31_drho + V23 * dr23_drho,
            V31 * dr31_dTheta + V23 * dr23_dTheta
        ]

    def V(self, r, theta, rho, Theta):
        '''
        Calculate potential energy
        '''
        r23 = np.sqrt(rho**2 -
                        2*(self.m[0]/self.mu)*rho*r*np.cos(Theta - theta) +
                        (self.m[0]/self.mu)**2*r**2)
        r31 = np.sqrt(rho**2+
                        2*(self.m[1]/self.mu)*rho*r*np.cos(Theta - theta) +
                        (self.m[1]/self.mu)**2*r**2)
        return (- self.G * self.m[0] * self.m[1] / r
                - self.G * self.m[1] * self.m[2] / r23 
                - self.G * self.m[2] * self.m[0] / r31)

    def g(self, r, xi3, rho, theta, Theta):
        return newton_raphson(r,
                              lambda r: self.V(r, theta, rho, Theta) - xi3,
                              lambda r: self.dV(r, theta, rho, Theta)[0],
                              1.0e-3,
                              1000)

    def dH(self):
        [r, theta, rho, Theta, p, l, P, L] = self.x
        [r_dot, theta_dot, rho_dot, Theta_dot, p_dot, l_dot, P_dot, L_dot] = self.dx()
        [Vr, Vtheta, Vrho, VTheta] = self.dV(r, theta, rho, Theta)
        return [
            p * p_dot / self.g1 + (l * r * r * l_dot - r * l * l * r_dot) / (self.g1 * r * r * r * r),
            P * P_dot / self.g2 + (L * rho * rho * L_dot - rho * L * L * rho_dot) / (self.g2 * rho * rho * rho * rho),
            Vr * r_dot + Vtheta * theta_dot + Vrho * rho_dot + VTheta * Theta_dot
        ]

    def inverse_jacobi(self):
        [r, theta, rho, Theta, p, l, P, L] = self.x
        rho_vector = rho * np.array([np.cos(Theta), np.sin(Theta)])
        r_vector = r * np.array([np.cos(theta), np.sin(theta)])
        r3 = (self.mu / self.M) * rho_vector
        r1 = r3 - rho_vector - (self.m[1] / self.mu) * r_vector
        r2 = r3 - rho_vector + (self.m[0] / self.mu) * r_vector       
        return (r1, r2, r3)

def parse_args():
    '''
    Parse command line arguments
    '''
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--figs', default='./figs', help=f'Path to plots')
    parser.add_argument('--show',default=False,action='store_true',help='Used to display figure')
    parser.add_argument('-N', default=80000, type=int)
    parser.add_argument('--step',default=1.0e-4,type=float)
    return parser.parse_args()

def get_angle(r):
    abs_theta = 0 if r[0] == 0 else np.atan(r[1] / r[0])
    return abs_theta + adjust_quadrant(r)

def get_r_velocity(velocity, theta):
    '''
    Get the radial component of velocity
    '''
    return np.dot(np.array([np.cos(theta),  np.sin(theta)]), velocity)


def get_theta_dot(zdot, theta, r):
    '''
    Get the angular component of velocity
    '''    
    return np.dot(np.array([-np.sin(theta),  np.cos(theta)]), zdot)/r

def adjust_quadrant(r):
    '''
    Determive the correct quadrant for a vector.
    
    Parameters:
        r       Vector
    Returns:
       Minumum angle for the quadrant
    '''
    if   r[0] >= 0 and r[1] >= 0:  return 0
    elif r[0] < 0 and r[1] >= 0: return np.pi / 2
    elif r[0] < 0 and r[1] < 0:  return np.pi
    else:
        return 3 * np.pi / 2 

def main():
    rcParams['text.usetex'] = True
    start  = time()    
    args = parse_args()
    hamiltonian = ThreeBody(                  #Equation (34)
        np.array([[0.97000436, -0.24308753],
                  [0, 0],
                  [-0.97000436, 0.24308753]]),
        np.array([[0.46620369, 0.43236573],
                  [-0.93240737, -0.86473146],
                  [0.46620369, 0.43236573]]),
        np.array([1.0, 1.0,1.0])
    )

    hamiltonian0 = hamiltonian.get_energy()
    uv = np.zeros((args.N,2))
    wz = np.zeros((args.N,2))
    rs = np.zeros((args.N,2))
    pq = np.zeros((args.N,2))
    integrator = KotovychBowman(args.step, hamiltonian)
    for i in range(args.N):
        integrator.integrate()
        (r1, r2, r3) = hamiltonian.inverse_jacobi()
        rs[i,:] = np.array(r1)
        pq[i,:] = np.array(r2)
        uv[i,:] = hamiltonian.x[0] * np.array([np.cos(hamiltonian.x[1]),np.sin(hamiltonian.x[1])])
        wz[i,:] = hamiltonian.x[2] * np.array([np.cos(hamiltonian.x[3]),np.sin(hamiltonian.x[3])])
  
    fig = figure(figsize=(8, 8))
    ax1 = fig.add_subplot(1, 1, 1)  
    ax1.plot(uv[:,0], uv[:,1], 'b',label='uv')
    ax1.plot(wz[:,0], wz[:,1], 'r',label='wz')
    ax1.plot(rs[:,0], rs[:,1], 'g',label='rs')
    ax1.plot(pq[:,0], pq[:,1], 'm',label='pq')
    ax1.legend()
    ax1.set_title(rf'N={args.N:,}, $\delta H=${hamiltonian.get_energy()-hamiltonian0:.4e}')
    fig.tight_layout(h_pad=2)
    fig.savefig(Path(args.figs)/Path(__file__).stem)    
    elapsed = time() - start
    minutes = int(elapsed/60)
    seconds = elapsed - 60*minutes
    print (f'Elapsed Time {minutes} m {seconds:.2f} s')
    if args.show:
        show()    

if __name__ == '__main__':
    main()