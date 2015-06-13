# Copyright (C) 2015 Greenweaves Software Pty Ltd

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

import math, utilities, integrators



class ThreeBody(integrators.Hamiltonian):
    
    def __init__(self,r1,r2,r3,r1dot, r2dot,r3dot,m1,m2,m3,G=1,clone=False):
        self.G=G
        self.M = m1 + m2 + m3
        self.mu = m1 + m2
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.g1 = m1*m2/self.mu
        self.g2 = m3 * self.mu /self.M
        if clone: return

        r = [r2_-r1_ for (r1_,r2_) in zip(r1,r2)]
        r_dot = [r2_-r1_ for (r1_,r2_) in zip(r1dot,r2dot)]        
        r_polar=utilities.get_r(r)
        theta=utilities.get_angle(r)
        p = self.g1 * utilities.get_r_velocity(r_dot,theta)
        l = self.g1 * r_polar * r_polar * utilities.get_theta_dot(r_dot,theta,r_polar)

        rho = [self.M*r3_/self.mu for r3_ in r3]
        rho_dot = [self.M*r3_/self.mu for r3_ in r3dot]        
        rho_polar=utilities.get_r(rho)
        Theta=utilities.get_angle(rho)
        P = self.g2 * utilities.get_r_velocity(rho_dot,Theta)
        L = self.g2 * rho_polar * rho_polar * utilities.get_theta_dot(rho_dot,Theta,rho_polar)

        self.x = [r_polar,theta,rho_polar,Theta,p,l,P,L]
        self.transform()
        
    def dx(self):
            [r,theta,rho,Theta,p,l,P,L] = self.x
            [Vr,Vtheta,Vrho,VTheta]=self.dV(r,theta,rho,Theta)
            return [
                p/self.g1,
                l/(self.g1*r*r),
                P/self.g2,
                L/(self.g2*rho*rho),
                l*l/(self.g1*r*r*r)-Vr,
                -Vtheta,
                L*L/(self.g2*rho*rho*rho)-Vrho,
                -VTheta
            ]
     
    def d_eta(self):
        [r_dot,theta_dot,rho_dot,Theta_dot,p_dot,l_dot,P_dot,L_dot] = self.dx()
        hamiltonian_components = self.hamiltonian_components()
        return [
             hamiltonian_components[0],
             hamiltonian_components[1],
             hamiltonian_components[2],
             rho_dot,
             l_dot,
             L_dot,
             theta_dot,
             Theta_dot
         ]
     
    def create(self,x):
        product = ThreeBody(0,0,0,0,0,0,self.m1,self.m2,self.m3,self.G,True)
        product.g1 = self.g1
        product.g2 = self.g2
        product.x  = [x0 for x0 in self.x]
        product.transform()        
        return product
  
    def transform(self):
        [r,theta,rho,Theta,p,l,P,L] = self.x
        self.eta = [
             (p*p + (l/r)*(l/r))/(2*self.g1), 
             (P*P + (L/rho)*(L/rho))/(2*self.g2),
             self.V(r,theta,rho,Theta),
             rho,
             l,
             L,
             theta,
             Theta,
         ]
         
    def invert(self,hamiltonian):
        [r,_,_,_,p,_,P,_]           = self.x
        [_,_,_,rho, l,L,theta,Theta] = self.eta
        r = self.g(r,self.eta[2],rho,theta,Theta)
        p = utilities.signum(p)*utilities.guarded_sqrt(2*self.g1*(self.eta[0]-l*l/(2.0*self.g1*r*r)))
        P = utilities.signum(P)*utilities.guarded_sqrt(2*self.g2*(self.eta[1]-L*L/(2.0*self.g2*rho*rho)))
         
        self.x = [r, theta, rho, Theta, p, l, P, L ]
         
    def hamiltonian(self):
        [r,theta,rho,Theta,p,l,P,L] = self.x
        return (p*p/self.g1 +P*P/self.g2 +l*l/(self.g1*r*r)+L*L/(self.g2*rho*rho))/2 + self.V(r,theta,rho,Theta)

    def dV(self,r,theta,rho,Theta):
        r23_sq      = rho*rho - 2*(self.m1/self.mu)*rho*r*math.cos(Theta-theta) + (self.m1/self.mu)*(self.m1/self.mu)*r*r  
        r23         = math.sqrt(r23_sq)       
        dr23_drho   = (rho - (self.m1/self.mu)*r*math.cos(Theta-theta))/r23
        dr23_dr     = (-(self.m1/self.mu)*rho*math.cos(Theta-theta) + (self.m1/self.mu)*(self.m1/self.mu)*r)/r23
        dr23_dTheta = 2 * (self.m1/self.mu)*r*rho * math.sin(Theta-theta)/r23
        dr23_dtheta = - dr23_dTheta        

        r31_sq      = rho*rho + 2*(self.m2/self.mu)*rho*r*math.cos(Theta-theta) + (self.m2/self.mu)*(self.m2/self.mu)*r*r
        r31         = math.sqrt(r31_sq)       
        dr31_drho   = (rho + (self.m2/self.mu)*r*math.cos(Theta-theta))/r31
        dr31_dr     = ((self.m2/self.mu)*rho*math.cos(Theta-theta) + (self.m2/self.mu)*(self.m2/self.mu)*r)/r31

        dr31_dTheta = -2 * (self.m2/self.mu)*r*rho * math.sin(Theta-theta)/r31
        dr31_dtheta = - dr31_dTheta
        
        V23 = self.G*self.m2*self.m3/r23_sq
        V31 = self.G*self.m3*self.m1/r31_sq
        V12 = self.G*self.m1*self.m2/(r*r)
                
        return [
                V12 + V31*dr31_dr + V23*dr23_dr,
                V31*dr31_dtheta +V23*dr23_dtheta,
                V31*dr31_drho +V23*dr23_drho,
                V31*dr31_dTheta +V23*dr23_dTheta
        ]
        
 
   
    def V(self,r,theta,rho,Theta):
        r23_sq      = rho*rho - 2*(self.m1/self.mu)*rho*r*math.cos(Theta-theta) + (self.m1/self.mu)*(self.m1/self.mu)*r*r  
        r23         = math.sqrt(r23_sq)         
        r31_sq      = rho*rho + 2*(self.m2/self.mu)*rho*r*math.cos(Theta-theta) + (self.m2/self.mu)*(self.m2/self.mu)*r*r
        r31         = math.sqrt(r31_sq)              
        return -self.G*self.m1*self.m2/r - self.G*self.m2*self.m3/r23  - self.G*self.m3*self.m1/r31
        


    def g(self,r,eta3,rho,theta,Theta):
        return utilities.newton_raphson(r,\
                                        lambda r: self.V(r,theta,rho,Theta)-eta3, \
                                        lambda r: self.dV(r,theta,rho,Theta)[0], \
                                        1.0e-3,
                                        1000)
    
    def hamiltonian_components(self):
        [r,theta,rho,Theta,p,l,P,L] = self.x
        [r_dot,theta_dot,rho_dot,Theta_dot,p_dot,l_dot,P_dot,L_dot] = self.dx()
        [Vr,Vtheta,Vrho,VTheta]=self.dV(r,theta,rho,Theta)        
        return [
            p*p_dot/self.g1 + (l*r*r*l_dot-r*l*l*r_dot)/(self.g1*r*r*r*r),
            P*P_dot/self.g2 + (L*rho*rho*L_dot-rho*L*L*rho_dot)/(self.g2*rho*rho*rho*rho),
            Vr*r_dot + Vtheta*theta_dot + Vrho*rho_dot + VTheta*Theta_dot
        ]

    def inverse_jacobi(self):
        r3 = [(self.mu/self.M)*rho for rho in self.rho]
        r1 = [rr3-r]
        
if __name__=='__main__':
    import matplotlib.pyplot as plt, sys
    
    u=[]
    v=[]
    w=[]
    z=[]
    hamiltonian=ThreeBody(
        [0.97000436,-0.24308753],
        [0,0],
        [-0.97000436,0.24308753],
        [0.46620369,0.43236573],
        [-0.93240737,-0.86473146],
        [0.46620369,0.43236573],
        1.0,
        1.0,
        1.0
    )
    
    integrator = integrators.Integrate2(1.0e-5,hamiltonian)
    
    nn = 2000000
    mm = 1000
    print hamiltonian.hamiltonian()
    try:
        for i in range(nn):
            integrator.integrate()
            if i>mm:
                u.append(hamiltonian.x[0]*math.cos(hamiltonian.x[1]))
                v.append(hamiltonian.x[0]*math.sin(hamiltonian.x[1]))
                w.append(hamiltonian.x[2]*math.cos(hamiltonian.x[3]))
                z.append(hamiltonian.x[2]*math.sin(hamiltonian.x[3]))        
        plt.plot(u,v,'b',w,z,'r')
        print hamiltonian.hamiltonian()
    except:
        plt.plot(u,v,'b',w,z,'r')
        print "Unexpected error:", sys.exc_info()[0]        