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

import math

def denom(x1,x2,y1,y2):
    return ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))**(3.0/2.0)

    
def dx(x,G,M1,M2):
    x1=x[0]
    y1=x[1]
    x2=x[2]
    y2=x[3]    
    x3=x[4]
    y3=x[5]
    denom12=denom(x1,x2,y1,y2)

    denom31=denom(x3,x1,y3,y1)
    denom23=denom(x2,x3,y2,y3)
    return [
        x[6], #x1'
        x[7], #y1'
        x[8], #x2'
        x[9],  
        x[10],
        x[11],
        -G*M2*(x1-x2)/denom12, #dx1'
        -G*M2*(y1-y2)/denom12, #dy1'
        -G*M1*(x2-x1)/denom12, #dx2'
        -G*M1*(y2-y1)/denom12, #dy2'
        -G*M1*(x3-x1)/denom31 - G*M2*(x3-x2)/denom23, #dx3'
        -G*M1*(y3-y1)/denom31 - G*M2*(y3-y2)/denom23, #dy3'
    ]

def T(x,G,M1,M2):
    return M1*(x[6]*x[6]+x[7]*x[7]) + M2*(x[8]*x[8]+x[9]*x[9])

if __name__=='__main__':
    import rki, matplotlib.pyplot as plt
    
    nn=100000
    h = 0.001
    M1 = 250000
    epsilon=0.1
    M2 = 1.0
    G = 1.0
    L = 25.0
    R2 = L*M1/(M1+M2)
    x2 = R2
    y2 = 0
    x1 = -L*M2/(M1+M2)
    y1 = 0.0
    x3 = 0.5*(x1+x2)*(1+epsilon)
    y3 = L* math.sqrt(3.0)/2.0
    omega = math.sqrt(G*(M1+M2)/(L*L*L))
    x=[
        x1, y1,
        x2, y2,
        x3, y3, 
        0, x1*omega,
        0, x2*omega,
        -0.5*math.sqrt(3)*L*omega, 0.5*L*omega
    ]
    
    rk = rki.ImplicitRungeKutta4(lambda (x): dx(x,1,M1,M2),200,1e-18)
    
    driver = rki.Driver(rk,0.000000001,0.5,1.0,0.000000001)
    
    x1s=[]
    y1s=[]
    x2s=[]
    y2s=[]    
    x3s=[]
    y3s=[]
    for i in range(nn):
        x= driver.step(x)
        x1s.append(x[0])
        y1s.append(x[1])
        x2s.append(x[2])
        y2s.append(x[3])        
        x3s.append(x[4])
        y3s.append(x[5])
    plt.plot(x1s,y1s,'b', x2s,y2s,'g', x3s,y3s,'r') 
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('M1={0:5.1f},M2={1:5.0f},R2={2:5.1f},L={3:5.1f},nn={4:5.1f},h={5:.1e}'.format(M1,M2,R2,L,nn,h))