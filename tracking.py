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

import os, re, tempfile
from array import array

class Reader(object):
    def __init__(self,partner):
        self.names=[]
        self.current=0
        for name in os.listdir(partner.path):
            m=re.match(partner.header+'(\\d+)',name)
            if m:
                self.names.append(name)
                
    def __iter__(self):
        return self
    def next(self):
        if self.current<len(self.names):
            self.current+=1
            return self.names[self.current-1]
        else:
            raise StopIteration
    
class Tracking(object):
    def __init__(self,header,suffix='dat',path=tempfile.gettempdir()):
        self.latest = ("",-1)
        self.header = header
        self.path   = path
        self.suffix = suffix
        
        for name in os.listdir(path):
            m=re.match(header+'(\\d+)',name)
            if m:
                seq=int(m.group(1))
                (nn,ii)=self.latest
                if seq>ii:
                    self.latest=(name,seq)
    
    def get_reader(self):
        return Reader(self)
            
    def output_status(self,status):
        (_,seq)=self.latest
        seq+=1
        name = '%(header)s%(seq)08d.%(suffix)s'%{
            'header' : self.header,
            'seq'    : seq,
            'suffix' : self.suffix
        }
        with open(os.path.join(self.path,name),'wb') as f:
            status_array = array('d', status)
            status_array.tofile(f)        
    
    def get_status(self):
        (name,_)=self.latest
        return self.get_one_status(name)
    
    def get_one_status(self,name):
        try:
            with open(os.path.join(self.path,name),'rb') as f:
                status_array = array('d')
                status_array.fromstring(f.read())        
                return status_array.tolist()
        except IOError:
            return []

if __name__=='__main__':
    import math, restricted,rki
    tracking=Tracking('run')
        
    old_trace= tracking.get_status()
    
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
    x=[]
    if len(old_trace)==0:
        x=[
            x1, y1,
            x2, y2,
            x3, y3, 
            0, x1*omega,
            0, x2*omega,
            -0.5*math.sqrt(3)*L*omega, 0.5*L*omega
        ]
    else:
        for i in range(12):
            x.append(old_trace[len(old_trace)-12+i])

    trace=[]
    
    rk = rki.ImplicitRungeKutta4(lambda (x): restricted.dx(x,1,M1,M2),200,1e-18)
    
    driver = rki.Driver(rk,1.e-8,0.5,1.0,1e-12)
    
    trace=trace+x
    for i in range(nn):
        x= driver.step(x)
        trace=trace+x
    tracking.output_status(trace)
    