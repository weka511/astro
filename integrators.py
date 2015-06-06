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

class Hamiltonian:
    def dx(self):
        raise NotImplementedError
    def d_eta(self):
        raise NotImplementedError
    def create(self,x):
        raise NotImplementedError
    def invert(self,hamiltonian):
        raise NotImplementedError
    def hamiltonian(self):
        raise NotImplementedError
    
class Integrate2:
    def __init__(self,h,hamiltonian):
        self.h=h
        self.hamiltonian=hamiltonian
        
    def predict(self):
        return [x0 + self.h * f0 for (x0,f0) in zip(self.hamiltonian.x,self.hamiltonian.dx())]    
    
    def correct(self,hamiltonian1):
        return [e0 + (self.h/2) * (f0 + f1) for (e0,f0,f1) in zip(self.hamiltonian.eta,self.hamiltonian.d_eta(),hamiltonian1.d_eta())]
    
    def integrate(self):
        hamiltonian1=self.hamiltonian.create(self.predict())
        self.hamiltonian.eta=self.correct(hamiltonian1)
        self.hamiltonian.invert(hamiltonian1)
        