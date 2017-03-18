# Copyright (C) 2015-2017 Greenweaves Software Pty Ltd

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

# Used to convert between CGS and SI
class Conversion:
    cm_per_metre = 100
    gm_per_Kg = 1000
    cm3_per_meter3 = cm_per_metre*cm_per_metre*cm_per_metre
    
class CO2:
    condensation_temperature = 145
    latent_heat = 574               # http://www.engineeringtoolbox.com/fluids-evaporation-latent-heat-d_147.html    
    albedo = 0.6
    
class Radiation:
    stefan_bolzmann = 5.670374e-8
    
    @staticmethod
    def bolzmann(t):
        t2=t*t    # we call this function often, so don't use exponentiation
        return Radiation.stefan_bolzmann*t2*t2 
    
    @staticmethod
    def reverse_bolzmann(radiation):
        return (radiation/Radiation.stefan_bolzmann)**0.25
    
if __name__=='__main__':
    print (Radiation.reverse_bolzmann(1350*(1-0.3)/4))