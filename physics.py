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

class CO2:
    condensation_temperature = 190
    latent_heat = 574               # http://www.engineeringtoolbox.com/fluids-evaporation-latent-heat-d_147.html    
    albedo = 0.6
    
class Radiation:
    stefan_bolzmann = 5.670374e-8
    
    @staticmethod
    def bolzmann(t):
        t2=t*t
        return Radiation.stefan_bolzmann*t2*t2    