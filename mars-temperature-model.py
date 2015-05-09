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

import thermalmodel, planet, solar, utilities

with open('output.txt', 'w') as f:
     mars = planet.Mars()
     solar = solar.Solar(mars)
     history = utilities.ExternalTemperatureLog(f)    
     thermal=thermalmodel.ThermalModel(22.3,0,[(9,0.015),(10,0.3)],solar,mars,history)
     
     thermal.runModel(0,720,10)     