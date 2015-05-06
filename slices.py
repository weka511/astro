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

import math, planet, solar

class Layer:
    def __init__(self,name,latitude,longitude,thickness,depth):
        self.name=name
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.thickness = thickness
        self.top_temperature =0
        self.bottom_temperature=0
    def __str__(self):
        return "{0}: Latitude={1:6.1f},"            +\
               "Longitude={2:6.1f}, "               +\
               "Depth={3:6.3f}, "                   +\
               "Thickness={4:6.3f}, "               +\
               "Temperatures = ({5:6.3f},{6:6.3f})". \
               format(self.name,self.latitude,self.longitude,self.depth,self.thickness,self.top_temperature,self.bottom_temperature)

class MedialLayer(Layer):
    def __init__(self,latitude,longitude,thickness,depth):
        Layer.__init__(self,"Medial",latitude,longitude,thickness,depth)
        
class Surface(Layer):
    def __init__(self,latitude,longitude):
        Layer.__init__(self,"Surface",latitude,longitude,0,0)
        
class Bottom(Layer):
    def __init__(self,layer):
        Layer.__init__(self,"Bottom",layer.latitude,layer.longitude,layer.thickness,layer.depth)
        
class ThermalModel:
    def __init__(self,latitude,longitude,spec):
        self.layers=[]
        z=0
        self.layers.append(Surface(latitude,longitude))
        for (n,dz)in spec:
            for i in range(n):
                self.layers.append(MedialLayer(latitude,longitude,dz,z))
                z+=dz
        bottom=self.layers.pop()
        self.layers.append(Bottom(bottom))
            
if __name__=="__main__":
    import matplotlib.pyplot as plt
    
    mars = planet.Mars()
    solar = solar.Solar(mars)
        
    thermal=ThermalModel(45,0,[(9,0.015),(10,0.3)])
    for layer in thermal.layers:
        print layer