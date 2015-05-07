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
    next_id=0
    def __init__(self,name,latitude,longitude,thickness,depth):
        self.id=Layer.next_id
        Layer.next_id+=1
        self.name=name
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.thickness = thickness
        self.top_temperature = float('NaN')
        self.bottom_temperature = float('NaN')
        
    def inititialize_temperatures(self,t0):
        self.top_temperature = t0
        self.bottom_temperature = t0    
    
    def propagate_temperature(self,above,areocentric_longitude,T):
        return self
    
    def __str__(self):
        return (                                  \
            "{0}: Latitude={1:6.1f},"            +\
            "Longitude={2:6.1f}, "               +\
            "Depth={3:6.3f}, "                   +\
            "Thickness={4:6.3f}, "               +\
            "Temperatures = ({5:6.1f},{6:6.1f})"  \
            ).format(self.name,
                     self.latitude,
                     self.longitude,
                     self.depth,
                     self.thickness,
                     self.top_temperature,
                     self.bottom_temperature)


class MedialLayer(Layer):
    def __init__(self,latitude,longitude,thickness,depth):
        Layer.__init__(self,"Medial",latitude,longitude,thickness,depth)
        
class Surface(Layer):
    stefan_bolzmann = 5.670374e-8
    def __init__(self,latitude,longitude,solar):
        Layer.__init__(self,"Surface",latitude,longitude,0,0)
        self.solar=solar
        
    def propagate_temperature(self,above,areocentric_longitude,T):
        previous= Layer.propagate_temperature(self,above,areocentric_longitude,T)
        irradiance=solar.surface_irradience(areocentric_longitude,self.latitude,T)
        t2=self.top_temperature*self.top_temperature
        outflow=Surface.stefan_bolzmann*t2*t2
        nett_gain = irradiance - outflow
        print "Irradiance ={0:6.1f}, outflow={1:6.1f}, nett_gain={2:6.1f}".format(irradiance,outflow,nett_gain)
        return previous
    
class Bottom(Layer):
    def __init__(self,layer):
        Layer.__init__(self,"Bottom",layer.latitude,layer.longitude,layer.thickness,layer.depth)
        
class ThermalModel:
    def __init__(self,latitude,longitude,spec,solar):
        self.layers=[]
        z=0
        self.layers.append(Surface(latitude,longitude,solar))
        for (n,dz)in spec:
            for i in range(n):
                self.layers.append(MedialLayer(latitude,longitude,dz,z))
                z+=dz
        bottom=self.layers.pop()
        self.layers.append(Bottom(bottom))
        
    def inititialize_temperatures(self,t0):
        for layer in self.layers:
            layer.inititialize_temperatures(t0)
            
    def propagate_temperature(self,areocentric_longitude,T):
        above=None
        for layer in self.layers:
            above=layer.propagate_temperature(above,areocentric_longitude,T)


             
if __name__=="__main__":
    import matplotlib.pyplot as plt
    
    mars = planet.Mars()
    solar = solar.Solar(mars)
        
    thermal=ThermalModel(22.3,0,[(9,0.015),(10,0.3)],solar)
    thermal.inititialize_temperatures(150)
    
    #for layer in thermal.layers:
        #print layer
     
    thermal.propagate_temperature(153,12)
    
    thermal.propagate_temperature(153,18)
    
    thermal.propagate_temperature(153,20)
    
    #for layer in thermal.layers:
        #print layer    