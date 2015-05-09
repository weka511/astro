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

import math, planet, solar, utilities

class Layer:
    next_id=0
    def __init__(self,name,latitude,longitude,thickness,depth,planet):
        self.id=Layer.next_id
        Layer.next_id+=1
        self.name=name
        self.latitude = latitude
        self.longitude = longitude
        self.depth = depth
        self.thickness = thickness
        self.temperature = planet.average_temperature
        self.planet = planet
    
    def propagate_temperature(self,above,below,areocentric_longitude,T,dT,record):
        raise NotImplementedError("propagate_temperature")
    
    # heat flow to me from neighbour
    def heat_flow(self,neighbour):
        temperature_gradient =                          \
            (self.temperature-neighbour.temperature) /  \
            (0.5*(self.thickness + neighbour.thickness))
        return - (self.planet.K * temperature_gradient)
            
    def update_temperature(self,nett_gain,dT,planet):    
        heat = nett_gain*dT*3600
        delta_temperature= heat / (self.planet.C * self.planet.rho * self.thickness)
        self.temperature += delta_temperature

    
    def __str__(self):
        return (                                  \
            "{0}: Latitude={1:6.1f},"            +\
            "Longitude={2:6.1f}, "               +\
            "Depth={3:6.3f}, "                   +\
            "Thickness={4:6.3f}, "               +\
            "Temperature = {5:6.1f}"  \
            ).format(self.name,
                     self.latitude,
                     self.longitude,
                     self.depth,
                     self.thickness,
                     self.top_temperature,
                     self.bottom_temperature)


class MedialLayer(Layer):
    def __init__(self,latitude,longitude,thickness,depth,planet):
        Layer.__init__(self,"Medial",latitude,longitude,thickness,depth,planet)
        
    def propagate_temperature(self,above,below,areocentric_longitude,T,dT,record):
        nett_gain = self.heat_flow(above) + self.heat_flow(below)
        self.update_temperature(nett_gain,dT,planet)
        record.add(self.temperature)

class Surface(Layer):
    stefan_bolzmann = 5.670374e-8
    
    def __init__(self,latitude,longitude,thickness,solar,planet):
        Layer.__init__(self,"Surface",latitude,longitude,thickness,0,planet)
        self.solar=solar
        
    def propagate_temperature(self,above,below,areocentric_longitude,T,dT,record):
        irradiance=self.planet.F*self.solar.surface_irradience(areocentric_longitude,self.latitude,T)
        outflow=self.bolzmann(self.temperature)
        nett_gain = irradiance - outflow + self.heat_flow(below)
        self.update_temperature(nett_gain,dT,planet)
        record.add(self.temperature)
    
    def bolzmann(self,t):
        t2=t*t
        return self.planet.E*Surface.stefan_bolzmann*t2*t2
        
class Bottom(Layer):
    def __init__(self,layer):
        Layer.__init__(self,"Bottom",layer.latitude,layer.longitude,layer.thickness,layer.depth,layer.planet)
        
    def propagate_temperature(self,above,below,areocentric_longitude,T,dT,record):
        nett_gain = self.heat_flow(above)
        self.update_temperature(nett_gain,dT,planet)
        record.add(self.temperature)
   
class ThermalModel:
    def __init__(self,latitude,longitude,spec,solar,planet,history):
        self.layers=[]
        self.planet=planet
        (n,dz)=spec[0]
        z=dz
        self.layers.append(Surface(latitude,longitude,dz,solar,planet))
        for (n,dz)in spec:
            for i in range(n):
                self.layers.append(MedialLayer(latitude,longitude,dz,z,planet))
                z+=dz
        bottom=self.layers.pop()
        self.layers.append(Bottom(bottom))
        self.history=history
        self.record=None
        self.zipper_layers = utilities.slip_zip(self.layers)

    
    def propagate_temperature(self,areocentric_longitude,T,dT):
        for above,layer,below in self.zipper_layers:
            layer.propagate_temperature(above,below,areocentric_longitude,T,dT,self.record)
            
    def runModel(self,start_day,number_of_days,number_of_steps_in_hour):
        step_size=1/float(number_of_steps_in_hour)
        for day in range(start_day,start_day+number_of_days):
            for hour in range(self.planet.hours_in_day):
                areocentric_longitude=self.planet.get_areocentric_longitude(day,hour)
                for step in range(number_of_steps_in_hour):
                    self.record = utilities.TemperatureRecord(day,hour,self.planet.hours_in_day)
                    self.propagate_temperature(areocentric_longitude,hour,step_size)
            self.history.add(self.record)


        
if __name__=="__main__":
    import matplotlib.pyplot as plt
    
    mars = planet.Mars()
    solar = solar.Solar(mars)
    history = utilities.InternalTemperatureLog()    
    thermal=ThermalModel(22.3,0,[(9,0.015),(10,0.3)],solar,mars,history)
    
    thermal.runModel(0,720,10)
    (days,surface_temp) = history.extract(0)
    (_,t1) = history.extract(1)
    (_,t2) = history.extract(2)
    (_,t3) = history.extract(3)
    (_,t4) = history.extract(4)
    (_,t5) = history.extract(5)
    (_,t6) = history.extract(10)
    (_,t7) = history.extract(19)
    plt.plot(days,surface_temp,days,t1,days,t2,days,t3,days,t4,days,t5,days,t6,days,t7)
    plt.show()
