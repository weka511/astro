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

import math, planet

class Solar:
    def __init__(self,planet):
        self.planet=planet
#   Beam Irradience in W/m2          (1)

    def beam_irradience(self,r):
        return self.planet.S/(r*r)
    
    def surface_irradience(self,areocentric_longitude,latitude,T):
        cos_zenith_angle = planet.cos_zenith_angle(areocentric_longitude,latitude,T)
        beam_irradience = self.beam_irradience(planet.instantaneous_distance(areocentric_longitude))
        return cos_zenith_angle*beam_irradience
    
if __name__=="__main__":
    import matplotlib.pyplot as plt
    def generate_points(planet,areocentric_longitude,latitude):
        x=[]
        y=[]
        for T in range(12,22):
            irradiance=solar.surface_irradience(areocentric_longitude,latitude,T)
            x.append(T)
            y.append(irradiance)
            if irradiance<0: break
        return (x,y)
     
    mars = planet.Mars()
    solar = Solar(mars)
    
    beam_irradience_top=solar.beam_irradience(mars.a)
    print "Mean beam irradience at top of atmosphere = {0:6.2f}".format(beam_irradience_top)
    
    plt.figure(3)
    plt.title("Mean beam irradience at top of Mars atmosphere")
    plt.xlabel("Areocentric longitude")
    plt.ylabel("Beam irradience")
    x=[]
    y=[]
    for i in range(360):
        x.append(i)
        y.append(solar.beam_irradience(mars.instantaneous_distance(i)))
    plt.plot(x,y)
    
    plt.figure(4)
    plt.title("Variation of solar declination angle")
    plt.xlabel("Areocentric longitude")
    plt.ylabel("Solar Declination Angle")
    x=[]
    y=[]
    for i in range(360):
        x.append(i)
        y.append(math.degrees(math.asin(mars.sin_declination(i))))
    plt.plot(x,y)

    plt.figure(6)
    plt.title("Diurnal Variation of Beam Irradience on a horizontal surface")
    plt.xlabel("Solar Time")
    plt.ylabel("Beam Irradiance")
    (x1,y1)=generate_points(mars,69,22.3)
    (x2,y2)=generate_points(mars,120,22.3)
    (x2,y2)=generate_points(mars,120,22.3)
    (x3,y3)=generate_points(mars,153,22.3)
    (x4,y4)=generate_points(mars,249,22.3)
    (x5,y5)=generate_points(mars,299,22.3)
    plt.plot(x1,y1,"r",x2,y2,"g",x3,y3,"b",x4,y4,"c",x5,y5,"m")
    plt.axis([12, 20, 0, 600])
    plt.grid(True)    
    plt.show()    