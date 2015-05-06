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

# TODO:

#The daily and yearly position of the sun

#Heat flow into and out of the subsurface

#The effect of condensation and evaporation on the energy balance

##################################################################
#The daily and yearly position of the sun

#The total amount of volatile at the surface at each point in time

# The Solar model is based on NASA Technical Memorandum 102299
# Solar Radiation on Mars
# Joseph Appelbaum & Dennis J Flood


import math, matplotlib.pyplot as plt



class Planet:
    def __init__(self):
        self.S = 1371 # Solar constant at the mean Sun-Earth distance of l AU, in N/m2
        self.a = 1.0 # the  semimajor axis in AU,
        self.e = 0.017 #  eccentricity
        self.obliquity = 23.4
        self.hours_in_day = 24

#   Beam Irradience in W/m2          (1)

    def beam_irradience(self,r):
        return self.S/(r*r)
    
#   Instantaneaous Distance from Sun in AU

    def instantaneous_distance(self,areocentric_longitude):
        theta = math.radians(areocentric_longitude - 248)
        return self.a*(1-self.e*self.e)/(1 + self.e * math.cos(theta))
    
    def sin_declination(self,areocentric_longitude):
        return math.sin(math.radians(self.obliquity)) * \
               math.sin(math.radians(areocentric_longitude))
        
    def cos_zenith_angle(self,areocentric_longitude,latitude,T):
        sin_declination=self.sin_declination(areocentric_longitude)
        cos_declination=math.sqrt(1-sin_declination*sin_declination)
        return math.sin(math.radians(latitude))*sin_declination +            \
            math.cos(math.radians(latitude))*cos_declination*math.cos(math.radians(self.hour_angle(T)))

    def hour_angle(self,T):
        return 360*T/self.hours_in_day-180
    
class Mars(Planet):
    def __init__(self):
        Planet.__init__(self)
        self.a = 1.5236915
        self.e = 0.093377
        self.obliquity = 24.936
        self.hours_in_day = 24 # should be 24.65
        
if __name__=="__main__":
    def generate_points(areocentric_longitude,latitude):
        x=[]
        y=[]
        for T in range(12,22):
            cos_zenith_angle = planet.cos_zenith_angle(areocentric_longitude,latitude,T)
            beam_irradience = planet.beam_irradience(planet.instantaneous_distance(areocentric_longitude))
            irradiance = cos_zenith_angle*beam_irradience
            x.append(T)
            y.append(irradiance)
            if irradiance<0: break
        return (x,y)
     
    planet = Mars()
    beam_irradience_top=planet.beam_irradience(planet.a)
    print "Mean beam irradience at top of atmosphere = {0:6.2f}".format(beam_irradience_top)
    
    plt.figure(3)
    plt.title("Mean beam irradience at top of Mars atmosphere")
    plt.xlabel("Areocentric longitude")
    plt.ylabel("Beam irradience")
    x=[]
    y=[]
    for i in range(360):
        x.append(i)
        y.append(planet.beam_irradience(planet.instantaneous_distance(i)))
    plt.plot(x,y)
    
    plt.figure(4)
    plt.title("Variation of solar declination angle")
    plt.xlabel("Areocentric longitude")
    plt.ylabel("Solar Declination Angle")
    x=[]
    y=[]
    for i in range(360):
        x.append(i)
        y.append(math.degrees(math.asin(planet.sin_declination(i))))
    plt.plot(x,y)

    plt.figure(6)
    plt.title("Diurnal Variation of Beam Irradience on a horizontal surface")
    plt.xlabel("Solar Time")
    plt.ylabel("Beam Irradiance")
    (x1,y1)=generate_points(69,22.3)
    (x2,y2)=generate_points(120,22.3)
    (x3,y3)=generate_points(153,22.3)
    (x4,y4)=generate_points(249,22.3)
    (x5,y5)=generate_points(299,22.3)
    plt.plot(x1,y1,"r",x2,y2,"g",x3,y3,"b",x4,y4,"c",x5,y5,"m")
    plt.axis([12, 20, 0, 600])
    plt.grid(True)    
    plt.show()