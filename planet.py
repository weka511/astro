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

import math

class Planet:
    def __init__(self,name):
        self.name = name
        self.S = 1371 # Solar constant at the mean Sun-Earth distance of l AU, in N/m2
        self.a = 1.0 # the  semimajor axis in AU,
        self.e = 0.017 #  eccentricity
        self.obliquity = 23.4
        self.hours_in_day = 24
        self.average_temperature = 300
        
    def __str__(self):
        return ("{0}\nSemimajor axis = {1:9.7f} AU\n" +\
               "eccentricity = {2:8.6f}\n" + \
               "obliquity = {3:6.3f}\n" + \
               "hours in day = {4:6.4f}\n" + \
               "absorption = {5:4.2f}\n"+ \
               "emissivity={6:4.2f}\n" + \
               "conductivity={7:5.1f}\n" + \
               "specific heat={8:6.1f}\n" + \
               "rho={9:6.1f}\n" \
               "average temperature={10:5.1f}" \
               ).format(\
            self.name,   \
            self.a, \
            self.e, \
            self.obliquity, \
            self.hours_in_day, \
            self.F , \
            self.E,  \
            self.K, \
            self.C, \
            self.rho, \
            self.average_temperature \
        )
  
      
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
        Planet.__init__(self,"Mars")
        self.a = 1.5236915
        self.e = 0.093377
        self.obliquity = 24.936
        self.hours_in_day = 24 # should be 24.65
        self.F = 0.85 # absorption fraction
        self.E = 0.85 # Emissivity
        self.K = 2.50e-2 # soil conductivity
        self.C = 3300 # specific heat
        self.rho = 1600 # density
        self.average_temperature = 150
        
if __name__=="__main__":
    mars = Mars()
    print mars