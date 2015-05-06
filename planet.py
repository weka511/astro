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
    def __init__(self):
        self.S = 1371 # Solar constant at the mean Sun-Earth distance of l AU, in N/m2
        self.a = 1.0 # the  semimajor axis in AU,
        self.e = 0.017 #  eccentricity
        self.obliquity = 23.4
        self.hours_in_day = 24
        
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
