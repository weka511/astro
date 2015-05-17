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

# Repository for basic data about planets

import math

# Used to convert between CGS and SI
class Conversion:
    cm_per_metre = 100
    gm_per_Kg = 1000
    cm3_per_meter3 = cm_per_metre*cm_per_metre*cm_per_metre
    
class Planet:     
    def __init__(self,name):
        self.name = name
        self.S = 1371   # Solar constant at the mean Sun-Earth distance of l AU, in N/m2
                        # Appelbaum & Flood
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
               "conductivity={7:6.2g}\n" + \
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
#   Appelbaum & Flood equations (2) & (3)
    def instantaneous_distance(self,areocentric_longitude):
        theta = areocentric_longitude - 248 # True anomaly
        return (self.a*(1-self.e*self.e)/
                (1 + self.e * math.cos(math.radians(theta))))

    #   Sine of declination
    #   Appelbaum & Flood equation (7)
    def sin_declination(self,areocentric_longitude):
        return math.sin(math.radians(self.obliquity)) * \
               math.sin(math.radians(areocentric_longitude))
     
    #   Cosine of zenith angle
    #   Appelbaum & Flood equation (6)        
    def cos_zenith_angle(self,areocentric_longitude,latitude,T):
        sin_declination=self.sin_declination(areocentric_longitude)
        cos_declination=math.sqrt(1-sin_declination*sin_declination)
        return math.sin(math.radians(latitude))*sin_declination +            \
            math.cos(math.radians(latitude))*cos_declination*math.cos(math.radians(self.hour_angle(T)))

    #   Hour angle
    #   Appelbaum & Flood equation (8) 
    def hour_angle(self,T):
        return 15*T-180

    def get_days_in_year(self):
        return 365.256363004*math.sqrt(self.a*self.a*self.a)
         
    def get_areocentric_longitude(self,day,hour):
        return 360*float(day)/get_days_in_year()
    
class Mars(Planet):
    def __init__(self):
        Planet.__init__(self,"Mars")
        self.a = 1.523679  # Wikipedia Mars page
        self.e = 0.093377  # Appelbaum & Flood
        self.obliquity = 24.936 # Appelbaum & Flood
        self.hours_in_day = 24 # should be 24.6597 http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
        self.F = 0.85 # absorption fraction - Leighton & Murray
        self.E = 0.85 # Emissivity - Leighton & Murray
        self.K = 6e-5 * Conversion.cm_per_metre # soil conductivity - Leighton & Murray
        self.C = 3.3 * Conversion.gm_per_Kg # specific heat
        self.rho = 1.6 * Conversion.cm3_per_meter3 / Conversion.gm_per_Kg # density
        self.average_temperature = 210 #http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html

        
if __name__=="__main__":
    import unittest
    
    class TestMarsMethods(unittest.TestCase):
        def setUp(self):
            self.mars = Mars()
        def test_get_days_in_year(self):
            self.assertAlmostEqual(687,self.mars.get_days_in_year(),places=1)
            
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise