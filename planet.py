# Copyright (C) 2015-2017 Greenweaves Software Pty Ltd

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

'''Repository for basic data about planets'''

import math, physics


    
class Planet: 
    '''Store information about a planet
    
    Attributes:
                name
                a             Semimajor axis
                e             eccentricity
                obliquity
                hours_in_day
                F                   absorption
                E                   emissivity
                K                   conductivity
                C                   specific heat
                rho                 Density
                average_temperature
    '''    
    def __init__(self,name):
        '''Create planet and initialize'''
        self.name = name
        self.S = 1371   # Solar constant at the mean Sun-Earth distance of l AU, in N/m2
                        # Appelbaum & Flood

    def __str__(self):
        '''Convert planet to string'''
        return ('{0}\n'
                'Semimajor axis      = {1:9.7f} AU\n'
                'eccentricity        = {2:8.6f}\n' 
                'obliquity           = {3:6.3f}\n' 
                'hours in day        = {4:6.4f}\n' 
                'absorption          = {5:4.2f}\n'
                'emissivity          = {6:4.2f}\n' 
                'conductivity        = {7:5.2g}\n'
                'specific heat       = {8:6.1f}\n'
                'rho                 = {9:6.1f}\n' 
                'average temperature = {10:5.1f}'
                ).format(
                    self.name,
                    self.a,
                    self.e,
                    self.obliquity,
                    self.hours_in_day,
                    self.F,
                    self.E,
                    self.K,
                    self.C,
                    self.rho,
                    self.average_temperature
        )
  
    def instantaneous_distance(self,areocentric_longitude):
        '''Instantaneaous Distance from Sun in AU
        Appelbaum & Flood equations (2) & (3)
        '''
        theta = areocentric_longitude - 248 # True anomaly
        return (self.a*(1-self.e*self.e)/
                (1 + self.e * math.cos(math.radians(theta))))

    def sin_declination(self,areocentric_longitude):
        '''Sine of declination
        Appelbaum & Flood equation (7)'''
        return math.sin(math.radians(self.obliquity)) * \
               math.sin(math.radians(areocentric_longitude))
     
    def cos_zenith_angle(self,areocentric_longitude,latitude,T):
        '''Cosine of zenith angle
        Appelbaum & Flood equation (6)
        '''
        sin_declination=self.sin_declination(areocentric_longitude)
        cos_declination=math.sqrt(1-sin_declination*sin_declination)
        return math.sin(math.radians(latitude))*sin_declination +            \
            math.cos(math.radians(latitude))*cos_declination *               \
            math.cos(math.radians(self.hour_angle(T)))

    def hour_angle(self,T):
        '''Hour angle
        Appelbaum & Flood equation (8)
        '''
        return 15*T-180

    def get_days_in_year(self):
        return Earth.get().earth.get_days_in_year()*math.sqrt(self.a*self.a*self.a)
         
    def get_areocentric_longitude(self,day,hour):
        return 360*float(day)/self.get_days_in_year()

class Mercury(Planet):
    '''Data for the planet Mercury'''
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Mercury')
        self.a = 0.387098  # Wikipedia Mercury page
        self.e = 0.205630  # Wikipedia Mercury page
        self.obliquity = 0 # Wikipedia Axial Tilt page
        
class Venus(Planet):
    '''Data for the planet Venus'''
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Venus')
        self.a = 0.723327  # Wikipedia Venus page
        self.e = 0.0067  # Wikipedia Venus page
        self.obliquity = 177.36 # Wikipedia Axial Tilt pagee
        
class Earth(Planet):
    '''Data for the planet Earth'''
    earth = None
    @classmethod
    def get(cls):
        if Earth.earth==None:
            Earth.earth=Earth()
        return Earth.earth
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Earth')
        self.a = 1.0 # the  semimajor axis in AU,
        self.e = 0.017 #  eccentricity
        self.obliquity = 23.4 # Wikipedia Axial Tilt page
        self.hours_in_day = 24
        self.average_temperature = 300        
    def get_days_in_year(self):
        return 365.256363004
    
class Mars(Planet):
    '''Data for the planet Mars'''
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Mars')
        self.a = 1.523679  # Wikipedia Mars page
        self.e = 0.093377  # Appelbaum & Flood
        self.obliquity = 24.936 # Appelbaum & Flood
        self.hours_in_day = 24 # should be 24.6597 http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
        self.F = 0.85 # absorption fraction - Leighton & Murray
        self.E = 0.85 # Emissivity - Leighton & Murray
        self.K = 6e-5 * physics.Conversion.cm_per_metre # soil conductivity - Leighton & Murray
        self.C = 3.3 * physics.Conversion.gm_per_Kg # specific heat
        self.rho = 1.6 * physics.Conversion.cm3_per_meter3 / physics.Conversion.gm_per_Kg # density
        self.average_temperature = 210 #http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
        
class Jupiter(Planet):
    '''Data for the planet Jupiter'''
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Jupiter')
        self.a = 5.204267  # Wikipedia Jupiter page
        self.e = 0.048775  # Wikipedia Jupiter page
        self.obliquity = 3.13 # Wikipedia Axial Tilt page

class Saturn(Planet):
    '''Data for the planet Saturn'''
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Saturn')
        self.a = 9.5820172   # Wikipedia Saturn page
        self.e = 0.055723219  # Wikipedia Saturn page
        self.obliquity = 26.73 # Wikipedia Axial Tilt page
        
class Uranus(Planet):
    '''Data for the planet Uranus'''
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Uranus')
        self.a = 19.189253   # Wikipedia Jupiter page
        self.e = 0.047220087  # Wikipedia Uranus page
        self.obliquity = 97.77 # Wikipedia Axial Tilt page
        
class Neptune(Planet):
    '''Data for the planet Neptune'''
    def __init__(self):
        '''Create data for planet'''
        Planet.__init__(self,'Neptune')
        self.a = 30.070900  # Wikipedia Jupiter page
        self.e = 0.00867797  # Wikipedia Neptune page
        self.obliquity = 28.32 # Wikipedia Axial Tilt page

def create(name):
    '''Create a named Planet'''
    planets=[Mercury(),
             Venus(),
             Earth(),
             Mars(),
             Jupiter(),
             Saturn(),
             Uranus(),
             Neptune()
    ]
    for planet in planets:
        if planet.name.upper()==name.upper(): return planet
    return None

if __name__=='__main__':
    import unittest
    
    class TestMarsMethods(unittest.TestCase):
        def setUp(self):
            self.mars = create('mars')
            print(self.mars)
        def test_get_days_in_year(self):
            self.assertAlmostEqual(687,self.mars.get_days_in_year(),places=1)
            
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise