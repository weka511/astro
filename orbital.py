# Copyright (C) 2019 Greenweaves Software Limited

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

from numpy import matmul,sign,arctan2,clip
from math import cos,sin,radians,isclose,degrees,pi,sqrt,floor,modf

# compose
#
# Rotate coordinates in 3D: MS 2.119, 2.120, & 2.121
#
# Parameters:
#    omega  Argument of pericentre
#    I      Inclination       
#    Omega  Longitude of ascending node

def compose(omega=0,I=0,Omega=0):
    cos_omega = cos(omega)
    sin_omega = sin(omega)
    cos_I     = cos(I)
    sin_I     = sin(I)
    cos_Omega = cos(Omega)
    sin_Omega = sin(Omega)
    
    P1        = [[cos_omega, -sin_omega, 0],\
                 [sin_omega,  cos_omega, 0],\
                 [0,          0,         1]]
   
    P2        = [[1,          0,      0],     \
                 [0,          cos_I, -sin_I], \
                 [0,          sin_I,  cos_I]]
    
    P3        =  [[cos_Omega, -sin_Omega, 0],\
                  [sin_Omega,  cos_Omega, 0],\
                  [0,          0,         1]]
   
    return matmul(P3, matmul(P2, P1))

# get_XYZ
#
# Calculate position of particle in 3D: MD 2.122
#
# Parameters:
#    omega  Argument of pericentre
#    I      Inclination       
#    Omega  Longitude of ascending node
#    f      True anomaly
#    r      radius
def get_XYZ(omega=0,I=0,Omega=0,f=0,r=1):
    cos_Omega   = cos(Omega)
    sin_Omega   = sin(Omega)
    cos_I       = cos(I)
    sin_I       = sin(I)
    cos_omega_f = cos(omega + f)
    sin_omega_f = sin(omega + f)
    return (r * (cos_Omega*cos_omega_f - sin_Omega*sin_omega_f*cos_I),
            r * (sin_Omega*cos_omega_f + cos_Omega*sin_omega_f*cos_I),
            r * sin_omega_f*sin_I)

# get_lambda
#
# Calculate mean longitude: MD A.19
#
# Parameters:
#       T          Number of Julian centuries since start of epoch
#       lambda0    mean longitude at T = 0 in degrees
#       lambda_dot Rate of change of mean longitude in arc seconds per century 
#       Nr
#
# Returns: mean longitude in radians

def get_mean_longitude(T=0,lambda0=34.40438,lambda_dot=557078.35,Nr=8):
    return radians(lambda0 + (lambda_dot/3600 + 360 * Nr) *T)

# get_eccentricity
#
# Calculate eccentricity: MD A.15
#
# Parameters:
#     T        Number of Julian centuries since start of epoch
#     e0       Eccentricity at t=0
#     e_dot
#
def get_eccentricity(T=0,e0=0.04839266,e_dot=-12280/1e8):
    return e0 + e_dot * T

#  get_semimajor_axis
#
# Calculate semi-major axis: MD A.14
#
# Parameters:
#     T        Number of Julian centuries since start of epoch
#     a0       Semi-major axis at start of epoch
#     a_dot

def get_semimajor_axis(T=0,a0=5.20336301,a_dot=60377/1e8):
    return a0 + a_dot * T

#  get_varpi
#
# Calculate longitude of perihelion:  MD A.17
#
# Parameters:
#     T          Number of Julian centuries since start of epoch
#     varpi0     Longitude at start of epoch
#     varpi_dot

def get_longitude_of_perihelion(T=0,varpi0=14.75385,varpi_dot=839.93):
    return radians(varpi0 + varpi_dot * T / 3600)

# get_true_anomaly
#
# Solve equation 2.43 by the Newton-Raphson method to estimate f.
#
# Parameters: E         Eccentric anomally
#             e         Eccentricity
#             N         Maximum number of iterations
#             tolerance Used to assess correction: throws Assetion Error 
#                       if correction still exceeds tolerance after N iterations

def get_true_anomaly(E,e,N=1000,tolerance=1e-12):
    cos_f = (cos(E)-e)/(1-e*cos(E))
    f     = E  # E and F should be in same quadrant, so use E as starting value
    for i in range(N):
        if abs(f) < tolerance: return f
        correction = (cos(f)-cos_f)/sin(f)
        f += correction
        if (abs(correction)<tolerance):
            return f
        
    assert False,'Correction {0} is still greater than tolerance {1} after {2} iterations.'.format(correction,tolerance,N)        
        
# get_xy
#
# Compute 2D position using 2.41
#
#    Parameters:
#        T               Number of Julian centuries since start of epoch
#        lambdaT         Mean longitude in radians
#        e               eccentricity
#        a               Semi-major axis
#        varpi           longitude of perihelion in radians

def get_xy(T=0,lambdaT=0,e=0.0484007,a=5.20332,varpi = 14.7392):
    M = lambdaT - varpi #MD: paragraph before 2.123
    E = get_eccentric_anomaly(e=e,M=M)
    return a * (cos(E) - e),a * sqrt(1 - e*e)*sin(E)

# get_eccentric_anomaly
#
# Calculate eccentric anomaly by solveing kepler's equation - 2.52
#
# Parameters:
#             e             Eccentricty of orbit
#             M  Mean anomaly
#             tolerance     Used to assess correction: throws Assetion Error 
#                           if correction still exceeds tolerance after N iterations
#             N             Maximum number of iterations
#             k             Paramter used in Danby's starting value - MD 2.64
def get_eccentric_anomaly(e=0,M=0,rel_tol=1e-7, abs_tol=1e-9,N=10000,k=0.85):
   
    E = M + sign(sin(M) * k * e)  # Danby's starting value - MD 2.64
    
    for i in range(N):
        correction = (E - e * sin(E) - M)/(1 - e * cos(E))
        if isclose(correction,0,rel_tol=rel_tol,abs_tol=abs_tol): return E
        E -= correction
        
    assert False, 'Correction {0} is still greater than tolerance {1} after {2} iterations.'.format(correction,
        max(rel_tol * abs(correction), abs_tol),N) 
    
# Times
#
# Generator used for iterating over times
#
#         From        Starting time (Julian days)
#         To          End time (Julian days)
#         Incr        Interval from one sample to the next

def Times(From=0,To=10,Incr=1,Epoch=2451545.0):
    t = From
    while t < To +Incr:
        yield t,(t-Epoch)/36525
        t+=Incr

# create_orbit
#
# Calculate positions in orbit
#
#     Parameters:
#         planet      Elements for planet
#         lambda_dot  Used to calculate mean longitude
#         Nr          Used to calculate mean longitude
#         From        Starting time (Julian centuries)
#         To          End time (Julian centuries)
#         Incr        Interval from one sample to the next
#         is2D        Used to force a 2D calculcation

def create_orbit(planet,
                 lambda_dot = 1293740.63,
                 Nr         = 99,
                 From       = 0,
                 To         = 10,
                 Incr       = 1,
                 is2D       = False):
    a,e,I,varpi,Omega,lambda0 = planet
    if is2D: I                = 0
    Xs                        = []
    Ys                        = []
    Zs                        = []
    ts                        = []
    Rotation                  = compose(omega = radians(varpi - Omega), #MD 2.118
                                         I     = radians(I),
                                         Omega = radians(Omega))

    for t,T in Times(From=From,To=To,Incr=Incr):
        x,y = get_xy(T       = T, 
                       lambdaT = get_mean_longitude(T,
                                                    lambda0    = lambda0,
                                                    lambda_dot = lambda_dot,
                                                    Nr         = Nr), 
                       e       = e,
                       a       = a,
                       varpi   = radians(varpi))

        W   = matmul(Rotation,[[x],[y],[0]])
        Xs.append(W[0][0])
        Ys.append(W[1][0])
        Zs.append(W[2][0])
        ts.append(t)

    return (Xs,Ys,Zs,ts)

# is_minimum
#
# Verify that value if a minimum
#
#    Parameters:
#       a
#       b
#       c
#
# Returns: True iff b is less than both a and c

def is_minimum(a,b,c):
    return a>b and b < c

# get_distance
#
# Get Euclidean distance between two points

def get_distance(x0,y0,z0,x1,y1,z1):
    return sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1))

# get_julian_date
#
# Calculate Julian Date MD A.1, A.2 and A-3.
#
# Parameters: Y  Year (NB there is no year 0 - we got from 10 to 1
#             M  Month:  1-12
#             D  Day:    1-31    
#             UT Universal Time:  
def get_julian_date(Y,M,D,UT=12):
    # is_gregorian
    #
    # Verify that date is within Gregorian Era
    def is_gregorian():
        if Y<1582: return False
        if Y>1582: return True
        assert Y==1582
        if M<10: return False
        if M>10: return True
        assert M==10
        if D <= 4: return False
        if D >= 15: return True
        assert False,'There is no such date as {0}-{1}-{2}'.format(Y,M,D)
    y,m = (Y-1,M+12) if M<=2 else (Y,M)
    B   = floor(y/400) - floor(y/100) if is_gregorian() else -2
    return floor(365.25 * y) + floor(30.6001*(m+1)) + B + 1720996.5 + D +UT/24

# get_calendar_date
#
# Convert Julian Date to Calendar Date, MD A.4-A.13

def get_calendar_date(JD):
    frac_JD,a = modf(JD + 0.5)
    c         = a + 1524
    if a>=2299161:
        b     = floor((a-1867216.25)/36524.25)
        c     = a + b-floor(b/4)+1525
    d         = floor((c-122.1)/365.25)
    e         = floor(365.25*d)
    f         = floor((c-e)/30.6001)
    D         = c - e -floor(30.6001*f) + frac_JD
    M         = f - 1 - 12 * floor(f/14)
    Y         = d -4715-floor((7+M)/10)
    return (Y,M,D)

if __name__=='__main__':
    import unittest
    
    class TestJulian(unittest.TestCase):
        
        def test_toJulian(self):
            self.assertAlmostEqual(2431855.933,get_julian_date(1946,2,4,10.4),places=3)
            
        def test_toJulian0(self):
            self.assertAlmostEqual(0,get_julian_date(-4712,1,1,12),places=3) #noon at start of 4713 BC
            
        def test_toCalendar(self):
            Y,M,D = get_calendar_date(2434903.75)
            self.assertEqual(1954,Y)
            self.assertEqual(6,M)
            self.assertEqual(10.25,D)
            
        def test_toCalendar0(self):
            Y,M,D = get_calendar_date(0)
            self.assertEqual(-4712,Y)
            self.assertEqual(1,M)
            self.assertEqual(1.5,D)        
            
    class TestJupiter(unittest.TestCase):
    
        # test_mult
        #
        # This test uses Murray & Dermott (2.124)
        def test_compose(self):
            rotation = compose(omega=radians(14.7392 -100.535 ),I=radians(1.30537),Omega=radians(100.535))
            self.assertAlmostEqual(0.966839,   rotation[0][0], places=5)
            self.assertAlmostEqual(-0.254401,  rotation[0][1], places=6)
            self.assertAlmostEqual(0.0223971,  rotation[0][2], places=6)
            self.assertAlmostEqual(0.254373,   rotation[1][0], places=5)
            self.assertAlmostEqual(0.967097,   rotation[1][1], places=6)
            self.assertAlmostEqual(0.00416519, rotation[1][2], places=6)
            self.assertAlmostEqual(-0.0227198, rotation[2][0], places=6)
            self.assertAlmostEqual(0.00167014, rotation[2][1], places=7)
            self.assertAlmostEqual(0.99974,    rotation[2][2], places=6)
        
        # test_mult
        #
        # This test uses Murray & Dermott (2.124)            
        def test_get_XYZ(self):
            x     = -5.39027
            y     = -0.818277
            r     = sqrt(x*x + y * y)
            f     = arctan2(y,x)
            X,Y,Z = get_XYZ(omega=radians(14.7392 -100.535 ),I=radians(1.30537),Omega=radians(100.535),r=r,f=f)
            self.assertAlmostEqual(-5.00336, X, places=5)
            self.assertAlmostEqual(-2.16249, Y, places=5)
            self.assertAlmostEqual(0.121099, Z, places=6)

        def test_get_T(self):
            self.assertAlmostEqual(-0.06266423,(2449256.189-2451545.0)/36525) # check Julian centuries
            
        def test_get_mean_longitude(self):
            self.assertAlmostEqual(radians(204.234),get_mean_longitude(T=-0.06266423)%(2*pi),places=3)
            
        def test_get_xy(self):
            x,y = get_xy(T            = -0.06266423,
                         lambdaT      = get_mean_longitude(T=-0.06266423),
                         e            = get_eccentricity(T=-0.06266423),
                         a            = get_semimajor_axis(-0.06266423),
                         varpi        = get_longitude_of_perihelion(-0.06266423))
            self.assertAlmostEqual(-5.39027,x,places=5)
            self.assertAlmostEqual(-0.818277,y,places=5)
            
        def test_eccentricty(self):
            self.assertAlmostEqual(0.0484007,
                                   get_eccentricity(T=-0.06266423),places=6)
            
        def test_semimajor_axis(self):
            self.assertAlmostEqual(5.20332,
                                   get_semimajor_axis(T=-0.06266423),places=4)
            
        def test_get_longitude_of_perihelion(self):
            self.assertAlmostEqual(radians(14.7392),
                                   get_longitude_of_perihelion(T=-0.06266423),
                                   places=4)        
               
    class TestKepler(unittest.TestCase):
        # test_kepler_inverse
        #
        # Verify solution to Kepler's equation for a range of parameters
        def test_kepler_inverse(self):
            N = 25
            for i in range(N):
                mean_anomaly = 2 * pi /N
                E            = get_eccentric_anomaly(e=0.205635,M=mean_anomaly)
                M            = E - 0.205635 * sin(E)
                self.assertAlmostEqual(mean_anomaly,M)
            
        # test_kepler_jupiter
        #
        # See paragraph before MD 2.123
        def test_kepler_jupiter(self):
            self.assertAlmostEqual(radians(189.059),
                                   get_eccentric_anomaly(
                                       M=radians(189.495),
                                       e=0.0484007),
                                   places=4)
            
    class TestOrbit(unittest.TestCase):
        
        def test_perihelion(self): # see https://www.timeanddate.com/astronomy/perihelion-aphelion-solstice.html 
            Xs,Ys,Zs,Ts     = create_orbit(( 1.00000011, 0.01671022, 0.00005, 102.94719, 348.93936, 100.46435),
                                       lambda_dot = 1293740.63,
                                       Nr         = 99,
                                       From       = get_julian_date(2018,12,31),
                                       To         = get_julian_date(2019,12,31),
                                       Incr       = 1)
            
            solar_distances = [get_distance(Xs[i],Ys[i],Zs[i],0,0,0) for i in range(len(Ts))]
            
            minima          = [i for i in range(1,len(solar_distances)-1)
                               if is_minimum(solar_distances[i-1],
                                             solar_distances[i],
                                             solar_distances[i+1])]
            Y,M,D           = get_calendar_date(Ts[minima[0]])
            
            self.assertEqual(1,len(minima))
            self.assertEqual(2019,Y)
            self.assertEqual(1,M)
            self.assertEqual(3.5,D)
    
    unittest.main()