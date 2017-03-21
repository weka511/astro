'''This module looks after orbital calculations. It is responsible for
determining the distance of a planet from the Sun at a particular time.
'''
import re
__version__ = re.sub('[^0-9]','','$LastChangedRevision: 1119 $')

# $HeadURL: https://server/svn/sandbox/trunk/climate/kepler.py $
# $LastChangedDate: 2016-11-27 17:40:14 +1300 (Sun, 27 Nov 2016) $
# $LastChangedRevision: 1119 $

# Copyright (C) 2016 Greenweaves Software Pty Ltd

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

#
# References:
# MD Murray & Dermott, Solar System Dynamics
# GA Goosse H., P.Y. Barriat, W. Lefebvre, M.F. Loutre and V. Zunz, (2008-2010).
#    Introduction to climate dynamics and climate modeling.
#    Online textbook available at http://www.climate.be/textbook.


import math,physics


def get_mean_anomaly(n,t,tau=0):
    '''Calculate True Anomaly using Murray & Dermott (2.39)
    Argumants: n     Mean Motion
               t     Time
               tau   Time of pericentre passage
    '''
    return n*(t-tau)

def get_true_anomaly(E,eccentricity):
    '''Calculate True Anomaly using Murray & Dermott (2.46)
    Arguments: E             Eccentric anomaly
               eccentricity  Eccentricty of orbit
    '''
    def reflect(theta):
        return 2*math.pi-theta
    def get_true_anomaly_upper_hemisphere(E):
        return 2.0*math.atan( math.sqrt((1+eccentricity)/(1-eccentricity))*math.tan(0.5*E) )
    if E<math.pi:
        return get_true_anomaly_upper_hemisphere(E)
    else:
        return reflect(get_true_anomaly_upper_hemisphere(reflect(E)))

def get_eccentric_anomaly(M,eccentricity,tolerance=1.0e-9,k=0.85):
    '''Calculate Eccentric Anomaly using Murray & Dermott (2.52).
       Solve using Newton Raphson
       Arguments: M              Mean anomaly
                  eccentricity   Eccentricity of orbit
                  tolerance      Maximum error allowed in Newton Raphson step
                  k              Correction for starting value Murray & Dermott 2.64
    '''
    E = M + math.copysign(k*eccentricity,math.sin(M)) # MD 2.64
    delta = 2*tolerance*E
    while abs(delta)>tolerance*E:
        delta = (E-eccentricity*math.sin(E)-M)/(1-eccentricity*math.cos(E)) # MD (2.58
        E-=delta
    return E


def get_distance_from_focus(true_anomaly,semi_major_axis,eccentricity=0.0167):
    '''Calculate distance from focus using Goosse et al 2.15 and Murray & Dermott 2.19
    
    Arguments: true_anomaly
               semi_major_axis
               eccentricity
    '''
    return semi_major_axis*(1-eccentricity*eccentricity)/(1+eccentricity*math.cos(true_anomaly))

def get_mean_distance_from_focus(semi_major_axis,eccentricity=0.0167):
    '''Calculate mean distance of planet from the focus that is is orbiting
    around using Goosse et al 2.16
    
    Arguments: semi_major_axis
               eccentricity
    '''
    return semi_major_axis*math.sqrt(1-eccentricity*eccentricity)

    
def true_longitude_from_true_anomaly(true_anomaly,PERH=102.04):
    '''
    Calculate True Longitude using Goosse et al equation 2.19
    http://www.climate.be/textbook/chapter2_node1.html
    
    Argumenbts: true_anomaly
                PERH
    '''
    def clip_angle(angle):
        while angle>2*math.pi:
            angle-=2*math.pi
        while angle<0:
            angle+=2*math.pi        
        return angle    
    return clip_angle(math.radians(180+PERH)+true_anomaly)

def true_anomaly_from_true_longitude(true_longitude,PERH=102.04):
    '''
    Calculate True Anomaly
    Arguments: true_longitude
               PERH
    '''
    return true_longitude-math.radians(180+PERH)    # FIXME


if __name__=='__main__':
    import pylab
    semi_major_axis = physics.Conversion.au2meters(39.2851)  #Pluto
    eccentricity = 0.246682  #Pluto    
    
    figure = 1
    pylab.figure(figure,figsize=(10,10))
    xs = []
    ys = []
    areas = []

    r0 = -1
    nu = -1
    for i in range(360):
        M = get_mean_anomaly(1,math.radians(i))
        E = get_eccentric_anomaly(M,eccentricity)
        nu = get_true_anomaly(E,eccentricity)
        r = get_distance_from_focus(nu,semi_major_axis,eccentricity)
        xs.append(r*math.cos(nu))
        ys.append(r*math.sin(nu))
        if r0>0:
            areas.append(0.5*math.sin(nu-nu0)*r*r0)
        r0 = r
        nu0 = nu
    print ('Relative difference between maximum and minimum area: {0:.2g}'.format(2.0*(max(areas)-min(areas))/(min(areas)+max(areas))))
    pylab.plot(xs,ys,'b.')
    pylab.plot([0],[0],'rs')
    pylab.title(r'$\mathrm{Show\ that}\ \mathit{get\_distance\_from\_focus}\ \mathrm{produces\ an\ ellipse}$')
    lim0=min(min(xs),min(ys))
    lim1=max(max(xs),max(ys))         
    pylab.xlim(lim0,lim1)
    pylab.ylim(lim0,lim1)
    pylab.show()
    