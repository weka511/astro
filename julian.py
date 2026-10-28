#!/usr/bin/env python

#   Copyright (C) 2019 Simon Crase

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''Methods for converting between Calendar dates and Julian dates'''

from unittest import TestCase,main
import numpy as np

__version__ = '1.0'
__author__ = 'Simon Crase'


class Calendar:
    '''
    This class contains methods for converting between Calendar dates and Julian dates.
    It is based on Appendex A.3 of Murray and Dermott
    '''
    @staticmethod
    def get_julian_date(Y,M,D,UT=12):
        '''
        Calculate Julian Date MD Appendex A
       
        Parameters: 
            Y  Year (NB there is no year 0 - we got from 10 to 1
            M  Month:  1-12
            D  Day:    1-31    
            UT Universal Time:  
        '''
        def is_gregorian():
            '''
            Verify that date is within Gregorian Era. The Julian calendar day Thursday,
            4 October 1582 was followed by the first day of the Gregorian calendar, 
            Friday, 15 October 1582 (Wikipedia).
            '''
            if Y < 1582: return False
            if Y > 1582: return True
            assert Y == 1582
            if M < 10: return False
            if M > 10: return True
            assert M == 10
            if D <= 4: return False
            if D >= 15: return True
            raise ValueError(f'There is no such date as {Y}-{M}-{D}')
            
        y,m = (Y-1,M+12) if M <= 2 else (Y,M)                                        # A.1
        B = np.floor(y/400) - np.floor(y/100) if is_gregorian() else -2                    # A.2
        return np.floor(365.25 * y) + np.floor(30.6001*(m+1)) + B + 1720996.5 + D + UT/24  # A.3
    
    @staticmethod
    def get_calendar_date(JD):
        '''
        Convert Julian Date to Calendar Date, MD Appendex A
        
        Parameters:
             JD     Julian Date
             
        Returns:
             Y,M,D
        '''
        frac_JD,a = np.modf(JD + 0.5)    # A.4
        
        # Deal with transition between calendars
        
        if a < 2299161:               # A.5
            c = a + 1524              # A.5
        else:
            b = np.floor((a-1867216.25)/36524.25) # A.6
            c = a + b - np.floor(b/4) + 1525      # A.7
            
        d = np.floor((c-122.1)/365.25)                # A.8
        e = np.floor(365.25*d)                        # A.9
        f = np.floor((c-e)/30.6001)                   # A.10
        D = c - e - np.floor(30.6001*f) + frac_JD     # A.11
        M = int(f - 1 - 12*np.floor(f/14))            # A.12
        Y = int(d - 4715 - np.floor((7+M)/10))        # A.13
        
        return (Y,M,D)
 
class TestJulian(TestCase):
    
    def test_toJulian(self):
        self.assertAlmostEqual(2431855.933,Calendar.get_julian_date(1946,2,4,10.4),places=3)
        
    def test_toJulian0(self):
        self.assertAlmostEqual(0,Calendar.get_julian_date(-4712,1,1,12),places=3) #noon at start of 4713 BC
        
    def test_toCalendar(self):
        Y,M,D = Calendar.get_calendar_date(2434903.75)
        self.assertEqual(1954,Y)
        self.assertEqual(6,M)
        self.assertEqual(10.25,D)
        
    def test_toCalendar0(self):
        Y,M,D = Calendar.get_calendar_date(0)
        self.assertEqual(-4712,Y)
        self.assertEqual(1,M)
        self.assertEqual(1.5,D)        
        
if __name__ == '__main__':
    main()

