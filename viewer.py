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

import utilities, io, os,matplotlib.pyplot as plt, sys, getopt, glob

def help():
      print 'viewer.py -i <outputfile>'



def display(inputfile,figure=1):
      with open(inputfile, 'r') as f:
            history = utilities.ExternalTemperatureLog(f)
            plt.figure(figure)
            plt.title("Diurnal variation in temperature for {0}".format(inputfile))
            plt.xlabel("Time")
            plt.ylabel("Temperature - Kelvin")
            plt.grid(True)
            x,y1=history.extract(1)
            _,y2=history.extract(2) 
            _,y3=history.extract(3)
            _,y4=history.extract(4)
            _,y5=history.extract(5)
            _,y10=history.extract(10)
            _,y20=history.extract(20)
            plt.plot(x,y1,'r-',x,y2,'g-',x,y3,'b-',x,y4,'c-',x,y5,'m-',x,y10,'y-',x,y20,'k-')
            plt.savefig(os.path.splitext(inputfile)[0]+"-all")

def display_maxmin(inputfile,figure=1):
      with open(inputfile, 'r') as f:
            history = utilities.ExternalTemperatureLog(f)
            (xxx,ymin,ymax)=history.get_max_min(1)   
            plt.figure(figure)
            plt.title("Diurnal variation in temperature for {0}".format(inputfile))
            plt.xlabel("Time")
            plt.ylabel("Temperature - Kelvin")
            plt.grid(True)            
            plt.plot(xxx,ymin,'b-',xxx,ymax,'r-')
            plt.savefig(os.path.splitext(inputfile)[0]+"-minmax")


def display_daily_minima(inputfile,figure,colour):
      with open(inputfile, 'r') as f:
            history = utilities.ExternalTemperatureLog(f)
            latitude=history.get('latitude')
            (xxx,ymin,ymax)=history.get_max_min(1)
            plt.figure(figure)
            (NS,latitude)=utilities.format_latitude(latitude)
          
            plt.plot(xxx,ymin,colour,label="{0:5.1f}{1}".format(abs(latitude),NS)) 
            

def display_daily_minima_all_latitudes(figure):
      index=0;
      for name in glob.glob('*.txt'):
            display_daily_minima(name,figure,utilities.get_colour(index))
            index+=1
      plt.title("Minimum temperature for each Latitude")
      plt.xlabel("Time")
      plt.ylabel("Temperature - Kelvin")            
      plt.legend(loc='upper right')
      plt.grid(True)        
      plt.savefig("SurfaceTemperatureAllLatitudes")
      
def main(argv):
      inputfile='output.txt'    
      figure=1
      all = False
      minmax= False
      daily_minima=False
      if len(argv)>0:
            try:
                  opts, args = getopt.getopt( \
                        argv,\
                        "hi:amd",\
                        ["help","ifile=","allpoints","maxmin","daily"])
                        
            except getopt.GetoptError:
                  help()
                  sys.exit(2)
            for opt, arg in opts:
                  if opt == '-h':
                        help()
                        sys.exit()
                  elif opt in ("-i", "--ifile"):
                        inputfile = arg
                  elif opt == '-a':
                        all=True
                  elif opt == "-m":
                        minmax=True
                  elif opt == "-d":
                        daily_minima=True
      
      if all:
            display(inputfile,figure)
            figure += 1
      if minmax:
            display_maxmin(inputfile,figure)
            figure +=1
      if daily_minima:
            display_daily_minima_all_latitudes(figure)
            figure +=1
            
      plt.show()

if __name__=="__main__":
      main(sys.argv[1:])
      