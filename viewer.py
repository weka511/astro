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

import utilities, io, os,matplotlib.pyplot as plt, sys, getopt

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
            count=0
            xs,ys=history.extract(1)
            x_previous=-1
            ys_for_period=[]
            xxx=[]
            ymin=[]
            ymax=[]
            for x,y in zip(xs,ys):
                  if x_previous<0:
                        x_previous=x
                        ys_for_period=[]
                  if x-x_previous>=1:
                        xxx.append(x_previous)
                        ymin.append(min(ys_for_period))
                        ymax.append(max(ys_for_period))
                        x_previous=x
                        ys_for_period=[]
                  ys_for_period.append(y)
            xxx.append(x_previous)
            ymax.append(max(ys_for_period))                  
            ymin.append(min(ys_for_period))
            plt.figure(figure)
            plt.title("Diurnal variation in temperature for {0}".format(inputfile))
            plt.xlabel("Time")
            plt.ylabel("Temperature - Kelvin")
            plt.grid(True)            
            plt.plot(xxx,ymin,'b-',xxx,ymax,'r-')
            plt.savefig(os.path.splitext(inputfile)[0]+"-minmax")
            
def main(argv):
      inputfile='output.txt'    
      figure=1
      all = False
      minmax= False
      
      if len(argv)>0:
            try:
                  opts, args = getopt.getopt( \
                        argv,\
                        "hi:am",\
                        ["help","ifile=","allpoints","maxmin"])
                        
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
                        print "-a"
                        all=True
                  elif opt == "-m":
                        print "-m"
                        minmax=True
      
      if all:
            display(inputfile,figure)
            figure += 1
      if minmax:
            display_maxmin(inputfile,figure)
            figure +=1
            
      plt.show()

if __name__=="__main__":
      main(sys.argv[1:])
      