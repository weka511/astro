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

import thermalmodel, planet, solar, utilities, sys,getopt, string, physics

def help():
      print 'leighton.py -o <outputfile> -f <hour from> -t <hour to> -l <latitude> -s <steps per hour>'
      
def main(argv):
      outputfile=''
      from_date=0
      to_date=720
      latitude = 0
      step = 10
      temperature = -1
      co2 = True
      spec=[(9,0.015),(10,0.3)]
      
      if len(argv)>0:
            try:
                  opts, args = getopt.getopt( \
                        argv,\
                        "ho:f:t:l:s:m:p:c",\
                        ["ofile=","from","to","latitude","step","temperature","co2","spec"])
            except getopt.GetoptError:
                  help()
                  sys.exit(2)
            for opt, arg in opts:
                  if opt == '-h':
                        help()
                        sys.exit()
                  elif opt in ("-o", "--ofile"):
                        outputfile = arg
                  elif opt in ("-f", "--from"):
                        from_date=int(arg)
                  elif opt in ("-t", "--to"):
                        to_date=int(arg)
                  elif opt in ("-l","--latitude"):
                        latitude=float(arg)
                  elif opt in ("-s","--step"):
                        step=int(arg)                  
                  elif opt in ("-m","--temperature"):
                        temperature=int(arg)
                  elif opt in ("-c","co2"):
                        co2=False
                  elif opt in ("-p","spec"):
                        spec=[]
                        for run in arg.strip('[]').split(';'):
                              try:
                                    couple=run.strip('()').split(',')
                                    spec.append((int(couple[0]),float(couple[1])))
                              except ValueError:
                                    print "Could not parse {0}".format(run)
                                    sys.exit(2)
                        
      if outputfile=='':
            (lat,ns)=utilities.format_latitude(latitude)
            outputfile="{1:3.0f}{0}.txt".format(lat,ns).strip()
           
      with open(outputfile, 'w') as f:
            mars = planet.Mars()
            solar_model = solar.Solar(mars)
            if temperature<0: temperature=mars.average_temperature
            
            history = utilities.ExternalTemperatureLog(f)
            history.write("Semimajor axes={0:10.7f} AU".format(mars.a))
            history.write("Eccentricty={0:10.7f}".format(mars.e))
            history.write("Obliquity={0:6.3f}".format(mars.obliquity))
            history.write("Hours in Day={0:10.7f}".format(mars.hours_in_day))
            history.write("Absorption Fraction={0:6.3f}".format(mars.F))
            history.write("Emissivity={0:6.3f}".format(mars.E))            
            history.write("Soil Conductivity={0:7.3f} W/M/K".format(mars.K))
            history.write("Specific Heat={0:7.3f} J/Kg/K".format(mars.C))
            history.write("Density={0:6.3f} Kg/M3".format(mars.rho))         
            history.write("Latitude={0:6.1f}".format(latitude))
            history.write("Step={0:6.1f}".format(step))
            history.write("Starting Temperature={0:6.1f} K".format(temperature))
            history.write("Layering (from top down)")
            for n,thickness in spec:
                  history.write("{0:d} layers, thickness {1:5.2f} metres each.".format(n,thickness))
            history.write("Albedo of snowcap = {0:5.2f}".format(physics.CO2.albedo))
            
            thermal=thermalmodel.ThermalModel(latitude,spec,solar_model,mars,history,temperature,co2)
            thermal.runModel(from_date,to_date,step)
     
if __name__ == "__main__":
      main(sys.argv[1:])