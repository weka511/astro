#!/bin/bash

#   Copyright (C) 2026 Simon Crase

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
 
# Generate figure 3.17 (b) fro Murray and Dermott
 
N=10000
SHOW=""
step=0.01
integrator=rki2
dx="-1.02745 0"
dv="0 0.04032"
mu2=0.000953875
name="none"

while getopts "N:s:wi:e:" opt; do
  case ${opt} in
    N)
      N="$OPTARG"
      ;;
	w)
      SHOW='--show'
      ;;
	s)
      step="$OPTARG"
      ;;
	i)
      integrator="$OPTARG"
	  ;;
	e)
	  name=$OPTARG
	  case $OPTARG in
	     3_16_a)
			dx="0.0065 0.0065"
			dv="0 0"
			mu2=0.001
			;;
		3_16_b)
			dx="0.008 0.008"
			dv="0 0"
			mu2=0.001
			;;
        3_17_a)
			dx="-0.97668 0"
			dv="0 -0.06118"
			mu2=0.000953875
			;;
		3_17_b)
			dx="-1.02745 0"
			dv="0 0.04032"
			mu2=0.000953875
			;; 
	  esac
	  ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  
    esac
  done
  
./restricted.py $SHOW  --mu2 $mu2 --N $N --dx $dx --dv $dv --name fig_3_17b --step $step --integrator $integrator --L5 --norepeat
