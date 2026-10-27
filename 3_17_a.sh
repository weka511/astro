#!/bin/bash

#  Copyright (C) 2026 Simon Crase

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

while getopts "N:s:wi:" opt; do
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
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  
    esac
  done
  
./restricted.py $SHOW --mu2 0.000953875 --N $N --dx -0.97668 0 --dv 0 -0.06118 --name fig_3_17a --step $step --integrator $integrator ==L5
