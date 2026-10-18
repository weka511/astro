#!/bin/bash

N=10000
SHOW=""
step=0.01
while getopts "N:s:w" opt; do
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
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  
    esac
  done
  
./restricted.py $SHOW --mu2 0.000953875 --N $N --dx -0.97668 0 --dv 0 -0.06118 --name fig_3_17a --step $step
