# Copyright (C) 2017 Greenweaves Software Pty Ltd

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

# snarfed from http://stackoverflow.com/questions/19086106/how-to-utilize-all-cores-with-python-multiprocessing

'''
Execute leighton.py on multiple cores.

Do not try to run this under the IDE, as it fights with multiprocessing!
'''
import multiprocessing as mp, leighton, time

def init_worker(mps, fps, cut):
    global memorizedPaths, filepaths, cutoff
    global DG

    print ('process initializing', mp.current_process())
    memorizedPaths, filepaths, cutoff = mps, fps, cut
    DG = 1##nx.read_gml('KeggComplete.gml', relabel = True)

def work(item):
    print (item)
    leighton.main(item.split())



if __name__ == '__main__':
    start_time = time.time()    
    m = mp.Manager()
    memorizedPaths = m.dict()
    filepaths = m.dict()
    cutoff = 1 ##
    # use all available CPUs
    p = mp.Pool(initializer=init_worker, 
                initargs=(memorizedPaths,
                          filepaths,
                          cutoff))
    parameters_list=[
        '-f 0 -t 2880 -s 1 -l -80 -c',
         '-f 0 -t 2880 -s 1 -l -70 -c',
         '-f 0 -t 2880 -s 1 -l -50 -c',
         '-f 0 -t 2880 -s 1 -l -30 -c',
         '-f 0 -t 2880 -s 1 -l -10 -c',
         '-f 0 -t 2880 -s 1 -l 0 -c',
         '-f 0 -t 2880 -s 1 -l 10 -c',
         '-f 0 -t 2880 -s 1 -l 30  -c',
         '-f 0 -t 2880 -s 1 -l 50 -c',
         '-f 0 -t 2880 -s 1 -l 70 -c',
         '-f 0 -t 2880 -s 1 -l 90 -c'       
    ] 
    for _ in p.imap_unordered(work, parameters_list, chunksize=1):
        pass
    p.close()
    p.join()
    print('--- %s seconds ---' % (time.time() - start_time))