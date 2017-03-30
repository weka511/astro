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
import multiprocessing as mp, leighton, time, getopt,sys

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
    input_file  = 'cores.dat'
    if len(sys.argv[1:])>0:
        try:
            opts, args = getopt.getopt(sys.argv[1:],'',[])
            for opt, arg in opts: 
                pass
            input_file = args[0]
            m = mp.Manager()
            memorizedPaths = m.dict()
            filepaths = m.dict()
            cutoff = 1 ##
            # use all available CPUs
            p = mp.Pool(initializer=init_worker, 
                        initargs=(memorizedPaths,
                                  filepaths,
                                  cutoff))
            print ('Opening {0}'.format(input_file))
            with open(input_file) as f:
                parameters_list = [x.strip() for x in  f.readlines() if len(x.strip())>0] 
                for _ in p.imap_unordered(work, parameters_list, chunksize=1):
                    pass
                p.close()
                p.join()
                print('--- %s seconds ---' % (time.time() - start_time))            
        except getopt.GetoptError:
            help()
            sys.exit(2)
   
