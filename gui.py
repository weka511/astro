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

'''
User interface for Viewer

Do not try to run this under the IDE, as it fights with multiprocessing

Allows user to select text files from previous runs of leighton.py and
diplay them graphically.

'''

import fnmatch,os,tkinter as tk,viewer,matplotlib.pyplot as plt,multiprocessing as mp,re
from tkinter import *

# The next few methods need to be declared before multiprocessing calls
def exec_display_ext(inputfile,figure): 
    print('Display {0} as figure {1}'.format(inputfile,figure))
    viewer.display(inputfile,figure=figure)         
    plt.show()
    
def exec_display_maxmin_ext(inputfile,figure):
    print('Display Max Min{0} as figure {1}'.format(inputfile,figure))
    viewer.display_maxmin(inputfile,figure=figure)         
    plt.show() 
    
def exec_display_daily_minima_ext(inputfile,figure):
    print('Display Daily {0} as figure {1}'.format(inputfile,figure))
    m=re.match('(^[0-9]+-[0-9]+-[0-9]+-.*co2-)([0-9]+[NS]?)(.txt)',inputfile)  
    pattern='{0}[0-9]+([NS])?{1}'.format(m.group(1),m.group(3))
    print (m.group(1),m.group(2),m.group(3),pattern)
    viewer.display_daily_minima_all_latitudes(figure=figure,pattern=pattern)         
    plt.show()     

'''
This class represents to Frame and UI elements
'''
class Viewer(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()
        self.figure=1
        self.tasks=[]
        
    def createWidgets(self):
        scrollbar = tk.Scrollbar(root, orient="vertical")
        self.file_list=Listbox(self, width=50, yscrollcommand=scrollbar.set)
        scrollbar.config(command=self.file_list.yview)
        
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, '*.txt'):
                self.file_list.insert('end',file)        
 
        self.file_list.grid(row=0,column=1)
        self.all_button = tk.Button(self)
        self.all_button['text'] = 'All'
        self.all_button['command'] = self.exec_all
        self.all_button.grid(row=1,column=0) 
        
        self.minmax_button = tk.Button(self)
        self.minmax_button['text'] = 'MinMax'
        self.minmax_button['command'] = self.exec_maxmin
        self.minmax_button.grid(row=1,column=1) 
        
        self.daily_button = tk.Button(self)
        self.daily_button['text'] = 'Daily'
        self.daily_button['command'] = self.exec_daily_minima
        self.daily_button.grid(row=1,column=2)         

        self.QUIT = tk.Button(self, text='QUIT', command=self.close_all)
        self.QUIT.grid(row=2,column=1) 

    def close_all(self):
        for task in self.tasks:
            if task.is_alive():
                task.terminate()
                task.join()
        root.destroy()
        
    def exec_all(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        self.tasks.append(mp.Process(target=exec_display_ext,
                                     args=(inputfile,self.figure)))
        self.tasks[-1].start()  
        self.figure+=1
 
        
    def exec_maxmin(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        self.tasks.append(mp.Process(target=exec_display_maxmin_ext,
                                     args=(inputfile,self.figure)))
        self.tasks[-1].start() 
        self.figure+=1
        
    def exec_daily_minima(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        self.tasks.append( mp.Process(target=exec_display_daily_minima_ext,
                                      args=(inputfile,self.figure)))
        self.tasks[-1].start()
        self.figure+=1  
        

if __name__=='__main__':   
    mp.freeze_support() # Need this for MS Windows
    root = tk.Tk()
    app = Viewer(master=root)
    app.mainloop()

