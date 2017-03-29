# Copyright (C) 2016-2017 Greenweaves Software Pty Ltd

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
'''

import fnmatch,os,tkinter as tk,viewer,matplotlib.pyplot as plt,multiprocessing as mp
from tkinter import *

def exec_display_ext(inputfile,figure):
    print(inputfile,figure)
    viewer.display(inputfile,figure=figure)         
    plt.show()
    
def exec_display_maxmin_ext(inputfile,figure):
    print(inputfile,figure)
    viewer.display_maxmin(inputfile,figure=figure)         
    plt.show() 
    
def exec_display_daily_minima_ext(inputfile,figure):
    print(inputfile,figure)
    viewer.display_daily_minima_all_latitudes(figure=figure)         
    plt.show()     
    
class Viewer(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()
        self.figure=1
        
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

        self.QUIT = tk.Button(self, text='QUIT', fg="red", command=root.destroy)
        self.QUIT.grid(row=2,column=1) 

    def exec_all(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        print (inputfile)
        job_for_another_core = mp.Process(target=exec_display_ext,args=(inputfile,self.figure))
        job_for_another_core.start()  
        self.figure+=1
 
        
    def exec_maxmin(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        print (inputfile)
        job_for_another_core = mp.Process(target=exec_display_maxmin_ext,args=(inputfile,self.figure))
        job_for_another_core.start() 
        self.figure+=1
        
    def exec_daily_minima(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        print (inputfile)
        job_for_another_core = mp.Process(target=exec_display_daily_minima_ext,args=(inputfile,self.figure))
        job_for_another_core.start() 
        self.figure+=1  
        

if __name__=='__main__':   
    mp.freeze_support()
    root = tk.Tk()
    app = Viewer(master=root)
    app.mainloop()

