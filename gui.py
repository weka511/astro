import fnmatch,os,tkinter as tk,viewer
from tkinter import *


        
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
        viewer.display(inputfile,figure=self.figure)
        self.figure+=1
        
    def exec_maxmin(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        print (inputfile)
        viewer.display_maxmin(inputfile,figure=self.figure)
        self.figure+=1
        
    def exec_daily_minima(self):
        inputfile=self.file_list.get(self.file_list.curselection())
        print (inputfile)
        viewer.display_daily_minima(inputfile,figure=self.figure)
        self.figure+=1
        
root = tk.Tk()
app = Viewer(master=root)
app.mainloop()

