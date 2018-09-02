# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 13:50:32 2018

@author: Paulo Augusto
"""

from Tkinter import *
 
# these two imports are important
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import time
import threading


class dataPlotter():



    
    def __init__(self,xlabel,ylabel,indata):
        self.data=indata
        self.root = Tk()
        self.continuePlotting = False
        self.xlabel=xlabel
        self.ylabel=ylabel
    def change_state(self):

        if self.continuePlotting == True:
            self.continuePlotting = False
        else:
            self.continuePlotting = True
     
    def app(self):
        # initialise a window.

        self.root.config(background='white')
        self.root.geometry("600x400")
    #    root.protocol("WM_DELETE_WINDOW", lambda: ser.close())
        

        self.lab = Label(self.root, text="Live Plotting", bg = 'white').pack()

        self.fig = Figure()
        
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.grid()
     
        self.graph = FigureCanvasTkAgg(self.fig, master=self.root)
        self.graph.get_tk_widget().pack(side="top",fill='both',expand=True)
     ##########################################################################

    def plotter(self):
        while self.continuePlotting:
            self.ax.cla()           
            self.ax.grid()
            dpts = self.data#data_points()
            self.ax.plot(dpts, marker='o', color='orange')
            self.graph.draw()
            time.sleep(1)
 
    def gui_handler(self):
        self.change_state()
        threading.Thread(target=self.plotter).start()
#            threading.Thread(target=reader).start()
        
    def start(self):
        self.app()
        self.b = Button(self.root, text="Start/Stop", command=self.gui_handler, bg="red", fg="white")
        self.b.pack()
        
        self.root.mainloop()
 
