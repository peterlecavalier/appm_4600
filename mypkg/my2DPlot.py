import matplotlib.pyplot as plt
import numpy as np
import math
class my2DPlot:
    def __init__(self,f,a,b):
        # This init initializes the plotting window
        x = np.arange(a,b,0.01)
        y = f(x)
        self.p = plt.plot(x,y)
    def show(self):
        # Shows the plot
        plt.show()
    def dotted(self):
        # Changes the most recent line to dotted
        self.p[-1].set_linestyle('dotted')
    def labels(self, x, y):
        # Adds custom x and y labels to the plot
        plt.xlabel(x)
        plt.ylabel(y)
    def addPlot(self, f):
        # Adds a new plot to the window
        x = self.p[0].get_data()[0]
        y = f(x)
        self.p = plt.plot(x, y)
    def color(self, colorName):
        # Changes the color of the most recent plot
        self.p[-1].set_color(colorName)
    def logy(self):
        # Sets the scale of the y-axis to log
        plt.yscale("log")
    def logx(self):
        # Sets the scale of the x-axis to log
        plt.xscale("log")
    def save(self, fileName):
        # Saves the plot to a custom filename
        plt.savefig(fileName)