"""
TOPAS cartesian scorer helpers
"""
import matplotlib
matplotlib.use('TkAgg')

from glob import glob
from plotting import TopasResults, displaySlice, getBin
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2
from matplotlib.backend_bases import key_press_handler

import matplotlib.pyplot as plt

import sys
if sys.version_info[0] < 3:
    from ttk import Tkinter as tk
    import tkFileDialog
    import ttk
else:
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog as tkFileDialog

DEFAULT_TOPAS_RESULTS_FOLDER = ''

class topasFrameImg:
    def __init__(self, master=None):
        self.fig = plt.Figure(figsize=(15,5))
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.axes = {   'x':self.fig.add_subplot(131),
                        'y':self.fig.add_subplot(132),
                        'z':self.fig.add_subplot(133)}
        self.fig.tight_layout()
        for view in self.axes:
            self.axes[view].img = None

class topasResultsGUI(tk.Frame,object):
    def __init__(self, topasFile ='', rootDir=DEFAULT_TOPAS_RESULTS_FOLDER, master=None):
        super(topasResultsGUI, self).__init__(master)
        self.grid()
        self.results = None
        self.master = master
        self.master.protocol('WM_DELETE_WINDOW', self.quit)
        self.master.geometry('1600x640')
        self.rootDir = rootDir
        self.title = 'Topas Results Viewer'
        self.master.title(self.title)
        self.topasFile = topasFile
        self.create_widgets()

    def quit(self, *args):
        self.master.destroy()
        raise SystemExit

    def create_widgets(self):
        self.createInputsFrame()
        self.createViewerFrame()

    def createInputsFrame(self):
        self.inputsFrame = tk.Frame(self, width=800, height=100)
        self.inputsFrame.pack(fill="both", expand=True)
        self.entryPath = tk.Entry(self.inputsFrame, width=100)
        self.entryPath.insert(tk.END, self.topasFile)
        self.entryPath.bind("<Return>", (lambda event: self.entryPathReturn()))
        self.lblLoad = tk.Label(self.inputsFrame)
        self.getFileButton = tk.Button(self.inputsFrame, text="Select topas csv", command=self.getFile)
        self.getFileButton.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        self.entryPath.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        self.lblLoad.grid(row=0, column=2, sticky="w", padx=5, pady=5)
        if self.topasFile:
            self.readFile()

    def createViewerFrame(self):
        self.viewerFrame = tk.Frame(self, width=1600, height=500)
        self.viewerFrame.pack(fill="both", expand=True)
        self.orthViews = topasFrameImg(master=self.viewerFrame)

        self.orthViews.canvas.get_tk_widget().grid(row=0, column=0, rowspan=2, columnspan=6)
        self.toolbar = NavigationToolbar2(self.orthViews.canvas)
        self.toolbar.update()
        self.orthViews.canvas.get_tk_widget().grid(row=0, column=0, rowspan=2, columnspan=6)
        # canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        self.climMax = tk.Scale(
            self.viewerFrame, from_=100,to=0, resolution= 1, state='disabled', command=self._adjContrast)
 
        self.climMin = tk.Scale(
            self.viewerFrame, from_=100,to=0, resolution= 1, state='disabled', command=self._adjContrast)
        
        self.climMinLabel = tk.Label(master=self.viewerFrame, text='Min').grid(row=0,column=6,sticky="s")
        self.climMaxLabel = tk.Label(master=self.viewerFrame, text='Max').grid(row=0,column=7,sticky="s")
        self.climMax.grid(row=1,column=7,sticky="n")
        self.climMin.grid(row=1,column=6,sticky="n")
        self.orthScales = {
            'x': tk.Spinbox(self.viewerFrame, values=[''], state='disabled', command=self._updatePlot),
            'y': tk.Spinbox(self.viewerFrame, values=[''], state='disabled', command=self._updatePlot),
            'z': tk.Spinbox(self.viewerFrame, values=[''], state='disabled', command=self._updatePlot)}
        self.orthScaleLabels = {
            'x':tk.Label(master=self.viewerFrame, text='x bin centre [cm]'),
            'y':tk.Label(master=self.viewerFrame, text='y bin centre [cm]'),
            'z':tk.Label(master=self.viewerFrame, text='z bin centre [cm]')
        }
        self.statsCbox = ttk.Combobox(master=self.viewerFrame, values = [''], state="readonly")
        self.statsCbox.current(0)
        self.statsCbox.bind("<<ComboboxSelected>>", self._statsSelect)
        self.statsCbox.grid(row=2,column=0,sticky='w')
        for ind,view in enumerate(sorted(self.orthScales)):
            self.orthScales[view].grid(row=2+ind,column=1,sticky="e")
            self.orthScaleLabels[view].grid(row=2+ind,column=1,sticky="w")

# ---------------------------- #
    def _statsSelect(self,*args):
        self._updatePlot()

    def entryPathReturn(self):
        self.topasFile = self.entryPath.get()
        self.readFile()

    def getFile(self):
        self.topasFile = tkFileDialog.askopenfilename(
            title= 'Select a topas csv file',
            initialdir=self.rootDir,
            filetypes=[('Comma-seperated values','.csv'),('All files','.*')])
        if self.topasFile == '':
            return 

        self.entryPath.delete(0, tk.END)
        self.entryPath.insert(tk.END, self.topasFile)
        self.readFile()

    def readFile(self):
        try:
            self.results = TopasResults(self.topasFile)
            self.lblLoad['text'] = 'File read successfully'
            self._initialisePlots()
        except:
            self.results = None
            self.lblLoad['text'] = 'Error reading file'

    def _initialisePlots(self):
        if self.results:
            for view in self.orthViews.axes:
                self.orthViews.axes[view].img = displaySlice( 
                    self.results, self.results.header.scoredQuantity.stats[0], 
                    (view,None), self.orthViews.axes[view])
            
            self.orthViews.fig.subplots_adjust(right=0.90)
            self.cbar_ax = self.orthViews.fig.add_axes([0.93, 0.05, 0.01, 0.90])
            self.orthViews.fig.colorbar(self.orthViews.axes['x'].img, cax=self.cbar_ax, format='%.1e')
            self.cbar_ax.tick_params(labelsize=8)
            self.orthViews.canvas.draw()
            self.climMin['state'] = 'normal'
            self.climMax['state'] = 'normal'
            self.climMax.set(100)
            self.statsCbox['values'] = self.results.header.scoredQuantity.stats
            self.statsCbox.current(0)
            orthScalesHelper = {
                'x':tuple(self.results.X_cm_cent),
                'y':tuple(self.results.Y_cm_cent),
                'z':tuple(self.results.Z_cm_cent)}
            for view in self.orthScales:
                
                self.orthScales[view]['values'] = orthScalesHelper[view]
                self.orthScales[view]['state'] = 'normal'
    
    def _updatePlot(self, *args):
        if self.statsCbox['values'] != ['']:
            for view in self.orthViews.axes:
                self.orthViews.axes[view].img = displaySlice( 
                    self.results, self.statsCbox.get(), 
                    (view,getBin(self.results, view, float(self.orthScales[view].get()))), self.orthViews.axes[view])
            self.orthViews.fig.colorbar(self.orthViews.axes['x'].img, cax=self.cbar_ax, format='%.1e')
            self.orthViews.canvas.draw()
    
    def _adjContrast(self, *args):
        if self.climMin.get() < self.climMax.get():
            for view in self.orthViews.axes:
                    imgMax = self.orthViews.axes[view].img.get_array().max()
                    self.orthViews.axes[view].img.set_clim(
                        [imgMax*self.climMin.get()/100.0, imgMax*self.climMax.get()/100.0]
                        )
            self.orthViews.canvas.draw()        


# run the GUI
def main():
    root =tk.Tk()
    view= topasResultsGUI(master=root)
    view.mainloop()

if __name__ == '__main__':
    main()