"""
TOPAS cartesian scorer helpers
"""

from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import matplotlib.gridspec as gridspec



import numpy as np
import csv
import re

TOPAS_SCORERS = [   'ProtonLET', 
                    'DoseToMedium',
                    'DoseToWater',
                    'DoseToMaterial',
                    'EnergyDeposit',
                    'Fluence',
                    'EnergyFluence',
                    'StepCount',
                    'OpticalPhotonCount',
                    'Charge',
                    'EffectiveCharge']
TOPAS_STATS = [     'Sum', 
                    'Mean', 
                    'Histories', 
                    'Count_In_Bin', 
                    'Second_Moment', 
                    'Variance', 
                    'Standard_Deviation', 
                    'Min', 
                    'Max']

class TopasParameterFile:
    def __init__(self, parameterFile):
        self.geometry = {}
        self.source = {}

        for line in open(parameterFile,'rt'):
            # first get rid of comments
            if '#' in line:
                line = line[:line.index('#')]
            
            try:
                component = line.split(':')[1][:2]
                if component == 'Ge':
                    self._parseGeometryLine(line)
                if component == 'So':
                    self._parseSourceLine(line)
            
            except:
                pass
    
    def _parseGeometryLine(self, line):
        name = re.search(':Ge/(.*?)/',line).group(1)
        if name not in self.geometry:
            self.geometry[name] = {}
        topasType = line.split(':')[0]
        key = re.search(name+'/(.*?)=',line).group(1).strip()
        value = self._getTopasValue(line.split('=')[1].strip(), topasType)
        self.geometry[name][key] = value
    
    def _parseSourceLine(self, line):
        name = re.search(':So/(.*?)/',line).group(1)
        if name not in self.source:
            self.source[name] = {}
        topasType = line.split(':')[0]
        key = re.search(name+'/(.*?)=',line).group(1).strip()
        value = self._getTopasValue(line.split('=')[1].strip(), topasType)
        self.source[name][key] = value

    def _getTopasValue(self, lineVal, topasType):
        if topasType == 's':
            # single string; remove the qupte marks
            return lineVal[1:-1]
        if topasType == 'i':
            # single integer
            return int(lineVal)
        if topasType == 'b':
            # boolean like "FALSE"
            return lineVal[1:-1].upper() == 'TRUE'
        if topasType == 'd':
            # double with unit. Return as tuple.
            return (float(lineVal.split()[0]), lineVal.split()[1])
        if topasType == 'sv':
            # string vector. return as list of strings
            return [val[1:-1] for val in lineVal.split()[1:]]

class TopasHeader:
    def __init__(self, csvFile):
        self.topasVersion = ''
        self.parameterFile = ''
        self.pInfo = None
        self.scorer = ''
        self.scoredComponent = ''
        self.X = readBins('')
        self.Y = readBins('')
        self.Z = readBins('')
        self.scoredQuantity = readQuantities('')
        for line in open(csvFile,'rt'):
            if not line.startswith('#'):
                break
            if line.startswith('# TOPAS Version'):
                self.topasVersion = line.rstrip().split(': ')[1]
            elif line.startswith('# Parameter File'):
                self.parameterFile = line.rstrip().split(': ')[1]
            elif line.startswith('# Results for scorer '):
                self.scorer = line[len('# Results for scorer '):].rstrip()
            elif line.startswith('# Scored in component'):
                self.scoredComponent = line.rstrip().split(': ')[1]
            elif line.startswith('# X in '):
                self.X = readBins(line)
            elif line.startswith('# Y in '):
                self.Y = readBins(line)
            elif line.startswith('# Z in '):
                self.Z = readBins(line)
            elif any([line.startswith('# '+scorer) for scorer in TOPAS_SCORERS]):
                self.scoredQuantity = readQuantities(line)

class readQuantities:
    def __init__(self, line):
        self.unit = ''
        self.name = ''
        self.stats = []
        
        try:
            self.unit = line.split()[3]
            self.name = line.split()[1]
            self.stats = line.split(': ')[1].split()
        except:
            pass

class readBins:
    def __init__(self, line):
        self.bins = 0
        self.size = 0.0
        self.unit = ''
        self.scale = np.linspace(0,self.size*self.bins,self.bins)
        try:
            self.bins = int(line.split(' bins of ')[0].split()[-1])
            self.size = float(line.split(' bins of ')[1].split()[0]) # always work in cm
            self.unit = 'cm'
            unit = line.split(' bins of ')[1].split()[1]
            if unit == 'mm':
                self.size = self.size/10.0
            if unit == 'm':
                self.size = self.size*10.0
            if unit == 'nm':
                self.size = self.size/1E8
            # Add more in future if necessary...
        except:
            pass
        self.scale = np.linspace(0,self.size*self.bins,self.bins)

def val_cm(val):
    """
    arg val is tuple (2.0,'mm').
    Return just the number in cm.
    """
    if val[1] == 'cm':
        return val[0]
    if val[1] == 'mm':
        return val[0]/10.0
    if val[1] == 'm':
        return val[0]*10
    if val[1] == 'nm':
        return val[0]/1E8

class TopasResults:
    """
    Class to parse a topas csv file written by a scoring component
    
    Args: 
        csvFile (string)
    """

    def __init__(self, csvFile):
        self.fileName = csvFile
        self.header = TopasHeader(csvFile)
        self.param = TopasParameterFile(self.header.parameterFile)
        self._read_csvFile()
        #ASSUMES NO ROTATION - REVISIT LATER
        self.X_cm = np.linspace(
            val_cm(self.param.geometry[self.header.scoredComponent]['TransX']) - val_cm(self.param.geometry[self.header.scoredComponent]['HLX']), 
            val_cm(self.param.geometry[self.header.scoredComponent]['TransX']) + val_cm(self.param.geometry[self.header.scoredComponent]['HLX']), 
            self.header.X.bins, endpoint=False)
        self.Y_cm = np.linspace(
            val_cm(self.param.geometry[self.header.scoredComponent]['TransY']) - val_cm(self.param.geometry[self.header.scoredComponent]['HLY']), 
            val_cm(self.param.geometry[self.header.scoredComponent]['TransY']) + val_cm(self.param.geometry[self.header.scoredComponent]['HLY']), 
            self.header.Y.bins, endpoint=False)
        self.Z_cm = np.linspace(
            val_cm(self.param.geometry[self.header.scoredComponent]['TransZ']) - val_cm(self.param.geometry[self.header.scoredComponent]['HLZ']), 
            val_cm(self.param.geometry[self.header.scoredComponent]['TransZ']) + val_cm(self.param.geometry[self.header.scoredComponent]['HLZ']), 
            self.header.Z.bins, endpoint=False)
        self.X_cm_cent = self.X_cm + 0.5*self.header.X.size
        self.Y_cm_cent = self.Y_cm + 0.5*self.header.Y.size
        self.Z_cm_cent = self.Z_cm + 0.5*self.header.Z.size
        self.X_cm_extent = (self.X_cm[0], self.X_cm[-1]+self.header.X.size)
        self.Y_cm_extent = (self.Y_cm[0], self.Y_cm[-1]+self.header.Y.size)
        self.Z_cm_extent = (self.Z_cm[0], self.Z_cm[-1]+self.header.Z.size)
      
    def _read_csvFile(self):
        self.data = {}
        
        for stat in self.header.scoredQuantity.stats:
            self.data[stat] = np.zeros((self.header.X.bins, self.header.Y.bins, self.header.Z.bins))

        for line in open(self.fileName, 'rt'):
            if not line.startswith('#'):
                [x,y,z] = [int(i) for i in line.split(', ')[0:3]]
                for ind,val in enumerate(line.split(', ')[3:]):
                    self.data[self.header.scoredQuantity.stats[ind]][x,y,z] = float(val)

        # This part I will need to be revisited. 
        # TOPAS seems to start binning from the rear to front face, 
        # not 100% sure what is happening in x and y
        for stat in self.header.scoredQuantity.stats:
            self.data[stat] = self.data[stat][:,:,::-1]


def displaySlice(topasResults, quantity, fixedDim, ax=None):
    """
    Use pyplot imshow to display a 2D slice from the topasResults.data[quantity] array
    fixedDim is tuple (dim, bin index) e.g. (x,3). 
    fixedDim[1] = None default behaviour is a central slice.
    
    ax is the image axis we update
    """
    if fixedDim[0] not in ['x','y','z']:
        return

    title = '{} {} [{}]\nScorer: {}; Scored component: {}\nParameterFile: {}\nResultsFile:{}'.format(
        topasResults.header.scoredQuantity.name, quantity, topasResults.header.scoredQuantity.unit,
        topasResults.header.scorer, 
        topasResults.header.scoredComponent,
        topasResults.header.parameterFile, 
        topasResults.fileName)

    fixSlice = fixedDim[1]

    if fixedDim[0] == 'x':
        left, right = topasResults.Z_cm_extent
        bottom, top = topasResults.Y_cm_extent

        (xlabel, ylabel) = ('Z (cm)', 'Y (cm)')
        if not fixSlice:
            fixSlice = int(0.5*topasResults.header.X.bins)
        
        title += '\n(Slice: {} = {} cm)'.format(
            fixedDim[0],
            topasResults.X_cm_cent[fixSlice])

        plotData = topasResults.data[quantity][fixSlice,:,:]
    
    elif fixedDim[0] == 'y':
        left, right = topasResults.Z_cm_extent
        bottom, top = topasResults.X_cm_extent

        (xlabel, ylabel) = ('Z (cm)', 'X (cm)')
        
        if not fixSlice:
            fixSlice = int(0.5*topasResults.header.Y.bins)
        
        title += '\n(Slice: {} = {} cm)'.format(
            fixedDim[0],
            topasResults.Y_cm_cent[fixSlice])

        plotData = topasResults.data[quantity][:,fixSlice,:]
    
    elif fixedDim[0] == 'z':
        left, right = topasResults.X_cm_extent
        bottom, top = topasResults.Y_cm_extent

        (xlabel, ylabel) = ('X (cm)', 'Y (cm)')

        if not fixSlice:
            fixSlice = int(0.5*topasResults.header.Z.bins)

        title += '\n(Slice: {} = {} cm)'.format(
            fixedDim[0],
            topasResults.Z_cm_cent[fixSlice])

        plotData = topasResults.data[quantity][:,:,fixSlice]
    
    img= None
    if ax:
        img = ax.imshow(plotData, extent = (left, right, bottom, top), interpolation = 'none')
        ax.set_title(title, fontsize=10); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    
    return img


def getBin(topasResults, dim, cmValue, discrete=True):
    """
    Return the bin number for a given cmValue  
    """
    cmCents = {'x':topasResults.X_cm_cent, 'y':topasResults.Y_cm_cent, 'z':topasResults.Z_cm_cent}
    
    if cmValue < cmCents[dim].min() or cmValue > cmCents[dim].max():
        return None
    else:
        for ind,val in enumerate(cmCents[dim][:-1]):
            if val-cmValue * cmCents[dim][ind+1]-cmValue <= 0:
                if discrete:
                    return int(round(ind + (cmValue - val)/float(cmCents[dim][ind+1]-val)))
                else:
                    return ind + (cmValue - val)/float(cmCents[dim][ind+1]-val)

def getMiddle_cm(topasResults, dim):
    cmCents = {'x':topasResults.X_cm_cent, 'y':topasResults.Y_cm_cent, 'z':topasResults.Z_cm_cent}
    return cmCents[dim][int(0.5*len(cmCents[dim]))]

def orthViews(topasResults):
    """
    matplotlib 
    """
    plt.ion()
    gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1.5, 5, 1])
    fig = plt.figure(figsize=(15,8))
    axRadBut = plt.subplot(gs[0])
    radButs = RadioButtons(axRadBut, labels=topasResults.header.scoredQuantity.stats, active=0)

    sposSetup = {
        'x': {
            'valmin':topasResults.X_cm_cent.min(), 
            'valmax':topasResults.X_cm_cent.max(),
            'valinit':getMiddle_cm(topasResults, 'x')}, 
            #'valstep':topasResults.header.X.size},
        'y': {
            'valmin':topasResults.Y_cm_cent.min(), 
            'valmax':topasResults.Y_cm_cent.max(),
            'valinit':getMiddle_cm(topasResults, 'y')}, 
            #'valstep':topasResults.header.Y.size},
        'z': {
            'valmin':topasResults.Z_cm_cent.min(), 
            'valmax':topasResults.Z_cm_cent.max(),
            'valinit':getMiddle_cm(topasResults, 'z')}, 
            #'valstep':topasResults.header.Z.size}
            }   

    spos, axes, axSlide = dict(), dict(), dict()

    for ind,ax in enumerate(['x','y','z']):
        axes[ax] = plt.subplot(gs[ind+3])

        axSlide[ax] = plt.subplot(gs[ind+6])
        spos[ax] = Slider(
            axSlide[ax], ax+' bin [cm]', **sposSetup[ax])
        startBin = getBin(topasResults,ax,sposSetup[ax]['valinit'])
        displaySlice(topasResults, radButs.value_selected, (ax, startBin), axes[ax])

    #plt.tight_layout()



    def updateAll(val):
        updateX(0.0)
        updateY(0.0)
        updateZ(0.0)

    def updateX(val):
        dim = 'x'
        selBin = getBin(topasResults, dim, spos[dim].val)
        displaySlice(topasResults, radButs.value_selected, (dim, selBin), axes[dim])
        fig.canvas.draw_idle()
    def updateY(val):
        dim = 'y'
        selBin = getBin(topasResults, dim, spos[dim].val)
        displaySlice(topasResults, radButs.value_selected, (dim, selBin), axes[dim])
        fig.canvas.draw_idle()
    def updateZ(val):
        dim = 'z'
        selBin = getBin(topasResults, dim, spos[dim].val)
        displaySlice(topasResults, radButs.value_selected, (dim, selBin), axes[dim])
        fig.canvas.draw_idle()
    spos['x'].on_changed(updateX)
    spos['y'].on_changed(updateY)
    spos['z'].on_changed(updateZ)
    radButs.on_clicked(updateAll)    


myFile = r'/Users/jonathanmohajer/2ProtonLET.csv'
