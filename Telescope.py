from IFU import IFU
import numpy as np

class Telescope(object):
    """
    Telescope class:

    Parameters:
    -----------

    name: str
    	Telescope name. Syntax is LVM[160,1000]-[SCI,CAL,SKY]-[N,S]. So for example, the spectrophotometric
    	calibration 0.16m telescope at LCO (i.e. South) would be "LVM160-CAL-S"
    """

    def __init__ (self, name):
        """
        Initialize for Telescope class
        Uses hardcoded dictionaries and atributes are set using telescope name as key
        """
        self.siteDict = {"LVM160-SCI-S":"LCO"}
        self.siteCoordinatesDict = {'LCO':[-29.0146, -17.6926]}
        self.apertureDict = {"LVM160-SCI-S":160}
        self.apertureAreaDict= {"LVM160-SCI-S":np.pi*160}
        self.obstructionAreaDict = {"LVM160-SCI-S":0.3*self.apertureADict[self.name]} # This number is absolutely a guess.
        self.fRatioDict   = {"LVM160-SCI-S":6.2}

        self.name = name
        self.site = self.siteDict[self.name]
        self.siteCoordinates = self.siteCoordinatesDict[self.name]
        self.aperture = self.siteDict[self.name]
        self.apertureArea = self.apertureDict[self.name]
        self.obstructionArea = self.obstructionDict[self.name]
        self.fRatio = self.fRatioDict[self.name]
        
        # IFU model object: ID, x, y, hexagon radius(center to corner)
        self.ifu = IFU(self.name)

    def platescale(self, x=0, y=0):
        """Returns the plate scale of a telescope with aperture diameter Ap, and f-ratio fRatio"""
        return(206265/self.apertureADict[self.name]/self.fRatioDict[self.name])
