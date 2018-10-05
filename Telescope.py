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
        self.ifuDict = {"LVM160-SCI-S":5}
        self.siteCoordinatesDict = {'LCO':[-29.0146, -17.6926]}
        self.apertureDict = {"LVM160-SCI-S":160}
        self.apertureAreaDict= {"LVM160-SCI-S":np.pi*160}
        self.obstructionAreaDict = {"LVM160-SCI-S":0.3*self.apertureAreaDict["LVM160-SCI-S"]} # This number is absolutely a guess.
        self.fRatioDict   = {"LVM160-SCI-S":6.2}

        self.name = name
        self.site = self.siteDict[self.name]
        self.siteCoordinates = self.siteCoordinatesDict[self.site]
        self.aperture = self.apertureDict[self.name]
        self.apertureArea = self.apertureAreaDict[self.name]
        self.obstructionArea = self.obstructionAreaDict[self.name]
        self.fRatio = self.fRatioDict[self.name]
        
        # IFU model object: ID, x, y, hexagon radius(center to corner)
        self.ifu = IFU(self.ifuDict[self.name])

    def platescale(self, x=0, y=0):
        """Returns the plate scale of a telescope with aperture diameter Ap, and f-ratio fRatio"""
        return(206265/self.aperture/self.fRatio)

    def ifu2sky(self, ra,dec,theta):
        thetarad=theta*np.pi/180. # position angle in radians
        rotlensx=np.cos(thetarad)*self.ifu.lensx+np.sin(thetarad)*self.ifu.lensy
        rotlensy=-np.sin(thetarad)*self.ifu.lensx+np.cos(thetarad)*self.ifu.lensy
        lensdec=dec+rotlensy*self.platescale(rotlensx, rotlensy)/3600.
        lensra=ra+rotlensx*self.platescale(rotlensx, rotlensy)/3600./np.cos(lensdec*np.pi/180.)
        return(lensra, lensdec)

