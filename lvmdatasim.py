# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Top-level manager for LVM spectroscopic simulation.
"""

from __future__ import print_function, division

import numpy as np
import pickle
import astropy.io.fits as fits
import specsim

class LVMSimulator(object):
    """
    Manages the simulation of LVM data
    """

    def __init__ (self, config, input, telescope, psfmodel, inputType='fitscube', savelenscube=False, savepsfcube=False ):
        """ 
        Initialize for Simulator
        """
        self.telescope= self.settelescope(telescope)
        self.psf= self.makepsf(psfmodel)
        self.data= self.readinput()
        self.procdata= self.processinput()

        self.input = input
        self.inputType = input
        self.savelenscube = savelenscube
        self.savepsfcube = savepsfcube
        
    def settelescope(self, telescope):
        return Telescope()

    def makepsf(self, psfmodel):
        if isinstance(psfmodel, (float or int)):
            """
            Make GAussian 2D PSF
            """
        elif isinstance(psfmodel, str):
            """
            Read 2D PSF from fits file
            """
        elif psfmodel is False:
            return psfmodel

    def readinput(self):
        if self.inputType == 'fitscube':
            """
            - Read self.input as fits cube, sample, save if requested, and return data
            """
            data = fits.open(self.input)
        elif self.inputType == 'sampledcube':
            """
            - Read self.input as fits cube, do nothing, save if requested, and return data
            """
        elif self.inputType == 'fitsrss':
            """ 
            - Read self.input as fits RSS file with spectra for each spaxel and return data
            """
        elif self.inputType == 'asciirss':
            """ 
            - Read self.input as ascii file with one spectrum per each spaxel (1st column = wavelength, each following column is one spectrum), and return data
            """
        return(data)
            
    def processinput(self):
        if self.inputType == 'fitscube':
            """ 
            - sample with lenslets, save if requested, and return processed data, and del(data)
             """
            procdata=self.convolvelenslet()
        elif self.inputType == 'sampledcube':
            procdata=self.data
        elif self.inputType == 'fitsrss':
            procdata=self.data
        elif self.inputType == 'asciirss':
            procdata=self.data
        return(procdata)
        
    def convolvelenslet(self):
        """
        - convolve with a hexagon of the right size given self.telescope.IFUmodel
        """

    def convolvepsf(self):
        if self.psfmodel is not False and self.inputType == ('fitscube' or 'sampledcube'):
            """
            - connvolve with a 2D PSF kernel, store it as convdata, save if requested, and return it
            """
        else:
           convdata=procdata
        return convdata


    def lvmsimulate(self):
        """
        Convolve with the time averaged focal plane PSF (if necessary)
        """
        self.convdata=self.convolvepsf()
        """
        Crate the simspec Simulator object
        """
        






    
class Telescope(object):

    def __init__ (self):
        """
        Initialize for Telescope class
        """
        self.site='LCO'
        self.ifu=lvmdatasim.IFUmodel('science')


class IFUmodel(object):

    def __init__ (self, ifuname):
        pickle.load(ifuname+'.pkl')


        

        

        
def main():
    telescope = lvmdatasim.Telescope()
    telescope.set_Telescope("lvm160-s")
    telescope.IFUModel.setNFibers(500)



        
