# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Top-level manager for LVM spectroscopic simulation.
"""

from __future__ import print_function, division

import numpy as np
import pickle
import astropy.io.fits as fits
import astropy.convolution as convolution
import specsim

class LVMSimulator(object):
    """
    Manages the simulation of LVM data

    Parameters:
    -----------

    inputType: str
    	Can be one of the following:
    	'fitscube' = native input datacube in fits format
    	'lenscube' = lenslet convolved datacube in fits format
    	'psfcube' = psf+lenslet convolved datacube in fits format
    	'fitsrss' = RSS file with one spectrum per lenslet
    	'asciirss' = ascii file with one spectrum per lenslet
    """

    def __init__ (self, config, input, telescopename, psfmodel, inputType='fitscube', savelenscube=False, savepsfcube=False):
        """ 
        Initialize for Simulator
        """
        self.telescope= self.settelescope()
        self.psf= self.makepsf()
        self.data= self.readinput()
        self.procdata= self.processinput()

        self.input = input
        self.inputType = input
        self.savelenscube = savelenscube
        self.savepsfcube = savepsfcube
        
    def settelescope(self):
        return Telescope(self.telescopename)

    def makepsf(self):
        if isinstance(self.psfmodel, (float or int)):
            """
            Make GAussian 2D PSF
            - will use astropy.convolution.Gaussian2DKernel
            - need to know pixel scale of input cube
            - therefore only makes sense to run this if inputType is fitscube, lenscube, or psfcube
            """

        elif isinstance(self.psfmodel, str):
            """
            Read 2D PSF from fits file
            """
        elif self.psfmodel is False:
            return self.psfmodel

    def readinput(self):
        if self.inputType == 'fitscube':
            """
            - Read self.input as fits cube, sample, save if requested, and return data
            """
            data = fits.open(self.input)
        elif self.inputType == 'lenscube':
            """
            - Read self.input as fits cube, do nothing, save if requested, and return data
            """
        elif self.inputType == 'psfcube':
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
        return(fftconvolve(self.data, self.telescope.ifu.lenslet_psf, mode="same"))

    def convolvepsf(self):
        if self.psf is not False and self.inputType == ('fitscube' or 'sampledcube'):
            """
            - connvolve with a 2D PSF kernel, store it as convdata, save if requested, and return it
            """
            convdata = fftconvolve(self.data, self.psf, mode="same")
        else:
           convdata=self.procdata
        return convdata

    def get_data_fluxes(self, data, center_x, center_y):
        """
        - Extract the spectral fluxes from data based on the IFU foot print using
        -- data - allows passing data, sampled data or convolved data
        -- self.ifu
        -- center_xy - x,y coordinates in the plane of the sky where the fluxes will be extracted. Default to center
        -- potentially self.skycor and telescope model, otherwise
        """

        pass

    def lvmsimulate(self):
        """
        Convolve with the time averaged focal plane PSF (if necessary)
        """
        self.convdata=self.convolvepsf()
        """
        Create the simspec Simulator object
        """
        self.fluxes = self.get_fluxes()


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
        """
        self.site='LCO'
        self.ifu = IFUmodel('science')


class IFUmodel(object):
    """Read an existing IFU model stored in the data directory as a pickle"""
    def __init__ (self, ifuname):
        (self.lenslet_psf, self.ifu_xy_positions) = pickle.load(ifuname+'.pkl')

        

def main():
    """
    main: main should do something if someone calls the function form the command line. Perhaps just report the contents of each class
    """
    pass

if __name__ == '__main__':
    main()



        
