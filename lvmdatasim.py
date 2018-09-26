# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Top-level manager for LVM spectroscopic simulation.
"""

from __future__ import print_function, division
import sys
import numpy as np
import pickle
import astropy.io.fits as fits
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
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



    psfModel: float or int or list or str or None
        If float or int generates a symmetric 2D Gaussian kernel of FWHM=psfModel
        If list=[FWHM_a,FWHM_b,theta] and len(list)==3 generates elliptical and rotated Gaussian kernel 
        If str reads a fits file image as the kernel
        If None or False no psf convolution is performed

    """

    def __init__ (self, config, input, telescopeName, psfModel, inputType='fitscube', saveLensCube=False, savePsfCube=False):
        """ 
        Initialize for Simulator
        """
        self.input = input
        self.inputType = inputType
        self.telescopeName = telescopeName
        self.psfModel = psfModel
        self.savelenscube = savelenscube
        self.savepsfcube = savepsfcube


        self.data= self.readInput()
        self.telescope= self.setTelescope()
        self.procdata= self.processInput()
        

        
    def settelescope(self):
        return Telescope(self.telescopename)

    def makeKernel(self, kernelModel):
        if isinstance(self.psfModel, (float or int)):
            """
            Make Symmetric Gaussian 2D PSF
            - Need to calculate the scaling between the plate scale and the PSF model
            """
            scale = 1.0
            return(Gaussian2DKernel(x_stddev=scale*self.psfModel/2.355, mode='integral'))

        elif isinstance(self.psfModel, list):
            """
            Make rotated Gaussian 2D elliptical PSF
            
            """
            if len(self.psfModel) == 3:
                # Need to calculate the scaling between the plate scale and the PSF model
                scale = 1.0 # Place holder.
                #Extract PSF model parameters from the list
                (a_stddev, b_stddev) = psfmodel[0:2]/2.355
                (a_stddev, b_stddev) = (a_stddev*scale, b_stddev*scale)
                theta = psfmodel[2]

                return(Gaussian2DKernel(x_stddev=a_stddev, y_stddev=b_stddev, theta=theta, mode='integral'))
            else:
                sys.exit("The provided PSF model is a list, but not of length three. It can not be interpreted as a_FWHM, b_FWHM, theta")

        elif isinstance(self.psfModel, str):
            """
            Read 2D PSF from fits file
            """
            psf = fits.open(self.psfModel)
            return(psf[0].data)

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
        elif self.inputType == 'lencube':
            procdata=self.data
        elif self.inputType == 'psfcube':
            procdata=self.data
        elif self.inputType == 'fitsrss':
            procdata=self.data
        elif self.inputType == 'asciirss':
            procdata=self.data
        return(procdata)
        
    def convolveinput(self):
        if self.psfModel is not (False or None):
            self.psfKernel= self.makeKernel(self.psfModel)
            if self.inputType == ('fitscube'):
                """
                Connvolve with a 2D PSF kernel, store it as convdata, save if requested, and return it
                """            
                self.kernel=convolve_ftt(self.telescope.IFUmodel.lensletKernel, self.psfKernel)
            
            elif self.inputType == ('lenscube'):
                self.kernel=self.psfKernel 

            convdata = convolve_fft(self.data, kernel, normalize_kernel=True)


        else:
           convdata=self.procdata
        return convdata

    def get_data_fluxes(self, data, center_x, center_y):
        """
        - Extract the spectral fluxes from data based on the IFU foot print
        Parameters
        ----------
        'data' = data, sampled data or convolved data
        'center_x' = image pixel coordinate x  where the FOV will be centered and fluxes will be extracted. Default to center
        'center_y' = image pixel coordinate y  where the FOV will be centered and fluxes will be extracted. Default to center
        potentially self.skycor and telescope model, otherwise use pixel
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
        self.fluxes = self.get_data_fluxes(self.convdata, x, y) #intentionally broken, x and y are not defined


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
        (self.lensletKernel, self.lensletPositions) = pickle.load(ifuname+'.pkl')

        

def main():
    """
    main: main should do something if someone calls the function form the command line. Perhaps just report the contents of each class
    """
    pass

if __name__ == '__main__':
    main()



        
