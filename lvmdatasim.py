# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Top-level manager for LVM spectroscopic simulation.
"""

from __future__ import print_function, division
import sys
import os
import numpy as np
import pickle
from astropy.io import fits as fits
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from reproject import reproject_interp
import astropy.wcs as wcs
import hexagonlib as hexlib
from PIL import Image, ImageDraw
from Telescope import Telescope
import scipy.interpolate as interpolate

try:
    sys.path.append(os.environ['SIMSPEC_DIR'])
except:
    sys.exit("No SIMSPEC_DIR environmental variable defined. Good luck!")
try:
    import specsim
except:
    sys.exit("The SIMSPEC_DIR environmental varible is defined, but the import of specsim failed: Wrong path or specsim is not installed")
    

class LVMSimulator(object):
    """
    Manages the simulation of LVM data

    Parameters:
    -----------
        
    input: str
        path to input data file

    telescopeName: str
        Telescope name. Syntax is LVM[160,1000]-[SCI,CAL,SKY]-[N,S]. So for example, the spectrophotometric
        calibration 0.16m telescope at LCO (i.e. South) would be "LVM160-CAL-S"

    psfModel: float or int or list or str or None
        If float or int generates a symmetric 2D Gaussian kernel of FWHM=psfModel
        If list=[FWHM_a,FWHM_b,theta] and len(list)==3 generates elliptical and rotated Gaussian kernel 
        If str reads a fits file image as the kernel
        If None or False no psf convolution is performed

    inputType: str
        Can be one of the following:
        'fitscube' = native input datacube in fits format
        'lenscube' = lenslet convolved datacube in fits format
        'psfcube' = psf+lenslet convolved datacube in fits format
        'fitsrss' = RSS file with one spectrum per lenslet, first row is wavelength in A
        'asciirss' = ascii file with one spectrum per lenslet, first column is wavelength in A, headers must be commented with "#"

    fluxType: str
        Can be one of the following:
        'intensity' = input is assumed to be in units of erg/s/cm2/arcsec2 (for both cubes and RSS files)
        'flux' = input is assumed to be in units of erg/s/cm2/pixel (for cubes) or erg/s/cm2 (for RSS files)

    """

    def __init__ (self, input, telescopeName, psfModel, inputType='fitscube', fluxType='intensity', saveConvCube=True, yamlfile='lvmdatasim.yaml'):
        """ 
        Initialize for Simulator
        """
        self.input = input
        self.inputType = inputType
        self.fluxType = fluxType
        self.telescopeName = telescopeName
        self.psfModel = psfModel
        self.saveConvCube = saveConvCube
        self.yamlfile=yamlfile

        self.data, self.hdr = self.readInput()
        
        self.telescope = Telescope(telescopeName)
        """ still need to define how user sets parameters of exposure and simulation"
        """

        import lvmdatasimdefaults
        self.simparam = lvmdatasimdefaults.param
        self.convdata = None

    """
    Lensed cube should store PA in header, and if imputType=lenscube or psfcube code should check that PA in header is consistent with PA of observation, otherwise raise error.
    """

    def readInput(self):
        if self.inputType == ('fitscube' or 'lenscube' or 'psfcube'):
            """
            - Read self.input as fits cube and return data and header
            - Add a PIXSCALE keyword to the header if not present before passing it on
            """
            data = fits.open(self.input)
            if ('PIXSCALE' not in data[0].header.keys()) and ('CDELT1' not in data[0].header.keys()) and ('CD1_1' not in data[0].header.keys()):
                sys.exit('No WCS Information in input FITS header')
            elif 'PIXSCALE' not in data[0].header.keys(): 
                mywcs=wcs.WCS(data[0].header)
                pixscale=wcs.utils.proj_plane_pixel_scales(mywcs).mean()
                data[0].header.set('PIXSCALE', pixscale, 'Pixel scale calculated from WCS by LVMSimulator')    
            return(data[0].data, data[0].header)

        elif self.inputType == 'fitsrss':
            """ 
            - Read self.input as fits RSS file with spectra for each spaxel and return data
            """
            data = fits.open(self.input)
            return(data.data, data.header)
        elif self.inputType == 'asciirss':
            """ 
            - Read self.input as ascii file with one spectrum per each spaxel (1st column = wavelength, each following column is one spectrum), and return data
            """
            data=np.genfromtxt(self.input, comments="#", unpack=True)
            return(data, '')
        else:
            sys.exit('Input Type \"'+self.inputType+'\" not recognized.')

   
    def makePsfKernel(self):
        try:
            pixscalecube = self.hdr['PIXSCALE']
        except:
            sys.exit("Something went wrong. You made it this far with no 'PIXSCALE' defined in the data header.")

        if isinstance(self.psfModel, (float or int)):
            """
            Make Symmetric Gaussian 2D PSF
            - Need to calculate the scaling between the plate scale and the PSF model
            """
            return(Gaussian2DKernel(pixscalecube*self.psfModel/2.355, mode='oversample', factor=25))

        elif isinstance(self.psfModel, list):
            """
            Make rotated Gaussian 2D elliptical PSF
            
            """
            if len(self.psfModel) == 3:
                #Extract PSF model parameters from the list
                (a_stddev, b_stddev) = self.psfModel[0:2]/2.355
                (a_stddev, b_stddev) = (a_stddev*pixscalecube, b_stddev*pixscalecube)
                theta = self.psfModel[2]

                return(Gaussian2DKernel(x_stddev=a_stddev, y_stddev=b_stddev, theta=theta, mode='oversample', factor=25))
            else:
                sys.exit("The provided PSF model is a list, but not of length three. It can not be interpreted as a_FWHM, b_FWHM, theta")

        elif isinstance(self.psfModel, str):
            """
            Read 2D PSF from fits file
            - Right now rebining fits kernel with interpolation, it would be better to do an actual integration over a surface fit
            """
            if self.psfModel.endswith() != ".fits": 
                sys.exit("psfModel is str but does not end in \".fits\" as expected")
            else:
                psf = fits.open(self.psfModel)
                pixscalepsf=psf.header['PIXSCALE']
                wcs0=wcs.WCS(naxis=2)
                wcs0.wcs.crpix=[0,0]
                wcs0.wcs.crval=[0,0]
                wcs0.wcs.cdelt=np.full(2, pixscalepsf)
                wcs1=wcs.WCS(naxis=2)
                wcs1.wcs.crpix=[0,0]
                wcs1.wcs.crval=[0,0]
                wcs1.wcs.cdelt=np.full(2, pixscalecube)

                psf[0].header=wcs0.to_header()
                psfrebin=reproject_interp(psf, wcs1)

            return(psfrebin[0].data)

        elif self.psfModel is False:
            return self.psfModel


    def updateyaml(self):
        with open("lvmdatasimTemplate.yaml", 'rb') as f:
            data = f.read()  # produces single string
        for key in self.simparam.keys():
            data.replace("__%s__placeholder"%key, self.simparam[key])
        with open(self.yamlfile, 'w') as f:
            f.writelines(data)


    def makeLensKernel(self):
        """
        Generate a 2D numpy array with the characteristic function of a hexagon of size 'radius' (center to corner)
        in units of pixels of the array. The image will have dimensions of imgsize x imgsize and the hexagon will
        be at the integer center.
        """
        rlensmm=np.mean(self.telescope.ifu.lensr)
        rlensarcsec=rlensmm*self.telescope.platescale()
        rlenspix=rlensarcsec/self.hdr['PIXSCALE']
        imgsize=2*rlenspix
        antialias=5
        imgsize=int(imgsize*antialias)
        center = imgsize//2+1
        hexLayout = hexlib.Layout(hexlib.layout_pointy, rlenspix*antialias, hexlib.Point(center,center))
        polygon = hexlib.polygon_corners(hexLayout,hexlib.Hex(0,0,0))
        img = Image.new('L', (imgsize, imgsize), 0)
        ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
        kernel = np.array(img,dtype=float)
        kernel /= np.sum(kernel)
        kernel = int_rebin(kernel, (imgsize//antialias,imgsize//antialias))
        return kernel

    
    def convolveInput(self):
        """ results of call:
        0,0 - no convolution
        1,0 - psf convolution, no lens convolution
        0,1 - no psf convolution. lens convolution
        1,1 - psf and lens convolution
        """
        if self.inputType == "psfcube":
            # if is a psfcube dont do anything (0,0)
            self.convdata=self.data
        else:
            # if not psfcube it needs some sort of convolution
            if self.psfModel is not (False or None):
                # if psfModel defined then make psf kernel
                self.psfKernel= self.makePsfKernel()
                if self.inputType == ('fitscube'):
                    # if its fitscube kernel is convolution of psf+lenslet (1,1)
                    self.telescope.ifu.lensKernel=self.makeLensKernel()                
                    self.kernel=convolve_fft(self.telescope.ifu.lensKernel, self.psfKernel)            
                elif self.inputType == ('lenscube'):
                    # if its lenscube kernel is only the psf (1,0)
                    self.kernel=self.psfKernel 
            else:
                # if psfModel is not defined then kernel is lenslet only (0,1)
                self.telescope.ifu.lensKernel=self.makeLensKernel()
                self.kernel=self.telescope.ifu.lensKernel
            self.convdata = convolve_fft(self.data, self.kernel, normalize_kernel=True)


    def getDataFluxes(self):
        """
        - Extract the spectral fluxes from data. 
        - If working on cube sample convolved cube at lenslet positions and scale fluxes by lenslet area with respect to average lenslet area
        - This version does not correct for DAR!!!!!
        - If working on an RSS file simply pass on the fluxes for each spaxel
        - Return a 2D np array of dimensions (Nspax, Nwave) in the format needed by specsim (do in one step, check scipy.interpolate.interp1d)
        -

        Main Dependent Parameters
        ----------
        'self.data' = data, sampled data or convolved data
        'self.telescope.ifu.lensx' = image pixel coordinate x  where the FOV will be centered and fluxes will be extracted. Default to center
        'self.telescope.ifu.lensy' = image pixel coordinate y  where the FOV will be centered and fluxes will be extracted. Default to center
        'self.telescope.ifu.lensr' = hexagon radius: This does not result in a convolution of the image with each indivdual PSF of each lens. That's too expensive. What is done is to scale the fluxes of an indivudal fiber by the area of the average lens size.
        potentially self.skycor and telescope model, otherwise use pixel
        """

        waveout=self.simulator.wavelength # get this wavelength from the simspec config member
        nlens=len(self.telescope.ifu.lensID)
        lensrsky=self.telescope.ifu.lensr*self.telescope.platescale(self.telescope.ifu.lensx, self.telescope.ifu.lensy)
        lensareasky=3*np.sqrt(3)*lensrsky**2/2 # lenslet area in arcsec2
        lensrpix=self.telescope.ifu.lensr*self.hdr['PIXSCALE']
        lensareapix=3*np.sqrt(3)*lensrpix**2/2 # lenslet area in number of pixels
                
        if self.inputType == ('fitsrss' or 'asciirss'):
            wavein=self.convdata[0,:]
            fluxesin=self.convdata[1:,:]            
            interp=interpolate.RectBivariateSpline(np.range(nlens), wavein, fluxesin)
            fluxesout=interp(np.range(nlens), waveout)

            if self.fluxType == 'intensity':
                """Multiply input spaxel area in arcsec2
                """
                fluxesout *= lensareasky 

        elif self.inputType == ('fitscube' or 'lenscube' or 'psfcube'):
            # compute lenslet coordinates, do mask, evaluate spectra
            # resample data to output wavelength sampling
            lensra, lensdec = self.telescope.ifu2sky(self.simparam['ra'], self.simparam['dec'], self.simparam['theta'])
            mywcs = wcs.WCS(self.hdr)
            lenscubex, lenscubey = np.array(mywcs.wcs_world2pix(lensra, lensdec, 1))

            # We will improve this with a cube of references
            fluxout=np.zeros((nlens, len(waveout)))
            for i in range(nlens):
                fluxout[i,:]=self.convdata[lenscubex[i], lenscubey[i],:]

            if self.fluxType == 'intensity':
                """Multiply input spaxel area in arcsec2
                """
                fluxout *= lensareasky

            elif self.fluxType == 'flux':
                """Multiply input by  spaxel area in pixels
                """
                fluxout *= lensareapix

        return fluxout

    def simulate(self, forceConv=False):
        
        """
        Measure fluxes for each spaxel, create the simspec Simulator object, update it with user defined parameters, and run simulation
        """

        if (self.convdata is None) or forceConv:
            # If convdata does not exist or the user wants to reconvolve the input (i.e. forceConv=True) then convolve the input
            self.convolveInput()
        self.fluxes = self.getDataFluxes() #intentionally broken, x and y are not defined
        self.updateyaml()
        self.simulator = specsim.simulator.Simulator(self.yamlfile, num_fibers=len(self.telescope.ifu.lensID))
        # TO DO: add right keywords to simultate so we pass on fluxes array
        self.simulator.simulate()


def int_rebin(a, new_shape):
    """
    Resizes a 2d array by averaging or repeating elements,
    new dimensions must be integral factors of original dimensions
    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are summed,
        if the new shape is bigger array elements are repeated and divided
        by the subsampling factor so that the total sum is conserved.
    See Also
    --------
    resize : Return a new array with the specified shape.
    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = int_rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])
    """
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,M//m,n,N//n)).sum(axis=(3,1))
    else:
        return np.repeat(np.repeat(a, m//M, axis=0), n//N, axis=1) / (m//M*n//N)
        


def main():
    """
    main: main should do something if someone calls the function form the command line. Perhaps just report the contents of each class
    """
    pass

if __name__ == '__main__':
    """Run the code"""
    main()
