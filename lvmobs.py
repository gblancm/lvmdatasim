import numpy as np
import scipy.ndimage
import scipy.signal
import scipy.interpolate
import scipy.fftpack
from astropy.io import fits
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

import hexagonlib as hexlib

def lvm_sample_image (fimg, platescale, lenslet_size=36.2, reconstruct=False, debug=False):
    """
    Sample an image to simulating observing with LVM.

    Inputs:
    fimg          fits image to observe with LVM
    platescale    fits image plate scale in arcsec/pix
    lenslet_size  size of the LVM lenslet elements. 36.2 by default (LVM160). Use 6.8 for LVM1000.
                  note that to properly simulate the 1m observations you will have to make sure 
                  the input images are properly sampled. 1" or better is required. The LVM160 simulations
                  can use a coarser input image scale.
    reconstruct   return list of samples (false) or reconstructed image (true)
    debug         show intermediate plots for illustration (false by default)

    Return:
    x,y,f         x and y positions and values f of the image sampled at the dither locations
                  OR
    img           reconstructed image from the samples if reconstruct==True (simple interpolation)
    """
    
    img, hdr = fits.getdata(fimg, 0, header=True)
    hex_diameter=lenslet_size

    # convolve image with lenslet/fiber footprint
    fiber = make_hex_kernel(hex_diameter/2.0/platescale, 21, debug=False)    
    imgf = convolve_with_kernel(img, fiber)
    if debug==True:
        plt.imshow(fiber)
        plt.show()

    # generate the IFU and the locations of the dithers
    ifu = mkifu(30)  # make one large IFU for now. 30 Rings. Later simulate proper tiling
    layout = hexlib.Layout(hexlib.layout_flat, hex_diameter/2.0,
                           hexlib.Point(img.shape[0]/2*platescale,img.shape[1]/2*platescale))
    dither = mkdither(layout,ifu,debug=False)

    # transform all the fiber centers to X,Y from hex coordinates
    ifuc = cubetoxy(ifu,layout)
    ifux = np.array([i.x for i in ifuc])
    ifuy = np.array([i.y for i in ifuc])
    if debug==True:
        plt.imshow(imgf,clim=(0.4, 7.0))
        for i in dither:
            plt.scatter((ifux+i[0])/platescale, (ifuy+i[1])/platescale, s=1)
        plt.show()

    # sample the convolved image in the IFU+dithers
    x=[]
    y=[]
    q = scipy.interpolate.interp2d( range(imgf.shape[1]), range(imgf.shape[0]), imgf, kind='cubic' )
    for i in dither:
        x = np.append(x, (ifux+i[0])/platescale)
        y = np.append(y, (ifuy+i[1])/platescale)
    f = []
    for i in range(len(x)):
        f = np.append(f, q(x[i], y[i]))

    if reconstruct==True:
        # attempt at reconstruction using the samples
        X, Y = np.meshgrid(range(imgf.shape[1]), range(imgf.shape[0]))
        imgr = scipy.interpolate.griddata( (x,y), f, (X,Y), method='cubic')
        inf = np.invert(np.isfinite(imgr))
        imgr[inf] = 0.
        return imgr
    else:
        return x,y,f


def mkdither(layout,ifu,debug=False):
    """
    Create and return a list of 2-vectors with the dither offsets for a 9-point dither.
    The dither positions are along the edges of a equilateral triangle spaced at 0.5*hexagon radius

    *** ASSUMES A FLAT-TOP LAYOUT for now ***
    """
    radius = layout.size  # center of hex to vertex
    zero = -np.array([0.0,0.0])

    # basis vectors for an equilateral triangle
    theta = np.radians(60)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    b1 = np.array([1.0,0.0])
    b2 = np.matmul(R,b1)
    c1 = zero

    # first 3 points: center and two neighboring vertices of the hexagon
    dither = [c1, c1+radius*b1, c1+radius*b2]

    # remaining 6 in smaller hexagon
    b3 = np.array([0.0,1.0])
    h = np.sqrt(3.0)/3.0*radius
    dither.append(h*b3)
    for i in range(5):
        dither.append(np.matmul(R,dither[-1]))

    if debug==True:
        plt.scatter(*zip(*dither))
        plt.show()
    return dither


def mkifu(N):
    """
    An IFU of N rings above the central fiber is returned as a list of hexlib.Hex objects.
    Coordinates are in Cube coordinates (3 numbers per hexagon).
    """
    ifu = []
    for dx in range(-N,N+1):
        for dy in range(max(-N, -dx-N), min(N, -dx+N)+1):
            dz = -dx-dy
            ifu.append(hexlib.Hex(dx, dy, dz))
    return ifu


def equal_hex(a,b):
    """
    predicate to test equality between two Hex structures
    """
    return (a.q == b.q and a.s == b.s and a.r == b.r)


def find_hex(ifu,hex):
    """
    return indices i of ifu where ifu[i]==hex
    """
    return [i for i, j in enumerate(ifu) if equal_hex(j, hex)]


def cubetoxy(ifu,layout):
    """
    Convert from Cube coordinates to XY coordinates of the centers of the hexagons in the plane.
    """
    return [hexlib.hex_to_pixel(layout, x)  for x in ifu]


def make_hex_kernel(radius,imgsize,antialias=5,debug=False):
    """
    Generate a 2D numpy array with the characteristic function of a hexagon of size 'radius' (center to corner)
    in units of pixels of the array. The image will have dimensions of imgsize x imgsize and the hexagon will
    be at the integer center.
    """
    antialias=5
    imgsize*=antialias
    center = imgsize//2+1
    pointy = hexlib.Layout(hexlib.layout_pointy, radius*antialias, hexlib.Point(center,center))
    polygon = hexlib.polygon_corners(pointy,hexlib.Hex(0,0,0))
    img = Image.new('L', (imgsize, imgsize), 0)
    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
    mask = np.array(img,dtype=float)
    mask /= np.sum(mask)
    mask = int_rebin(mask, (imgsize//antialias,imgsize//antialias))
    if debug==True:
        plt.imshow(mask)
        plt.show()
    return mask


def make_gaussian_kernel(sig,imagesize,debug=False):
    """
    return normalized Gaussian kernel with side length imagesize and a sigma of sig
    """
    l = imagesize
    ax = np.arange(-l // 2 + 1., l // 2 + 1.)
    xx, yy = np.meshgrid(ax, ax)
    kernel = np.exp(-(xx**2 + yy**2) / (2. * sig**2))
    if debug==True:
        plt.imshow(kernel)
        plt.show()
    return kernel / np.sum(kernel)


def convolve_with_kernel(img, kernel, subsample=1):
    """
    Convolve an image with a kernel representing the fiber/lenslet footprint on the
    sky. The image and kernel are optionally subsampled by a given factor for improved
    accuracy before the convolution.
    """
    if subsample>1:
        img = (img.repeat(subsample, axis=0).repeat(subsample, axis=1)) / (subsample*subsample)
        kernel = (kernel.repeat(subsample, axis=0).repeat(subsample, axis=1)) / (subsample*subsample)

    out = scipy.signal.convolve2d(img, kernel, mode='same', boundary='fill', fillvalue=0)
    if subsample>1:
        out.reshape(a.shape[0]//subsample, subsample, a.shape[1]//subsample, subsample).sum(axis=(0,1))

    return out


def integrate_flux(ifu, fibersize, image, platescale):
    """
    iterate over pixels in the 2d numpy array image and identify the hexagon the pixel is
    in in the list of Hex structures ifu. Return a list of floats with the summed fluxes
    in each hexagon.
    """
    pointy = hexlib.Layout(hexlib.layout_pointy, fibersize, hexlib.Point(0,0))
    fluxes = np.zeros(len(ifu), dtype=np.float)
    for (x,y), flux in np.ndenumerate(image):
        p = Point(x*platescale,y*platescale)
        hex = pixel_to_hex(pointy, p)
        i = find_hex(ifu, hex)
        if i>=0:
            flux[i] += image(x,y)

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
