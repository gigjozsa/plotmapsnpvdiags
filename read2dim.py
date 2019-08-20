#! /usr/bin/env python
# Please use python3 and the aplPY from https://github.com/gigjozsa/aplpy
#from astropy import log
#log.setLevel('ERROR')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from astropy.utils.data import download_file
#from astropy.io import fits
#from astropy import wcs
#from astropy import coordinates
#from astropy import units as u
#import os
#from pvextractor import extract_pv_slice
#from pvextractor.geometry import Path
#import matplotlib
#matplotlib.use('WxAgg')
#import aplpy
#import sys
#import numpy as np
#from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
#import math
#import copy

def read2dim(hdu,plane):
    """
    Read a plane from 3-D data, assuming it is the third axis, strip all non-3d info and return a newly created hdu.

    Input:
    hdu: astropy.io.fits HDU
    plane: number of plane, starting with 0
    """
    #hdu = fits.open(cube)[0]
    theshape = hdu.data.shape
    length = len(theshape)
    
    if length == 2:
        return hdu
    
    # We assume that the data are in that order
    if length > 3:
        newdata = hdu.data.reshape(theshape[-3],theshape[-2],theshape[-1])

    # Find out if channel number is valid, if not choose plane #0
    if hdu.header['NAXIS3'] <= plane:
        plane = 0
        
    # Assume first index to be velocity
    newdata = hdu.data[plane,:,:]

    # Now rearrange header
    hdu.header['naxis'] = 2
    while length > 2:
        try:
            hdu.header.pop('NAXIS{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CDELT{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CRVAL{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CRPIX{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CTYPE{:d}'.format(length))
        except:
            pass
        length -= 1

    # Re-generate the thing
    return(fits.PrimaryHDU(data=newdata, header=hdu.header))

