#! /usr/bin/env python
from astropy.io import fits
import numpy as np
import scipy.signal as signal
import sys
#import matplotlib.pyplot as plt
# Convolves with kernel (HPBW, arcsec) symmetrical Gaussian and then trims by trim arcsec, can also return only trimmed image
# example: convfits.py infile='eso149_optical_dss_r.fits' kernel=9 trim=45 outconv='con_3.fits' outtrim='tri_3.fits'

def getkey(thedict, key):
    for i in sys.argv[1:]:
        splitted = i.split('=')
        #print(splitted, key, splitted[0])
        if key == splitted[0]:
            thedict[key] = splitted[1]
    try:
        thedict[key]
    except:
        raise(BaseException)

kvdict = {}
    
try:
    getkey(kvdict,'infile')
    infile = kvdict['infile']
except:
    print('Provide infile')
    raise(BaseException)

try:
    getkey(kvdict,'kernel')
    kernel = float(kvdict['kernel'])/np.sqrt(np.log(256))
except:
    print('Kernel is 3 arcsec')
    kernel = 3./np.sqrt(np.log(256))

try:
    getkey(kvdict,'trim')
    trim = float(kvdict['trim'])
except:
    print('Image(s) will be trimmed by 3 arcsec')
    trim = 15.

try:
    getkey(kvdict,'outconv')
    outcfile = kvdict['outconv']
except:
    outcfile = ''

    #'.'.join((sys.argv[1].split('.')[:-1]))+'_conv_{:d}_{:d}'.format(int(kernel_x),int(kernel_y))+'.fits'
    print('No output convolved file')

try:
    getkey(kvdict,'outtrim')
    outtfile = kvdict['outtrim']
except:
    outtfile = ''
    #'.'.join((sys.argv[1].split('.')[:-1]))+'_conv_{:d}_{:d}'.format(int(kernel_x),int(kernel_y))+'.fits'
    print('No output trimmed file')

hdu = fits.open(infile)[0]

theshape=hdu.data.shape
sizey = theshape[-2]
sizex = theshape[-1]

hdu.data.reshape(theshape[-2],theshape[-1])

ker_x = kernel/np.fabs(hdu.header['CDELT1']*3600.)
ker_y = kernel/np.fabs(hdu.header['CDELT2']*3600.)
tri_x = int(np.ceil(trim/np.fabs(hdu.header['CDELT1']*3600.)))
tri_y = int(np.ceil(trim/np.fabs(hdu.header['CDELT2']*3600.)))

kernel = np.outer(signal.gaussian(int(20.*ker_y), ker_y), signal.gaussian(int(20.*ker_x), ker_x))
convolved = (signal.fftconvolve(hdu.data,kernel,mode='same')/kernel.sum())[tri_y:sizey-tri_y,tri_x:sizex-tri_x]
trimmed = hdu.data[tri_y:sizey-tri_y,tri_x:sizex-tri_x]
#plt.imshow(convolved)
#plt.show()

conhdu = fits.PrimaryHDU(convolved)
trimhdu = fits.PrimaryHDU(trimmed)

for thehdu in [conhdu, trimhdu]:
    for i in hdu.header.keys():
        try:
            thehdu.header[i]  = hdu.header[i]
        except:
            pass

    thehdu.header['CRPIX1'] =  thehdu.header['CRPIX1']-tri_x
    thehdu.header['CRPIX2'] =  thehdu.header['CRPIX2']-tri_y
    thehdu.header['DATAMAX'] = thehdu.data.min()
    thehdu.header['DATAMIN'] = thehdu.data.max()

if outcfile != '':
    conhdu.writeto(outcfile, overwrite=True, output_verify = 'ignore')
if outtfile != '':
    trimhdu.writeto(outtfile, overwrite=True, output_verify = 'ignore')
