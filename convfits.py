#! /usr/bin/env python
from astropy.io import fits
import numpy as np
import scipy.signal as signal
import sys
#import matplotlib.pyplot as plt

try:
    infile = sys.argv[1]
except:
    raise('Hey, you need to give at least one argument')
try:
    kernel_x = float(sys.argv[2]) # arcsec
except:
    print('Kernel 1 is 3 arcsec')
    kernel_x = 3.
    
try:
    kernel_y = float(sys.argv[3]) # arcsec
except:
    print('Kernel 2 is 3 arcsec')
    kernel_y = 3.
    
try:
    outfile = sys.argv[4]
except:
    outfile =  '.'.join((sys.argv[1].split('.')[:-1]))+'_conv_{:d}_{:d}'.format(int(kernel_x),int(kernel_y))+'.fits'
    print('Output file is {}'.format(outfile))

hdu = fits.open(infile)[0]

theshape=hdu.data.shape
hdu.data.reshape(theshape[-2],theshape[-1])

ker_x = kernel_x/np.fabs(hdu.header['CDELT1']*3600.)
ker_y = kernel_y/np.fabs(hdu.header['CDELT2']*3600.)

kernel = np.outer(signal.gaussian(int(20.*ker_x), ker_x), signal.gaussian(int(20.*ker_y), ker_y))
convolved = signal.fftconvolve(hdu.data,kernel,mode='same')/kernel.sum()
#plt.imshow(convolved)
#plt.show()

conhdu = fits.PrimaryHDU(convolved)
for i in hdu.header.keys():
    try:
        conhdu.header[i]  = hdu.header[i]
    except:
        pass
    
conhdu.header['DATAMAX'] = convolved.min()
conhdu.header['DATAMIN'] = convolved.max()

conhdu.writeto(outfile, overwrite=True, output_verify = 'ignore')
