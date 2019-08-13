#! /usr/bin/env python
from astropy.io import fits
import numpy as np
import scipy.signal as signal
import sys
#import matplotlib.pyplot as plt

# Calculate total flux and HI mass from a fitsfile input and a distance from the command line

try:
    infile = sys.argv[1]
except:
    raise('Hey, you need to give at least one argument')

try:
    distance = float(sys.argv[2])
except:
    raise('Hey, you need to give at least one argument')

hdu = fits.open(infile)[0]

bmaj = float(hdu.header['BMAJ'])
bmin = float(hdu.header['BMIN'])
psize = float(hdu.header['CDELT2'])
psum = hdu.data.sum()
flux = psum*psize*psize/(1.13309*bmaj*bmin)
mass = 2.36E5*distance*distance*flux
logmass = np.log(mass)/np.log(10.)
print('{:41s}:    Total flux: {:.2f} Jy km/s    HI mass:{:.2e} solar Masses    Logmass: {:.1f}'.format(infile,flux,mass,logmass))
