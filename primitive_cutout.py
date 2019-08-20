#! /usr/bin/env python
from astropy.io import fits
import astropy.wcs as wcs
import numpy as np
import scipy.signal as signal
import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
import kapteyn.wcs as kapwcs

"""
Read in infile, cut to size either given by tempfile or ramax, ramin, demax, demin (degrees, fk5, J2000). Example:
primitive_cutout.py infile=yo.fits tempfile=yot.fits outfile=yout.fits
primitive_cutout.py infile=yo.fits ramax=358.1022, ramin=357.9340 demin=-52.6305 demax=-52.5283 outfile=yout.fits
"""


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

def getkey(thedict, key):
    for i in sys.argv[1:]:
        splitted = i.split('=')
        if key == splitted[0]:
            thedict[key] = splitted[1]
    try:
        thedict[key]
    except:
        raise(BaseException)

kvdict = {}
    
try:
    getkey(kvdict,'infile')
except:
    print('Provide infile')
    raise(BaseException)

try:
    getkey(kvdict,'tempfile')

except:
    try:
        getkey(kvdict,'ramax')
        kvdict['ramax'] = float(kvdict['ramax'])
    except:
        print('Provide higher right ascension ramax')
        raise(BaseException)

    try:
        getkey(kvdict,'ramin')
        kvdict['ramin'] = float(kvdict['ramin'])
    except:
        print('Provide lower right ascension ramin')
        raise(BaseException)

    try:
        getkey(kvdict,'demax')
        kvdict['demax'] = float(kvdict['demax'])
    except:
        print('Provide higher declination demax')
        raise(BaseException)

    try:
        getkey(kvdict,'demin')
        kvdict['demin'] = float(kvdict['demin'])
    except:
        print('Provide lower declination demin')
        raise(BaseException)

try:
    getkey(kvdict,'outfile')
except:
    print('Provide output file outfile')
    raise(BaseException)

# Open hdu
hduc = fits.open(kvdict['infile'])[0]
w = wcs.WCS(hduc.header, naxis=2)
hdu = read2dim(hduc,0)
proj = kapwcs.Projection(hdu.header)
proj.skyout = (kapwcs.equatorial, kapwcs.fk5, 'J2000')

if 'tempfile' in kvdict.keys():
    # Get info
    hduinfoc = fits.open(kvdict['tempfile'])[0]
    
    # This makes the whole thing slightly unprecise. Write into the doc that it's always the first plane as a reference plane
    hduinfo = read2dim(hduinfoc,0)
    nypsinfo = hduinfo.shape[-2]
    nxpsinfo = hduinfo.shape[-1]
    #wnew = wcs.WCS(hduinfo.header, naxis=2)
    #solution = wnew.wcs_pix2world(nxps,nyps,0)
    projinfo = kapwcs.Projection(hduinfo.header)
    projinfo.skyout = (kapwcs.equatorial, kapwcs.fk5, 'J2000')
    pixelinfo = ([1, nxpsinfo],[1, nypsinfo])
    worldinfo = projinfo.toworld(pixelinfo)
    kvdict['ramin'] = float(worldinfo[0][1])
    kvdict['demax'] = float(worldinfo[1][1])
    kvdict['ramax'] = float(worldinfo[0][0])
    kvdict['demin'] = float(worldinfo[1][0])

    #print(solution)
    #solution = w.wcs_pix2world(50.,600.,0)
    #print(solution)
    #c = SkyCoord(ra=solution[0]*u.degree, dec=solution[1]*u.degree)
    #print(c.to_string('hmsdms'))

world = ([kvdict['ramax'],kvdict['ramin']],[kvdict['demin'],kvdict['demax']])
pixel = proj.topixel(world)

theshape = hdu.data.shape
sizey = theshape[-2]
sizex = theshape[-1]

if int(np.rint(pixel[0][0]))-1 < 0.:
    print('aBox out of image')
    raise(BaseException)
if int(np.rint(pixel[1][0]))-1 < 0.:
    print('bBox out of image')
    raise(BaseException)
if int(np.rint(pixel[0][1])) > sizex:
    print('cBox out of image')
    raise(BaseException)
if int(np.rint(pixel[1][1])) > sizey:
    print('dBox out of image')
    raise(BaseException)

datanew = hdu.data[int(np.rint(pixel[1][0]))-1:int(np.rint(pixel[1][1])), int(np.rint(pixel[0][0]))-1:int(np.rint(pixel[0][1]))]

newhdu = fits.PrimaryHDU(datanew)
for i in hdu.header.keys():
        try:
            newhdu.header[i]  = hdu.header[i]
        except:
            pass
        
newhdu.header['CRPIX1'] =  hdu.header['CRPIX1']-pixel[0][0]+1.
newhdu.header['CRPIX2'] =  hdu.header['CRPIX2']-pixel[1][0]+1.
newhdu.header['DATAMAX'] = datanew.min()
newhdu.header['DATAMIN'] = datanew.max()

newhdu.writeto(kvdict['outfile'], overwrite=True, output_verify = 'ignore')

#print(pixcoords_xmin, pixcoords_ymin, pixcoords_xmax, pixcoords_ymax)
#solution = w.wcs_pix2world(10.,400.,0)
#print(solution)
#solution = w.wcs_pix2world(50.,600.,0)
#print(solution)
#c = SkyCoord(ra=solution[0]*u.degree, dec=solution[1]*u.degree)
#print(c.to_string('hmsdms'))
