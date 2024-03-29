#! /usr/bin/env python
# Please use python3 and the aplPY from https://github.com/gigjozsa/aplpy
from astropy import log
log.setLevel('ERROR')

import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.utils.data import download_file
from astropy.io import fits
from astropy import wcs
from astropy import coordinates
from astropy import units as u
import os
from pvextractor import extract_pv_slice
from pvextractor.geometry import Path
import matplotlib
#matplotlib.use('WxAgg')
import aplpy
import sys
import numpy as np
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
import math
import copy

# A cm in inch
cm_in_inch = 0.3937008

# A4width_in_cm
A4widht_in_cm = np.power(2,1./4.)*100*np.power(2,-(4.+1)/2.)

# A4widht_in_inch
A4widht_in_inch = cm_in_inch*A4widht_in_cm

class Slastroplot:
    """
    Serial astronomical plotting
    """
    def __init__(self):
        """
        Initialize any potential variable to either plotmaps or plotpvs
        """
        # Layout re lengths
        self.figsize = (A4widht_in_inch,A4widht_in_inch*np.sqrt(2.)) # only for plotpvdiagrams
        self.width = A4widht_in_inch
        self.plotmargin = 0.5*cm_in_inch
        self.labelmargin_left = 3.5*cm_in_inch
        self.labelmargin_bottom = 1.75*cm_in_inch
        self.labelmargin_right = 0.

        self.frametickcolour = 'black'
        self.frameticklinewidth = 2
        self.frameticklength = 5

        self.annotationfontsize = 'small'
        
        self.suppress_xlab = False
        self.suppress_ylab = False
        self.x_major_formatter = None
        self.x_ticksn = None
        self.hide_ticl = False
        self.hide_ticl_left = True
        self.hide_ticl_min = 0.
        self.hide_ticl_max = 0.
        self.hide_ticl_top = 0.
        self.hide_ticl_bottom = 0.
        
        # Colourbar and lengths
        self.colourbar_width = 0.2
        self.colourbar_pad = 0.05
        self.colourbar_label = None

        # Number of plots
        self.nx = 0
        self.ny = 0
        
        # Background images and their design
        self.basemap = None
        self.vmin = None
        self.vmax= None
        self.vmid = None
        self.exponent = None
        self.pmin = 0.25
        self.pmax = 99.75
        self.stretch = 'linear'
        self.cmap = None
        self.invert = False
        self.colourbar = False
        self.colourbar_label = None

        self.chans = None

        self.individual_maps = False # Indicating if contours and pixelmap formats are for each map/plane or for all
        
        self.contoursets = None
        self.contourlevs = None
        self.contourcols = None
        self.contourstyles = None
        self.contouralphas = None
        self.contourlinewidths = None

        self.olays = None
        self.individual_olays = False
        
        self.showscale = None # 4000.*u.pc
        self.distance = 10
        self.scale_borderpad = 1
        self.scale_sep = 10
        self.scalebarlinewidth = 1
        self.showscalelab = True
        self.reducescalelab = True
        self.scaleunits = 'arcsecs' # alternatively arcmin
        self.scalebarfontsize = self.annotationfontsize
        self.scalebarcolor = 'black'
        self.scale_ffc = None
        self.scale_ffa = 1.0
        self.scale_fec = None
        self.scale_fea = 1.0

        self.showbeam = True
        self.beamfc = '0.7'
        self.beamfa = 0.5
        self.beamec = 'black'
        self.beamea = 1.

        self.plotvelo = None
        self.velofontsize = self.annotationfontsize
        self.veloposx = 0.25
        self.veloposy = 0.9
        self.velbbfc = None
        self.velbbec = None
        self.velbbfa = 1
        self.velbbea = 1

        self.aspect = 'equal' # axis have same unit lengths, alternative is 'auto', meaning fill the are, see matplotlib ax.set_aspect


    def set_standard(self, border = None, background = None):
        """
        Delivers standard setup for given parameters
        """
        
        # General border
        if border != None:
            if border == True:
                self.plotmargin = 0.35*cm_in_inch
            else:
                self.plotmargin = 0.

        # General background
        if background == 'bright':
            self.scale_ffc = None
            self.scale_ffa = 0.7
            self.scale_fec = '0.7'
            self.scale_fea = 1.

            self.velbbfc = self.scale_ffc
            self.velbbec = self.scale_fec
            self.velbbfa = self.scale_ffa
            self.velbbea = self.scale_fea
        elif background == 'dark':
            self.scale_ffc = 'white'
            self.scale_ffa = 0.7
            self.scale_fec = '0.7'
            self.scale_fea = 1.

            self.velbbfc = self.scale_ffc
            self.velbbec = self.scale_fec
            self.velbbfa = self.scale_ffa
            self.velbbea = self.scale_fea

                
        if self.nx == 1:
            self.labelmargin_left = 3.0*cm_in_inch # default/2
            self.labelmargin_bottom = 1.7*cm_in_inch # default/2

            self.annotationfontsize = 'xx-large'

            self.frameticklinewidth = 2
            self.frameticklength = 7.5
            
            # Colourbar and lengths
            if self.colourbar:
                self.colourbar_width = 0.25
                self.colourbar_pad = 0.05
                self.veloposx = 0.145
                self.veloposy = 0.95
                self.labelmargin_top = 0.25*cm_in_inch

                if self.colourbar_label != None:
                    self.labelmargin_right = 2.35*cm_in_inch
                else:
                    self.labelmargin_right = 1.75*cm_in_inch
            else:
                self.veloposx = 0.1325
                self.veloposy = 0.95
                self.labelmargin_top = 0.
                self.labelmargin_right = 0.

            self.scale_borderpad = 1.25
            self.scale_sep = 7.5
            self.scalebarlinewidth = 2
            self.scalebarfontsize = self.annotationfontsize

            self.velofontsize = self.annotationfontsize
        elif self.nx == 2:
            self.labelmargin_left = 2.2*cm_in_inch # default/2
            self.labelmargin_bottom = 1.2*cm_in_inch # default/2

            self.annotationfontsize = 'large'

            self.frameticklinewidth = 2
            self.frameticklength = 6.5
            
            # Colourbar and lengths
            if self.colourbar:
                self.colourbar_width = 0.175
                self.colourbar_pad = 0.05
                self.veloposx = 0.2
                self.veloposy = 0.925
                self.labelmargin_top = 0.125*cm_in_inch

                if self.colourbar_label != None:
                    self.labelmargin_right = 2.3*cm_in_inch
                else:
                    self.labelmargin_right = 1.75*cm_in_inch
            else:
                self.veloposx = 0.175
                self.veloposy = 0.935
                self.labelmargin_top = 0.
                self.labelmargin_right = 0.

            self.scale_borderpad = 0.85
            self.scale_sep = 6.5
            self.scalebarlinewidth = 2
            self.scalebarfontsize = self.annotationfontsize

            self.velofontsize = self.annotationfontsize
        elif self.nx == 3:
            self.labelmargin_left = 1.9*cm_in_inch # default/2
            self.labelmargin_bottom = 1.025*cm_in_inch # default/2

            self.annotationfontsize = 'medium'

            self.frameticklinewidth = 2
            self.frameticklength = 6
            
            # Colourbar and lengths
            if self.colourbar:
                self.colourbar_width = 0.15
                self.colourbar_pad = 0.05
                self.veloposx = 0.26
                self.veloposy = 0.9
                self.labelmargin_top = 0.135*cm_in_inch

                if self.colourbar_label != None:
                    self.labelmargin_right = 1.9*cm_in_inch
                else:
                    self.labelmargin_right = 1.3*cm_in_inch
            else:
                self.veloposx = 0.235
                self.veloposy = 0.91
                self.labelmargin_top = 0.
                self.labelmargin_right = 0.

            self.scale_borderpad = 0.75
            self.scale_sep = 6
            self.scalebarlinewidth = 2
            self.scalebarfontsize = self.annotationfontsize

            self.velofontsize = self.annotationfontsize
        elif self.nx == 4:
            self.labelmargin_left = 1.75*cm_in_inch # default/2
            self.labelmargin_bottom = 0.875*cm_in_inch # default/2

            self.annotationfontsize = 'small'

            self.frameticklinewidth = 2
            self.frameticklength = 5
            
            # Colourbar and lengths
            if self.colourbar:
                self.colourbar_width = 0.1
                self.colourbar_pad = 0.05
                self.veloposx = 0.3
                self.veloposy = 0.89
                self.labelmargin_top = 0.125*cm_in_inch

                if self.colourbar_label != None:
                    self.labelmargin_right = 1.65*cm_in_inch
                else:
                    self.labelmargin_right = 1.25*cm_in_inch
            else:
                self.veloposx = 0.3
                self.veloposy = 0.9
                self.labelmargin_top = 0.
                self.labelmargin_right = 0.

            self.scale_borderpad = 0.7
            self.scale_sep = 5
            self.scalebarlinewidth = 2
            self.scalebarfontsize = self.annotationfontsize

            self.velofontsize = self.annotationfontsize
        else:
            self.labelmargin_left = 1.75*cm_in_inch # default/2
            self.labelmargin_bottom = 0.875*cm_in_inch # default/2

            self.annotationfontsize = 'small'

            self.frameticklinewidth = 2
            self.frameticklength = 5
            self.hide_ticl_top = -0.02
            self.hide_ticl_bottom = -0.1
           
            # Colourbar and lengths
            if self.colourbar:
                self.colourbar_width = 0.1
                self.colourbar_pad = 0.05
                self.veloposx = 0.3
                self.veloposy = 0.89
                self.labelmargin_top = 0.125*cm_in_inch

                if self.colourbar_label != None:
                    self.labelmargin_right = 1.65*cm_in_inch
                else:
                    self.labelmargin_right = 1.25*cm_in_inch
            else:
                self.veloposx = 0.35
                self.veloposy = 0.87
                self.labelmargin_top = 0.
                self.labelmargin_right = 0.

            self.scale_borderpad = 0.5
            self.scale_sep = 3
            self.scalebarlinewidth = 2
            self.scalebarfontsize = self.annotationfontsize

            self.velofontsize = self.annotationfontsize
        
        return
    
def setifnotthere(theinstance, thevariable, thevalue):
    if thevariable in vars(theinstance).keys():
        pass
    else:
        setattr(theinstance, thevariable, thevalue)
        
def getdefaults(nrows = 1, colorbar = False, showbeam = None, showscale = None, brightbackground = True):
    """
    Return a default kwarg dict for nrows rows
    """
    retdict = {}
    figsize = (A4widht_in_inch, A4widht_in_inch)
    # asf
    return {}

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

def putcontours(F, hdu, i, contourlevs, contourcols = None, contouralphas = None, contourstyles = None, contourlinewidths = None, usepixel = False):
    """
    Helper to conturcalls sorting out contourlevels
    F is a aplpy graph, i is the contour-number, contourlevs the contour level (can be a list for each i one or a list of lists for each i a set of aontourlevels), contourcols the colors (can be None, a list, or a list of lists), contourstyles the contour-styles (again, None, a list, or a list of lists)). If usepixel == True, draws contours only using pixel coordinates, otherwise maps of different sizes can be overlaid.
    """

    if type(contourlevs[0]) == type([]):
        if i < len(contourlevs):
            levelshere = contourlevs[i]
        else:
            levelshere = contourlevs[-1]
    else:
        levelshere = contourlevs


    if type(contourcols) == type([]):
        if i < len(contourcols):
            colorshere = contourcols[i]
        else:
            colorshere = contourcols[-1]
    else:
        colorshere = contourcols

    if type(contouralphas) == type([]):
        if i < len(contouralphas):
            alphashere = contouralphas[i]
        else:
            alphashere = contouralphas[-1]
    else:
        alphashere = contouralphas

    if type(contourstyles) == type([]):
        if i < len(contourstyles):
            styleshere = contourstyles[i]
        else:
            styleshere = contourstyles[-1]
    else:
        styleshere = contourstyles

    if type(contourlinewidths) == type([]):
        if i < len(contourlinewidths):
            contourlwshere = contourlinewidths[i]
        else:
            contourlwshere = contourlinewidths[-1]
    else:
        contourlwshere = contourlinewidths
        
    F.show_contour(data=hdu, levels = levelshere, colors = colorshere, alpha = alphashere, linewidths = contourlwshere, linestyles = styleshere, usepixel = usepixel)
    return


def pvdiagrams(filename, outname, centra, centdec, pa, majlen, minlen, mindist):
    """
    Produce four PV diagrams from filename

    Input:
    filename (str)     : Input file name
    outname (str)      : Prefix to output file names
    centra (float)     : Right Ascension of the centre in degrees
    centdec (float)    : Declination of the centre in degrees
    pa (float)         : Position angle, N over E, of the receding side major axis in degrees
    majlen (float)     : Length of cut along major axis in arcsec
    minlen (float)     : Length of cut along minor axis and along parallel lines in arcsec
    mindist (float)    : Distance of cuts along parallel lines from cut along minor axis

    Will produce four pv diagrams from a cube filename, one along the
    "major" axis, defined by the central coordinates centra (deg),
    centdec (deg), position angle (pa in degrees), and length majlen
    in arcsec, three perpendicular with respect to the major axis, one
    through the centre, two (at each side) at a distance of mindist
    (in arcsec) from the centre, of the length minlen. Output PV
    diagrams will be called outname+'_pvmaj.fits',
    outname+'_pvmin_cent.fits', outname+'_pvmin_left.fits',
    outname+'_pvmin_right.fits'. Position angle should be defined as
    the receding side major axis.

    """
    
    print('Creating Pvdiagrams')
#    return
    fitsfile = fits.open(filename)
    hdu = fitsfile[0]
    w = wcs.WCS(hdu.header, fitsfile)

#    fitsfile_mom = fits.open(overplotname)
#    hdu_mom = fitsfile_mom[0]

#hdulist = fits.open(filename)
#w = wcs.WCS(hdulist[('sci',1)].header, hdulist)

    #    print hdu.data.shape
    
    # find a list of endpoints
    cenco = np.array([[centra,centdec]])
    central_pix_x, central_pix_y = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix(centra,centdec, 0)
    scale_xy = hdu.header['cdelt2']*3600. # convert arcsec to pixel

    # Make a grid
    endpoints_y_b = np.array([majlen/2., -majlen/2.,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_b = np.array([       0.,         0., -minlen/2., minlen/2., -minlen/2., minlen/2., -minlen/2.-0.0001, minlen/2.+0.0001])
    endpoints_y_bt = np.array([majlen/2.+2.5, -majlen/2.+2.5,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_bt = np.array([       0.,         0., -minlen/2.-2.5, minlen/2.-2.5, -minlen/2.-2.5, minlen/2.-2.5, -minlen/2.-0.0001-2.5, minlen/2.+0.0001-2.5])

    # Rotate
    endpoints_x = endpoints_x_b*np.cos(np.pi*pa/180.)-endpoints_y_b*np.sin(np.pi*pa/180.)+central_pix_x
    endpoints_y = endpoints_x_b*np.sin(np.pi*pa/180.)+endpoints_y_b*np.cos(np.pi*pa/180.)+central_pix_y
    endpoints_xt = endpoints_x_bt*np.cos(np.pi*pa/180.)-endpoints_y_bt*np.sin(np.pi*pa/180.)+central_pix_x
    endpoints_yt = endpoints_x_bt*np.sin(np.pi*pa/180.)+endpoints_y_bt*np.cos(np.pi*pa/180.)+central_pix_y
    
    i = -1
    for names in ['_pvmaj.fits', '_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']:
        i = i + 1
        endpoints = [(endpoints_x[2*i],endpoints_y[2*i]),(endpoints_x[2*i+1],endpoints_y[2*i+1])]
        endpointst = [(endpoints_xt[2*i],endpoints_yt[2*i]),(endpoints_xt[2*i+1],endpoints_yt[2*i+1])]

        xy = Path(endpoints)
        pv = extract_pv_slice(hdu, xy)
    
        header = pv.header
        #    print header
        pixels =  header['NAXIS1']
        pv.header['CRPIX1'] = pixels/2
        pv.header['CDELT1'] = pv.header['CDELT1']
        pv.writeto(outname+names, clobber = True)
    
    return

def plotpvdiagrampositions(filename, centra, centdec, pa, majlen, minlen, mindist, arrowcolour, plotname = None, **kwargs):
    """
    OLD STUFF, CAN BE PHASED OUT
    Will draw lines onto an image of overplotname to indicate position of PV-diagrams

    Input:
    filename (str)     : Input file name
    centra (float)     : Right Ascension of the centre in degrees
    centdec (float)    : Declination of the centre in degrees
    pa (float)         : Position angle, N over E, of the receding side major axis in degrees
    majlen (float)     : Length of cut along major axis in arcsec
    minlen (float)     : Length of cut along minor axis and along parallel lines in arcsec
    mindist (float)    : Distance of cuts along parallel lines from cut along minor axis
    arrowcolour (float): Colour of the arrows
    plotname (str)     : Name of plot
    **kwargs           : Additional arguments, see plotmaps_prep

    Will draw lines onto an image of overplotname one along the
    "major" axis, defined by the central coordinates centra (deg),
    centdec (deg), position angle (pa in degrees), and length majlen
    in arcsec, three perpendicular with respect to the major axis, one
    through the centre, two (at each side) at a distance of mindist
    (in arcsec) from the centre, of the length minlen. Position angle
    should be defined as the receding side major axis.

    """
    print('Plotpvdiagrampositions')
    plt.rc('text', usetex=True)

    fitsfile = fits.open(filename)
    hdu = fitsfile[0]
    w = wcs.WCS(hdu.header, fitsfile)

    # find a list of endpoints
    cenco = np.array([[centra,centdec]])
    central_pix_x, central_pix_y = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix(centra,centdec, 0)
    scale_xy = hdu.header['cdelt2']*3600. # convert arcsec to pixel

    # Make a grid
    endpoints_y_b = np.array([majlen/2., -majlen/2.,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_b = np.array([       0.,         0., -minlen/2., minlen/2., -minlen/2., minlen/2., -minlen/2.-0.0001, minlen/2.+0.0001])
    endpoints_y_bt = np.array([majlen/2.+2.5, -majlen/2.+2.5,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_bt = np.array([       0.,         0., -minlen/2.-2.5, minlen/2.-2.5, -minlen/2.-2.5, minlen/2.-2.5, -minlen/2.-0.0001-2.5, minlen/2.+0.0001-2.5])

    # Rotate
    endpoints_x = endpoints_x_b*np.cos(np.pi*pa/180.)-endpoints_y_b*np.sin(np.pi*pa/180.)+central_pix_x
    endpoints_y = endpoints_x_b*np.sin(np.pi*pa/180.)+endpoints_y_b*np.cos(np.pi*pa/180.)+central_pix_y
    endpoints_xt = endpoints_x_bt*np.cos(np.pi*pa/180.)-endpoints_y_bt*np.sin(np.pi*pa/180.)+central_pix_x
    endpoints_yt = endpoints_x_bt*np.sin(np.pi*pa/180.)+endpoints_y_bt*np.cos(np.pi*pa/180.)+central_pix_y
    
    # Start mom0 plot
    F = plotmaps_prep(**kwargs)

    i = -1
    letters = ['A','B','C','D']
    for names in ['_pvmaj.fits', '_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']:
        i = i + 1
        endpoints = [(endpoints_x[2*i],endpoints_y[2*i]),(endpoints_x[2*i+1],endpoints_y[2*i+1])]
        endpointst = [(endpoints_xt[2*i],endpoints_yt[2*i]),(endpoints_xt[2*i+1],endpoints_yt[2*i+1])]
    
        endpoints_wcs  = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[x,y] for x,y in endpoints], 0)
        endpoints_wcst = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[x,y] for x,y in endpointst], 0)

        xs = endpoints_wcs.T[0,0]
        ys = endpoints_wcs.T[1,0]
        exs = endpoints_wcs.T[0,1]-endpoints_wcs.T[0,0]
        eys = endpoints_wcs.T[1,1]-endpoints_wcs.T[1,0]
        F.show_arrows(xs,ys,exs,eys, width = 0.2, head_width = 1.5, head_length = 1.5, length_includes_head = True, color = arrowcolour)
        
        F.add_label(endpoints_wcst.T[0,0], endpoints_wcst.T[1,0], letters[i], relative = False, size=24, weight=400, horizontalalignment='center', verticalalignment='center', color = arrowcolour)
    if plotname != None:
        F.save(plotname)
    F.close()

def getpvdiagramolay(filename, centra, centdec, pa, majlen, minlen, mindist, arrow_colour = 'red', arrow_width = 0.2, arrow_head_width = 1.5, arrow_head_length = 1.5, arrow_length_includes_head = True, labels = ['A','B','C','D'], label_size=24, label_weight=400, label_horizontalalignment='center', label_verticalalignment='center', label_dist = 2.5):
    """
    Will draw lines onto an image of overplotname to indicate position of PV-diagrams

    Input:
    filename (str)                   : Input file name
    centra (float)                   : Right Ascension of the centre in degrees
    centdec (float)                  : Declination of the centre in degrees
    pa (float)                       : Position angle, N over E, of the receding side major axis in degrees
    majlen (float)                   : Length of cut along major axis in arcsec
    minlen (float)                   : Length of cut along minor axis and along parallel lines in arcsec
    mindist (float)                  : Distance of cuts along parallel lines from cut along minor axis
    arrow_colour (float)             : Colour of the arrows
    arrow_width (float)              : Width of arrows (see parameter width for aplpy.show_arrows)
    arrow_head_width (float)         : Head width of arrows (see parameter head_width for aplpy.show_arrows)
    arrow_head_length (float):       : Head length of arrows (see parameter head_length for aplpy.show_arrows)
    arrow_length_includes_head (bool): See parameter length_includes_head for aplpy.show_arrows)
    labels (list of str)             : Names of labels to be put on arrows
    label_size (float)               : Size of labels
    label_weight (float)             : Weight of labels 
    label_horizontalalignment (str)  : Horizontal alignment of labels
    label_verticalalignment (str)    : Vertical alignment of labels

    Will draw lines onto an image of overplotname one along the
    "major" axis, defined by the central coordinates centra (deg),
    centdec (deg), position angle (pa in degrees), and length majlen
    in arcsec, three perpendicular with respect to the major axis, one
    through the centre, two (at each side) at a distance of mindist
    (in arcsec) from the centre, of the length minlen. Position angle
    should be defined as the receding side major axis.

    """
    print('Getpvdiagramolay')
    #plt.rc('text', usetex=True)

    fitsfile = fits.open(filename)
    hdu = fitsfile[0]
    w = wcs.WCS(hdu.header, fitsfile)

    # find a list of endpoints
    cenco = np.array([[centra,centdec]])
    central_pix_x, central_pix_y = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix(centra,centdec, 0)
    scale_xy = hdu.header['cdelt2']*3600. # convert arcsec to pixel

    # Make a grid
    endpoints_y_b = np.array([majlen/2., -majlen/2.,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_b = np.array([       0.,         0., -minlen/2., minlen/2., -minlen/2., minlen/2., -minlen/2.-0.0001, minlen/2.+0.0001])
    endpoints_y_bt = np.array([majlen/2.+label_dist, -majlen/2.+label_dist,  mindist  ,  mindist ,         0.,        0.,   -mindist,  -mindist])
    endpoints_x_bt = np.array([       0.,         0., -minlen/2.-label_dist, minlen/2.-label_dist, -minlen/2.-label_dist, minlen/2.-label_dist, -minlen/2.-0.0001-label_dist, minlen/2.+0.0001-label_dist])

    # Rotate
    endpoints_x = endpoints_x_b*np.cos(np.pi*pa/180.)-endpoints_y_b*np.sin(np.pi*pa/180.)+central_pix_x
    endpoints_y = endpoints_x_b*np.sin(np.pi*pa/180.)+endpoints_y_b*np.cos(np.pi*pa/180.)+central_pix_y
    endpoints_xt = endpoints_x_bt*np.cos(np.pi*pa/180.)-endpoints_y_bt*np.sin(np.pi*pa/180.)+central_pix_x
    endpoints_yt = endpoints_x_bt*np.sin(np.pi*pa/180.)+endpoints_y_bt*np.cos(np.pi*pa/180.)+central_pix_y
    
    # Start mom0 plot
    #F = plotmaps_prep(**kwargs)

    i = -1
#    letters = labels
    returnolays = []
    for names in ['_pvmaj.fits', '_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']:
        i = i + 1
        endpoints = [(endpoints_x[2*i],endpoints_y[2*i]),(endpoints_x[2*i+1],endpoints_y[2*i+1])]
        endpointst = [(endpoints_xt[2*i],endpoints_yt[2*i]),(endpoints_xt[2*i+1],endpoints_yt[2*i+1])]
    
        endpoints_wcs  = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[x,y] for x,y in endpoints], 0)
        endpoints_wcst = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[x,y] for x,y in endpointst], 0)

        xs = endpoints_wcs.T[0,0]
        ys = endpoints_wcs.T[1,0]
        exs = endpoints_wcs.T[0,1]-endpoints_wcs.T[0,0]
        eys = endpoints_wcs.T[1,1]-endpoints_wcs.T[1,0]
        returnolays.append(['arrows', {'pv_x_world': [xs], 'pv_y_world': [ys], 'pv_dx': exs, 'pv_dy': eys, 'width' : arrow_width, 'head_width' : arrow_head_width, 'head_length' :  arrow_head_length, 'length_includes_head' : arrow_length_includes_head, 'color' : arrow_colour}])
        returnolays.append(['label',{'pv_x_world': endpoints_wcst.T[0,0], 'pv_y_world': endpoints_wcst.T[1,0], 'pv_label' : labels[i], 'relative': False, 'size': label_size, 'weight': label_weight, 'horizontalalignment': label_horizontalalignment, 'verticalalignment': label_verticalalignment, 'color': arrow_colour}])
    return returnolays

    
def plotprep_general(F, vmin = None, vmax= None, pmin = None, pmax = None, stretch = 'linear', cmap = None, vmid = None, exponent = None, invert = False, colourbar = False, colourbar_label = None, colourbar_width = 0.2, colourbar_pad = 0.05, frametickcolour = 'black', frameticklinewidth = None, frameticklength = None, showscale = None, showscalelab = True, distance = 10, scaleunits = 'arcsec',  scale_borderpad = 1, scale_ffc = None, scale_ffa = 1.0, scale_fec = None, scale_fea = 1.0, scale_sep = 10, annotationfontsize = 'xx-large', suppress_xlab = False, suppress_ylab = False, scalebarlinewidth = 1, scalebarfontsize = None, scalebarcolor = 'black', aspect='equal', offsetunits = 'arcsec'):
    """
    Prepare a plot, called by any plotting routine

    Input:
    vmin (float)            : Minimum for colourscale, see aplpy
    vmax (float)            : Maximum for colourscale, see aplpy
    pmin (float)            : Percentile minimum for colourscale, see aplpy
    pmax (float)            : Percentile maxiimum for colourscale, see aplpy
    stretch (str)           : Stretch scale, see aplpy
    cmap (str)              : Colour map, greyscale if None 
    vmid (float)            : Baseline value used for the log and arcsinh stretches
    exponent (float)        : If stretch is set to ‘power’, this is the exponent to use
    invert (bool)           : Invert scaling, only valid if cmap == None
    colourbar (bool)         : Add colour bar if True 
    colourbar_label (str)    : Label for colourbar
    colourbar_width (str)    : Colourbar width in inches
    colourbar_pad (str)      : Colourbar padding in inches
    frametickcolour (str)   : Colour of frames and ticks
    frameticklinewidth (float)  : Linewidth of frames and ticks
    frameticklength (float) : Length of ticks
    showscale (None, astropy quantity): Show a physical and angular scale bar if either an astropy distance (e.g. 1000.*u.pc) or an angle (e.g. 1*u.arcmin) is given
    showscalelab (bool)     : Show a label for the scale bar if True 
    distance = 10,          : Distance of object to calculate scalebar length
    scaleunits = 'arcmin'   : Units used for the angle in the scale bar label
    scale_borderpad (float) : Distance of scalebar from the border
    scale_sep (float)       : Distance of scalebar from scale bar labelling
    scale_ffc (same as colour)   : scalebar frame face colour (None: do not plot)
    scale_ffa = (float)          : scalebar frame face alpha (1: visible, 0: invisible)
    scale_fec = (same as colour) : scalebar frame edge colour (None: do not plot)
    scale_fea = (float)          : scalebar frame edge alpha (1: visible, 0: invisible)
    annotationfontsize (str): General font size
    suppress_xlab (bool)    : Suppress labeling of x-axis (axis- and tick labels)
    suppress_ylab = False   : Suppress labeling of y-axis (axis- and tick labels)
    aspect ('str')          : See matplotlib ax.set_aspect
    """
    
    
#    F.set_system_latex(True)
    if suppress_xlab:
        F.axis_labels.hide_x()
        F.tick_labels.hide_x()
    if suppress_ylab:
        F.axis_labels.hide_y()
        F.tick_labels.hide_y()

    F.axis_labels.set_font(size=annotationfontsize, weight='medium', \
                      stretch='normal', family='serif', \
                           style='normal', variant='normal', color = 'black')
    F.tick_labels.set_font(size=annotationfontsize, weight='medium', \
                         stretch='normal', family='serif', \
                           style='normal', variant='normal', color = 'black')

    F.ticks.set_length(frameticklength)  # points
    F.ticks.set_color(frametickcolour)  # points
    F.ticks.set_linewidth(frameticklinewidth)  # points
    F.ticks._ax.coords[F.ticks.x].ticks.set_tick_out(False)

    kwargs = {'pmin': pmin, 'pmax': pmax, 'stretch': stretch, 'vmid' : vmid, 'exponent': exponent, 'aspect': aspect}

#    print(vmin)
    if vmin != None:
        kwargs['vmin'] = vmin
    if vmax != None:
        kwargs['vmax'] = vmax

    if invert == True:
        kwargs['invert'] = True

    if cmap == None:
        F.show_grayscale(**kwargs)
    else:
        kwargs.pop('invert')
        kwargs['cmap'] = cmap
        F.show_colorscale(**kwargs)

    if colourbar:
        F.add_colorbar()
        #F.set_system_latex(False)
        F.colorbar.set_width(colourbar_width)
        F.colorbar.set_pad(colourbar_pad)
        F.colorbar.set_frame_linewidth(frameticklinewidth)
        F.colorbar.set_frame_color(frametickcolour)
        F.colorbar.set_tick_length(frameticklength)
        F.colorbar.set_font(size=annotationfontsize, weight='medium', \
                            stretch='normal', family='serif', \
                            style='normal', variant='normal')
        if colourbar_label != None:
            F.colorbar.set_axis_label_font(size=annotationfontsize, weight='medium', \
                            stretch='normal', family='serif', \
                            style='normal', variant='normal')
            F.colorbar.set_axis_label_text(colourbar_label)
            #print('box', F.colorbar._colorbar_axes)
            #for i in dir(F.colorbar._colorbar.outline):
            #    if '' in i:
            #        print(i)
    F.ticks._ax.coords[F.ticks.x].ticks.set_tick_out(False)
    F.ticks._ax.coords[F.ticks.y].ticks.set_tick_out(False)
    

    ####
    #### old
    ####
    if False:
#    if showscale != None:
        arcsectokpc = distance.to(u.pc).value*np.pi/(180.*3600.*1000.)
        if scaleunits == 'arcmin':
            scalehere = 1./60.
            symbolhere = '^{{\prime}}'
        else:
            scalehere = 1.
            symbolhere = '^{{\prime\prime}}'
            
        try:
            physicalkpc = showscale.to(u.pc).value/1000.
            angular = u.arcsec*physicalkpc/arcsectokpc
            if scaleunits == 'arcmin':
                text = r'holo$\bf{{{0:.1f}\, \mathrm{{\bf kpc}} \,\bf{{\widehat{{\bf{{=}}}}}} \,{1:.1f}{2:s}}}$'.format(physicalkpc, angular.value*scalehere, symbolhere)
            else:
                text = r'$\bf{{{0:.1f}\, \mathrm{{\bf kpc}} \,\bf{{\widehat{{\bf{{=}}}}}} \,{1:.f}{2:s}}}$'.format(physicalkpc, angular.value*scalehere, symbolhere)
        except:
            angular = showscale.to(u.arcsec)
            physical = (angular.value*arcsectokpc)
            if scaleunits == 'arcmin':
                text = r'$\bf{{{1:.1f}{2:s}\,\widehat{{=}}\, {0:.1f}\, \mathrm{{kpc}}}}$'.format(physical, angular.value*scalehere, symbolhere)
            else:
                text = r'$\bf{{{1:.1f}{2:s}\,\widehat{{=}}\, {0:.f}\, \mathrm{{kpc}}}}$'.format(physical, angular.value*scalehere, symbolhere)
        F.add_scalebar(angular)
        F.scalebar.show(angular, borderpad = scale_borderpad, sep = scale_sep)
        F.scalebar.set_color(scalebarcolor)
        if True:
            F.scalebar.set_frame(True)
            F.scalebar.set_alpha(0.7)
        else:
            F.scalebar.set_frame(False)
            
        F.scalebar.set_corner('bottom right')
        F.scalebar.set_linewidth(scalebarlinewidth)
        if scalebarfontsize == None:
            scalebarfontsize = annotationfontsize
        F.scalebar.set_font(size = scalebarfontsize, weight = 'bold', stretch = 'normal', family = 'serif', style = 'normal', variant = 'normal')
        if showscalelab:
            F.scalebar.set_label(text)
    

    ####
    ####
    ####

#    if False:
    if showscale != None:
        arcsectokpc = distance.to(u.pc).value*np.pi/(180.*3600.*1000.)
            
        if scaleunits == 'arcmin':
            scalehere = 1./60.
            symbolhere = '^{{\prime}}'
        else:
            scalehere = 1.
            symbolhere = '^{{\prime\prime}}'
            
        try:
            physicalkpc = showscale.to(u.pc).value/1000.
            angular = u.arcsec*physicalkpc/arcsectokpc

            if scaleunits == 'arcmin':
                text = r'$\bf{{{0:.1f}\, \mathrm{{\bf kpc}} \,\bf{{\widehat{{\bf{{=}}}}}} \,{1:.1f}{2:s}}}$'.format(physicalkpc, angular.value*scalehere, symbolhere)
            else:
                text = r'$\bf{{{0:.1f}\, \mathrm{{\bf kpc}} \,\bf{{\widehat{{\bf{{=}}}}}} \,{1:.0f}{2:s}}}$'.format(physicalkpc, angular.value*scalehere, symbolhere)
        except:
            angular = showscale.to(u.arcsec)
            physical = (angular.value*arcsectokpc)
            if scaleunits == 'arcmin':
                text = r'$\bf{{{1:.1f}{2:s}\,\widehat{{=}}\, {0:.1f}\, \mathrm{{kpc}}}}$'.format(physical, angular.value*scalehere, symbolhere)
            else:
                text = r'$\bf{{{1:.0f}{2:s}\,\widehat{{=}}\, {0:.1f}\, \mathrm{{kpc}}}}$'.format(physical, angular.value*scalehere, symbolhere)

        if offsetunits != 'arcsec':
            angular = angular/60.
        F.add_scalebar(angular)
        F.scalebar.show(angular, borderpad = scale_borderpad, sep = scale_sep)
        F.scalebar.set_color(scalebarcolor)
        F.scalebar.set_corner('bottom right')
        F.scalebar.set_linewidth(scalebarlinewidth)
        
        if showscalelab:
            if scale_ffc != None or scale_fec != None:
                F.scalebar.set_frame(True)
                if scale_ffc == None:
                    scale_ffc = (1.,1.,1.,0.)
                    scale_ffa = 0.
                if scale_fec == None:
                    scale_fec = (1.,1.,1.,0.)
                    scale_fea = 0.
                grofo = mpl.colors.to_rgba(scale_ffc)
                scale_ffc = (grofo[0], grofo[1], grofo[2], scale_ffa)
                grofo = mpl.colors.to_rgba(scale_fec)
                scale_fec = (grofo[0], grofo[1], grofo[2], scale_fea)
                
#                sigh = (1.,1.,1.,0.7)
                F.scalebar._scalebar.patch.set_facecolor(scale_ffc)
#                sigh = (0,0,0,0)
                F.scalebar._scalebar.patch.set_edgecolor(scale_fec)
#                for hg in dir(F.scalebar._scalebar.patch):
#                    if '' in hg:
#                        print('found it', hg)
#                a = F.scalebar._scalebar.get_children()
#                print('ohn',a[0])
                #F.scalebar.set_alpha(0.7)
            else:
                F.scalebar.set_frame(False)
                
            if scalebarfontsize == None:
                scalebarfontsize = annotationfontsize
            F.scalebar.set_font(size = scalebarfontsize, weight = 'bold', stretch = 'normal', family = 'serif', style = 'normal', variant = 'normal')
            F.scalebar.set_label(text)
    return

def plotmaps_prep(figure = None, figsize = None, subplot=[0.0,0.0,1.,1.], basemap = None, vmin = None, vmax= None, pmin = 0, pmax = 1, stretch = 'linear', vmid = None, exponent = None, cmap = None, invert = False, colourbar = False, colourbar_label = None, colourbar_width = 0.2, colourbar_pad = 0.05, annotationfontsize = 'x-large', suppress_xlab = False, suppress_ylab = False, contoursets = [], contourcols = None,  contourstyles = None, contouralphas = None, contourlevs = [], contourlinewidths = None, frametickcolour = 'black', frameticklinewidth = None, frameticklength = None, showbeam = False, beamfc = '0.7', beamfa = 0.5, beamec = 'black', beamea = 1., showscale = None, showscalelab = True, distance = 10, scaleunits = 'arcsecs', scale_borderpad = 1, scale_sep = 10, scale_ffc = None, scale_ffa = 1.0, scale_fec = None, scale_fea = 1.0, plane = 0, plotvelo = None, velofontsize = None, veloposx = 0.25, veloposy = 0.9, velbbfc = None, velbbec = None, velbbfa = 1, velbbea = 1, scalebarlinewidth = 1, scalebarfontsize = None, scalebarcolor = 'black', aspect = 'equal', x_major_formatter = None, x_ticksn = None, hide_ticl = False, hide_ticl_left = True, hide_ticl_min = 0., hide_ticl_max = 0., hide_ticl_top = 0., hide_ticl_bottom = 0.):
    """
    Plot one map

    Input:
    figure (mpl.figure)                          : Matplotlib figure to add plot onto. Required if used in sub-plot mode, None if sub-plot mode is not wanted
    figsize (float tuple)                        : Figure size in inches, only if figure == None
    subplot (list of floats)                     : Position of sub-plot on figure, first two numbers lower left corner, last two numbers sub-plot width and height, in units of figsize (owned either by figure or specified in this function)
    basemap (str)                                : Name of the file containing the background image
    vmin (float)                                 : Minimum for colourscale, see aplpy
    vmax (float)                                 : Maximum for colourscale, see aplpy
    pmin (float)                                 : Percentile minimum for colourscale, see aplpy
    pmax (float)                                 : Percentile maxiimum for colourscale, see aplpy
    stretch (str)                                : Stretch scale, see aplpy
    cmap (str)                                   : Colour map, greyscale if None 
    invert (bool)                                : Invert scaling, only valid if cmap == None
    colourbar (bool)                              : Add colour bar if True 
    colourbar_label (str)                         : Label for colourbar
    colourbar_width (str)                         : Colourbar width in inches
    colourbar_pad (str)                           : Colourbar padding in inches
    annotationfontsize (str or number)           : General type size
    suppress_xlab (bool)                         : Suppress labeling of x-axis (axis- and tick labels)
    suppress_ylab = False                        : Suppress labeling of y-axis (axis- and tick labels)
    x_major_formatter (None or str)               : (dummy, parameter is used by plotmaps) Formatter of major xticklabels, e.g. 'hh:mm'
    x_ticksn (None or int)                          : (dummy, parameter is used by plotmaps) Numer of xticklabels
    hide_ticl (bool or list of bool)              : (dummy, parameter is used by plotmaps) Create a white patch to hide a portion of the x-axis ticklabels, either one bool for all sub-plots or a list of nx bools for nx individual plots
    hide_ticl_left (bool)                         : (dummy, parameter is used by plotmaps) Indication if the patch is on the left (True) or on the right (False)
    hide_ticl_min (float)                         : (dummy, parameter is used by plotmaps) Horizontal start position of the patch for each sub-plot in units of width of the subplot
    hide_ticl_max (float)                         : (dummy, parameter is used by plotmaps) Horizontal end position of the patch for each sub-plot in units of width of the subplot
    hide_ticl_top (float)                         : (dummy, parameter is used by plotmaps) Vertical start position of the patch for each sub-plot in units of height of the subplot
    hide_ticl_bottom (float)                      : (dummy, parameter is used by plotmaps) Vertical end position of the patch for each sub-plot in units of height of the plot
    contoursets (list of str)                    : List of data sets to generate contours from, which will be overlaid on background image
    contourlevs (list of float/lists)            : List of contours valid for all data sets in contoursets or list of lists of contours valid for each member of contoursets. If fewer controur level lists are given than contoursets, the last given contours are repeated for consecutive contoursets.
    contourcols (None, list, or list of lists)   : Colours corresponding to contourlevs
    contourstyles (None, lists, or list of lists): Contour styles corresponding to contourlevs (see aplpy description, e.g. 'solid' or 'dashed')
    contouralphas (None, lists, or list of lists): Contour alphas corresponding to contourlevs
    contourlinewidths (None, lists, or list of lists): Contour linewidths corresponding to contourlevs (see aplpy description, e.g. 'solid' or 'dashed')
    frametickcolour (str)                        : Colour of frames and ticks
    frameticklinewidth (float)                   : Linewidth of frames and ticks
    frameticklength (float)                      : Length of ticks
    showbeam (bool)                              : Show beam if True
    beamfc (same as colour)                      : beam face colour
    beamfa = (float)                             : beam face alpha (1: visible, 0: invisible)
    beamec = (same as colour)                    : beam edge colour
    beamea = (float)                             : beam edge alpha (1: visible, 0: invisible)
    showscale (bool)                             : Show a physical and angular scale bar
    showscalelab (bool)                          : Show a label for the scale bar if True 
    distance = 10,                               : Distance of object to calculate scalebar length
    scaleunits (str)                             : Units used for the angle in the scale bar label, arcsec or arcmin
    scale_ffc (same as colour)                   : scalebar frame face colour (None: do not plot)
    scale_ffa = (float)                          : scalebar frame face alpha (1: visible, 0: invisible)
    scale_fec = (same as colour)                 : scalebar frame edge colour (None: do not plot)
    scale_fea = (float)                          : scalebar frame edge alpha (1: visible, 0: invisible)
    plane                                        : Number of plane to be plotted/overlaid if input is a cube, starting at 0. Defaults to 0 if plane does not exist or file does not have a third dimension 
    plotvelo (None or str)                       : If a file is give, use the file to plot the velocity information corresponding to plane on the panel
    velofontsize (float)                         : Font size for velocity information, defaulting to annotationfontsize
    veloposx (float)                             : Relative position of velocity label on the panel, x
    veloposy (float)                             : Relative position of velocity label on the panel, y
    velbbfc (None or str)                        : Face colour of bounding box around velocity information, None means transparent
    velbbec (None or str)                        : Line colour of bounding box around velocity information, None means no edge
    velbbfa (float)                              : Alpha of face of bounding box arund velocity information 1: opaque, 0: transparent
    velbbea (float)                              : Alpha of edge of bounding box arund velocity information 1: opaque, 0: transparent
    
    All parameters are self-explaining. If contourlevs is a list of
    floats, these are the contour levels valid for all contoursets. If
    it is a list of lists of floats, they are the contour levels for
    the single contour sets (if fewer elements than contoursets are
    present, then the last element in the list will be used). If
    contourcols is a list of str, these are the contour colors for all
    contoursets. If it is a list of lists of str, they are the contour
    colors for the single contour sets (if fewer elements than
    contoursets are present, then the last element in the list will be
    used). Mixed is possible, like [['red','red','green'], 'red',
    ['green','blue']]. For a list of colours with fewer elements than
    corresponding contours, the contour color is repeated from the
    beginning of the list (matplotlib convention). If contourstyles is
    a list of str, these are the contour styles for all
    contoursets. If it is a list of lists of str, they are the contour
    styles for the single contour sets (if fewer elements than
    contoursets are present, then the last element in the list will be
    used). Mixed is possible, like [['dashed','solid','solid'],
    'dashed', ['solid','solid']]. For a list of styles with fewer
    elements than corresponding contours, the contour style is
    repeated from the beginning of the list (matplotlib convention).

    """

    if basemap == None:
        return
    
    fitsfile_mom = fits.open(basemap)

    if plotvelo != None:
        velotesthead = fits.open(plotvelo)[0].header
        try:
            refpixvelo = float(velotesthead['CRPIX3'])
            refvalvelo = float(velotesthead['CRVAL3'])
            cdeltvelo = float(velotesthead['CDELT3'])

            # If cdeltvelo is larger than some large value it's in units of m/s
            if math.fabs(cdeltvelo) > 1000.:
                cdeltvelo = cdeltvelo/1000.
                refvalvelo = refvalvelo/1000.
        except:
            plotvelo = None
            
    hdu_mom = read2dim(fitsfile_mom[0], plane)

    # Start mom0 plot
    #F = aplpy.FITSFigure(hdu_mom, figure = figure, subplot=subplot)

    #F.set_system_latex(True)
    #F.close()
    if figure == None:
        # This should actually never happen any more
        F = aplpy.FITSFigure(hdu_mom, figsize = figsize)
    else:
        F = aplpy.FITSFigure(hdu_mom, figure=figure, subplot=subplot)

    # Subregion, this does not work because of aplpy!!!
    # Originally I had a a keyword of plotmaps and this function 
    #    subregion (None or list of 4 floats)          : All plot areas will cut to the subregion with the coordinates subimage = [xmin, ymin, xmax, ymax], where xmin, ymin is the lower left corner and  xmax, ymax the upper right one
    #if subregion == None:
    #    print('got haha here')
    #    pass
    #else:
    #    print('got here')
    #    srcentx = (subregion[0]+subregion[2])/2.
    #    srcenty = (subregion[1]+subregion[3])/2. 
    #    srcentw = float(np.fabs(subregion[0]-subregion[2]))/2
    #    srcenth = float(np.fabs(subregion[1]-subregion[3]))/2
    #    F.recenter(srcentx, srcenty, width=srcentw, height=srcenth)


    kwargs = {'vmin' : vmin, 'vmax': vmax, 'pmin' : pmin, 'pmax' : pmax, 'stretch' : stretch, 'vmid' : vmid, 'exponent': exponent, 'cmap' : cmap, 'invert' : invert, 'colourbar' : colourbar, 'colourbar_label': colourbar_label, 'colourbar_width':  colourbar_width, 'colourbar_pad': colourbar_pad, 'frametickcolour' : frametickcolour, 'frameticklinewidth' : frameticklinewidth, 'frameticklength': frameticklength, 'annotationfontsize' : annotationfontsize, 'showscale' : showscale, 'distance' : distance, 'showscalelab' : showscalelab, 'scaleunits' : scaleunits, 'scale_borderpad': scale_borderpad, 'scale_sep': scale_sep, 'scale_ffc': scale_ffc, 'scale_ffa': scale_ffa , 'scale_fec': scale_fec, 'scale_fea': scale_fea, 'suppress_xlab' : suppress_xlab, 'suppress_ylab' : suppress_ylab, 'scalebarlinewidth' : scalebarlinewidth, 'scalebarfontsize': scalebarfontsize, 'scalebarcolor': scalebarcolor, 'aspect': aspect}

    
    plotprep_general(F, **kwargs)

    if showbeam:
#        F.add_beam(major=2, minor = 0.008)
        F.add_beam()
        F.beam.show()
        F.beam.set_corner('bottom left')
        F.beam.set_frame(False)
        F.beam.set_facecolor(beamfc)
        F.beam.set_edgecolor(beamec)
        # If alpha should be set individually to facecolor or edgecolor, then this should be None
        F.beam.set_alpha(None)
        # This is not required
        #        F.beam._beam.ellipse.set_fill(False)
        # This is how to add alpha to a colour, it is the fourth argument
        yoman = F.beam._beam.ellipse.get_facecolor()
        sigh = (yoman[0],yoman[1],yoman[2],beamfa)
        F.beam.set_facecolor(sigh)
        yoman = F.beam._beam.ellipse.get_edgecolor()
        sigh = (yoman[0],yoman[1],yoman[2],beamea)
        F.beam.set_edgecolor(sigh)
#        F.beam._beam.ellipse.set_facecolor(sigh)
#        F.beam.set(facecolor = sigh)
#        for hg in dir(F.beam._beam.ellipse):
#            if '' in hg:
#                print(hg)

    if plotvelo != None:
        if velofontsize == None:
            velofontsize = annotationfontsize
        velo = r'${{\bf {0:.1f}}}\,\mathrm{{\bf{{km}}}}\,\mathrm{{\bf{{s}}}}^{{\bf{{-1}}}}$'.format(refvalvelo+cdeltvelo*(plane-refpixvelo+1))
        F.add_label(veloposx, veloposy, velo, relative = True, size=velofontsize, weight = 'bold', stretch = 'normal', family = 'serif', style = 'normal', variant = 'normal', horizontalalignment='center', verticalalignment='center', color = 'black')
#        velbblw = None
#        velbbfill = True
#        if velbbec == None:
#            velbblw = 0
#        if velbbfc == None:
#            velbbfill = False
        if velbbfc == None:
            velbbfc = (0,0,0,0)
            velbbfa = 0
        if velbbec == None:
            velbbec = (0,0,0,0)
            velbbea = 0
            
        grofo = mpl.colors.to_rgba(velbbfc)
        velbbfc = (grofo[0], grofo[1], grofo[2], velbbfa)
        grofo = mpl.colors.to_rgba(velbbec)
        velbbec = (grofo[0], grofo[1], grofo[2], velbbea)
        F._layers['label_'+str(F._label_counter)].set_bbox(dict(facecolor=velbbfc, edgecolor = velbbec))
#        F._layers['label_'+str(F._label_counter)].set_bbox(dict(facecolor=velbbfc, edgecolor = velbbec, alpha=velbbfa, linewidth = velbblw, fill = velbbfill))
            # Don't delete
            #print(dir(F._layers['label_'+str(F._label_counter)]._bbox_patch))
            #for hy in dir(F._layers['label_'+str(F._label_counter)]._bbox_patch._alpha):
            #    if 'alpha' in hy:
            #        print(hy)

            #F._layers['label_'+str(F._label_counter)]._bbox_patch.set_facecolor('red')
            #F._layers['label_'+str(F._label_counter)]._bbox_patch.set_edgecolor('red')
    for i in range(len(contoursets)):
        fitsfile_here = fits.open(contoursets[i])
        hdu_here = read2dim(fitsfile_here[0], plane)

        # This is obsolete: OK, since this has taken me a day, I should describe it a bit. So in the method show_contour in aplpy there is a call to the method ax.get_transform(frame), where frame is the wcs frame of the data set. If this is a moment-0 map generated by Miriad, this has 3 dimensions which gets get_transform to fail. The local hack is to force usepixel to only use the pixel map, which means that the maps should have the same dimensions. Intrinsically one uses ax.contour(transform=ax.get_transform('pixel')) So aplpy has been modified for this...
#        F.show_contour(hdu_here, levels = contourlevs, colors = contourcols[i], usepixel = True)
        putcontours(F, hdu_here, i, contourlevs = contourlevs, contourcols = contourcols, contourstyles = contourstyles, contouralphas = contouralphas, contourlinewidths = contourlinewidths, usepixel = False)
        #F.show_contour(hdu_here, levels = levelshere, colors = colorshere, linestyles = styleshere)
        
    F.ticks._ax.coords[F.ticks.x].ticks.set_tick_out(False)
    F.ticks._ax.coords[F.ticks.y].ticks.set_tick_out(False)
    F.frame.set_color(frametickcolour)
    if frameticklinewidth != None:
        F.frame.set_linewidth(frameticklinewidth)
        F.ticks.set_linewidth(frameticklinewidth)
    if frameticklength != None:
        F.ticks.set_length(frameticklength)
    return F

def putolays(F, olhere):
    """
    Put overlays onto aplpy figure F

    Input:
    F (aplpy FITSfigure): Input figure to put overlays on
    olhere (list)       : List of overlays

    Uses aplpy to overlay figures on aplpy fits figure F. olhere is a
    list of pairs [figurename (str), kwargs (dict)], where kwargs is a
    dictionary of argument-value pairs. Some arguments in kwargs are required,
    depending on the figure. Others are the kwargs passed to the
    corresponding aplpy overlay function. putoverlays goes through
    all members of the input list, adding the corresponding overlays on
    F. The following table shows compulsory kwargs, and the aplpy
    method that will be called. The exception is label with four additional arguments:
    label_bbfc: label bounding box face coulour
    label_bbfa: label bounding box face alpha
    label_bbec: label bounding box edge colour
    label_bbea: label bounding box edge alpha

    Example:
    olhere = [['circle', {'pv_x_world' : 120, 'pv_y_world' : -30, 'radius': 0.001}]]
    
    For the corresponding description of
    the parameters, see aplpy description
    https://aplpy.readthedocs.io/en/stable/api/aplpy.FITSFigure.html

    figurename   compulsory kwargs  corresponding aplpy method

    'label'      pv_x_world (float or list) F.add_label(pv_x_world, pv_y_world, pv_label, **kwargs)
                 pv_y_world (float or list)
                 pv_label (str)

    'markers'    pv_x_world (float or list) F.show_markers(pv_x_world, pv_y_world, **kwargs)
                 pv_y_world (float or list)

    'circles'    pv_x_world (float or list) F.show_circles(pv_x_world, pv_y_world, pv_radius, **kwargs)
                 pv_y_world (float or list)
                 pv_radius (float or list)

    'ellipses'   pv_x_world (float or list) F.show_circles(pv_x_world, pv_y_world, pv_width, pv_height, **kwargs)
                 pv_y_world (float or list)
                 pv_width (float or list)
                 pv_height (float or list)
    
    'rectangles' pv_x_world (float or list) F.show_rectangles(pv_x_world, pv_y_world, pv_width, pv_height, **kwargs)
                 pv_y_world (float or list)
                 pv_width (float or list)
                 pv_height (float or list)

    'arrows'     pv_x_world (float or list) F.show_arrows(pv_x_world, pv_y_world, pv_dx, pv_dy, **kwargs)
                 pv_y_world (float or list)
                 pv_dx (float or list)
                 pv_dy (float or list)

    """

    # Lazy variant to get all the parameters for whichever overlay.
    for overlay in olhere:
        try:
            label = overlay[1]['pv_label']
            overlay[1].pop('pv_label')
        except:
            pass
        try:
            x_world = overlay[1]['pv_x_world']
            overlay[1].pop('pv_x_world')
        except:
            pass
        try:
            y_world = overlay[1]['pv_y_world']
            overlay[1].pop('pv_y_world')
        except:
            pass
        try:
            width   = overlay[1]['pv_width']
            overlay[1].pop('pv_width')
        except:
            pass
        try:
            height  = overlay[1]['pv_height']
            overlay[1].pop('pv_height')
        except:
            pass
        try:
            dx      = overlay[1]['pv_dx']
            overlay[1].pop('pv_dx')
        except:
            pass
        try:
            dy      = overlay[1]['pv_dy']
            overlay[1].pop('pv_dy')
        except:
            pass

        # Now make the right overlay
        if overlay[0] == 'label':
            try:
                label_bbfa = overlay[1]['label_bbfa']
                retlabel_bbfa = overlay[1]['label_bbfa']
                overlay[1].pop('label_bbfa')
            except:
                label_bbfa = 1.
                #retlabel_bbfa = None
            try:
                label_bbfc = overlay[1]['label_bbfc']
                retlabel_bbfc = overlay[1]['label_bbfc']
                overlay[1].pop('label_bbfc')
            except:
                label_bbfc = 'black'
                label_bbfa = 0.
                #retlabel_bbfc = None
            try:
                label_bbea = overlay[1]['label_bbea']
                retlabel_bbea = overlay[1]['label_bbea']
                overlay[1].pop('label_bbea')
            except:
                label_bbea = 1.
                #retlabel_bbea = None
            try:
                label_bbec = overlay[1]['label_bbec']
                retlabel_bbec = overlay[1]['label_bbec']
                overlay[1].pop('label_bbec')
            except:
                label_bbec = 'black'
                label_bbea = 0.
                #retlabel_bbec = None
            grofo = mpl.colors.to_rgba(label_bbfc)
            label_bbfc = (grofo[0], grofo[1], grofo[2], label_bbfa)
            grofo = mpl.colors.to_rgba(label_bbec)
            label_bbec = (grofo[0], grofo[1], grofo[2], label_bbea)
            F.add_label(x_world, y_world, label, **overlay[1])
            F._layers['label_'+str(F._label_counter)].set_bbox(dict(facecolor=label_bbfc, edgecolor = label_bbec))
            try:        
                overlay[1]['label_bbfc']   = retlabel_bbfc   
            except:     
                pass    
            try:        
                overlay[1]['label_bbfa']   = retlabel_bbfa   
            except:     
                pass    
            try:        
                overlay[1]['label_bbec']   = retlabel_bbec   
            except:     
                pass    
            try:        
                overlay[1]['label_bbea']   = retlabel_bbea   
            except:     
                pass    

        if overlay[0] == 'markers':
            F.show_markers(x_world, y_world, **overlay[1])

        if overlay[0] == 'circles':
            F.show_circles(x_world, y_world, radius, **overlay[1])

        if overlay[0] == 'ellipses':
            F.show_ellipses(x_world, y_world, width, height, **overlay[1])
    
        if overlay[0] == 'rectangles':
            F.show_rectangles(x_world, y_world, width, height, **overlay[1])

        if overlay[0] == 'arrows':
            F.show_arrows(x_world, y_world, dx, dy, **overlay[1])

        # Then put them back
        try:
            overlay[1]['pv_label'] = label  
        except:     
            pass    
        try:
            overlay[1]['pv_x_world'] = x_world  
        except:     
            pass    
        try:        
            overlay[1]['pv_y_world'] = y_world 
        except:     
            pass    
        try:        
            overlay[1]['pv_width']   = width   
        except:     
            pass    
        try:        
            overlay[1]['pv_height']  = height  
        except:     
            pass    
        try:        
            overlay[1]['pv_dx']      = dx      
        except:     
            pass    
        try:        
            overlay[1]['pv_dy']      = dy      
        except:
            pass
    return

def plotmaps(width = A4widht_in_inch, plotmargin = 0.5*cm_in_inch, labelmargin_left = 3.5*cm_in_inch, labelmargin_bottom = 1.75*cm_in_inch, labelmargin_right = 0., labelmargin_top = 0., vmin = None, vmax= None, pmin = None, pmax = None, stretch = None, cmap = None, colourbar = False, colourbar_width = 0.2, colourbar_pad = 0.05, invert = False, chans = None, nx = 0, ny =0, showscale = None, showbeam = True, showscalelab = True, reducescalelab = True, plotname = None, individual_maps = False, contoursets = None, contourlevs = None, contourcols = None, contourstyles = None, contouralphas = None, individual_olays = False, olays = None, x_major_formatter = None, x_ticksn = None, hide_ticl = False, hide_ticl_left = True, hide_ticl_min = 0., hide_ticl_max = 0., hide_ticl_top = 0., hide_ticl_bottom = 0., **kwargs):
    """

    Plot a grid of maps
    width (float)             : Width in inches, 8.27 is A4 width
    plotmargin (float)        : A margin around the whole plot, 0.1968504 is 0.5 cm
    labelmargin_right (float)  : An additional margin to the right, needed for colour bars
    labelmargin_left (float)  : An additional margin to the left, 1.377953 is 3.5 cm
    labelmargin_bottom (float): An additional margin to the bottom for labels, 0.6889764 is 1.75
    chans (None or list)      : List of channels to plot.
    nx (int)                  : Number of vertical panels
    ny (int)                  : Number of horizontal panels
    showscalelab (bool)       : Show label to the scalebar if True
    reducescalelab (bool)     : Show only label to the scalebar for left panels if True

    basemap (str or list of str)                 : (List of) name(s) of the file(s) containing the background image
    vmin (float or list of float)                : (List of) minimum(s) for colourscale(s), see aplpy
    vmax (float or list of float)                : (List of) Maximum for colourscale, see aplpy
    pmin (float or list of float)                : (List of) Percentile minimum for colourscale, see aplpy
    pmax (float or list of float)                : (List of) Percentile maxiimum for colourscale, see aplpy
    stretch (str or list of str)                 : (List of) Stretch scale, see aplpy
    cmap (str or list of str)                    : (List of) Colour map, greyscale if None 
    invert (bool or list of bool)                : (List of) Invert scaling, only valid if cmap == None
    colourbar (bool)                              : Add colour bar to the right side of the plot raser if True. This is never a list.
    colourbar_label (str)                         : Label for colourbar
    colourbar_width (str)                         : Colourbar width in inches
    colourbar_pad (str)                           : Colourbar padding in inches
    annotationfontsize (str or number)            : General type size
    suppress_xlab (bool)                          : Generally suppress labeling of x-axis (axis- and tick labels)
    suppress_ylab (bool)                          : Generally suppress labeling of y-axis (axis- and tick labels)
    x_major_formatter (None or str)               : Formatter of major xticklabels, e.g. 'hh:mm'
    x_ticksn (None or int)                          : Numer of xticklabels
    hide_ticl (bool or list of bool)              : Create a white patch to hide a portion of the x-axis ticklabels, either one bool for all sub-plots or a list of nx bools for nx individual plots
    hide_ticl_left (bool)                         : Indication if the patch is on the left (True) or on the right (False)
    hide_ticl_min (float)                         : Horizontal start position of the patch for each sub-plot in units of width of the subplot
    hide_ticl_max (float)                         : Horizontal end position of the patch for each sub-plot in units of width of the subplot
    hide_ticl_top (float)                         : Vertical start position of the patch for each sub-plot in units of height of the subplot
    hide_ticl_bottom (float)                      : Vertical end position of the patch for each sub-plot in units of height of the plot

    individual_maps (bool)                   : Instead of applying contoursets, contourlevs, contourcols, contourstyles for all data sets given, the corresponding quantities are given for each input cube in a list
    contoursets (list of str or list of lists)   : List of data sets to generate contours from, which will be overlaid on background images
    contourlevs (list of float/lists)            : List of contours valid for all data sets in contoursets or list of lists of contours valid for each member of contoursets. If fewer controur level lists are given than contoursets, the last given contours are repeated for consecutive contoursets.
    contourcols (None, list, or list of lists)   : Colours corresponding to contourlevs
    contourstyles (None, lists, or list of lists): Contour styles corresponding to contourlevs (see aplpy description, e.g. 'solid' or 'dashed')
    contouralphas (None, lists, or list of lists): Contour alphas corresponding to contourlevs
    frametickcolour (str)                        : Colour of frames and ticks
    frameticklinewidth (float)                   : Linewidth of frames and ticks
    frameticklength (float)                      : Length of ticks
    showbeam (bool)                              : Show beam if True 
    beamfc (same as colour)                      : beam face colour
    beamfa = (float)                             : beam face alpha (1: visible, 0: invisible)
    beamec = (same as colour)                    : beam edge colour
    beamea = (float)                             : beam edge alpha (1: visible, 0: invisible)
    showscale (bool)                             : Show a physical and angular scale bar
    distance (float)                             : Distance of object to calculate scalebar length
    scaleunits (str)                             : Units used for the angle in the scale bar label, arcsec or arcmin
#    plane (None, or list of int)                 : Number of plane to be plotted/overlaid if input is a cube, starting at 0. Defaults to 0 if plane does not exist or file does not have a third dimension 
    plotvelo (None or str)                       : If a file is give, use the file to plot the velocity information corresponding to plane on the panel
    veloposx (float)                             : Relative position of velocity label on the panel, x
    veloposy (float)                             : Relative position of velocity label on the panel, y
    velbbfc (None or str)                        : Face colour of bounding box around velocity information, None means transparent
    velbbec (None or str)                        : Line colour of bounding box around velocity information, None means no edge
    velbbalpha (float)                           : Alpha of face of bounding box arund velocity information 1: opaque, 0: transparent
    indiviual_olays (bool)                       : Overlays to be interpreted as list, for each plane, or general
    olays (None, list, or list of lists)         : Overlays to each plot, None, list, or list of pairs
    plotname (str)                               : Name of output plot

    Will create an nx x ny raster of panels with each panel being a
    map, sharing the axes. The graph will be filled from top left over
    top right and so on to bottom right (like reading). Basemap is the
    background cube, or a list of cubes. If a list of cubes, the first
    panel will be generated from the first background cube in basemap,
    the second from the second background cube, and so on. Chans is a
    list of channels (starting at 0) to be plotted, starting from the
    first panel, extrapolated from the last given number. So, for
    example, if basemap is ['a.fits','b.fits'] and chans = [1,2,3],
    the first panel will show the second channel of the background
    cube a.fits, the second the third channel of the background cube
    b.fits, and the third panel will show the fourth channel of the
    background cube b.fits. The parameters vmin, vmax, pmin, pmax,
    stretch, cmap, invert, colourbar can be given as lists and their
    nth elements will apply to nth panel, extrapolting from the last
    element.  If a plane does not exist, the last plane of the cube is
    shown, if a third axis does not exist, the image will be show. The
    list of data sets to overplot contours contoursets is valid for
    all panels unless the parameter individual_maps is set to True. In that
    case contoursets, contourlevs, contourcols, contourstyles are
    lists of the single quantities (or further lists), for each
    channel or basemap individually, with extrapolation being
    applied. Also here, the single planes are shown as listed in
    chans. The
    list of scaling variables for the background image is valid for
    all panels unless the parameter individual_maps is set to True. In that
    case vmin, vmax, pmin, pmax, stretch, cmap, invert are
    lists of the single quantities (or further lists), for each
    channel or basemap individually, with extrapolation being
    applied. Also here, the single planes are shown as listed in
    chans. The same is true for olays with the switch individual_olays, which is either None, a list of
    overlays or a list of lists of overlays. The latter is expected if
    individual_olays == True. For the description of overlays see function
    putolays.

    """
    print('Plotmaps {}'.format(plotname))
    plt.rc('text', usetex=True)

    # Check if basemap is a list
    maps = kwargs['basemap']
    if type(maps) == type([]):
        kwargs['basemap'] = maps[0]

    # Get intel
    hdu = fits.open(kwargs['basemap'])[0]
    try:
        nchans = hdu.header['NAXIS3']
    except:
        nchans = 0

    sizex = hdu.header['NAXIS1']
    sizey = hdu.header['NAXIS2']
    x_over_y = float(sizex)/(sizey)

    if chans == None:
        chans = [0]
    if type(maps) == type([]):
        while len(chans) < len(maps):
            chans.append(0)
    
    if nx == 0:
        if ny == 0:
            nx = int(np.ceil((np.sqrt(float(len(chans))))))
        else:
            nx = int(np.ceil(float(len(chans))/float(ny)))
    if ny == 0:
        ny = int(np.ceil(float(len(chans))/float(nx)))
        
    npanels = nx*ny
    
    # Width of a single panel, in inches
    # swidth = (1.0-labwidth_over_width)/float(nx)
    # Caution! if a colourbar is expected this needs to be expanded by the width of the colour bar and the padding, which is given in inches
    if colourbar:
        adwidth = colourbar_width+colourbar_pad
    else:
        adwidth = 0.
    swidth = (width-labelmargin_left-2.*plotmargin-labelmargin_right)/float(nx)+adwidth
    shiftwidth = (width-labelmargin_left-2.*plotmargin-labelmargin_right)/float(nx)
    
    # Height of a single panel, in inches
    sheight = (swidth-adwidth)/x_over_y
                 
    # Plot height in inches
    plotheight = labelmargin_bottom+2.*plotmargin+ny*sheight+labelmargin_top

    fig = plt.figure(figsize=(width, plotheight))

    # y-spaces, in units of plot height
    labmarg_over_height = labelmargin_bottom/plotheight
    pnoh = sheight/plotheight
    plotmarg_over_height = plotmargin/plotheight
    btoh = labmarg_over_height+plotmarg_over_height
  
    # x-spaces, in units of plot width
    labmarg_over_width = labelmargin_left/width
    pnow = swidth/width
    pnowshift = shiftwidth/width
    plotmarg_over_width = plotmargin/width
    lfow = labmarg_over_width+plotmarg_over_width

    #subplot=[lfow+posnx*pnow, btoh+posny*pnoh, pnow, pnoh]

    # For purposes of feeding into FITSFigure, the plot heights and widths
    # sheight = (1.0-labwidth_over_width)/float(ny)

    # Put common axis labels
    # To make sure that the background is white, just as the labels
    fig.patch.set_facecolor('white')

    # Do all the labeling using
    contoursets_intro = contoursets
    contourlevs_intro = contourlevs
    contourcols_intro = contourcols
    contourstyles_intro = contourstyles
    contouralphas_intro = contouralphas
    vmin_intro    = vmin   
    vmax_intro    = vmax   
    pmin_intro    = pmin   
    pmax_intro    = pmax   
    stretch_intro = stretch
    invert_intro  = invert 

    a = showscale
    b = showbeam
    c = colourbar
    kwargs['showscale'] = None
    kwargs['showbeam'] = False
    kwargs['colourbar'] = False
    kwargs['vmin'] = -10000002.
    kwargs['vmax'] = -10000001.
    d = kwargs['plotvelo']
    kwargs['plotvelo'] = None
    
    F = plotmaps_prep(figure = fig, subplot=[lfow, btoh, pnow, pnoh*float(ny)], plane = 0, **kwargs)
    kwargs['showscale'] = a
    kwargs['showbeam'] = b
    kwargs['colourbar'] = c
    
    #F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow, btoh, pnow, pnoh*float(ny)])
    F.axis_labels.hide_x()
    F.tick_labels.hide_x()
    F.ticks.hide()
    #F.colourbar.hide()
    F.axis_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'black')
    F.tick_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'white')
    F.frame.set_color('white')
    F.close()

    kwargs['showscale'] = None
    kwargs['showbeam'] = False
    kwargs['colourbar'] = False
    kwargs['vmin'] = -10000002.
    kwargs['vmax'] = -10000001.
    F = plotmaps_prep(figure = fig, subplot=[lfow, btoh, pnow*float(nx), pnoh], plane = 0, **kwargs)
    kwargs['showscale'] = a
    kwargs['showbeam'] = b
    kwargs['colourbar'] = c
    
    #    F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow, btoh, pnow*float(nx), pnoh])
    F.axis_labels.hide_y()
    F.tick_labels.hide_y()
    F.ticks.hide()
    #F.colourbar.hide()
    F.axis_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'black')
    F.tick_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'white')
    F.frame.set_color('white')
    F.close()

    kwargs['colourbar_width'] = colourbar_width
    kwargs['colourbar_pad'] = colourbar_pad

    contoursets_intro = contoursets
    contourlevs_intro = contourlevs
    contourcols_intro = contourcols
    contourstyles_intro = contourstyles
    contouralphas_intro = contouralphas
    vmin_intro    = vmin   
    vmax_intro    = vmax   
    pmin_intro    = pmin   
    pmax_intro    = pmax   
    stretch_intro = stretch
    cmap_intro    = cmap   
    invert_intro  = invert 
    kwargs['plotvelo'] = d

    # Now loop over the panels
    chanpospre = []
    for i in range(len(chans)):
        chanpospre += [[i, i % nx, ny - i // nx - 1]]
    if hide_ticl_left:
        chanpos = []
        for iy in range(ny):
            for ix in range(nx):
                if nx*iy+(nx-ix-1) < len(chans):
                    chanpos += [chanpospre[nx*iy+(nx-ix-1)]]
    else:
        chanpos = chanpospre
                    
        
#        posnx = i % nx
#        posny = ny - i // nx - 1

    for jojo in range(len(chanpos)):
        i, posnx, posny = chanpos[jojo]
        
        if posnx == 0:
            showscalelab_here = showscalelab
        else:
            showscalelab_here = False

        if reducescalelab:
            if i != 0:
                showscalelab_here = False
                
        # Get the graph
        if type(maps) == type([]):
            if i < len(maps):
                kwargs['basemap'] = maps[i]
            else:
                kwargs['basemap'] = maps[-1]
                
        # Same for contours, contour colors, and contour types
        for j in [['contoursets',   contoursets_intro],
                  ['contourlevs',   contourlevs_intro],
                  ['contourcols',   contourcols_intro],
                  ['contourstyles', contourstyles_intro],
                  ['contouralphas', contouralphas_intro],
                  ['vmin',          vmin_intro],
                  ['vmax',          vmax_intro],
                  ['pmin',          pmin_intro],
                  ['pmax',          pmax_intro],
                  ['stretch',       stretch_intro],
                  ['cmap',          cmap_intro],
                  ['invert',        invert_intro]
                  ]:
            if individual_maps and type(j[1]) == type([]):
                if i < len(j[1]):
                    kwargs[j[0]] = j[1][i]
                else:
                    kwargs[j[0]] = j[1][-1]
            else:
                kwargs[j[0]] = j[1]

        F = plotmaps_prep(figure = fig, subplot=[lfow+posnx*pnowshift, btoh+posny*pnoh, pnow, pnoh], plane = chans[i], showscalelab = showscalelab_here, **kwargs)
        # F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow+posnx*pnow, btoh+posny*pnoh, pnow, pnoh])

        # Now care about the olays
        if olays != None:
            if individual_olays:
                if i < len(olays):
                    olhere = olays[i]
                else:
                    olhere = olays[-1]
            else:
                olhere = olays
            putolays(F, olhere)

        F.axis_labels.hide_x()
        F.axis_labels.hide_y()
#        print('yoman: {}'.format(F.tick_labels._ax.coords[0].get_xticklabels()))
        if posnx != 0:
            F.tick_labels.hide_y()
        if colourbar and posnx != nx-1:
            F.colorbar.hide()
        # Leave tick labels if there is no panel below
        if (posny != 0) and nx*(i//nx+1)+posnx < len(chans):
            F.tick_labels.hide_x()
        #    for i in dir(F.ax.coords[F.x]):
#            print('jo {}'.format(i))
#        print('hji |{}|'.format(F.ax.coords[F.ticks.x].ticklabels.get_label()))
#        print('hji |{}|'.format(F.ax.get_xticklabels(which='major')))
#        thelist = dir(F.ax.coords[F.ticks.x].ticklabels)
        #    F.ticks._ax.coords[F.ticks.x].set_ticklabel(exclude_overlapping=True)
        #F.ax.coords[F.x].set_ticklabel(exclude_overlapping=True)
        if x_major_formatter != None:
            F.ax.coords[F.x].set_major_formatter(x_major_formatter)
        if x_ticksn != None:
            F.ax.coords[F.x].set_ticks(number = x_ticksn)
#        thelist = F.ax.get_xticklabels(which='major')
#        for hji in thelist:
#            hji.set_text('fu')
#            print('ji {}'.format(hji))
        F.close()
    
        # [lfow+posnx*pnowshift+pnow, btoh+posny*pnoh-0.05, pnow, 0.05]
        # Yes, these are redundant calculations, but I want to keep this together in one progamming block
        if type(hide_ticl) != type([]):
            hide_ticle = [hide_ticl]
        else:
            hide_ticle = hide_ticl        
        hidelabs = [i for i in hide_ticle]
        while len(hidelabs) < nx:
            hidelabs += [hidelabs[-1]]
        if hidelabs[posnx] == True:
            fhead = read2dim(fits.open(kwargs['basemap'])[0],0)
            fhead.data = 0.*fhead.data
            F = aplpy.FITSFigure(fhead, figure=fig, subplot=[lfow+posnx*pnowshift+hide_ticl_min*pnow, btoh+posny*pnoh+hide_ticl_bottom*pnoh, hide_ticl_max*pnow-(hide_ticl_min*pnow), (hide_ticl_top-hide_ticl_bottom)*pnoh])
            F.axis_labels.hide_x()
            F.tick_labels.hide_x()
            F.axis_labels.hide_y()
            F.tick_labels.hide_y()
            F.ticks.hide()
            F.frame.set_color('white')
            F.close()


    #    fig.canvas.draw()
    if plotname != None:
        fig.savefig(plotname)
    return    

def plotpvdiagrams(bgname_prefix = '', vmin = None, vmax = None, pmin = None, pmax= None, stretch = 'default', cmap=None, invert = True, colourbar = False, annotationfontsize = 'x-large', contoursets = [], postfixes = ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits'], letters = ['A', 'B', 'C', 'D'], lettercol = 'black', contourcols = None, contourstyles = None, contouralphas = None, contourlevs = [], frametickcolour = '#555555', frameticklinewidth = None, frameticklength = None, showbeam = True, beaminfoset = '', vres = 2, showscale = 1*u.arcmin, scaleunits = 'arcsec', distance = 7010000*u.pc, scale_borderpad = 1, scale_sep = 10, scale_ffc = None, scale_ffa = 1.0, scale_fec = None, scale_fea = 1.0, plotprefix = '', plotpostfix = '.pdf', figsize = None, aspect = 'auto', contourlinewidths = 1, scalebarlinewidth = 1, scalebarcolor = 'black', letterfontsize = 32, letterposx = 0.95, letterposy = 0.95, letter_ffc = None, letter_ffa = 0., letter_fec = None, letter_fea = 1., scalebarfontsize = 'medium', beamfc = '0.7', beamfa = 0.5, beamec = 'black', beamea = 1., showscalelab = True, reducescalelab = True, borderleft = 1.25, borderbottom = 1.25, borderright = 0., bordertop = 0.,x_major_formatter = None, x_ticksn = None, hide_ticl = False, hide_ticl_left = True, hide_ticl_min = 0., hide_ticl_max = 0., hide_ticl_top = 0., hide_ticl_bottom = 0., offsetunits = 'arcsec'):
    """
    Plot a number of PV diagrams

    plotmargin (float)                           : A margin around the whole plots, 0.1968504 is 0.5 cm
    bgname_prefix (str)                          : Name prefix of background fits files
    vmin (float)                                 : Minimum for colourscale, see aplpy
    vmax (float)                                 : Maximum for colourscale, see aplpy
    pmin (float)                                 : Percentile minimum for colourscale, see aplpy
    pmax (float)                                 : Percentile maxiimum for colourscale, see aplpy
    stretch (str)                                : Stretch scale, see aplpy
    cmap (str)                                   : Colour map, greyscale if None 
    invert (bool)                                : Invert scaling, only valid if cmap == None
    colourbar (bool)                              : Add colour bar if True 
    annotationfontsize (str or number)           : General type size
    contoursets (list of str)                    : List of data sets to generate contours from, which will be overlaid on background image
    contourlevs (list of float/lists)            : List of contours valid for all data sets in contoursets or list of lists of contours valid for each member of contoursets. If fewer controur level lists are given than contoursets, the last given contours are repeated for consecutive contoursets.
    contourcols (None, list, or list of lists)   : Colours corresponding to contourlevs
    contourstyles (None, lists, or list of lists): Contour styles corresponding to contourlevs (see aplpy description, e.g. 'solid' or 'dashed')
    contouralphas (None, lists, or list of lists): Contour alphas corresponding to contourlevs
    postfixes (list of str)                      : Name postfix of background fits files
    letters (list of str)                        : List of strings to annotate output plots with
    lettercol (str)                              : Colour of strings to annotate output plots with
    frametickcolour (str)                        : Colour of frames and ticks
    frameticklinewidth (str)                     : Linewidth of frames and ticks
    frameticklength (float)                      : Length of ticks
    showbeam (bool)                              : Show beam if True 
    beaminfoset (str)                            : Name of data set to generate beam from 
    vres = 2                                     : Velocity resolutio in pixels
    showscale (bool)                             : Show a physical and angular scale bar
    scaleunits (str)                             : Units used for the angle in the scale bar label, arcsec or arcmin
    distance (float)                             : Distance of object to calculate scalebar length
    plotprefix (str)                             : Prefix to output plots
    plotpostfix (str)                            : Postfix (ending) of output plots, determines output file type (e.g. '.pdf', '.png')
    figsize (float tuple)                        : Size of figure in inches, x and y
x_major_formatter (None or str)               : Formatter of major xticklabels, e.g. 'hh:mm'
    x_ticksn (None or int)                          : Numer of xticklabels
    hide_ticl (bool or list of bool)              : Create a white patch to hide a portion of the x-axis ticklabels, either one bool for all sub-plots or a list of nx bools for nx individual plots
    hide_ticl_left (bool)                         : Indication if the patch is on the left (True) or on the right (False)
    hide_ticl_min (float)                         : Horizontal start position of the patch for each sub-plot in units of width of the subplot
    hide_ticl_max (float)                         : Horizontal end position of the patch for each sub-plot in units of width of the subplot
    hide_ticl_top (float)                         : Vertical start position of the patch for each sub-plot in units of height of the subplot
    hide_ticl_bottom (float)                      : Vertical end position of the patch for each sub-plot in units of height of the plot
    offsetunits (str)                             : Units used for the angle in the x-axis, arcsec or arcmin

    Produces a number of overlays reading in files
    bgname_prefix+postfixes[i] as background greyscale or colour
    images and contoursets+postfixes[i] to generate contours from. For
    each image, a letter from letters is put on the graph as an
    annotation. It is assumed that the input files are PV diagrams,
    with an angular x-axis and a velocity y-axis, as generated by
    pvdiagrams. The beam is then generated with the angular size of
    sqrt(BMAJ*BMIN) and the velocity size of npixels*CDELT3 in km/s,
    where BMAJ, BMIN, and CDELT3 are taken from beaminfoset. After the
    plot is generated, it is saved as
    plotprefix+'_'+postfixes[i]+bgname_prefix+contoursets[0]+contoursets[1]+...+plotpostfix,
    where postfixes are stripped from the file descriptors at the end
    (anything after the last '.'). As an example, bgname_prefix = 'o',
    contoursets = ['o','f'], postfixes = ['maj.fits', 'min.fits'],
    plotprefix = 'phew', plotpostfix = '.png' will generate the output
    files phew_maj_o_o_f.png and phew_min_o_o_f.png.

    """
    
    print('Plotpvdiagrams')
    plt.rc('text', usetex=True)

    levels = contourlevs

    if figsize == None:
        figsize = (A4widht_in_inch, A4widht_in_inch)
    borderleft = (borderleft*cm_in_inch)/figsize[0]
    borderright = (borderright*cm_in_inch)/figsize[0]
    bordertop = (bordertop*cm_in_inch)/figsize[1]
    borderbottom = (borderbottom*cm_in_inch)/figsize[1]

    for posnum in range(len(postfixes)):
        name = postfixes[posnum]
        letter = letters[posnum]
        hdulist = fits.open(bgname_prefix+name)
        header = hdulist[0].header
        image = hdulist[0].data
        hdulist.close()
#        sys.exit()
        
        
        # region of image we want to display
        if offsetunits == 'arcsec':
            header['CDELT1'] = header['CDELT1']*3600.
            header['CRVAL1'] = header['CRVAL1']*3600.
        else:
            header['CDELT1'] = header['CDELT1']*60.
            header['CRVAL1'] = header['CRVAL1']*60.
            
        header['CDELT2'] = header['CDELT2']/1000.
        header['CRVAL2'] = header['CRVAL2']/1000.
        header.set('BMAJ', 2)
        header.set('BMIN', 0.008)
        header.set('BPA', 0.)
        hdu = fits.PrimaryHDU(data=image,header=header)
        #fig = aplpy.FITSFigure(hdu, figsize = figsize)
        figure = plt.figure(figsize=figsize)
        fig = aplpy.FITSFigure(hdu, figure=figure, subplot=[borderleft,borderbottom,1.-borderleft+borderright,1.-bordertop-borderbottom])
        # To make sure that the background is white, just as the labels
        # fig.patch.set_facecolor('white')

#        fig.set_system_latex(True)
        if reducescalelab:
            if posnum > 0:
                scale_ffc = None
                scale_fec = None
                showscalelab = False
        kwargs = {'vmin' : vmin, 'vmax': vmax, 'pmin' : pmin, 'pmax' : pmax, 'stretch' : stretch, 'cmap' : cmap, 'invert' : invert, 'colourbar' : colourbar, 'frametickcolour' : frametickcolour, 'frameticklinewidth' : frameticklinewidth, 'frameticklength' : frameticklength,'showscale' : showscale, 'distance' : distance, 'scaleunits' : scaleunits, 'annotationfontsize' : annotationfontsize, 'aspect': aspect, 'scalebarlinewidth': scalebarlinewidth, 'scalebarfontsize': scalebarfontsize, 'scalebarcolor': scalebarcolor, 'scale_borderpad': scale_borderpad, 'scale_sep': scale_sep, 'scale_ffc': scale_ffc, 'scale_ffa': scale_ffa , 'scale_fec': scale_fec, 'scale_fea': scale_fea, 'showscalelab': showscalelab, 'offsetunits': offsetunits}
        plotprep_general(fig, **kwargs)

        fig.set_xaxis_coord_type('scalar')
        if offsetunits == 'arcsec':
            fig.axis_labels.set_xtext(r'OFFSET (arcsec)')
        else:
            fig.axis_labels.set_xtext(r'OFFSET (arcmin)')
            
        fig.axis_labels.set_ytext(r'VELOCITY (km$\,$s$^{-1}$)')

        if x_major_formatter != None:
            fig.ax.coords[fig.x].set_major_formatter(x_major_formatter)
        if x_ticksn != None:
            fig.ax.coords[fig.x].set_ticks(number = x_ticksn)

        name = plotprefix+'.'.join(postfixes[posnum].split('.')[:-1])+'_'+bgname_prefix
        for contournumber in range(len(contoursets)):
            name = name+'_'+contoursets[contournumber]
            hdulist = fits.open(contoursets[contournumber]+postfixes[posnum])
            header = hdulist[0].header
            image = hdulist[0].data
            hdulist.close()
            
            # region of image we want to display
            if offsetunits == 'arcsec':
                header['CDELT1'] = header['CDELT1']*3600.
                header['CRVAL1'] = header['CRVAL1']*3600.
            else:
                header['CDELT1'] = header['CDELT1']*60.
                header['CRVAL1'] = header['CRVAL1']*60.
                
            header['CDELT2'] = header['CDELT2']/1000.
            header['CRVAL2'] = header['CRVAL2']/1000.

            hdu = fits.PrimaryHDU(data=image,header=header)
            putcontours(fig, hdu, contournumber, contourlevs, contourlinewidths = contourlinewidths, contouralphas = contouralphas, contourcols = contourcols, contourstyles = contourstyles, usepixel = True)
            
            #fig.show_contour(data=hdu, colors=contourcols[contournumber], levels=contourlevs, usepixel = True)
        if letter_ffc == None:
            letter_ffa_ret = letter_ffa
            letter_ffa = 0.
        else:
            letter_ffc_ret = letter_ffc
            letter_ffa_ret = letter_ffa
        if letter_fec == None:
            letter_fea_ret = letter_fea
            letter_fea = 0.
        else:
            letter_fec_ret = letter_fec
            letter_fea_ret = letter_fea
            
        grofo = mpl.colors.to_rgba(letter_ffc)
        letter_ffc = (grofo[0], grofo[1], grofo[2], letter_ffa)
        grofo = mpl.colors.to_rgba(letter_fec)
        letter_fec = (grofo[0], grofo[1], grofo[2], letter_fea)
        fig.add_label(letterposx, letterposy, letters[posnum], relative = True, size=letterfontsize, weight='bold', horizontalalignment='center', verticalalignment='center', color=lettercol)
        fig._layers['label_'+str(fig._label_counter)].set_bbox(dict(facecolor=letter_ffc, edgecolor = letter_fec))
        if showbeam:
            hheader = fits.open(beaminfoset)[0].header
            hb = math.sqrt(hheader['BMAJ']*hheader['BMIN'])
            hcdelt = header['CDELT1']
            if offsetunits == 'arcsec':
                hp = hb*3600./hcdelt
            else:
                hp = hb*60./hcdelt
            vb = hheader['CDELT3']
            vp = vb*vres/(1000.*header['CDELT2'])
            ax = matplotlib.pyplot.gca()
            #ax.set_aspect('equal')
            theartist = AnchoredEllipse(ax.transData, width=hp,
                                     height=vp, angle=0, loc='lower left',
                                     pad=0.5, borderpad=0.4,
                                        frameon=False)
            #            theartist.set_alpha(0.3)
            if beamfc == None:
                beamfa = 0.
            if beamfa == None:
                beamfa = 1.
            if beamec == None:
                beamea = 0.
            if beamea == None:
                beamea = 1.
            
            
            grofo = mpl.colors.to_rgba(beamfc)
            beamfc = (grofo[0], grofo[1], grofo[2], beamfa)
            grofo = mpl.colors.to_rgba(beamec)
            beamec = (grofo[0], grofo[1], grofo[2], beamea)

            theartist.ellipse.set(facecolor=beamfc)
            theartist.ellipse.set(edgecolor=beamec)
            #theartist.ellipse.set(alpha=0.5)

            fig.ax.add_artist(theartist)

        fig.close()

        if hide_ticl == True:
            hdu.data = 0.*hdu.data
            pnow = 1-borderleft-borderright
            pnoh = 1-bordertop-borderbottom
            fig = aplpy.FITSFigure(hdu, figure=figure, subplot=[borderleft+hide_ticl_min*pnow,borderbottom+hide_ticl_bottom*pnoh,hide_ticl_max*pnow-(hide_ticl_min*pnow), (hide_ticl_top-hide_ticl_bottom)*pnoh])
            fig.axis_labels.hide_x()
            fig.tick_labels.hide_x()
            fig.axis_labels.hide_y()
            fig.tick_labels.hide_y()
            fig.ticks.hide()
            fig.frame.set_color('white')
            fig.close()

        figure.savefig(name+plotpostfix)
        plt.close()
        
if __name__ == "__main__":
    import copy

    a = Slastroplot()

    # General setting to show the scale
    a.showscale = 4000.*u.pc
    a.distance = 7010000*u.pc

    # Common overlays
    eh = ['ellipses', {'pv_x_world': [358.01], 'pv_y_world': [-52.5550], 'pv_height': 0.0275, 'pv_width': 0.0375, 'angle': 30., 'edgecolor': 'black', 'linewidth': 2}]
    # Mark the centre
    mh = ['markers',{'pv_x_world': [358.0128034], 'pv_y_world': [-5.2576952E+01], 'marker':'x', 'linewidth': 2, 's': 75, 'c': '0.9'}]
    mhbl = copy.deepcopy(mh)
    mhbl[1]['c'] = '0.1'
    #['markers',{'pv_x_world': [358.0128034], 'pv_y_world': [-5.2576952E+01], 'marker':'x', 'linewidth': 2, 's': 75, 'c': '0.1'}]

    # Generally only one colourbar
    a.colourbar = False
    a.colourbar_label = None

    ###
    ###
    # Cubes
    ###
    ###

    # Velocity information and plotting
    a.plotvelo = 'eso149_original.fits'

    # The maps to plot
    a.basemap=['eso149_original.fits']
    a.individual_maps = False # Indicating if contours and pixelmap formats are for each map/plane or for all

    a.vmin = -0.004
    a.vmax = 0.08
    a.stretch = 'sqrt'
    a.cmap = None
    a.invert = True

    # The contours to plot
    a.contoursets = ['eso149_original.fits', 'eso149_model.fits']
    a.contourcols = [['#9494c3','DarkBlue','DarkBlue','DarkBlue','DarkBlue','DarkBlue'],['#f2e4e8','DeepPink','DeepPink','DeepPink','DeepPink','DeepPink']]
    rms = 0.001
    a.contourlevs = [[-2.*rms,2.*rms,4.*rms,8.*rms,16.*rms,32.*rms]]
    a.contourstyles = [['dashed','solid','solid','solid','solid','solid']]
    a.contourlinewidths = [2,1.4]

    ###
    ###
    # Cube 12 channels
    ###
    ###

    # Things private to this plot
    a.plotname = 'ESO149-G003_data_cube_12.pdf'
    a.chans = [2+2*i for i in range(12)]
    a.individual_olays = True
    a.olays=[[mh],[mh],[mh],[mh],[mh,eh],[mh,eh],[mh,eh],[mh],[mh],[mh],[mh],[mh]]
    a.nx = 4

    # Always do this
    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Cube 24 channels
    ###
    ###

    a.chans = [2+i for i in range(24)]
    a.olays=[[mh],[mh],[mh],[mh],[mh],[mh],[mh],[mh,eh],[mh,eh],[mh,eh],[mh,eh],[mh,eh],[mh,eh],[mh,eh],[mh]]
    a.plotname = 'ESO149-G003_data_cube_24.pdf'

    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Cube 24 channels, only observations
    ###
    ###

    a.chans = [2+i for i in range(24)]
    a.olays=[[mh]]
    # The contours to plot
    a.contoursets = ['eso149_original.fits']
    a.contourcols = [['#9494c3','DarkBlue','DarkBlue','DarkBlue','DarkBlue','DarkBlue']]
    rms = 0.001
    a.contourlevs = [[-2.*rms,2.*rms,4.*rms,8.*rms,16.*rms,32.*rms]]
    a.contourstyles = [['dashed','solid','solid','solid','solid','solid']]
    a.contourlinewidths = 2

    a.plotname = 'ESO149-G003_data_cube_24_obs.pdf'
    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Moment-0
    ###
    ###
    a = Slastroplot()

    # General setting to show the scale
    a.showscale = 4000.*u.pc
    a.distance = 7010000*u.pc

    # Generally only one colourbar
    a.colourbar = False
    a.colourbar_label = None

    a.basemap='eso149_or_mask_mom0_sol.fits'
    a.vmin = -0.1
    a.vmax = 14.
    a.stretch = 'sqrt'
    a.cmap = None
    a.invert = True

    a.contoursets = ['eso149_or_mask_mom0_sol.fits', 'eso149_model_or_mask_mom0_sol.fits']
    a.annotationfontsize = 'xx-large'
    a.contourcols = ['DarkBlue','DeepPink']
    a.contourlevs = [0.25,0.5,1.0,2.0,4.0,8.0]
    a.contourlinewidths = [2,1.4]

    centre_ra= -1.9871734
    centre_dec= -5.2576952E+01
    pa = +3.32631E+02
    len_maj = 40
    len_min = 30
    mindist = 8
    pvpos = getpvdiagramolay('eso149_in.fits', centre_ra, centre_dec, pa, len_maj, len_min, mindist, arrow_colour = '0.4', label_size = 28)
    a.olays = [eh]+pvpos

    a.plotname = 'ESO149-G003_mom0.pdf'

    a.nx = 1

    a.set_standard(border = True, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Top view
    ###
    ###
    a.basemap='eso149_model_topview_mom0_sol.fits'
    a.plotname = 'ESO149-G003_topview_mom0.pdf'
    a.olays = [mh]
    a.vmax = 7.
    a.contoursets = ['eso149_model_topview_mom0_sol.fits']
    a.contourcols = ['DeepPink']
    a.contourlevs = [0.25,0.5,1.0,2.0,4.0,8.0]
    a.contourlinewidths = 1.4

    a.nx = 1
    a.ny = 0
    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Velocity field on mom0
    ###
    ###
    # This is to create contours around vsys
    nvelconts = 3
    svelconts = 5.
    vsys = +5.78158E+02
    # These need to be real python float, not numpy, and a real list
    velcontourlevs = np.arange(vsys-nvelconts*svelconts,vsys+(nvelconts+1)*svelconts,svelconts).tolist()

    a.basemap='eso149_or_mask_mom0_sol.fits'
    a.vmin = -0.1
    a.vmax = 14.
    a.contoursets = ['eso149_or_mask_mom1.fits', 'eso149_model_or_mask_mom1.fits']
    a.contourlevs = [velcontourlevs]
    a.contourcols = ['DarkBlue','DeepPink']
    a.contourlinewidths = [2,1.4]
    a.plotname = 'ESO149-G003_mom0_mom1cont.pdf'
    a.nx = 1
    a.ny = 0
    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Velocity field on mom1
    ###
    ###
    # This is to create contours around vsys

    a.basemap='eso149_or_mask_mom1_2_v0.fits'
    a.vmin = -25.
    a.vmax = 25.
    a.contoursets = ['eso149_or_mask_mom1_2.fits', 'eso149_model_or_mask_mom1_2.fits']
    a.contourlevs = [velcontourlevs]
    a.contourcols = ['DarkBlue','DeepPink']
    a.plotname = 'ESO149-G003_mom1_vf.pdf'
    a.cmap = 'rainbow'
    a.olays = [eh,mh]
    a.colourbar = True
    a.stretch = 'linear'
    a.colourbar_label = '$v_\mathrm{{r}}-v_\mathrm{sys}\,(\mathrm{km}\,\mathrm{s}^{-1})$'
    a.nx = 1
    a.ny = 0
    a.set_standard(border = False, background = 'dark')
    a.labelmargin_right = 3.35*cm_in_inch
    plotmaps(**vars(a))

    ###
    ###
    # Combination Velocity moment 0 and velocity field on mom1
    ###
    ###
    # This is to create contours around vsys

    a.basemap=['eso149_or_mask_mom0_sol.fits', 'eso149_or_mask_mom1_2_v0.fits']
    a.individual_maps = True
    a.vmin = [-0.1, -25.]
    a.vmax = [14., 25.]
    a.contoursets = [['eso149_or_mask_mom0_sol.fits', 'eso149_model_or_mask_mom0_sol.fits'],['eso149_or_mask_mom1_2.fits', 'eso149_model_or_mask_mom1_2.fits']]
    a.contourlevs = [[0.25,0.5,1.0,2.0,4.0,8.0],[velcontourlevs]]
    a.contourcols = [['DarkBlue','DeepPink']]
    a.plotname = 'ESO149-G003_mom0_mom1_vf.pdf'
    a.cmap = [None, 'rainbow']
    pvpos = getpvdiagramolay('eso149_in.fits', centre_ra, centre_dec, pa, len_maj, len_min, mindist, arrow_colour = '0.4', label_size = 16)
    a.olays = [[eh]+pvpos,[eh,mh]]
    a.individual_olays = True
    a.colourbar = True
    a.stretch = ['sqrt', 'linear']
    a.colourbar_label = r'$v_\mathrm{{r}}-v_\mathrm{sys}\,(\mathrm{km}\,\mathrm{s}^{-1})$'
    a.nx = 2
    a.ny = 0

    #a.plotvelo = 'eso149_original.fits'
    #a.colourbar_label = None
    #a.colourbar = False
    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Combination Velocity moment 0 and velocity field on mom1, only observations
    ###
    ###

    a.basemap=['eso149_or_mask_mom0_sol.fits', 'eso149_or_mask_mom1_2_v0.fits']
    a.individual_maps = True
    a.vmin = [-0.1, -25.]
    a.vmax = [14., 25.]
    a.contoursets = [['eso149_or_mask_mom0_sol.fits'],['eso149_or_mask_mom1_2.fits']]
    a.contourlevs = [[0.25,0.5,1.0,2.0,4.0,8.0],[velcontourlevs]]
    a.contourcols = [['DarkBlue']]
    a.cmap = [None, 'rainbow']
    pvpos = getpvdiagramolay('eso149_in.fits', centre_ra, centre_dec, pa, len_maj, len_min, mindist, arrow_colour = '0.4', label_size = 16)
    a.olays = [[mh],[mhbl]]
    a.individual_olays = True
    a.colourbar = True
    a.stretch = ['sqrt', 'linear']
    a.colourbar_label = r'$v_\mathrm{{r}}-v_\mathrm{sys}\,(\mathrm{km}\,\mathrm{s}^{-1})$'
    a.nx = 2
    a.ny = 0

    #a.plotvelo = 'eso149_original.fits'
    #a.colourbar_label = None
    #a.colourbar = False
    a.plotname = 'ESO149-G003_mom0_mom1_vf_obs.pdf'
    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Difference maps
    ###
    ###

    a.basemap=['eso149_cdiff_mom0_sol.fits','eso149_or_mdiff_mom0_sol.fits','eso149_not_mask_mom0_sol.fits']
    a.vmin = [-0.01,-0.01,-0.01]
    a.vmax = [1.4,1.4,0.7]
    a.stretch = ['sqrt']
    a.cmap= [None]
    a.invert = [True]
    a.colourbar = False
    a.contoursets = [['eso149_cdiff_mom0_sol.fits'],['eso149_or_mdiff_mom0_sol.fits'],['eso149_not_mask_mom0_sol.fits']]
    a.contourcols = [[['#9494c3','#9494c3','#9494c3','DarkBlue','DarkBlue']]]
    a.contourlevs = [[-0.5,-0.25, 0., 0.25,0.5]]
    a.contourstyles=[[['dashed', 'dashed', 'solid', 'solid','solid']]]
    a.contourlinewidths = 2
    a.plotname = 'ESO149-G003_mom0diffs.pdf'
    a.nx = 3
    a.ny = 0
    a.individual_olays = False
    #print(vars(a).keys())

    mh = ['markers',{'pv_x_world': [358.0128034], 'pv_y_world': [-5.2576952E+01], 'marker':'x', 'linewidth': 2, 's': 75, 'c': '0.1'}]
    a.olays=[mh,eh]
    a.set_standard(border = False, background = 'dark')
    #a.plotmargin = 0.02
    plotmaps(**vars(a))

    ###
    ###
    # Optical images unconvolved
    ###
    ###

    eh = ['ellipses', {'pv_x_world': [358.01], 'pv_y_world': [-52.5550], 'pv_height': 0.0275, 'pv_width': 0.0375, 'angle': 30., 'edgecolor': 'white', 'linewidth': 2, 'alpha': 0.7}]
    mh = ['markers',{'pv_x_world': [358.0128034], 'pv_y_world': [-5.2576952E+01], 'marker':'x', 'linewidth': 1, 's': 12, 'c': '0.9'}]

    a.basemap = []
    a.olays = []
    for letter in ['u','g','q','r','i','z']:
        a.basemap.append('eso149_opt_me_'+letter+'_tri.fits')
        a.olays.append([mh,eh,['label', {'pv_x_world': 0.9, 'pv_y_world': 0.9, 'pv_label': letter, 'relative' : True, 'size' : 'x-large', 'weight' : 400, 'color': 'black', 'label_bbfc': 'white', 'label_bbfa': 0.8, 'label_bbec': '0.7', 'label_bbea': 1.}]])

    a.vmin = None
    a.vmax = None
    a.pmin = [8.,2.]
    a.pmax= [95]
    a.stretch = 'linear'

    a.contoursets = [['eso149_or_mask_mom0_sol.fits']]
    a.contourcols = [['DarkBlue']]
    a.contourlevs = [[1.0,4.0]]
    a.contourstyles = [['solid']]
    a.contouralphas = 0.6
    a.contourlinewidths =2

    a.showbeam = False

    a.individual_olays = True
    a.nx = 3
    a.ny= 0

    #a.plotvelo = 'eso149_original.fits'
    a.colourbar_label = None
    a.colourbar = False

    a.plotname = 'ESO149-G003_optical_meerlicht.pdf'

    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    #    subregion = [358.01+0.16,-52.5550-0.16,358.01-0.16,-52.5550+0.16]

    ###
    ###
    # Optical images convolved
    ###
    ###

    a.basemap = []
    for letter in ['u','g','q','r','i','z']:
        a.basemap.append('eso149_opt_mc_'+letter+'_tri.fits')
    a.plotname = 'ESO149-G003_optical_meerlicht_con_6.pdf'
    a.pmin = [5.,3.]
    a.pmax= [95]
    plotmaps(**vars(a))

    ###
    ###
    # Optical DSS 
    ###
    ###
    a.basemap='eso149_opt_de_r_tri.fits'
    a.pmin = [5.]
    a.pmax= [95.]
    a.plotname = 'ESO149-G003_optical_dss_r.pdf'
    a.olays = [[mh,eh]]

    a.nx = 1
    a.ny = 0
    a.set_standard(border = False, background = 'dark')
    plotmaps(**vars(a))

    ###
    ###
    # Optical DSS convolved
    ###
    ###
    a.basemap='eso149_opt_dc_r_tri.fits'
    a.plotname = 'ESO149-G003_optical_dss_r_con_6.pdf'
    plotmaps(**vars(a))


    ###
    ###
    # Optical DSS normal and convolved
    ###
    ###
    a.basemap=['eso149_opt_de_r_tri.fits', 'eso149_opt_dc_r_tri.fits']
    a.plotname = 'ESO149-G003_optical_dss_r_norm_con_6.pdf'
    a.nx = 2
    a.ny = 0
    a.set_standard(border = False, background = 'dark')
    a.nx = 1
    a.ny = 0
    a.width = a.width/2.
    plotmaps(**vars(a))

    ###
    ###
    # 
    ###
    ###

    pvdiagrams('eso149_in.fits', 'original', centre_ra, centre_dec, pa, len_maj, len_min, mindist)
    pvdiagrams('final_eso149_out.fits', 'finalmodel', centre_ra, centre_dec, pa, len_maj, len_min, mindist)

    ###
    ###
    # 
    ###

    # Now this is not perfect, but we first create a standard, to then change it to our needs. PVdiagrams are different from other plots

    a =  Slastroplot()
    a.nx = 3
    a.ny = 0
    a.showscale = 4000.*u.pc
    a.distance = 7010000*u.pc
    a.set_standard(border = False, background = 'dark')

    #keywordlist = vars(a)
    #for i in ['annotationfontsize', 'frameticklength','scalebarfontsize', 'scale_sep', 'scalebarcolor', 'scale_borderpad', 'figsize', 'letters', 'lettercol', 'plotprefix', 'plotpostfix', 'aspect', 'vres', 'showbeam', 'contoursets', 'contourcols', 'contourlevs', 'contourstyles', 'contoursets', 'postfixes', 'beaminfoset', 'letterfontsize', 'letterposx', 'letterposy', 'letter_ffc', 'letter_ffa', 'letter_fec', 'letter_fea', 'bgname_prefix']:
    #    if i in keywordlist:
    #        print('{} does exist, value: {}'.format(i, keywordlist[i]))
    #    else:
    # print('Hey {} does not exist'.format(i))

    a.scale_sep = 3.5
    figsize = a.figsize[0]
    a.figsize = (figsize/2, figsize/2)
    a.aspect = 'auto'

    a.bgname_prefix = 'original'
    a.postfixes = ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']

    a.vmin = -0.004
    a.vmax = 0.08
    a.stretch = 'sqrt'
    a.cmap = None
    a.invert = True
    a.colourbar = False

    a.contoursets = ['original', 'finalmodel']
    a.contourcols = [['#9494c3','DarkBlue','DarkBlue','DarkBlue','DarkBlue','DarkBlue'],['#f2e4e8','DeepPink','DeepPink','DeepPink','DeepPink','DeepPink']]
    rms = 0.001
    a.contourlevs = [[-2.*rms,2.*rms,4.*rms,8.*rms,16.*rms, 32.*rms],[-2.*rms,2.*rms,4.*rms,8.*rms,16.*rms, 32.*rms],]
    a.contourstyles = [['dashed','solid','solid','solid','solid','solid'],['dashed','solid','solid','solid','solid','solid']]
    a.contoursets = ['original', 'finalmodel']
    a.contourlinewidths = [2,1.4]

    a.letters = ['A', 'B', 'C', 'D']
    a.lettercol = 'black'

    a.vres = 2 # Velocity resolution for 'beam'
    a.beaminfoset = 'eso149_in.fits'

    a.letterfontsize = 24
    a.letterposx = 0.9
    a.letterposy = 0.9
    a.letter_ffc = a.scale_ffc
    a.letter_ffa = a.scale_ffa
    a.letter_fec = a.scale_fec
    a.letter_fea = a.scale_fea

    a.plotprefix = 'ESO149-G003'
    a.plotpostfix = '.pdf'

    a.bordertop = 0.
    a.borderbottom = 1.25
    a.borderleft = 1.25
    a.borderright = 0

    keywordlist = vars(a)
    # Clean up, pop out parameters that plotmaps does not know but plotpvdiagrams doesn't
    [keywordlist.pop(i) for i in ['colourbar_label','exponent','suppress_xlab','suppress_ylab', 'velofontsize', 'nx', 'velbbea', 'velbbfa', 'labelmargin_left', 'veloposx', 'width', 'plotpostfix', 'chans', 'labelmargin_bottom', 'individual_maps', 'ny', 'labelmargin_right', 'individual_olays', 'colourbar_pad', 'velbbec', 'vmid', 'plotmargin', 'velbbfc', 'basemap', 'olays', 'plotvelo', 'labelmargin_top', 'veloposy', 'colourbar_width']]
    plotpvdiagrams(**keywordlist)
