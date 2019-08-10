#! /usr/bin/env python
# Please use python3 and the aplPY from https://github.com/gigjozsa/aplpy
from astropy import log
log.setLevel('ERROR')

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

def putcontours(F, hdu, i, contourlevs, contourcols = None, contourstyles = None, usepixel = False):
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

    if type(contourstyles) == type([]):
        if i < len(contourstyles):
            styleshere = contourstyles[i]
        else:
            styleshere = contourstyles[-1]
    else:
        styleshere = contourstyles

    F.show_contour(data=hdu, levels = levelshere, colors = colorshere, linestyles = styleshere, usepixel = usepixel)
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

def plotprep_general(F, vmin = None, vmax= None, pmin = None, pmax = None, stretch = 'linear', cmap = None, invert = False, colorbar = False, frametickcolour = 'black', showscale = None, showscalelab = True, distance = 10, scaleunits = 'arcmin',  scale_borderpad = 1, scale_sep = 10, annotationfontsize = 'xx-large', suppress_xlab = False, suppress_ylab = False):
    """
    Prepare a plot, called by any plotting routine

    Input:
    vmin (float)            : Minimum for colourscale, see aplpy
    vmax (float)            : Maximum for colourscale, see aplpy
    pmin (float)            : Percentile minimum for colourscale, see aplpy
    pmax (float)            : Percentile maxiimum for colourscale, see aplpy
    stretch (str)           : Stretch scale, see aplpy
    cmap (str)              : Colour map, greyscale if None 
    invert (bool)           : Invert scaling, only valid if cmap == None
    colorbar (bool)         : Add colour bar if True 
    frametickcolour (str)   : Colour of frames and ticks
    showscale (None, astropy quantity): Show a physical and angular scale bar if either an astropy distance (e.g. 1000.*u.pc) or an angle (e.g. 1*u.arcmin) is given
    showscalelab (bool)     : Show a label for the scale bar if True 
    distance = 10,          : Distance of object to calculate scalebar length
    scaleunits = 'arcmin'   : Units used for the angle in the scale bar label
    scale_borderpad (float) : Distance of scalebar from the border
    scale_sep (float)       : Distance of scalebar from scale bar labelling
    annotationfontsize (str): General font size
    suppress_xlab (bool)    : Suppress labeling of x-axis (axis- and tick labels)
    suppress_ylab = False   : Suppress labeling of y-axis (axis- and tick labels)
    """
    
    
    F.set_system_latex(True)
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

    F.ticks.set_length(10)  # points
    F.ticks.set_color(frametickcolour)  # points


    kwargs = {'pmin': pmin, 'pmax': pmax, 'stretch': stretch}

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
        kwargs['cmap'] = cmap
        kwargs.pop('invert')
        F.show_colorscale(**kwargs)
        
#    print(kwargs)
    
    if colorbar:
        F.add_colorbar()
        F.colorbar.set_font(size=annotationfontsize, weight='medium', \
                            stretch='normal', family='serif', \
                            style='normal', variant='normal')

    F.ticks._ax.coords[F.ticks.x].ticks.set_tick_out(False)
    F.ticks._ax.coords[F.ticks.y].ticks.set_tick_out(False)

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
            text = '${0:.1f}\, \mathrm{{kpc}} \,\widehat{{=}} \,{1:.0f}{2:s}$'.format(physicalkpc, angular.value*scalehere, symbolhere)
        except:
            angular = showscale.to(u.arcsec)
            physical = (angular.value*arcsectokpc)
            text = '${1:.0f}{2:s}\,\widehat{{=}}\, {0:.1f}\, \mathrm{{kpc}}$'.format(physical, angular.value*scalehere, symbolhere)
        F.add_scalebar(angular)
        F.scalebar.show(angular, borderpad = scale_borderpad, sep = scale_sep)
#        F.scalebar.set_frame(True)
        F.scalebar.set_corner('bottom right')
        F.scalebar.set_linewidth(1)
        F.scalebar.set_font(size = annotationfontsize, weight = 'medium', stretch = 'normal', family = 'serif', style = 'normal', variant = 'normal')
        if showscalelab:
            F.scalebar.set_label(text)
    return

def plotmaps_prep(figure = None, figsize = None, subplot=[0.0,0.0,1.,1.], basemap = None, vmin = None, vmax= None, pmin = 0, pmax = 1, stretch = 'linear', cmap = None, invert = False, colorbar = False, annotationfontsize = 'x-large', suppress_xlab = False, suppress_ylab = False, contoursets = [], contourcols = None,  contourstyles = None, contourlevs = [], frametickcolour = 'black', showbeam = False, showscale = None, showscalelab = True, distance = 10, scaleunits = 'arcmin', scale_borderpad = 1, scale_sep = 10, plane = 0, plotvelo = None, veloposx = 0.25, veloposy = 0.9):
    """Plot one map

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
    colorbar (bool)                              : Add colour bar if True 
    annotationfontsize (str or number)           : General type size
    suppress_xlab (bool)                         : Suppress labeling of x-axis (axis- and tick labels)
    suppress_ylab = False                        : Suppress labeling of y-axis (axis- and tick labels)
    contoursets (list of str)                    : List of data sets to generate contours from, which will be overlaid on background image
    contourlevs (list of float/lists)            : List of contours valid for all data sets in contoursets or list of lists of contours valid for each member of contoursets. If fewer controur level lists are given than contoursets, the last given contours are repeated for consecutive contoursets.
    contourcols (None, list, or list of lists)   : Colours corresponding to contourlevs
    contourstyles (None, lists, or list of lists): Contour styles corresponding to contourlevs (see aplpy description, e.g. 'solid' or 'dashed')
    frametickcolour (str)                        : Colour of frames and ticks
    showbeam (bool)                              : Show beam if True 
    showscale (bool)                             : Show a physical and angular scale bar
    showscalelab (bool)                          : Show a label for the scale bar if True 
    distance = 10,                               : Distance of object to calculate scalebar length
    scaleunits (str)                             : Units used for the angle in the scale bar label, arcsec or arcmin
    plane                                        : Number of plane to be plotted/overlaid if input is a cube, starting at 0. Defaults to 0 if plane does not exist or file does not have a third dimension 
    plotvelo (None or str)                       : If a file is give, use the file to plot the velocity information corresponding to plane on the panel
    veloposx (float)                             : Relative position of velocity label on the panel, x
    veloposy (float)                             : Relative position of velocity label on the panel, y

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
            if cdeltvelo > 1000.:
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
        F = aplpy.FITSFigure(hdu_mom, figsize = figsize)
    else:
        F = aplpy.FITSFigure(hdu_mom, figure=figure, subplot=subplot)
    kwargs = {'vmin' : vmin, 'vmax': vmax, 'pmin' : pmin, 'pmax' : pmax, 'stretch' : stretch, 'cmap' : cmap, 'invert' : invert, 'colorbar' : colorbar, 'frametickcolour' : frametickcolour, 'annotationfontsize' : annotationfontsize, 'showscale' : showscale, 'distance' : distance, 'showscalelab' : showscalelab, 'scaleunits' : scaleunits, 'scale_borderpad': scale_borderpad, 'scale_sep': scale_sep, 'suppress_xlab' : suppress_xlab, 'suppress_ylab' : suppress_ylab}
    plotprep_general(F, **kwargs)

    if showbeam:
        F.add_beam(major=2, minor = 0.008)
        F.beam.show()
        F.beam.set_corner('bottom left')
        F.beam.set_frame(False)
        F.beam.set_facecolor('0.7')
        F.beam.set_edgecolor('black')
        F.beam.set_alpha(0.5)

    if plotvelo != None:
        velo = '${0:.1f}\,\mathrm{{km}}\,\mathrm{{s}}^{{-1}}$'.format(refvalvelo+cdeltvelo*(plane-refpixvelo+1))
        F.add_label(veloposx, veloposy, velo, relative = True, size=annotationfontsize, weight = 'medium', stretch = 'normal', family = 'serif', style = 'normal', variant = 'normal', horizontalalignment='center', verticalalignment='center', color = 'black')
        
    for i in range(len(contoursets)):
        fitsfile_here = fits.open(contoursets[i])
        hdu_here = read2dim(fitsfile_here[0], plane)

        # This is obsolete: OK, since this has taken me a day, I should describe it a bit. So in the method show_contour in aplpy there is a call to the method ax.get_transform(frame), where frame is the wcs frame of the data set. If this is a moment-0 map generated by Miriad, this has 3 dimensions which gets get_transform to fail. The local hack is to force usepixel to only use the pixel map, which means that the maps should have the same dimensions. Intrinsically one uses ax.contour(transform=ax.get_transform('pixel')) So aplpy has been modified for this...
#        F.show_contour(hdu_here, levels = contourlevs, colors = contourcols[i], usepixel = True)
        putcontours(F, hdu_here, i, contourlevs = contourlevs, contourcols = contourcols, contourstyles = contourstyles, usepixel = False)
        #F.show_contour(hdu_here, levels = levelshere, colors = colorshere, linestyles = styleshere)
        
    F.ticks._ax.coords[F.ticks.x].ticks.set_tick_out(False)
    F.ticks._ax.coords[F.ticks.y].ticks.set_tick_out(False)
    F.frame.set_color(frametickcolour)
    return F

def plotmaps(width = 8.27, plotmargin = 0.1968504, labelmargin_left = 1.377953, labelmargin_bottom = 0.6889764, vmin = None, vmax= None, pmin = None, pmax = None, stretch = 'linear', cmap = None, invert = False, chans = None, nx = 0, ny =0, showscale = None, showbeam = True, showscalelab = True, reducescalelab = True, plotname = None, annotationplotsize = 'medium', individual = False, contoursets = None, contourlevs = None, contourcols = None, contourstyles = None, **kwargs):
    """
    Plot a grid of maps
    width (float)             : Width in inches, 8.27 is A4 width
    plotmargin (float)        : A margin around the whole plot, 0.1968504 is 0.5 cm
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
    colorbar (bool or list of bool)              : (List of) Add colour bar if True 
    annotationfontsize (str or number)           : General type size
    suppress_xlab (bool)                         : Generally suppress labeling of x-axis (axis- and tick labels)
    suppress_ylab (bool)                         : Generally suppress labeling of y-axis (axis- and tick labels)
    individual (bool)                            : Instead of applying contoursets, contourlevs, contourcols, contourstyles for all data sets given, the corresponding quantities are given for each input cube in a list
    contoursets (list of str or list of lists)   : List of data sets to generate contours from, which will be overlaid on background images
    contourlevs (list of float/lists)            : List of contours valid for all data sets in contoursets or list of lists of contours valid for each member of contoursets. If fewer controur level lists are given than contoursets, the last given contours are repeated for consecutive contoursets.
    contourcols (None, list, or list of lists)   : Colours corresponding to contourlevs
    contourstyles (None, lists, or list of lists): Contour styles corresponding to contourlevs (see aplpy description, e.g. 'solid' or 'dashed')
    frametickcolour (str)                        : Colour of frames and ticks
    showbeam (bool)                              : Show beam if True 
    showscale (bool)                             : Show a physical and angular scale bar
    distance (float)                             : Distance of object to calculate scalebar length
    scaleunits (str)                             : Units used for the angle in the scale bar label, arcsec or arcmin
#    plane (None, or list of int)                 : Number of plane to be plotted/overlaid if input is a cube, starting at 0. Defaults to 0 if plane does not exist or file does not have a third dimension 
    plotvelo (None or str)                       : If a file is give, use the file to plot the velocity information corresponding to plane on the panel
    veloposx (float)                             : Relative position of velocity label on the panel, x
    veloposy (float)                             : Relative position of velocity label on the panel, y
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
    stretch, cmap, invert, colorbar can be given as lists and their
    nth elements will apply to nth panel, extrapolting from the last
    element.  If a plane does not exist, the last plane of the cube is
    shown, if a third axis does not exist, the image will be show. The
    list of data sets to overplot contours contoursets is valid for
    all panels unless the parameter individual is set to True. In that case contoursets, contourlevs, contourcols, contourstyles are lists of the single quantities (or further lists), for each channel or basemap individually, with extrapolation being applied. Also here, the single planes are shown as listed in
    chans.

    """
    print('Plotmaps')

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
    swidth = (width-labelmargin_left-2.*plotmargin)/float(nx)

    # Height of a single panel, in inches
    sheight = swidth/x_over_y
                 
    # Plot height in inches
    plotheight = labelmargin_bottom+2.*plotmargin+ny*sheight

    fig = plt.figure(figsize=(width, plotheight))

    # y-spaces, in units of plot height
    labmarg_over_height = labelmargin_bottom/plotheight
    pnoh = sheight/plotheight
    plotmarg_over_height = plotmargin/plotheight
    btoh = labmarg_over_height+plotmarg_over_height
  
    # x-spaces, in units of plot width
    labmarg_over_width = labelmargin_left/width
    pnow = swidth/width
    plotmarg_over_width = plotmargin/width
    lfow = labmarg_over_width+plotmarg_over_width
  

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
    vmin_intro    = vmin   
    vmax_intro    = vmax   
    pmin_intro    = pmin   
    pmax_intro    = pmax   
    stretch_intro = stretch
    cmap_intro    = cmap   
    invert_intro  = invert 

    a = showscale
    b = showbeam
    kwargs['showscale'] = None
    kwargs['showbeam'] = False
    kwargs['vmin'] = -10000002.
    kwargs['vmax'] = -10000001.
    F = plotmaps_prep(figure = fig, subplot=[lfow, btoh, pnow, pnoh*float(ny)], plane = 0, **kwargs)
    kwargs['showscale'] = a
    kwargs['showbeam'] = b
    
    #F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow, btoh, pnow, pnoh*float(ny)])
    F.axis_labels.hide_x()
    F.tick_labels.hide_x()
    F.ticks.hide()
    F.axis_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'black')
    F.tick_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'white')
    F.frame.set_color('white')
    F.close()

    a = showscale
    b = showbeam
    kwargs['showscale'] = None
    kwargs['showbeam'] = False
    kwargs['vmin'] = -10000002.
    kwargs['vmax'] = -10000001.
    F = plotmaps_prep(figure = fig, subplot=[lfow, btoh, pnow*float(nx), pnoh], plane = 0, **kwargs)
    kwargs['showscale'] = a
    kwargs['showbeam'] = b
    #    F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow, btoh, pnow*float(nx), pnoh])
    F.axis_labels.hide_y()
    F.tick_labels.hide_y()
    F.ticks.hide()
    F.axis_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'black')
    F.tick_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'white')
    F.frame.set_color('white')
    F.close()

    contoursets_intro = contoursets
    contourlevs_intro = contourlevs
    contourcols_intro = contourcols
    contourstyles_intro = contourstyles
    vmin_intro    = vmin   
    vmax_intro    = vmax   
    pmin_intro    = pmin   
    pmax_intro    = pmax   
    stretch_intro = stretch
    cmap_intro    = cmap   
    invert_intro  = invert 

    # Now loop over the panels
    for i in range(len(chans)):
        posnx = i % nx
        posny = ny - i // nx - 1
        
        if posnx == 0:
            showscalelab_here = showscalelab
        else:
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
                  ['vmin',          vmin_intro],
                  ['vmax',          vmax_intro],
                  ['pmin',          pmin_intro],
                  ['pmax',          pmax_intro],
                  ['stretch',       stretch_intro],
                  ['cmap',          cmap_intro],
                  ['invert',        invert_intro]
                  ]:
            if individual:
                if i < len(j[1]):
                    kwargs[j[0]] = j[1][i]
                else:
                    kwargs[j[0]] = j[1][-1]
            else:
                kwargs[j[0]] = j[1]

        F = plotmaps_prep(figure = fig, subplot=[lfow+posnx*pnow, btoh+posny*pnoh, pnow, pnoh], plane = chans[i], showscalelab = showscalelab_here, **kwargs)
        # F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow+posnx*pnow, btoh+posny*pnoh, pnow, pnoh])
        F.axis_labels.hide_x()
        F.axis_labels.hide_y()
        if posnx != 0:
            F.tick_labels.hide_y()
        # Leave tick labels if there is no panel below
        if (posny != 0) and nx*(i//nx+1)+posnx < len(chans):
            F.tick_labels.hide_x()
        F.close()

    #    fig.canvas.draw()
    if plotname != None:
        fig.savefig(plotname)
    return    

def plotpvdiagrams(bgname_prefix = '', vmin = None, vmax = None, pmin = None, pmax= None, stretch = 'default', cmap=None, invert = True, colorbar = False, annotationfontsize = 'x-large', contoursets = [], postfixes = ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits'], letters = ['A', 'B', 'C', 'D'], lettercol = 'black', contourcols = None, contourstyles = None, contourlevs = [], frametickcolour = '#555555', showbeam = True, beaminfoset = '', vres = 2, showscale = 1*u.arcmin, scaleunits = 'arcmin', distance = 7010000*u.pc, plotprefix = '', plotpostfix = '.pdf', figsize = None):
    """
    Plot a number of PV diagrams

    bgname_prefix (str)                          : Name prefix of background fits files
    vmin (float)                                 : Minimum for colourscale, see aplpy
    vmax (float)                                 : Maximum for colourscale, see aplpy
    pmin (float)                                 : Percentile minimum for colourscale, see aplpy
    pmax (float)                                 : Percentile maxiimum for colourscale, see aplpy
    stretch (str)                                : Stretch scale, see aplpy
    cmap (str)                                   : Colour map, greyscale if None 
    invert (bool)                                : Invert scaling, only valid if cmap == None
    colorbar (bool)                              : Add colour bar if True 
    annotationfontsize (str or number)           : General type size
    contoursets (list of str)                    : List of data sets to generate contours from, which will be overlaid on background image
    contourlevs (list of float/lists)            : List of contours valid for all data sets in contoursets or list of lists of contours valid for each member of contoursets. If fewer controur level lists are given than contoursets, the last given contours are repeated for consecutive contoursets.
    contourcols (None, list, or list of lists)   : Colours corresponding to contourlevs
    contourstyles (None, lists, or list of lists): Contour styles corresponding to contourlevs (see aplpy description, e.g. 'solid' or 'dashed')
    postfixes (list of str)                      : Name postfix of background fits files
    letters (list of str)                        : List of strings to annotate output plots with
    lettercol (str)                              : Colour of strings to annotate output plots with
    frametickcolour (str)                        : Colour of frames and ticks
    showbeam (bool)                              : Show beam if True 
    beaminfoset (str)                            : Name of data set to generate beam from 
    vres = 2                                     : Velocity resolutio in pixels
    showscale (bool)                             : Show a physical and angular scale bar
    scaleunits (str)                             : Units used for the angle in the scale bar label, arcsec or arcmin
    distance (float)                             : Distance of object to calculate scalebar length
    plotprefix (str)                             : Prefix to output plots
    plotpostfix (str)                            : Postfix (ending) of output plots, determines output file type (e.g. '.pdf', '.png')
    figsize (float tuple)                        : Size of figure in inches, x and y

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
    
    for posnum in range(len(postfixes)):
        name = postfixes[posnum]
        letter = letters[posnum]
        hdulist = fits.open(bgname_prefix+name)
        header = hdulist[0].header
        image = hdulist[0].data
        hdulist.close()
#        sys.exit()
        
        
        # region of image we want to display
        header['CDELT1'] = header['CDELT1']*3600.
        header['CDELT2'] = header['CDELT2']/1000.
        header['CRVAL1'] = header['CRVAL1']*3600.
        header['CRVAL2'] = header['CRVAL2']/1000.
        header.set('BMAJ', 2)
        header.set('BMIN', 0.008)
        header.set('BPA', 0.)
        hdu = fits.PrimaryHDU(data=image,header=header)
        fig = aplpy.FITSFigure(hdu, figsize = figsize)
        fig.set_system_latex(True)

        kwargs = {'vmin' : vmin, 'vmax': vmax, 'pmin' : pmin, 'pmax' : pmax, 'stretch' : stretch, 'cmap' : cmap, 'invert' : invert, 'colorbar' : colorbar, 'frametickcolour' : frametickcolour, 'showscale' : showscale, 'distance' : distance, 'scaleunits' : scaleunits, 'annotationfontsize' : annotationfontsize}
        plotprep_general(fig, **kwargs)

        fig.set_xaxis_coord_type('scalar')
        
        fig.axis_labels.set_xtext(r'OFFSET (arcsec)')
        fig.axis_labels.set_ytext(r'VELOCITY (km$\,$s$^{-1}$)')

        name = plotprefix+'.'.join(postfixes[posnum].split('.')[:-1])+'_'+bgname_prefix
        for contournumber in range(len(contoursets)):
            name = name+'_'+contoursets[contournumber]
            hdulist = fits.open(contoursets[contournumber]+postfixes[posnum])
            header = hdulist[0].header
            image = hdulist[0].data
            hdulist.close()
            # region of image we want to display
            header['CDELT1'] = header['CDELT1']*3600.
            header['CDELT2'] = header['CDELT2']/1000.
            header['CRVAL1'] = header['CRVAL1']*3600.
            header['CRVAL2'] = header['CRVAL2']/1000.

            hdu = fits.PrimaryHDU(data=image,header=header)
            putcontours(fig, hdu, contournumber, contourlevs, contourcols = contourcols, contourstyles = contourstyles, usepixel = True)
            #fig.show_contour(data=hdu, colors=contourcols[contournumber], levels=contourlevs, usepixel = True)
            
        fig.add_label(0.95, 0.95, letters[posnum], relative = True, size=32, weight='bold', horizontalalignment='center', verticalalignment='center', color=lettercol)
        if showbeam:
            hheader = fits.open(beaminfoset)[0].header
            hb = math.sqrt(hheader['BMAJ']*hheader['BMIN'])
            hcdelt = header['CDELT1']
            hp = hb*3600./hcdelt
            vb = hheader['CDELT3']
            vp = vb*vres/(1000.*header['CDELT2'])
            ax = matplotlib.pyplot.gca()
            ax.set_aspect('equal')
            theartist = AnchoredEllipse(ax.transData, width=hp,
                                     height=vp, angle=0, loc='lower left',
                                     pad=0.5, borderpad=0.4,
                                        frameon=False)
#            theartist.set_alpha(0.3)
            theartist.ellipse.set(facecolor='0.7')
            theartist.ellipse.set(edgecolor='black')
            theartist.ellipse.set(alpha=0.5)

            fig.ax.add_artist(theartist)

        fig.save(name+plotpostfix)
        fig.close()
        
if __name__ == "__main__":

    # Plot data cube
    rms = 0.001
    plotmaps(basemap=['eso149_original.fits'], vmin = -0.004, vmax = 0.08, stretch = 'sqrt', cmap=None, invert = True, colorbar = False, contoursets = ['eso149_original.fits', 'eso149_model.fits'], contourcols = [['#9494c3','DarkBlue','DarkBlue','DarkBlue','DarkBlue'],['#f2e4e8', 'DeepPink', 'DeepPink', 'DeepPink', 'DeepPink']], contourlevs = [[-2.*rms,2.*rms,4.*rms,8.*rms,16.*rms]], contourstyles = [['dashed','solid','solid','solid','solid']], frametickcolour = '#555555', showbeam = True, showscale = 3000.*u.pc, showscalelab = True, reducescalelab = True, scaleunits = 'arcsec', distance = 7010000*u.pc, scale_borderpad = 0.5, scale_sep = 5, plotname = 'ESO149-G003_data_cube_12.pdf', chans = [2+2*i for i in range(12)], annotationfontsize = 'small', plotvelo = 'eso149_original.fits', nx = 3, labelmargin_left = 0.6889765, labelmargin_bottom = 0.3444883)

    plotmaps(basemap=['eso149_original.fits'], vmin = -0.004, vmax = 0.08, stretch = 'sqrt', cmap=None, invert = True, colorbar = False, contoursets = ['eso149_original.fits', 'eso149_model.fits'], contourcols = [['#9494c3','DarkBlue','DarkBlue','DarkBlue','DarkBlue'],['#f2e4e8', 'DeepPink', 'DeepPink', 'DeepPink', 'DeepPink']], contourlevs = [[-2.*rms,2.*rms,4.*rms,8.*rms,16.*rms]], contourstyles = [['dashed','solid','solid','solid','solid']], frametickcolour = '#555555', showbeam = True, showscale = 3000.*u.pc, showscalelab = True, reducescalelab = True, scaleunits = 'arcsec', distance = 7010000*u.pc, scale_borderpad = 0.5, scale_sep = 5, plotname = 'ESO149-G003_data_cube_24.pdf', chans = [2+i for i in range(24)], annotationfontsize = 'medium', plotvelo = 'eso149_original.fits', veloposx = 0.3, veloposy = 0.85, nx = 4, labelmargin_left = 0.6889765, labelmargin_bottom = 0.3444883)

    plotmaps(basemap='eso149_or_mask_mom0_sol.fits', vmin = -0.1, vmax = 14., pmin = 0.25, pmax= 99.75, stretch = 'sqrt', cmap= None, invert = True, colorbar = False, contoursets = ['eso149_or_mask_mom0_sol.fits', 'eso149_model_or_mask_mom0_sol.fits'], annotationfontsize = 'xx-large', contourcols = ['DarkBlue','DeepPink'], contourlevs = [0.25,0.5,1.0,2.0,4.0,8.0], frametickcolour = '#555555', showbeam = True, showscale = 4000.*u.pc, distance = 7010000*u.pc, plotname = 'ESO149-G003_mom0.pdf')

    # This is to create contours around vsys
    nvelconts = 3
    svelconts = 5.
    vsys = +5.78158E+02
    # These need to be real python float, not numpy, and a real list
    velcontourlevs = np.arange(vsys-nvelconts*svelconts,vsys+(nvelconts+1)*svelconts,svelconts).tolist()

    plotmaps(basemap='eso149_or_mask_mom0_sol.fits', vmin = -0.1, vmax = 14., pmin = 0.25, pmax= 99.75, stretch = 'sqrt', cmap= None, invert = True, colorbar = False, contoursets = ['eso149_or_mask_mom1.fits', 'eso149_model_or_mask_mom1.fits'], annotationfontsize = 'xx-large', contourcols = ['DarkBlue','DeepPink'], contourlevs = [velcontourlevs], frametickcolour = '#555555', showbeam = True, showscale = 4000.*u.pc, distance = 7010000*u.pc, plotname = 'ESO149-G003_mom1.pdf')
    
    plotmaps(basemap=['eso149_cdiff_mom0_sol.fits','eso149_or_mdiff_mom0_sol.fits','eso149_not_mask_mom0_sol.fits'], vmin = [-0.01,-0.01,-0.01], vmax = [1.4,1.4,0.7], pmin = [0.25], pmax= [99.75], stretch = ['sqrt'], cmap= [None], invert = [True], colorbar = False, individual = True, contoursets = [['eso149_cdiff_mom0_sol.fits'],['eso149_or_mdiff_mom0_sol.fits'],['eso149_not_mask_mom0_sol.fits']], annotationfontsize = 'medium', contourcols = [[['#9494c3','#9494c3','#9494c3','DarkBlue','DarkBlue']]], contourlevs = [[-0.5,-0.25, 0., 0.25,0.5]], contourstyles=[[['dashed', 'dashed', 'dashed', 'solid','solid']]], frametickcolour = '#555555', showbeam = True, showscale = 4000.*u.pc, distance = 7010000*u.pc, plotname = 'ESO149-G003_mom0diffs.pdf', nx = 3, labelmargin_left = 0.6889765, labelmargin_bottom = 0.3444883, scale_borderpad = 0.5, scale_sep = 5)
    
    
#    F = aplpy.FITSFigure('eso149_in.fits',slices=[29])
#    F.show_grayscale()
    # matplotlib.pyplot.show()
    
    mom0 = 'bla_0.fits'
    centre_ra= -1.98719660E+00
    centre_dec= -5.25764076E+01
    pa = +3.33388E+02
    len_maj = 40
    len_min = 30
    mindist = 8
    
    plotpvdiagrampositions('eso149_in.fits', centre_ra, centre_dec, pa, len_maj, len_min, mindist, '0.4', basemap='eso149_or_mask_mom0_sol.fits', plane = 0, vmin = -0.1, vmax = 14., pmin = 0.25, pmax= 99.75, stretch = 'sqrt', cmap=None, invert = True, colorbar = False, annotationfontsize = 'xx-large', contoursets = ['eso149_or_mask_mom0_sol.fits', 'eso149_model_or_mask_mom0_sol.fits'], contourcols = ['DarkBlue','DeepPink'], contourlevs = [0.25,0.5,1.0,2.0,4.0,8.0], frametickcolour = '#555555', showbeam = True, showscale = 4000*u.pc, distance = 7010000*u.pc, plotname = 'ESO149-G003_slicepositions.pdf', figsize = (8.27,8))

    pvdiagrams('eso149_in.fits', 'original', centre_ra, centre_dec, pa, len_maj, len_min, mindist)
    pvdiagrams('final_eso149_out.fits', 'finalmodel', centre_ra, centre_dec, pa, len_maj, len_min, mindist)

    # Most explains itself, but the ellipse height is set to vres times the channel width in the beaminfoset
    plotpvdiagrams(bgname_prefix = 'original', figsize = (8,8), vmin = None, vmax = 0.03, pmin = 0.25, pmax= 99.75, stretch = 'arcsinh', cmap=None, invert = True, colorbar = False, annotationfontsize = 'xx-large', contoursets = ['original','finalmodel'], postfixes = ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits'], letters = ['A', 'B', 'C', 'D'], lettercol = 'red', contourcols = ['DarkBlue','DeepPink'], contourlevs = [0.002,0.004,0.008,0.016,0.032], contourstyles=['solid','solid'], frametickcolour = '#555555', showbeam = True, beaminfoset = 'eso149_in.fits', vres = 2, showscale = 2000.*u.pc, scaleunits = 'arcsec', distance = 7010000*u.pc, plotprefix = 'ESO149-G003', plotpostfix = '.pdf')

    
