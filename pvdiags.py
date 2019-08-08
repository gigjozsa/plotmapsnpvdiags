#! /usr/bin/env python
# Please use python3
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

def read2dim(hdu,plane):
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
    
    if type(contourlevs[0]) == type(1.) or type(contourlevs[0]) == type(1):
            levelshere = contourlevs
    else:
        if i < len(contourlevs):
            levelshere = contourlevs[i]
        else:
            levelshere = contourlevs[-1]


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
    Will produce four pv diagrams, one along the "major" axis, defined by the central coordinates centra (deg), centdec (deg), position angle (pa in degrees), and length majlen in arcsec, three perpendicular with respect to the major axis, one through the centre, two (at each side) at a distance of mindist (in arcsec) from the centre, of the length minlen. Output PV diagrams will be called outname+'_pvmaj.fits', outname+'_pvmin_cent.fits', outname+'_pvmin_left.fits', outname+'_pvmin_right.fits'. Position angle should be defined as the receding side major axis.
    """
    print('Pvdiagrams')
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
    
    # Start mom0 plot
#    F = aplpy.FITSFigure(hdu_mom,slices=[0])
#    F.show_grayscale()
#    F.axis_labels.set_font(size='x-large', weight='medium', \
#                      stretch='normal', family='serif', \
#                      style='normal', variant='normal')
#    F.tick_labels.set_font(size='x-large', weight='medium', \
#                         stretch='normal', family='serif', \
#                         style='normal', variant='normal')
#    F.ticks.set_length(10)  # points

    i = -1
#    letters = ['A','B','C','D']
    for names in ['_pvmaj.fits', '_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']:
        i = i + 1
        endpoints = [(endpoints_x[2*i],endpoints_y[2*i]),(endpoints_x[2*i+1],endpoints_y[2*i+1])]
        endpointst = [(endpoints_xt[2*i],endpoints_yt[2*i]),(endpoints_xt[2*i+1],endpoints_yt[2*i+1])]

        xy = Path(endpoints)
#        print 'hoho'
#        print xy
        pv = extract_pv_slice(hdu, xy)
    
        header = pv.header
        #    print header
        pixels =  header['NAXIS1']
        pv.header['CRPIX1'] = pixels/2
        pv.header['CDELT1'] = pv.header['CDELT1']
        #print header
        #F2 = aplpy.FITSFigure(pv)
        #F2.show_grayscale(aspect='auto', stretch='arcsinh')
        #F2.save("testpv.png")
        pv.writeto(outname+names, clobber = True)
    
#        endpoints_wcs = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[x,y] for x,y in endpoints], 0)
#        endpoints_wcst = w.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[x,y] for x,y in endpointst], 0)
#        print(endpoints_wcs.T)
    
        # F.show_lines([endpoints_wcs.T],color='w')

#        xs = endpoints_wcs.T[0,0]
#        ys = endpoints_wcs.T[1,0]
#        exs = endpoints_wcs.T[0,1]-endpoints_wcs.T[0,0]
#        eys = endpoints_wcs.T[1,1]-endpoints_wcs.T[1,0]
#        F.show_arrows(xs,ys,exs,eys,color='w', width = 0.2, head_width = 1.5, head_length = 1.5, length_includes_head = True)
#        for i in range(len(xs)):
#        F.add_label(endpoints_wcst.T[0,0], endpoints_wcst.T[1,0], letters[i], relative = False, size=24, weight=400, horizontalalignment='center', verticalalignment='center', color='w')
#    print plotoutname
#    F.save(plotoutname)
#    F.close()

def plotpvdiagrampositions(filename, centra, centdec, pa, majlen, minlen, mindist, arrowcolour, **kwargs):
    """
    Will draw lines onto an image of overplotname one along the "major" axis, defined by the central coordinates centra (deg), centdec (deg), position angle (pa in degrees), and length majlen in arcsec, three perpendicular with respect to the major axis, one through the centre, two (at each side) at a distance of mindist (in arcsec) from the centre, of the length minlen. Plotoutname is the name of the plot of file overplotname together with the slice positions. Position angle should be defined as the receding side major axis.
    """
#    return
    print('Plotpvdiagrampositions')
    plt.rc('text', usetex=True)

    fitsfile = fits.open(filename)
    hdu = fitsfile[0]
    w = wcs.WCS(hdu.header, fitsfile)

#    fitsfile_mom = fits.open(basemap)
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
    
    # Start mom0 plot
    F = plotmaps_prep(**kwargs)

    i = -1
    letters = ['A','B','C','D']
    for names in ['_pvmaj.fits', '_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits']:
        i = i + 1
#        print('{:s}: {:s}'.format(letters[i],names))
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
    F.save(kwargs['plotname'])
    F.close()

def plotprep_general(F, vmin = None, vmax= None, pmin = 0, pmax = 1, stretch = 'linear', cmap = '', invert = False, colorbar = False, frametickcolour = 'black', showscale = None, showscalelab = True, distance = 10, scaleunits = 'arcmin', aspect = None, annotationfontsize = 'xx-large', suppress_xlab = False, suppress_ylab = False):
    
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
    if aspect != None:
        kwargs['aspect'] = aspect

    if invert == True:
        kwargs['invert'] = True
        
    if cmap == '':
        F.show_grayscale(**kwargs)
    else:
        kwargs['cmap'] = cmap
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
        F.scalebar.show(angular, borderpad = 1, sep = 10)
#        F.scalebar.set_frame(True)
        F.scalebar.set_corner('bottom right')
        F.scalebar.set_linewidth(1)
        F.scalebar.set_font(size = annotationfontsize, weight = 'medium', stretch = 'normal', family = 'serif', style = 'normal', variant = 'normal')
        if showscalelab:
            F.scalebar.set_label(text)
    return

def plotmaps_prep(figure = None, subplot=[0.0,0.0,1.,1.], basemap = '', vmin = None, vmax= None, pmin = 0, pmax = 1, stretch = 'linear', cmap = '', invert = False, colorbar = False, annotationfontsize = 'x-large', suppress_xlab = False, suppress_ylab = False, contoursets = [], contours = [], contourcols = None,  contourstyles = None, contourlevs = [], frametickcolour = 'black', showbeam = False, showscale = None, showscalelab = True, distance = 10, scaleunits = 'arcmin', plotname = '', plane = 0, plotvelo = None):
    """
    Plot moment maps on a base map
    """
    print('Plotmaps_prep')

    if basemap == '':
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
        F = aplpy.FITSFigure(hdu_mom)
    else:
        F = aplpy.FITSFigure(hdu_mom, figure=figure, subplot=subplot)
    kwargs = {'vmin' : vmin, 'vmax': vmax, 'pmin' : pmin, 'pmax' : pmax, 'stretch' : stretch, 'cmap' : cmap, 'invert' : invert, 'colorbar' : colorbar, 'frametickcolour' : frametickcolour, 'annotationfontsize' : annotationfontsize, 'showscale' : showscale, 'distance' : distance, 'showscalelab' : showscalelab, 'scaleunits' : scaleunits, 'suppress_xlab' : suppress_xlab, 'suppress_ylab' : suppress_ylab}
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
        F.add_label(0.2, 0.9, velo, relative = True, size=annotationfontsize, weight = 'medium', stretch = 'normal', family = 'serif', style = 'normal', variant = 'normal', horizontalalignment='center', verticalalignment='center', color = 'black')
        
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

def plotmaps(width = 8.27, plotmargin = 0.1968504, labelmargin_left = 1.377953, labelmargin_bottom = 0.6889764, chans = None, nx = 0, ny =0, showscalelab = True, reducescalelab = True, **kwargs):
    # width is width in inches, 8.27 is A4 width
    # plotmargin: a margin around the whole plot, 0.1968504 is 0.5 cm
    # labelmargin_bottom/left: an additional margin to the left and the bottom for labels, 0.6889764 is 1.75 cm, 1.377953 is 3.5 cm
    # chans None: all channels, provide list otherwise
    # nx = 0: auto-generate
    # nx = 1: auto-generate

    # labwidth is the room for labels, we fix it to 0.1 for now
    # labwidth_over_width = 0.1

    # Do this by hand, for some reason this needs to be done
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
    F = plotmaps_prep(figure = fig, subplot=[lfow, btoh, pnow, pnoh*float(ny)], plane = 0, **kwargs)
    #F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow, btoh, pnow, pnoh*float(ny)])
    F.axis_labels.hide_x()
    F.tick_labels.hide_x()
    F.ticks.hide()
    F.axis_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'black')
    F.tick_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'white')
    F.frame.set_color('white')
    F.close()

    F = plotmaps_prep(figure = fig, subplot=[lfow, btoh, pnow*float(nx), pnoh], plane = 0, **kwargs)
#    F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow, btoh, pnow*float(nx), pnoh])
    F.axis_labels.hide_y()
    F.tick_labels.hide_y()
    F.ticks.hide()
    F.axis_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'black')
    F.tick_labels.set_font(size=kwargs['annotationfontsize'], weight='medium', stretch='normal', family='serif', style='normal', variant='normal', color = 'white')
    F.frame.set_color('white')
    F.close()

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
                
        F = plotmaps_prep(figure = fig, subplot=[lfow+posnx*pnow, btoh+posny*pnoh, pnow, pnoh], plane = chans[i], showscalelab = showscalelab_here, **kwargs)
        # F = aplpy.FITSFigure('bla_0.fits', figure=fig, subplot=[lfow+posnx*pnow, btoh+posny*pnoh, pnow, pnoh])
        F.axis_labels.hide_x()
        F.axis_labels.hide_y()
        if posnx != 0:
            F.tick_labels.hide_y()
        # Leave tick labels if there is no panel below
        if (posny != 0) and nx*(i//nx+1)+posnx < len(chans):
            F.tick_labels.hide_x()
#        F.show_grayscale()
#        F.add_label(0.95, 0.95, '{:d}'.format(i), relative = True, size=32, weight='bold', horizontalalignment='center', verticalalignment='center', color='w')
        F.close()

    #    fig.canvas.draw()
#    fig.savefig('blatest.pdf')
    fig.savefig(kwargs['plotname'])
    return    

def plotpvdiagrams(greyname_prefix = '', vmin = -0.1, vmax = 14., pmin = 0.25, pmax= 99.75, stretch = 'sqrt', aspect = 'auto', cmap='', invert = True, colorbar = False, annotationfontsize = 'x-large', contoursets = ['original','finalmodel'], postfixes = ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits'], letters = ['A', 'B', 'C', 'D'], lettercol = 'black', contourcols = None, contourstyles = None, contourlevs = [0.002,0.004,0.008], frametickcolour = '#555555', showbeam = True, beaminfoset = '', vres = 2, showscale = 2*u.arcmin, scaleunits = 'arcmin', distance = 7010000*u.pc, plotprefix = 'ESO149-G003_', plotpostfix = '.pdf'):
    print('Plotpvdiagrams')
    plt.rc('text', usetex=True)

    levels = contourlevs
    
    for posnum in range(len(postfixes)):
        name = postfixes[posnum]
        letter = letters[posnum]
        hdulist = fits.open(greyname_prefix+name)
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
        fig = aplpy.FITSFigure(hdu)
        fig.set_system_latex(True)

        kwargs = {'vmin' : vmin, 'vmax': vmax, 'pmin' : pmin, 'pmax' : pmax, 'stretch' : stretch, 'cmap' : cmap, 'invert' : invert, 'colorbar' : colorbar, 'frametickcolour' : frametickcolour, 'showscale' : showscale, 'distance' : distance, 'scaleunits' : scaleunits, 'aspect' : aspect, 'annotationfontsize' : annotationfontsize}
        plotprep_general(fig, **kwargs)

        fig.set_xaxis_coord_type('scalar')
        
        fig.axis_labels.set_xtext(r'OFFSET (arcsec)')
        fig.axis_labels.set_ytext(r'VELOCITY (km$\,$s$^{-1}$)')

        name = plotprefix+'.'.join(postfixes[posnum].split('.')[:-1])+'_'+greyname_prefix
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
            ax.set_aspect(1.)
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

    #plotmaps(basemap='eso149_or_mask_mom0_sol.fits', vmin = -0.1, vmax = 14., pmin = 0.25, pmax= 99.75, stretch = 'sqrt', cmap='', invert = True, colorbar = False, contoursets = ['eso149_or_mask_mom0_sol.fits', 'eso149_model_or_mask_mom0_sol.fits'], contourcols = ['DarkBlue','DeepPink'], contourlevs = [0.25,0.5,1.0,2.0,4.0,8.0], frametickcolour = '#555555', showbeam = True, showscale = 4000.*u.pc, distance = 7010000*u.pc, plotname = 'ESO149-G003_mom0.pdf')
#    plotmaps(basemap='eso149_or_mask_mom0_sol.fits', vmin = -0.1, vmax = 14., pmin = 0.25, pmax= 99.75, stretch = 'sqrt', cmap='', invert = True, colorbar = False, contoursets = ['eso149_or_mask_mom0_sol.fits', 'eso149_model_or_mask_mom0_sol.fits'], contourcols = ['DarkBlue','DeepPink'], contourlevs = [0.25,0.5,1.0,2.0,4.0,8.0], frametickcolour = '#555555', showbeam = True, showscale = 2*u.arcmin, distance = 7010000*u.pc, plotname = 'ESO149-G003_mom0.pdf')

    # This is to create contours around vsys
    nvelconts = 3
    svelconts = 5.
    vsys = +5.78158E+02
    velcontourlevs = np.arange(vsys-nvelconts*svelconts,vsys+(nvelconts+1)*svelconts,svelconts)
#    print(velcontourlevs)
    #n = 2, nx = 1, ny = 2, startchan = 0, 
    plotmaps(basemap=['eso149_original.fits','eso149_model.fits'], suppress_xlab = False, suppress_ylab = False, vmin = None, vmax = 0.03, pmin = 0.25, pmax= 99.75, stretch = 'sqrt', cmap='', invert = True, colorbar = False, contoursets = ['eso149_original.fits', 'eso149_model.fits'], contourcols = [['grey','DarkBlue','DarkBlue','DarkBlue','DarkBlue'],'DeepPink'], contourlevs = [[-0.002,0.002,0.004,0.008,0.016],[0.004]], contourstyles = [['dashed','solid','solid','solid','solid'],'solid'], frametickcolour = '#555555', showbeam = True, showscale = 2*u.arcmin, showscalelab = True, reducescalelab = True, scaleunits = 'arcmin', distance = 7010000*u.pc, plotname = 'blatest.pdf', width = 10, chans = [0,2,4,6,8,10,12,13], nx=2, annotationfontsize = 'xx-large', plotvelo = 'eso149_original.fits')
    
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
    
    plotpvdiagrampositions('eso149_in.fits', centre_ra, centre_dec, pa, len_maj, len_min, mindist, '0.7', basemap='eso149_or_mask_mom0_sol.fits', plane = 0, vmin = -0.1, vmax = 14., pmin = 0.25, pmax= 99.75, stretch = 'sqrt', cmap='', invert = True, colorbar = False, annotationfontsize = 'xx-large', contoursets = ['eso149_or_mask_mom0_sol.fits', 'eso149_model_or_mask_mom0_sol.fits'], contourcols = ['DarkBlue','DeepPink'], contourlevs = [0.25,0.5,1.0,2.0,4.0,8.0], frametickcolour = '#555555', showbeam = True, showscale = 4000*u.pc, distance = 7010000*u.pc, plotname = 'ESO149-G003_slicepositions.pdf')

    pvdiagrams('eso149_in.fits', 'original', centre_ra, centre_dec, pa, len_maj, len_min, mindist)
    pvdiagrams('final_eso149_out.fits', 'finalmodel', centre_ra, centre_dec, pa, len_maj, len_min, mindist)
    #pvdiagrams('eso149_out_19.fits', 'finalfreerotation', centre_ra, centre_dec, pa, len_maj, len_min, mindist)
    #pvdiagrams('eso149_out_20.fits', 'finalsymmetricdisk', centre_ra, centre_dec, pa, len_maj, len_min, mindist)
    #pvdiagrams('eso149_out_17.fits', 'finalfixeddisp', centre_ra, centre_dec, pa, len_maj, len_min, mindist)

    # Most explains itself, but the ellipse height is set to vres times the channel width in the beaminfoset
    plotpvdiagrams(greyname_prefix = 'original', vmin = None, vmax = 0.03, pmin = 0.25, pmax= 99.75, stretch = 'arcsinh', aspect = 'auto', cmap='', invert = True, colorbar = False, annotationfontsize = 'xx-large', contoursets = ['original','finalmodel'], postfixes = ['_pvmaj.fits','_pvmin_left.fits', '_pvmin_cent.fits', '_pvmin_right.fits'], letters = ['A', 'B', 'C', 'D'], lettercol = 'red', contourcols = ['DarkBlue','DeepPink'], contourlevs = [0.002,0.004,0.008,0.016,0.032], contourstyles=['solid','dashed'], frametickcolour = '#555555', showbeam = True, beaminfoset = 'eso149_in.fits', vres = 2, showscale = 0.75*u.arcmin, scaleunits = 'arcsec', distance = 7010000*u.pc, plotprefix = 'ESO149-G003', plotpostfix = '.pdf')
    #plotpvdiagrams('original','original','finalfreerotation')
    #plotpvdiagrams('original','original','finalsymmetricdisk')
    #plotpvdiagrams('original','original','finalfixeddisp')
