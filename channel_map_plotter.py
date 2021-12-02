# channel_map_plotter.py

########################################################################
# Copyright (C) 2021  Ramlal Unnikrishnan

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version. Some parts of the program are
# also licensed under the Creative Commons Attribution-ShareAlike 4.0 
# license. See <https://creativecommons.org/licenses/by-sa/4.0/>.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# contact: ramlal.unnikrishnan@chalmers.se
########################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy.io import fits
import os

from modifier_functions import ChangeAxisColor

###function to generate tiled channel maps (PDF)
###----------------------------------------------------------------------------###
def plot_channel_maps(cube, imsize, tick_size, v_sys, sig):
    """docstring"""

    # creating necessary directories (if not already present)
    dirs_to_create = ['channel_maps/']
    for dir in dirs_to_create:
        if not os.path.exists('%s' %dir):
            os.makedirs('%s' %dir)

    ###finding central channel of line cube (channel at systemic velocity)
    center_channel = np.argmin(np.abs(vels-v_sys))

    ###reading header, data and beams from fits cube
    hdulist = fits.open('%s' %cube)
    header = fits.getheader('%s' %cube)
    data = hdulist[0].data
    beams   = hdulist[1].data

    ###pixel references
    CRPIX1 = header['CRPIX1'] ###RA reference pixel
    CDELT1 = header['CDELT1']  ###RA pixel width (degrees)
    CRPIX2 = header['CRPIX2'] ###DEC reference pixel
    CDELT2 = header['CDELT2']  ###DEC pixel width (degrees)

    ###finding center pixel
    center = [int(CRPIX1 - 1), int(CRPIX2 - 1)]

    ###calculating pixel size in arcseconds
    pixel_size = 3600 * abs(CDELT1)

    ###intensity unit reference
    BUNIT = header['BUNIT']

    ##defining and adjusting figure and subplots
    nrows = 5
    ncols = 5
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))#, sharex=True, sharey=True)
    fig.text(0.215, 0.0775, 'RA offset (")', ha='center')
    fig.text(0.09, 0.19, 'DEC offset (")', va='center', rotation='vertical')
    fig.subplots_adjust(wspace=-0.1)
    fig.subplots_adjust(hspace=0.0)

    ###choosing tile size and location of ticks
    tick_size_pix = tick_size/(3600 * abs(CDELT1)) ###pixels
    tick_locs = [-tick_size_pix, 0, tick_size_pix]
    x_tick_labels = np.round([3600 * i * CDELT1 for i in tick_locs],decimals=1)
    y_tick_labels = np.round([3600 * i * CDELT2 for i in tick_locs],decimals=1)
    x_tick_labels[1] = np.abs(x_tick_labels[1])
    y_tick_labels[1] = np.abs(y_tick_labels[1])

    ###defining plotting extent and upper and lower limits
    extent = int(0.5*imsize/(3600*abs(CDELT1)))
    low_lim = int(CRPIX1)-extent
    upp_lim = int(CRPIX1)+extent

    ###defining colourmap for imshow
    cmap = plt.get_cmap('cubehelix')

    ###printing calculated cube details to terminal
    print('Channel velocities are %s' %vels)
    print('Center pixel of channels is %s' %center)
    print('Pixel size is %s arcseconds\n' %'{0:.3f}'.format(pixel_size))

    ###choosing number of channels to skip at the beginning of cube
    skip_chans = center_channel - int((nrows*ncols - 1)/2)
    last_chan = center_channel + int((nrows*ncols - 1)/2)

    ###finding maximum intensity to plot for setting colourbar contrast
    max_int = 0.0
    for chan in range(skip_chans, last_chan+1):
        if np.nanmax(data[chan]) > max_int:
            max_int = np.nanmax(data[chan])

    ###finding minimum intensity to plot for setting colourbar contrast
    min_int = 0.0
    boundary = 25 ###arcsec
    y, x = np.indices(data[center_channel].shape)
    r = np.hypot(x - center[0], y - center[1])
    boundary_pixel = boundary/pixel_size - center[0]
    begin_chan = data[0].copy()
    end_chan = data[-1].copy()
    for i in range(len(r)):
        for j in range(len(r)):
            if(r[i,j] >= boundary_pixel):
                begin_chan[i,j] = 0
                end_chan[i,j] = 0
    rms_begin_chan = np.sqrt((begin_chan**2).mean())
    rms_end_chan = np.sqrt((end_chan**2).mean())
    av_rms = (rms_begin_chan+rms_end_chan)/2
    min_int = -av_rms

    print("Unit of pixel intensity is '%s'" %BUNIT)
    print('Average RMS noise of channels is %s %s\n' %('{0:.3e}'.format(av_rms),BUNIT))
    print('Minimum and maximum colourbar levels are %s and %s %s respectively' %('{0:.3e}'.format(min_int),max_int,BUNIT))
    print("Plotted channel size is %dx%d arcseconds (~ %sx%s pixels)" %(imsize,imsize, extent*2, extent*2))
    print('Channels %d to %d will be plotted\n' %(skip_chans,last_chan))

    ###plotting channel maps
    print('Plotting tiled channel maps of %s...' %cube)
    print('Tile \t Channel \t Velocity(km/s)')
    for i in range(0,5):
        for j in range(0,5):
            channel = i*5+j
            vel = vels[skip_chans+channel]
            print('(%d,%d) \t\t %d \t\t %s' %(i, j, skip_chans+channel, '{0:7.1f}'.format(vel)))

            ###writing velocity values to tiles
            axes[4-i][4-j].text(0.15, 0.875, '{0:7.1f}'.format(vel),
                            horizontalalignment='center',
                            verticalalignment='center',
                            fontsize=9, color='white',
                            #fontweight='bold',
                            transform=axes[4-i][4-j].transAxes)

            ###plotting channel maps
            imgplot = axes[4-i][4-j].imshow(data[skip_chans+channel][low_lim:upp_lim,low_lim:upp_lim],cmap=cmap,
                                            extent=[-extent,extent,-extent,extent],origin="lower",vmax=max_int,vmin=min_int)

            ###ploting contour at input specified sigma level
            axes[4-i][4-j].contour(data[skip_chans+channel][low_lim:upp_lim,low_lim:upp_lim], [sig*av_rms], colors='white', linewidths=0.5, extent=[-extent,extent,-extent,extent], origin="lower")

            ###setting ticks and empty ticklabels
            axes[4-i][4-j].set_xticks(tick_locs)
            axes[4-i][4-j].set_yticks(tick_locs)
            axes[4-i][4-j].set_xticklabels([])
            axes[4-i][4-j].set_yticklabels([])

            ###function call to change axes colour to white
            axes[4-i][4-j] = ChangeAxisColor(axes[4-i][4-j],'white')

    ###setting ticklables for bottom-left corner tile
    axes[4][0].set_xticklabels(x_tick_labels)
    axes[4][0].set_yticklabels(y_tick_labels)

    ###adding colourbar to plot
    cbar_ax = fig.add_axes([0.9, 0.125, 0.025, 0.755])
    fig.colorbar(imgplot, cax=cbar_ax,orientation='vertical',ticklocation='right')
    cbar_ax.set_ylabel('Intensity (Jy/beam)', rotation=270, labelpad=15)

    ###plotting beam patch
    av_bmaj = np.median(beams['BMAJ'])/(3600*abs(CDELT1))
    av_bmin = np.median(beams['BMIN'])/(3600*abs(CDELT2))
    av_bpa  = np.median(beams['BPA'])

    ###adding beam patch and bounding rectangle
    axes[4][0].add_patch(matplotlib.patches.Ellipse(xy=[-40,-40], width=av_bmaj, height=av_bmin, angle=(av_bpa-90), color='white', fill=True))
    axes[4][0].add_patch(matplotlib.patches.Rectangle(xy=[-40-10,-40-10], width=20, height=20, color='red', fill=False, linewidth=0.5))
    print('\nMedian synthesised beam size is %sx%s arcseconds\n' %('{0:.3f}'.format(np.median(beams['BMAJ'])), '{0:.3f}'.format(np.median(beams['BMIN']))))

    ###adding legends for contours
    line_contour = [Line2D([0],[0], color='white', lw=0.5)]
    label_contour = ['%s$\sigma$' %sig]
    axes[0][4].legend(line_contour, label_contour, loc='best', bbox_to_anchor=(0.95, 0.95))

    ###adding plot title text
    fig.text(0.515, 0.92, 'Channel maps with %s$\sigma$ contours: %s' %(sig, cube[:-5]),
                            horizontalalignment='center',
                            verticalalignment='center',
                            fontsize=10, color='k',
                            fontweight='bold')

    ###saving image to pdf file
    plt.savefig('./channel_maps/%s_channel_map.pdf' %cube[:-5], format='pdf', bbox_inches='tight')
    plt.close()

    return center_channel, center
###----------------------------------------------------------------------------###
