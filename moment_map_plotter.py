# moment_map_plotter.py

########################################################################
# Copyright (C) 2021  Ramlal Unnikrishnan

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

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
from matplotlib.ticker import AutoMinorLocator
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import os

from modifier_functions import ChangeAxisColor, truncate

# function to generate moment (0 or 8) maps (FITS)
###----------------------------------------------------------------------------###
def make_moment_maps(cube, phasecenter, v_sys, v_exp, moment):
    """makes moment 0 or 8 maps of input line cube in FITS format
       Inputs: 1. cube (str): FITS filename of input cube, Eg: 'IRAS_07454-7112_100.076392GHz_HC3N_a3.fits'
               2. phasecenter (str): source cooridnate in the exact format and number of characters as follows: 'J2000 07h45m02.411040s -71d19m45.72840s'
               3. v_sys (float): systemic velocity (km/s) of the source, Eg: -38.70
               4. v_exp (float): expansion velocity (km/s) of the source, Eg: -13.0
               5. moment (int): moment(0 or 8) to generate, Eg: 0"""
  
    # creating necessary directories (if not already present)
    dirs_to_create = ['fits_moment_%s_maps' %moment]
    for dir in dirs_to_create:
        if not os.path.exists('%s' %dir):
            os.makedirs('%s' %dir)
    
    # reading header from fits cube
    header = fits.getheader('%s' %cube)
    
    # velocity references
    NAXIS3 = header['NAXIS3']
    CTYPE3 = header['CTYPE3']
    CRVAL3 = header['CRVAL3']
    CDELT3 = header['CDELT3']
    CRPIX3 = header['CRPIX3']
    
    # finding channel velocities
    if CTYPE3 == 'VRAD':
        vels = np.zeros(NAXIS3)
        for i in range(NAXIS3):
            vels[i] = CRVAL3 + CDELT3*(i - CRPIX3 + 1)
    else:
        raise TypeError('AXIS3 is not in velocity units')

    # calculating channel numbers corresponding to v_sys +/- 1.2*v_exp (20% extra on both sides)
    plus_v_exp_chan  = np.argmin(np.abs(vels-(v_sys+1.2*v_exp)))
    minus_v_exp_chan = np.argmin(np.abs(vels-(v_sys-1.2*v_exp)))

    # calculating upper and lower limiting channel numbers for moment calculations
    upper_cut_chan = max(plus_v_exp_chan, minus_v_exp_chan)
    lower_cut_chan = min(plus_v_exp_chan, minus_v_exp_chan)

    # generating moment map and converting to FITS
    os.system("rm -rf ./fits_moment_%s_maps/'%s'.mom%s" %(moment,cube[:-5],moment)) # deleting pre-existing .momx file of same line, if any
    immoments(imagename='%s' %cube, moments=moment, chans='%d~%d' %(lower_cut_chan, upper_cut_chan), outfile='./fits_moment_%s_maps/%s.mom%s' %(moment,cube[:-5],moment))
    exportfits(imagename='./fits_moment_%s_maps/%s.mom%s' %(moment,cube[:-5],moment), fitsimage='./fits_moment_%s_maps/%s.mom%s.fits' %(moment,cube[:-5],moment), overwrite=True) # saving moment map as FITS file
    os.system("rm -rf ./fits_moment_%s_maps/'%s'.mom%s" %(moment,cube[:-5],moment)) # deleting CASA moment image file
###----------------------------------------------------------------------------###


# function to plot moment maps (PDF).
###----------------------------------------------------------------------------###
def plot_moment_maps(cube, moment, sig_mom):
    """plots moment maps in PDF format
       Inputs: 1. cube (str): FITS filename of input cube, Eg: 'IRAS_07454-7112_100.076392GHz_HC3N_a3.fits'
               (NOTE: This is not the filename of the FITS moment map. but rather that of the original FITS image cube. However, function call will fail if 
               FITS moment map is not present in the ./fits_moment_x_maps/ directory where 'x' is the moment number under consideration.)
               2. moment (int): moment(0 or 8) to generate, Eg: 0
               3. sig_mom (list of ints): sigma levels at which to plot contours, Eg: [3,5]"""
    
    # creating necessary directories (if not already present)
    dirs_to_create = ['pdf_moment_%s_maps' %moment]
    for dir in dirs_to_create:
        if not os.path.exists('%s' %dir):
            os.makedirs('%s' %dir)
    
    # reading data from moment fits file
    mom_filename = './fits_moment_%s_maps/%s.mom%s.fits' %(moment, cube[:-5], moment)
    mom_fits_file = fits.open(mom_filename)
    mom_header = fits.getheader(mom_filename)
    mom_data = mom_fits_file[0].data

    # reading header from moment fits file
    CRPIX1_mom = mom_header['CRPIX1'] - 1 # RA reference pixel
    CRVAL1_mom = mom_header['CRVAL1']     # RA value (degree) at RA reference pixel
    CDELT1_mom = mom_header['CDELT1']     # pixel size (degrees) along RA axis
    CRPIX2_mom = mom_header['CRPIX2'] - 1 # DEC reference pixel
    CRVAL2_mom = mom_header['CRVAL2']     # DEC value (degrees) at DEC reference pixel
    CDELT2_mom = mom_header['CDELT2']     # pixel size (degrees) along DEC axis
    BUNIT_mom  = mom_header['BUNIT']      # unit of moment map

    # readign beams from moment fits file
    BMAJ_mom = mom_header['BMAJ']/(abs(np.sqrt(CDELT1_mom**2 + CDELT2_mom**2)))
    BMIN_mom = mom_header['BMIN']/(abs(np.sqrt(CDELT1_mom**2 + CDELT2_mom**2)))
    BPA_mom = mom_header['BPA']

    # plotting moment image
    print('Plotting moment %s map of %s\n' %(moment, cube))
    fig_mom, ax_mom = plt.subplots()
    ax_mom = ChangeAxisColor(ax_mom,'white')
    cmap = 'cubehelix'
    im_mom = ax_mom.imshow(mom_data[0,:,:], cmap=cmap, origin="lower")

    # calculating center pixel of moment image
    center_mom = [int(CRPIX1_mom), int(CRPIX2_mom)]

    # calculating pixel size in arcesconds
    pixelsize_mom = 3600 * abs(CDELT1_mom)

    # calculating rms of moment image
    boundary_mom = 25 # arcsec
    y_mom, x_mom = np.indices(mom_data[0].shape)
    r_mom = np.hypot(x_mom - center_mom[0], y_mom - center_mom[1])
    boundary_pixel_mom = boundary_mom/pixelsize_mom - center_mom[0]
    chan_copy_mom = mom_data.copy()
    for i in range(len(r_mom)):
        for j in range(len(r_mom)):
            if(r_mom[i,j] <= boundary_pixel_mom):
                chan_copy_mom[0][i,j] = 0
    rms_mom = np.sqrt((chan_copy_mom**2).mean())

    # plotting contours at input specified sigma level
    ax_mom.contour(mom_data[0,:,:], [sig_mom[0]*rms_mom], colors='blue', linewidths=0.5, origin="lower")
    ax_mom.contour(mom_data[0,:,:], [sig_mom[1]*rms_mom], colors='white', linewidths=0.5, origin="lower")

    # adding beam patch and bounding rectangle
    ax_mom.add_patch(matplotlib.patches.Ellipse(xy=[15,15], width=BMAJ_mom, height=BMIN_mom, angle=(BPA_mom-90), color='white', fill=True))
    ax_mom.add_patch(matplotlib.patches.Rectangle(xy=[15-10,15-10], width=20, height=20, color='red', fill=False, linewidth=0.5))

    # setting minor ticks
    ax_mom.xaxis.set_minor_locator(AutoMinorLocator(n=5))
    ax_mom.yaxis.set_minor_locator(AutoMinorLocator(n=5))

    # setting major tick locations
    ten_arcsec_interval = 8/(3600*np.abs(CDELT1_mom))
    x_tick_locs = np.zeros(4)
    y_tick_locs = np.zeros(4)
    for i in range(len(x_tick_locs)):
        x_tick_locs[i] = ten_arcsec_interval*i
    for i in range(len(y_tick_locs)):
        y_tick_locs[i] = ten_arcsec_interval*i
    ax_mom.set_xticks(x_tick_locs)
    ax_mom.set_yticks(y_tick_locs)

    # finding tick labels in degrees
    x_ticks_pix = ax_mom.get_xticks()
    x_ticks_deg = np.zeros(len(x_ticks_pix))
    for i in range(len(x_ticks_pix)):
        x_ticks_deg[i] = CRVAL1_mom + (x_ticks_pix[i]-CRPIX1_mom)*CDELT1_mom/np.cos(CRVAL2_mom*np.pi/180)
    y_ticks_pix = ax_mom.get_yticks()
    y_ticks_deg = np.zeros(len(y_ticks_pix))
    for j in range(len(y_ticks_pix)):
        y_ticks_deg[j] = CRVAL2_mom + (y_ticks_pix[j]-CRPIX2_mom)*CDELT2_mom

    # converting ticklabels to hms-dms
    ras  = []
    decs = []
    for i in range(len(x_ticks_deg)):
        c = SkyCoord(ra=x_ticks_deg[i]*u.degree, dec=y_ticks_deg[i]*u.degree, frame='fk5')
        ras.append(str(truncate(c.ra.hms.s, 1))+'$^\mathrm{s}$')
        decs.append(str(truncate(np.abs(c.dec.dms.s), 1))+'"')
    ref_ra_hms = SkyCoord(ra=x_ticks_deg[1]*u.degree, dec=y_ticks_deg[1]*u.degree, frame='fk5')
    ref_ra_hms = str(int(ref_ra_hms.ra.hms.h))+'$^\mathrm{h}$'+str(int(ref_ra_hms.ra.hms.m))+'$^\mathrm{m}$'+str(truncate(ref_ra_hms.ra.hms.s,1))+'$^\mathrm{s}$'
    ref_dec_dms = SkyCoord(ra=x_ticks_deg[-1]*u.degree, dec=y_ticks_deg[-1]*u.degree, frame='fk5')
    ref_dec_dms = str(int(ref_dec_dms.dec.dms.d))+'$^\mathrm{d}$'+str(np.abs(int(ref_dec_dms.dec.dms.m)))+'$^\mathrm{m}$'+str(np.abs(truncate(ref_dec_dms.dec.dms.s,1)))+'"'
    ras[1] = ref_ra_hms
    decs[-1] = ref_dec_dms
    ras[0] = ''
    #decs[0] = ''

    # setting ticklabels
    ax_mom.set_xticklabels(ras)
    ax_mom.set_yticklabels(decs)

    # setting axis labels
    ax_mom.set_xlabel('Right Ascension (J2000)', labelpad=7.5)
    ax_mom.set_ylabel('Declination (J2000)', labelpad=-1.5)

    # adding plot title text
    fig_mom.text(0.525, 0.92, 'Moment %s: %s' %(moment, cube[:-5]),
                            horizontalalignment='center',
                            verticalalignment='center',
                            color='k',
                            fontweight='bold')

    # adding colourbar to plot
    cbar_ax_mom = fig_mom.add_axes([0.815, 0.125, 0.025, 0.755])
    fig_mom.colorbar(im_mom, cax=cbar_ax_mom,orientation='vertical')
    cbar_ax_mom.set_ylabel('%s' %BUNIT_mom, rotation=270, labelpad=15)
    #cbar_ax_mom.set_ylabel('Jy beam$^{-1}$ km s$^{-1}$', rotation=270, labelpad=15)

    # adding legends for contours
    lines_contour = [Line2D([0],[0], color='blue', lw=0.5), Line2D([0],[0], color='white', lw=0.5)]
    labels_contour = ['%s$\sigma$' %sig_mom[0], '%s$\sigma$' %sig_mom[1]]
    ax_mom.legend(lines_contour, labels_contour, loc='upper right')

    # saving image to pdf file
    plt.savefig('pdf_moment_%s_maps/%s_moment_%s.pdf' %(moment, cube[:-5], moment), format='pdf', bbox_inches='tight')
    plt.close()
###----------------------------------------------------------------------------###
