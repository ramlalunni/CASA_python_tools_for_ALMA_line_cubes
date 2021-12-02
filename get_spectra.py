import numpy as np
import matplotlib.pyplot as plt
import os

from modify_plot_axes import align_yaxis_origins

# function to extract and plot line spectra from image cubes
###----------------------------------------------------------------------------###
def get_spectrum_from_cube(cube, phasecenter, mid_region, large_region, v_sys):
    """extracts and plots spectra from three apertures (central pixel, mid_region, and large_region)
       Inputs: 1. cube (str): FITS filename, Eg: 'IRAS_07454-7112_100.076392GHz_HC3N_a3.fits'
               2. phasecenter (str): source cooridnate in the exact format and number of characters as follows: 'J2000 07h45m02.411040s -71d19m45.72840s'
               3. mid_region, large_region (str): diameter of circular regions (centerted on phasecenter) to extract spectra from. Eg: '9', '25'
               4. v_sys (float): systemic velocity of the source. Eg: -38.70"""
    
    # creating necessary directories (if not already present)
    dirs_to_create = ['ascii_spectra/', 'pdf_spectra/']
    for dir in dirs_to_create:
        if not os.path.exists('%s' %dir):
            os.makedirs('%s' %dir)

    # calling CASA function specflux() to calculate flux densities from specified apertures
    specflux(imagename=cube, region='circle[[%s, %s], %sarcsec]' %(phasecenter[6:22],phasecenter[23:],large_region), function='flux density', unit='km/s', logfile='./ascii_spectra/'+cube[:-4]+'spectrum_%s.txt' %large_region, overwrite=True)
    specflux(imagename=cube, region='circle[[%s, %s], %sarcsec]' %(phasecenter[6:22],phasecenter[23:],mid_region), function='flux density', unit='km/s', logfile='./ascii_spectra/'+cube[:-4]+'spectrum_%s.txt' %mid_region, overwrite=True)
    specflux(imagename=cube, region='circle[[%s, %s], 1pix]' %(phasecenter[6:22],phasecenter[23:]), function='flux density', unit='km/s', logfile='./ascii_spectra/'+cube[:-4]+'spectrum_cp.txt', overwrite=True)

    # retrieving spectral data from text files
    spec_large = np.genfromtxt('./ascii_spectra/'+cube[:-4]+'spectrum_%s.txt' %large_region, skip_header=4, usecols=[2,3,4])
    spec_mid = np.genfromtxt('./ascii_spectra/'+cube[:-4]+'spectrum_%s.txt' %mid_region, skip_header=4, usecols=[2,3,4])
    spec_cp = np.genfromtxt('./ascii_spectra/'+cube[:-4]+'spectrum_cp.txt', skip_header=4, usecols=[2,3,4])

    # plotting spectra
    fig = plt.figure()
    ax = fig.add_subplot(111)
    twin_ax = ax.twinx()

    # using step plot for histogram effect
    twin_ax.step(spec_cp[:,1], spec_cp[:,2]*1e3, where='mid', color='g', linewidth=0.75, label='central pixel (mJy)')
    ax.step(spec_mid[:,1], spec_mid[:,2], where='mid', color='b', linewidth=0.75, label='%s" region (Jy)' %mid_region)
    ax.step(spec_large[:,1], spec_large[:,2], where='mid', color='r', linewidth=0.75, label='%s" region (Jy)' %large_region)
    #ax.fill_between(spec_large[:,1], spec_large[:,2], step='mid', color='gray', alpha=0.35)

    # adding horizontal line at zero and vertical line at systemic velocity
    ax.axhline(0.0, color='k', linewidth=0.5)
    ax.axvline(v_sys, color='k', linestyle='--', linewidth=0.5, label='v$_{sys}$ = %s km/s' %v_sys)

    # setting axes labels
    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('Circular Region Flux Density (Jy)')
    twin_ax.set_ylabel('Central Pixel Flux Density (mJy)', rotation=270, labelpad=12.5)

    # function call to align the circular and cp flux density axes
    axes_to_align = np.array([ax, twin_ax])
    align_yaxis_origins(axes=axes_to_align)

    # setting title and legends
    ax.set_title('Spectra: %s' %cube[:-5], pad=20, fontweight='bold')
    ax.legend(loc='upper left')
    twin_ax.legend(loc='upper right')

    fig.align_ylabels()

    # saving figure to pdf file
    plt.savefig('./pdf_spectra/%s_spectra.pdf' %cube[:-5])
    plt.close()

    # extracting channel velocities
    vels  = spec_large[:,1] # velocity of each channel (km/s)
    return vels
###----------------------------------------------------------------------------###
