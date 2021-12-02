# radial_profile_plotter.py

########################################################################
# Simplified adpatation of 'radialprofile.py' by Adam Ginsburg.
# See <https://github.com/keflavich/image_tools/blob/793e93065afe2754a818da8b58f9b222a3acf59f/image_tools/radialprofile.py>"""
# License(s) of the original (see above) script apply.

# contact: ramlal.unnikrishnan@chalmers.se
########################################################################

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# function to generate azimuthally averaged radial profiles (ASCII and PDF)
###----------------------------------------------------------------------------###
def plot_az_av_rad_prof(cube, min_radius, max_radius, center, center_channel, binsize):
    """function to plot azimuthally averaged radial profile of a channel from the line cube.
       simplified adaptation of 'radialprofile.py' by Adam Ginsburg.
       see <https://github.com/keflavich/image_tools/blob/793e93065afe2754a818da8b58f9b222a3acf59f/image_tools/radialprofile.py>"""

    # creating necessary directories (if not already present)
    dirs_to_create = ['ascii_az_av_radial_profiles/', 'pdf_az_av_radial_profiles/']
    for dir in dirs_to_create:
        if not os.path.exists('%s' %dir):
            os.makedirs('%s' %dir)
    
    # reading header and data from fits cube
    hdulist = fits.open('%s' %cube)
    data = hdulist[0].data
    header = fits.getheader('%s' %cube)

    # defining pixel coordinates of image center, if not already input to function
    CRPIX1 = header['CRPIX1']
    CRPIX2 = header['CRPIX2']
    if center is None:
        center = [int(CRPIX1) - 1, int(CRPIX2) - 1]

    # calculating pixel size in arcesconds
    CDELT1 = header['CDELT1']
    pixel_size = 3600 * abs(CDELT1)

    # finding central channel image
    center_image = data[center_channel]

    # calculating image indices
    y, x = np.indices(center_image.shape)
    r = np.hypot(x - center[0], y - center[1])

    # defining binning parameters
    nbins = int(np.round(r.max() / binsize)+1)
    maxbin = nbins * binsize
    bins = np.linspace(0,maxbin,nbins+1)
    bin_centers = (bins[1:]+bins[:-1])/2.0

    # calculating azimuthal profile and radii
    azimuthal_profile = np.histogram(r, bins, weights=center_image)[0]/np.histogram(r, bins)[0]
    radii = bin_centers*pixel_size

    # removing trailing NaNs
    azimuthal_profile = azimuthal_profile[np.logical_not(np.isnan(azimuthal_profile))]
    radii = radii[:len(azimuthal_profile)]

    # saving azimuthal profile and radius values to text file
    azimuthal_profile_data = np.column_stack([radii, azimuthal_profile])
    profile_saving_path = 'ascii_az_av_radial_profiles/%s_az_av_rad_profile.txt' %cube[:-5]
    np.savetxt(profile_saving_path , azimuthal_profile_data, fmt=['%6.3f','%10.8e'])

    # selecting only pixels within the given input range of minimum and maximum radii
    lower_boundary = 0
    flipped_radii = np.flip(radii)
    for i in range(len(flipped_radii)):
        if flipped_radii[i] < min_radius:
            lower_boundary = len(flipped_radii) - i
            break
    upper_boundary = len(radii)
    for j in range(len(radii)):
        if radii[j] > max_radius:
            upper_boundary = j
            break
    selected_radii = radii[lower_boundary:upper_boundary]
    selected_azimuthal_profile = azimuthal_profile[lower_boundary:upper_boundary]
    print('Radial profile from %s to %s arcseconds will be plotted\n' %(min_radius, max_radius))

    # plotting azimuthally averaged radial profile
    print('Plotting azimuthally averaged radial profile of %s\n' %cube)
    fig1, ax1 = plt.subplots()
    ax1.plot(selected_radii, selected_azimuthal_profile)
    ax1.set_xlabel('Radius (")')
    ax1.set_ylabel('Azimuthally averaged intensity (Jy/beam)')
    fig1.text(0.515, 0.95, 'Radial Profile: %s' %cube[:-5], horizontalalignment='center', verticalalignment='center', fontsize=10, color='k', fontweight='bold')
    plt.savefig('pdf_az_av_radial_profiles/%s_az_av_rad_profile.pdf' %cube[:-5], format='pdf', bbox_inches='tight')
    plt.close()

    return radii, azimuthal_profile
###----------------------------------------------------------------------------###
