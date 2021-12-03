# demo.py

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
import subprocess
import sys

# listing all FITS files in $PWD
directory = "."
extension = ".fits"
files = [file for file in os.listdir(directory) if
file.lower().endswith(extension)]

# user inputs to script
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
# source name
source = 'IRAS_07454-7112'

# input parameters for function moment_map_plotter
sig_mom = [3,5] # sigma value at which to draw contours in moment maps

# input parameters for function channel_map_plotter
imsize = 25 # arcseconds
tick_size = 10 # arcseconds
sig = 5 # sigma value at which to draw contours in channel maps

# input parameters for function plot_az_av_rad_prof
min_radius = 0 # arcseconds
max_radius = 15 # arcseconds
binsize = 1 # pixels

# choosing functions to run
get_presets    = True
get_spectrum   = True
make_moments   = True
plot_moments   = True
plot_chan_maps = True
plot_rad_prof  = True
combine_plots  = True

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# function calls
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
for cube in files:

    print('Processing %s...\n' %cube)

    if get_presets:
        # function call to define source specific parameters
        v_sys, v_exp, phasecenter, mid_region = get_source_presets(source=source)

    if get_spectrum:
        # function call to extract and plot spectral line from image cubes
        vels = get_spectrum_from_cube(cube=cube, phasecenter=phasecenter, mid_region=mid_region, v_sys=v_sys)

    if make_moments:
        # function call to produce FITS moment maps
        make_moment_maps(cube=cube, vels=vels, v_sys=v_sys, v_exp=v_exp, moment=0)
        make_moment_maps(cube=cube, vels=vels, v_sys=v_sys, v_exp=v_exp, moment=8)

    if plot_moments:
        # function calls to plot moment maps in PDF
        plot_moment_maps(cube=cube, moment=0, sig_mom=sig_mom)
        plot_moment_maps(cube=cube, moment=8, sig_mom=sig_mom)

    if plot_chan_maps:
        # function call to plot channel maps
        center_channel, center = plot_channel_maps(cube=cube, vels=vels, imsize=imsize, tick_size=tick_size, v_sys=v_sys, sig=sig)

    if plot_rad_prof:
        # function call to generate azimuthally averaged radial profiles
        radii, azimuthal_profile = plot_az_av_rad_prof(cube=cube, min_radius=min_radius, max_radius=max_radius, center=center, center_channel=center_channel, binsize=binsize)

if combine_plots:
    # function call to combine all plots to single pdf file
    combine_plots_in_latex(source=source)
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###

# removing extra files
os.system('rm -rf *.last')
