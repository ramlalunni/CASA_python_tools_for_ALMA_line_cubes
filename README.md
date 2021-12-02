# CASA_python_tools_for_ALMA_cubes
*Common Astronomy Software Applications ([CASA](https://casa.nrao.edu))* python scripts for visualizing interferometric ALMA spectral line cubes in various ways: channel maps, radial profiles, flux density spectra, moment maps, etc. Useful to loop through in case of spectral survey observations with hundreds of lines. In case of observations of a single line, or very few lines, some of these tasks can be easily performed manually via the associated CASA GUI.

**NOTE**: Can be run only from within a CASA terminal, and not from standalone Python or IPython consoles.

**PREREQUISITES**: 
- NRAO CASA (tested only on versions 5.6.1-8 and later, but expected to work on previous versions also).
- Numpy, Matplotlib, Astropy and Pandas installed into the PYTHONENV used by CASA.

**INPUT**: ALMA spectral line image FITS cube, velocity corrected and not having Stokes axes, created using CASA.

Use **dropstokes=True** in [CASA.exportfits()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.data.exportfits.html#exportfits) to remove redundant Stokes axes, if any. Use **restfreq='xx.xxxGHz'** in [CASA.tclean()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.imaging.tclean.html?highlight=tclean#tclean) or [CASA.imreframe()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.imreframe.html?highlight=imreframe#imreframe) to assign correct channel velocities, and use **velocity=True** [CASA.exportfits()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.data.exportfits.html#exportfits) to ensure that the spectral axis of the input cube is in velocity units rather than frequency.

Current functionality includes the following:

1. Extract (ASCII) and plot (PDF) spectra from three apertures (central pixel of image, and two other user specified circles).
2. Plot 5x5 tiled channel maps (PDF) with the systemic velocity channel at the center.
3. Generate moment (0 and 8) maps (FITS) from the line cubes and plot them (PDF).
4. Extract (ASCII) and plot (PDF) azimuthally averaged radial flux density profiles of the line center (systemic velocity channel).

CASA is being developed by scientists based at the National Radio Astronomical Observatory (NRAO), the European Southern Observatory (ESO), and the National Astronomical Observatory of Japan (NAOJ), under the guidance of NRAO, and is licensed under version 2 or later of the [GNU LGPL](http://www.gnu.org/licenses/lgpl.html). Refer to: [McMullin, J. P.; et al. 2007](https://ui.adsabs.harvard.edu/abs/2007ASPC..376..127M/abstract).

CASA documentation is available [here](https://casadocs.readthedocs.io/en/stable/).

## 1. Extract and Plot Spectrum

Uses CASA.specflux() to extract and write out flux densities from two circular apertures of user specified radii and also the central pixel of the image, and plots these using matplotlib.pyplot.step(), with the circular apertures and central pixel intensities in two separate y-axes, and a shared x-axis of velocity. The origins of the twin y-axes are forced to align.

<img width="814" alt="Screenshot 2021-12-02 at 17 18 13" src="https://user-images.githubusercontent.com/87668393/144460708-11fca880-8f66-45a8-b833-79ff9328fabc.png">

## 2. Plot Tiled Channel Maps

Plots 5x5 tiled channel maps (**NOTE**: *cube should have atleast 12 channels on either side of central (systemic velocity) channel in case of 5x5*). Each plotted channel will also have a single contour (of user specified sigma level) overplotted on it. Region (default: 25 arcseconds) to plot, number and arrangement of tiles (default: 25, 5x5), etc. are easily customisable.

<img width="774" alt="Screenshot 2021-12-02 at 17 25 07" src="https://user-images.githubusercontent.com/87668393/144462006-b5e8a2be-25f1-4883-ad8b-fe190df4cf4e.png">

## 3. Generate and Plot Moment Maps

Generates moment 0 (integrated value of the spectrum) and moment 8 (maximum value of the spectrum) maps (FITS) using CASA.immoments() and plots them in python. The RA and DEC axes ticklabels are printed as hms-dms as in the standard CASA.viewer() instance. Astropy.coordinates.SkyCoord() is used for this purpose, and must be installed. Two contour levels (user specified) are overplotted on the moment maps.

Moment 0             |  Moment 8
:-------------------------:|:-------------------------:
 <img width="699" alt="Screenshot 2021-12-02 at 17 44 47" src="https://user-images.githubusercontent.com/87668393/144465541-9fed4a39-b185-44d6-8ac6-2a0e23f9bf18.png"> |  <img width="710" alt="Screenshot 2021-12-02 at 17 45 39" src="https://user-images.githubusercontent.com/87668393/144465663-5b4d3164-4053-46fe-b426-1484eaa84fd1.png">

## 4. Generate and Plot Azimuthally Averaged Radial Profiles

Azimuthally Averaged Radial Profiles provide a circularly averaged radial intensity distribution of the image. There is currently no direct CASA task that implements the same. These can be used to determine the extents of the emission morphology in the case of cricumstellar envelopes or other spherical, disk or ring like emitting regions.

This is a simplified and slightly modified version of the [*radialprofile.py*](https://github.com/keflavich/image_tools/blob/793e93065afe2754a818da8b58f9b222a3acf59f/image_tools/radialprofile.py) script by [Adam Ginsburg](https://github.com/keflavich). The original [*radialprofile.py*](https://github.com/keflavich/image_tools/blob/793e93065afe2754a818da8b58f9b222a3acf59f/image_tools/radialprofile.py) has more advanced functionality (weighting, masking, etc.). It must be attributed in all reproductions.

<img width="737" alt="Screenshot 2021-12-02 at 18 06 41" src="https://user-images.githubusercontent.com/87668393/144469214-aac4aae0-c83b-4c4f-813c-ad83628c47fa.png">

## Licensing

All scripts provided here are licensed under the **GNU General Public License v3.0**, except where code is explicitely mentioned to be adapted/taken from other sources, in which case the respective source licenses continue to apply. The GNU General Public License is free and copyleft. THE PROGRAMS ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
