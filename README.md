# CASA_python_tools_for_ALMA_cubes
*Common Astronomy Software Applications ([CASA](https://casa.nrao.edu))* python scripts for visualizing interferometric ALMA spectral line cubes in various ways: channel maps, radial profiles, flux density spectra, moment maps, etc.

**NOTE**: Can be run only from within a CASA terminal, and not from standalone Python or IPython consoles.

**INPUT**: ALMA spectral line image FITS cube, velocity corrected and not having Stokes axes, created using CASA.

Use **dropstokes=True** in [CASA.exportfits()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.data.exportfits.html#exportfits) to remove redundant Stokes axes, if any. Use **restfreq='xx.xxxGHz'** in [CASA.tclean()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.imaging.tclean.html?highlight=tclean#tclean) or [CASA.imreframe()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.imreframe.html?highlight=imreframe#imreframe) to assign correct channel velocities, and use **velocity=True** [CASA.exportfits()](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.data.exportfits.html#exportfits) to ensure that the spectral axis of the input cube is in velocity units rather than frequency.

Current functionality includes the following:

1. Extract (ASCII) and plot (PDF) spectra from three apertures (central pixel of image, and two other user specified circles).
2. Plot 5x5 tiled channel maps (PDF) with the systemic velocity channel at the center.
3. Generate moment (0 and 8) maps (FITS) from the line cubes and plot them (PDF).
4. Extract (ASCII) and plot (PDF) azimuthally averaged radial flux density profiles of the line center (systemic velocity channel).

CASA is being developed by scientists based at the National Radio Astronomical Observatory (NRAO), the European Southern Observatory (ESO), and the National Astronomical Observatory of Japan (NAOJ), under the guidance of NRAO, and is licensed under version 2 or later of the [GNU LGPL](http://www.gnu.org/licenses/lgpl.html). Refer to: [McMullin, J. P.; et al. 2007](https://ui.adsabs.harvard.edu/abs/2007ASPC..376..127M/abstract).

CASA documentation is available [here](https://casadocs.readthedocs.io/en/stable/).
