import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy import units as u
import pyspeckit
import pylab
import matplotlib.pyplot as plt


#def read_fits(Filename):
'''
This function reads the given fits file and converts the datafrom angstroms to microns.

Input: directory/to/fitsfile.fits
Output: Simple plot of flux vs wavelength. Interactive mode can be entered by pressing 'b' for baseline,
and following command prompt. Press 'f' for line fitting, this must be done AFTER baseline, follow command prompt.
File and information must be saved manually at this time.

Flux: assigned unit of erg/s/cm^2/Angstrom
Wavelength: Input should be in angstroms, converts angstroms to microns, plots in microns
'''
plt_lin = [['[ArIII] 7136',7136.97], #wavelengths are in angstroms
               ['OI',11287],#*
               ['[OII] 7319',7319.0],
               ['[SII] 6718',6718.95],
               ['[SIII]',9069],#*
               ['SiI',15880],
               ['CaI',22630],
               ['MgI 15750',15750],
               ['MgI 17110',17110],
               ['HeI 10830',10830],#*
               ['HeI 20520',20520],#*
               ['HeII 8239',8239.26],
               ['HeII 20580',20580],
               ['CaII 8544',8544.44],
               ['CaII 8498',8498],
               ['CaII 8662',8662],
               ['[FeII] 9202',9202],#*
               ['[FeII]12600',12600],
               ['[FeII]16440',16440],
               ['H$_{2}$ 19570',19570],
               ['H$_{2}$ 21210',21210],
               [r'Pa$_{\alpha}$',18750.1],
               [r'Pa$_{\beta}$',12818.1],
               [r'Pa$_{\gamma}$',10938],
               [r'Pa$\delta$',10049.8],
               [r'Pa$\epsilon$',9545],#9546.2
               ['NaI',22080]]

hdulist = fits.open('Filename.fits') # Opens fits file
data= hdulist[0].data # This is the flux

# Following lines look at header and extract the wavelength information
header = hdulist[0].header
gal_name = header['OBJECT']
wcs = WCS(header)
index = np.arange(header['NAXIS1']) # Make index array
wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
wavelength = wavelength.flatten() # Makes sure the wavelength has correct dimensions

hdulist.close()

# Giving units to flux
flux_unit = u.erg / (u.cm**2 * u.s * u.AA)
flux = data * flux_unit * 1e16 # Need larger numbers for pyspeckit.Spectrum to be able to read
flux_values = flux.value # Just the numbers

#Giving units to wavelength and converting to microns
wv_unit = u.AA
wavelength = wavelength * wv_unit
wavelength = wavelength.to(u.micron) # Converting Angstrom to microns
wv_values = wavelength.value # Just the numbers, in microns

# print('Object name:', gal_name)
# print('Flux:', flux_prime * flux_unit)
# print('Wavelength:', wavelength)

# Create whole spectreum
spec = pyspeckit.Spectrum(data=flux_values, xarr=wv_values, header=header, unit='erg/s/cm^2/AA')


x_min = wv_values[0] - .01
x_max = wv_values[-1] + .01
#Look at last fifth of spectrum, excluding last 50 elements, and gets min flux value for whole spectrum plot
fifth = np.floor(len(wv_values)*.2)
y_min = np.min(flux_values[-fifth:-50]) - .1
#Look at first fourth of spectrum and gets max flux value for whole spectrum plot
fourth = np.floor(len(wv_values)*.25) 
y_max = np.max(flux_values[0:fourth]) + .2

# print(wv_values[0:fourth])
# print('x-range:', x_min, x_max)
# print('y-range:', y_min, y_max)

# Plot whole spectrum 
spec.plotter(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max) 
# Enter baseline fitter by pressing 'b', '1' around basline area, '3' to plot baseline
# Enter line fit by pressing 'f', '1'at each end of line, '2' at peak then at FWHM, '3' to fit


for line in range(0,len(plt_lin)):
  lin_name = plt_lin[line][0]
  lin_value = plt_lin[line][1]*u.AA
  lin_value = lin_value.to(u.micron).value #Converting plt_lin units to microns
  base_min = lin_value-.04
  base_max = lin_value+.04
  lin_min = lin_value-.01
  lin_max = lin_value+.01
  print(lin_name, lin_value)
  print('Baseline Range:', base_min, base_max)
  print('Line Widith:', lin_min, lin_max)
  # spec.plotter(xmin=base_min, xmax=base_max, ymin=y_min, ymax=y_max) # Plot only baseline range, not working
  if (lin_value > x_min) & (lin_value < x_max): # Plotting my lines on the graph
    plt.axvline(lin_value,color='b',linestyle='--')
    lbl_pos = (y_max-y_min)*0.85 # at 85% up plot
    plt.text(lin_value,lbl_pos,lin_name,rotation=45) 
    plt.axvline(lin_value,color='b',linestyle='--')
    plt.text(lin_value,lbl_pos,lin_name,rotation=45)

pylab.show()

# *********************************
# Nothing below here works but should be the bulk of this script
# *********************************

# Fit a first order polynomial for the continuum non-interactively
# spec.baseline(clear_all_connections=True, order=1, highlight_fitregion=True, reset_selection=True, exclude=(lin_min,lin_max), plot_baseline=True, linewidth=2,baseline_fit_color='r',fit_plotted_area=True, fit_original=True)

#Fit gaussians to the spectral lines non-interactively, needs to be done AFTER contin fit.
# amplitude_guess = lin_value.value
# center_guess = lin_value.value
# width_guess = lin_max - lin_min
# guesses = [amplitude_guess, center_guess, width_guess]
# spec.specfit(guesses=guesses, fittype='gaussian', plot=True, continuum_as_baseline=True, components=True, x_min=lin_min, x_max=lin_max, xunits='microns', midpt_location='plt-center')










