#!/usr/bin/env python

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy import units as u
import pyspeckit
import pylab
from astropy.io import fits
from astropy.modeling import models, fitting


def read_fits(filename):
    """Read FITS file and return data, header, galaxy name and wavelengths"""
    hdulist = fits.open(filename)  # Opens fits file
    data = hdulist[0].data  # This is the flux
    # Following lines look at header and extract the wavelength information
    header = hdulist[0].header
    gal_name = header['OBJECT']
    wcs = WCS(header)
    index = np.arange(header['NAXIS1'])  # Make index array
    wavelengths = wcs.wcs_pix2world(index[:, np.newaxis], 0)
    # Makes sure the wavelength has correct dimensions:
    wavelengths = wavelengths.flatten()
    hdulist.close()
    return data, header, gal_name, wavelengths


def get_flux_values(data):
    """Extract flux values from FITS data; return unitless values"""
    flux_unit = u.erg / (u.cm**2 * u.s * u.AA)
    flux = data * flux_unit * 1e16  # Need larger numbers for pyspeckit
    flux_values = flux.value  # Just the numbers
    return flux_values


def get_wave_values(wavelengths):
    """Convert wavelength units to microns, return values only"""
    wv_unit = u.AA
    wavelengths = wavelengths * wv_unit
    wavelengths = wavelengths.to(u.micron)  # Converting Angstrom to microns
    wave_values = wavelengths.value  # Just the numbers, in microns
    return wave_values


def main():
    '''
    This function parses command line arguments, reads the given fits file
    and converts the data from angstroms to microns.

    Input: directory/to/fitsfile.fits
    Output: Simple plot of flux vs wavelength.
    Interactive mode can be entered by pressing 'b' for baseline 
    and following the command prompt.
    ('1' around basline area, '3' to plot baseline)
    Press 'f' for line fitting, this must be done AFTER baseline 
    and follow command prompt.
    ('1'at each end of line, '2' at peak then at FWHM, '3' to fit)
    Information must be saved manually at this time,
    although images are saved automatically.

    Flux: assigned unit of 10E-16 erg/s/cm^2/Angstrom
    Wavelength: Input should be in angstroms, converts angstroms to microns,
    plots in microns
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fits_file')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--directory', '--dir', default='../Images/',
                        help="Will be prepended to each .eps file written")
    parser.add_argument('--plot_type', '-t', default='lines',
                        help="Plots entire galaxy spectrum")

    args = parser.parse_args()
    fits_file = args.fits_file

    if not os.path.exists(fits_file):
        sys.stderr.write("Couldn't find %s. Exiting...\n".format(fits_file))
        sys.exit()

    if args.debug:
        # Wavelengths are in angstroms
        reference_wavelengths = [['SiI', 15880]]
    else:
        reference_wavelengths = [['[ArIII] 7136', 7136.97],
                                 ['OI', 11287],  #
                                 ['[OII] 7319', 7319.0],
                                 ['[SII] 6718', 6718.95],
                                 ['[SIII]', 9069],  #
                                 ['SiI', 15880],
                                 ['CaI', 22630],
                                 ['MgI 15750', 15750],
                                 ['MgI 17110', 17110],
                                 ['HeI 10830', 10830],  #
                                 ['HeI 20520', 20520],  #
                                 ['HeII 8239', 8239.26],
                                 ['HeII 20580', 20580],
                                 ['CaII 8544', 8544.44],
                                 ['CaII 8498', 8498],
                                 ['CaII 8662', 8662],
                                 ['[FeII] 9202', 9202],  #
                                 ['[FeII]12600', 12600],
                                 ['[FeII]16440', 16440],
                                 ['H$_{2}$ 19570', 19570],
                                 ['H$_{2}$ 21210', 21210],
                                 [r'Pa$_{\alpha}$', 18750.1],
                                 [r'Pa$_{\beta}$', 12818.1],
                                 [r'Pa$_{\gamma}$', 10938],
                                 [r'Pa$\delta$', 10049.8],
                                 [r'Pa$\epsilon$', 9545],  # 9546.2
                                 ['NaI', 22080]]

    data, header, gal_name, wavelengths = read_fits(fits_file)
    # Giving units to flux
    flux_values = get_flux_values(data)
    # Giving units to wavelength and converting to microns
    wave_values = get_wave_values(wavelengths)

    # Create whole spectreum
    spec = pyspeckit.Spectrum(data=flux_values, xarr=wave_values,
                              header=header, unit='erg/s/cm^2/AA')
    # Get plotting range
    x_min = wave_values[0] - .01
    x_max = wave_values[-1] + .01
    ''' Look at last fifth of spectrum, excluding last 50 elements
             and gets min flux value for whole spectrum plot'''
    fifth = np.floor(len(wave_values) * .2)
    y_min = np.min(flux_values[-fifth:-50]) - .1
    ''' Look at first fourth of spectrum
            and get max flux value for whole spectrum plot'''
    fourth = np.floor(len(wave_values) * .25)
    y_max = np.max(flux_values[0:fourth]) + .2

    if args.plot_type == 'whole':
        new_spec = spec.copy()
        new_spec.plotter(
            xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max)
    elif args.plot_type == 'none':
        print('You have chosen not to plot or fit anything.')
        sys.exit(0)
    # Plot each reference wavelength, one at a time.
    for line in range(0, len(reference_wavelengths)):
        line_name = reference_wavelengths[line][0]
        line_value = reference_wavelengths[line][1] * u.AA
        # Converting reference_wavelengths units to microns
        line_value = line_value.to(u.micron).value
        base_min = line_value - .03
        base_max = line_value + .03
        line_min = line_value - .01
        line_max = line_value + .01
        # Not working, not sure why, looks the same as above
        # fwhm_guess = .5 * peak_guess #Need corresponding x value to this calculation
        # Plotting my lines on the graph
        if (line_value > x_min) & (line_value < x_max):
            if args.plot_type == 'lines':
                line_index = np.argwhere((wave_values>line_min) & (wave_values<line_max)).flatten()
                flux_data = flux_values[line_index]
                wave_data = wave_values[line_index]
                x = wave_data
                y = flux_data
                peak = np.max(y) 
                mean = np.mean(y)
                std = np.std(y)
                print(peak, line_value, std)
                gg_init = models.Gaussian1D(peak, line_value, 0.001)  #+ models.Gaussian1D(peak, mean, std)
                fitter = fitting.SLSQPLSQFitter()
                gg_fit = fitter(gg_init, x, y)

                # Plot the data with the best-fit model
                plt.figure(figsize=(8,5))
                plt.plot(x, y, 'ko')
                plt.plot(x, gg_fit(x), 'r-', lw=2)
                plt.xlabel('Position')
                plt.ylabel('Flux')

                # plt.axvline(line_value, color='b', linestyle='--')
                # lbl_pos = (y_max - y_min) * 0.85  # at 85% up plot
                # plt.text(line_value, lbl_pos, line_name, rotation=45)
                plt.show()
            #     new_spec = spec.copy()
            #     new_spec.plotter(
            #         xmin=base_min, xmax=base_max, ymin=y_min, ymax=y_max)
            # If plotting lines, show plot and then save.
            # if args.plot_type == 'lines':
            #     plt.show()
            #     output_filename = args.directory
            #     output_filename += gal_name + '_' + line_name + '.jpeg'
            #     plt.savefig(output_filename, format='jpeg', dpi=300)
    plt.show()
    output_filename = args.directory + gal_name + '.jpeg'
    plt.savefig(output_filename, format='jpeg', dpi=300)

##########################


if __name__ == '__main__':
    main()
