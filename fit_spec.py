#!/usr/bin/env python

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import pyspeckit
import specinfo


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

    spec_info = specinfo.SpecInfo(fits_file)
    flux_values = spec_info.get_flux_values()
    wave_values = spec_info.get_wave_values()

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

    # Create whole spectreum
    spec = pyspeckit.Spectrum(data=flux_values, xarr=wave_values,
                              header=spec_info.header, unit='erg/s/cm^2/AA')
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
    for line in reference_wavelengths:
        line_name = line[0]
        line_value = line[1] * u.AA
        # Converting reference_wavelengths units to microns
        line_value = line_value.to(u.micron).value
        base_min = line_value - .03
        base_max = line_value + .03
        line_min = line_value - .01
        line_max = line_value + .01
        peak_guess, peak_location = spec_info.find_max_y_in_x_range(
            wave_values, flux_values, line_min, line_max)
        print(line_name)
        print("peak: ", peak_guess)
        print("peak location: ", peak_location)
        print("line value: ", line_value)
        # Plotting my lines on the graph
        if (line_value > x_min) & (line_value < x_max):
            if args.plot_type == 'lines':
                new_spec = spec.copy()
                new_spec.plotter(
                    xmin=base_min, xmax=base_max, ymin=y_min, ymax=y_max)
                new_spec.baseline(order=1, exclude=(line_min, line_max), fit_plotted_area=True,)
            #  , fit_original=True, fit_plotted_area=True,)
            #  linewidth=2, baseline_fit_color='r'
            #  highlight_fitregion=True, reset_selection=True,

            plt.axvline(line_value, color='b', linestyle='--')
            lbl_pos = (y_max - y_min) * 0.85  # at 85% up plot
            plt.text(line_value, lbl_pos, line_name, rotation=45)
            # If plotting lines, show plot and then save.
            if args.plot_type == 'lines':
                plt.show()
                output_filename = args.directory
                output_filename += spec_info.gal_name + \
                    '_' + line_name + '.jpeg'
                plt.savefig(output_filename, format='jpeg', dpi=300)
    plt.show()
    output_filename = args.directory + spec_info.gal_name + '.jpeg'
    plt.savefig(output_filename, format='jpeg', dpi=300)

##########################


if __name__ == '__main__':
    main()
