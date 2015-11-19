# Spec_Fitter
Python code to fit Gaussians to spectral lines. Lines listed are in the infrared range.

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
