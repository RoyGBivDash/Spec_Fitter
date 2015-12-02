import numpy as np
from astropy.wcs import WCS
from astropy import units as u
from astropy.io import fits



class SpecInfo(object):
    """Create all information needed to fit galaxy"""
    def __init__(self, filename):
        super(SpecInfo, self).__init__()
        self.filename = filename
        data, header, gal_name, wavelengths = self.read_fits()
        self.data = data
        self.header = header
        self.gal_name = gal_name
        self.wavelengths = wavelengths

    def read_fits(self):
        """Read FITS file and 
            returns data, header, galaxy name and wavelengths"""
        hdulist = fits.open(self.filename)  # Opens fits file
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

    def get_flux_values(self):
        """Extract flux values from FITS data; return unitless values"""
        flux_angstrom = self.data
        flux_micron = flux_angstrom * 10000
        flux_unit = u.erg / (u.cm**2 * u.s * u.micron)
        flux = flux_micron * flux_unit * 1e10
        flux_values = flux.value  # Just the numbers
        self.flux = flux
        return flux_values

    def get_wave_values(self):
        """Convert wavelength units to microns, return values only"""
        wv_unit = u.AA
        wavelengths = self.wavelengths * wv_unit
        # Converting Angstrom to microns
        wavelengths = wavelengths.to(u.micron)
        wave_values = wavelengths.value  # Just the numbers, in microns
        return wave_values

    def find_max_y_in_x_range(self, x_list, y_list, x_min, x_max):
        spec_line = x_list[x_min:x_max]  # Limits whole array to wanted range
        max_y = None
        loc_max_y = 0
        for i in range(len(x_list)):
            if x_list[i] >= x_min and x_list[i] <= x_max:
                if max_y is None or y_list[i] > max_y:
                    max_y = y_list[i]
                    loc_max_y = x_list[i]
        return max_y, loc_max_y
