"""Submodule conatining main classes for wavefront phase analysis.

This module contains the main classes that are used to define wavefront phase objects and perform wavefront analysis. """

import numpy as np
import sys
import copy
import numpy.lib.index_tricks as itricks
from scipy import optimize, ndimage
from scipy.signal import periodogram
from astropy.io import fits

#from telescope import Telescope, TELESCOPE_DICT
#from waveplot import implot, Zplot, PSDplot

###############################################################################
##
## Image Class
##
###############################################################################

class Image(object):
    """Class representing an image. Images have the following properties:
    
    Attributes: 
        data array: A Numpy array of the image.
        pixel scale: A float tracking the pixel scale.
    """
    
    def __init__(self, data, pixscale = 7.77/43):
        """Return an Image object whose data is a Numpy array and pixel scale is *pixscale*"""
        self.data = data
        self.pixscale = pixscale
    
    @classmethod
    def import_fits(self, file_path, pixscale = 7.77/43):
        """Return a single 2-D image array or a 3-D image data cube."""
        hdulist = fits.open(file_path, memmap=True)
        data = hdulist[0].data
        
        shape = data.shape
        
        ## Create Image objects
        if len(shape) == 2:
            return cls(data,pixscale)
        elif len(shape) == 3:
            image_list = []
            
            ## Iterate over data cube and intianlize Image objects
            for i in range(data.shape[0]):
                single_image_data = data[i,:,:]
                image_list.append(cls(single_image_data,pixscale))
            return image_list
        else:
            print shape
            sys.exit("FITs Read Error: Must be 2-D or 3-D Image datacube")
            
        def export_fits(self, mask=None, **kwargs):
            """Export Image as a NumPy array to a FITS file"""
            
            ## Check key word arguments
            save_file = kwargs.pop('save_file', 'image.fits')
            fill_value = kwargs.pop('fill_value',0.)
            
            ##  Check if mask provided matches data shape
            if self.is_valid_mask(mask):
                masked_data = np.ma.MasedArray()
    