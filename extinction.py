# ----------------------------------------------------------------------------
#
# TITLE - extinction.py
# AUTHOR - James Lane
# PROJECT - archipelago
# CONTENTS:
#	1.
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Example'''
__author__ = "James Lane"


#Imports
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import pdb


def GetDensityStats(hdu, region, filename=None):
    '''
    GetDensityStats:

    Use NICEST star density maps to get information about the median, mean,
    and extrema of the stars used to create an extinction map. Outputs
    region+'_density_stats.txt' with star count statistics.

    Args:
        hdu (astropy HDU object): The HDU of the star density map.
        region (str): The name of the region.

    Returns:
        None
    '''

    # If a filename is specified then load the hdu from there.
    if filename != None:
        hdu = fits.open(filename)
    ##fi

    data = hdu[0].data
    no_zero = np.where(data > 0)
    data = data[no_zero]
    med_data = np.median(data)
    mean_data = np.average(data)
    min_data = np.amin(data)
    max_data = np.amax(data)

    outfile = open(region+'_density_stats.txt', 'w')
    outfile.write('Only non-zero data contributes.\n')
    outfile.write('median, mean, min, max\n')
    outfile.write('{:.4f}\n{:.4f}\n{:.4e}\n{:.4f}'.format(med_data,mean_data,min_data,max_data))
    outfile.close()
#def

def GetMaplbExtents(hdu, region, filename=None):
    '''
    GetMaplbExtents:

    Get the maximum and minimum extents for a .fits map in galactic longitude
    and latitude.

    Args:
        hdu (Astropy HDU object): The HDU of the map.
        region (string): The name of the region.
        filename (string): A filename to open if no hdu is provided.

    Returns:
        None

    Outputs:
        ...extents.txt: A list of the extents.
    '''

    # If a filename is specified then load the hdu from there.
    if filename != None:
        hdu = fits.open(filename)
    ##fi

    gbs_data = hdu[0].data
    gbs_wcs = wcs.WCS(hdu[0].header)
    gbs_x_len = hdu[0].header['NAXIS1']
    gbs_y_len = hdu[0].header['NAXIS2']
    xm, ym = np.meshgrid(   np.arange(1,gbs_x_len+1).astype(int),
                            np.arange(1,gbs_y_len+1).astype(int)   )
    gbs_xpix = xm.flatten()
    gbs_ypix = ym.flatten()
    gbs_zpix = np.ones_like(gbs_xpix)

    #pdb.set_trace()
    where_no_nan = np.isfinite(gbs_data.flatten())
    gbs_xpix = gbs_xpix[where_no_nan]
    gbs_ypix = gbs_ypix[where_no_nan]
    gbs_zpix = gbs_zpix[where_no_nan]

    _ra, _dec, _ = gbs_wcs.all_pix2world( gbs_xpix, gbs_ypix, gbs_zpix, 1 )
    coords = SkyCoord( _ra, _dec, unit=u.deg )
    coords_gal = coords.transform_to('galactic')
    gbs_l = coords_gal.l.value
    gbs_b = coords_gal.b.value

    #pdb.set_trace()

    lmin = gbs_l.min()
    lmax = gbs_l.max()
    bmin = gbs_b.min()
    bmax = gbs_b.max()

    outfile = open(region+'_extent.txt', 'w')
    outfile.write('l min, l max, b min, b max\n')
    outfile.write('{:.4f}, {:.4f}, {:.4f}, {:.4f}'.format(lmin,lmax,bmin,bmax))
    outfile.close()
#def
