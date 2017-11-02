# ----------------------------------------------------------------------------
#
# TITLE - proto.py
# AUTHOR - James Lane
# PROJECT - archipelago
# CONTENTS:
#	1. AssociateFragments
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''
Functions for dealing with protostars in the archipelago project.
'''

__author__ = "James Lane"

#Imports
import numpy as np
from astropy.io import fits
from astropy import wcs as worldcoords
import glob
import pdb
import archipelago.misc
import archipelago.tasks

def AssociateFragments(region):
    '''
    AssociateFragments:

    Associate a protostellar catalog with a catalog of fragments. Each island
    contains at least one fragment. Each of these fragments will be checked
    to see whether or not a protostar is associated with it. In order to be
    associated, a protostar must lie within the footprint of emmission of the
    fragment. A close association occurs when the protostar is also within
    one beams width of the fragment peak.

    A catalog with the name region+'_frag_protos.txt' is the output product.
    It contains as its columns: the unque fragment identifier, the protostar
    identity number, the protostar class, whether or not it is close (1 or 0),
    and the slope of the SED. Non-applicable values will be nan.

    Args:
        region (string) - The region name.

    Returns:
        None
    '''

    # Check if the region is in Orion.
    orion = archipelago.tasks.GetInOrion(region)
    proto_catalog = 'Dunham_Catalog.FIT'
    if orion == True: proto_catalog = 'Megeath_Catalog.FIT'

    # Get fragment cutout files.
    sdf_files = glob.glob('fragment_sdfs/*.sdf')
    fits_files = glob.glob('fragment_sdfs/*.fits')
    archipelago.misc.sort_nicely(sdf_files)
    archipelago.misc.sort_nicely(fits_files)

    # Open fragment catalog.
    frag_data = fits.open('fragment_info/'+region+'_peak-cat.fits')[1].data
    n_frags = len(frag_data['ID'])
    frag_id = frag_data['ID']
    frag_ra = frag_data['RA']
    frag_dec = frag_data['DEC']

    # Open protostar catalog.
    proto_data = fits.open(proto_catalog)[1].data
    proto_ra = proto_data['RA']
    proto_dec = proto_data['Dec']
    proto_class = proto_data['Class']
    proto_id = proto_data['ID']
    proto_slope = proto_data['Slope']

    # Prepare arrays.
    frag_pid = np.empty(n_frags,dtype='S6')
    frag_class = np.empty(n_frags,dtype='S6')
    frag_slope = np.empty(n_frags,dtype='S6')
    frag_pid[:] = 'nan'
    frag_slope[:] = 'nan'
    frag_class[:] = 'None'
    frag_close = np.zeros(n_frags,dtype=int)
    frag_index = -1

    # Loop over each of the island extents
    for i in range(len(fits_files)):
        print(str(i+1))
        # Read in the island extent file
        isl_hdu = fits.open(fits_files[i])
        # Construct WCS object for the file.
        wcs = worldcoords.WCS(isl_hdu[0].header, naxis=2)
        #Find the number of fragments to deal with on this island.
        n_isl_frags = int(np.nanmax( isl_hdu[0].data ))

        # Loop over these fragments
        for j in range(n_isl_frags):

            frag_index += 1
            print(frag_index)

            # Find the closest protostar.
            dists = np.sqrt(np.square(np.multiply((proto_ra-frag_ra[frag_index]),
                                    np.cos((proto_dec+frag_ra[frag_index])/2)))+
                                    np.square(proto_dec-frag_dec[frag_index]))
            proto_close = np.argmin(dists)
            proto_deg = np.array(   [proto_ra[proto_close],proto_dec[proto_close]],
                                    ndmin=2)
            proto_x, proto_y = wcs.wcs_world2pix(proto_deg, 1).astype(int)[0,:]

            if  proto_x < 0 or proto_x >= isl_hdu[0].header['NAXIS1']: continue
            if  proto_y < 0 or proto_y >= isl_hdu[0].header['NAXIS2']: continue
            if isl_hdu[0].data[proto_y,proto_x] == j+1:
                frag_class[frag_index] = proto_class[proto_close]
                frag_slope[frag_index] = proto_slope[proto_close]
                frag_pid[frag_index] = proto_id[proto_close]
                if dists[proto_close]*3600.0 < 15.0:
                    frag_close[frag_index] = 1
                ##fi
            ##fi
        ###j
    ###i

    # Output the results - should output in .FIT format for permenance.
    frag_slope[ np.where(frag_slope == '')[0] ] = 'nan'
    data_out = np.array([frag_id,frag_pid,frag_class,frag_close,frag_slope],
                        dtype='object')
    data_fmt = ('%s', '%s', '%s', '%i', '%s')

    np.savetxt( 'fragment_info/'+region+'_frag_protos.txt', data_out.T,
                delimiter='\t', fmt=data_fmt)
