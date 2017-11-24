# ----------------------------------------------------------------------------
#
# TITLE - analysis.py
# AUTHOR - James Lane
# PROJECT - Archipelago
# CONTENTS:
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Main analysis program for the Archipelago project.'''

__author__ = "James Lane"

# Imports
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import patches as mpatches
from matplotlib import lines as mlines
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import table
from scipy.stats import linregress
import jlastro.core
import aplpy
import sys
import archipelago.tasks
import pdb
import os

class RegionProperties:
    '''
    RegionProperties:

    Class that handles the properties of the different regions in the
    archipelago analysis.

    Args:
        region (str) - Name of the region.
        noise (float) - The average noise across the region in mJy/As.Sq.
        distance (float) - Distance in parsecs [450, Orion default]
        beamsize (float) - The FWHM of the beam in arcseconds [15.0]
        pixel_resolution (float) - Map pixel size in arcseconds [3, IR3 default]
    '''
    def __init__(   self,
                    region,
                    noise = 0.045,
                    distance = 450.0,
                    beamsize = 15.0,
                    pixel_resolution = 3.0,
                ):
        self.region = str(region)
        self.noise = float(noise)
        self.orion = archipelago.tasks.GetInOrion(region)
        self.distance = float(distance)
        self.beamsize = float(beamsize)
        self.pixel_resolution = float(pixel_resolution)
    #def
#cls

def PrepareData(region,
                noise,
                distance,
                extinction_map,
                img_pixelres=3.0,
                extinction_pixelres=3.0,
                beamsize=15.0,
                extinction_test=False
                ):

    '''
    PrepareData:

    Do the primary data preparation for the Archipelago analysis. Read in the
    data for islands and fragments and derive all quantities for both. Examine
    island boundaries. Determine association of protostars with flux from
    Islands and Fragments. Output information in a catalog.

    Args:
        region (string) - GBS name of the region.
        noise (float) - Noise in the main map in mJy/As.Sq.
        distance (float) - Distance to the molecular cloud.
        extinction_map (string) - Path to extinction map.
        img_pixelres (float) - Pixel resolution for the map in arcseconds (3)
        extinction_pixelres (float) - Pixel resolution for the extinction map
                                        in arcseconds (3)
        beamsize (float) - Size of the beam in arcseconds (15.0)

    Returns:
        None

    Outputs:

    '''

    ############################################################################

    #############
    # PREPARATION
    #############

    # Determine if the region is in Orion.
    orion = archipelago.tasks.GetInOrion(region)

    # Write distance directories to delineate results:
    distance_dir = 'D'+str(int(distance))
    os.system('mkdir -p {Test,Catalogs,Plots}/'+distance_dir)


    ###########
    # READ DATA
    ###########
    print('Reading data...')

    # Read in the regional fits image.
    try:
        img_hdu = fits.open('../'+region+'.fits')
        img_name = '../'+region+'.fits'
    except FileNotFoundError:
        try:
            img_hdu = fits.open('../'+region+'_clipped.fits')
            img_name = '../'+region+'_clipped.fits'
        except FileNotFoundError:
            sys.exit('Exited... No file found.\n\n')
        #try
    #try

    # Define the image data.
    img_data = img_hdu[0].data
    img_hdr = img_hdu[0].header
    img_wcs = wcs.WCS(img_hdr)

    # Read in the fits catalog data: island data, island derived properties,
    # and fragment data.
    data_isl= fits.getdata('ClumpFind/'+region+'_islands_outcat.FIT', 0)
    # data_isl_der = fits.getdata('ClumpFind/'+region+'_islands_outcat_derived.fits', 0)
    data_frag = fits.getdata('ClumpFind/'+region+'_peak-cat.fits', 0)

    # Island measured properties (From Clumpfind)
    n_islands = len(data_isl['PIDENT'])
    island_id = data_isl['PIDENT']                  # Unitless
    island_rapeak = data_isl['Peak1']               # Degrees RA
    island_decpeak = data_isl['Peak2']              # Degrees Dec
    island_volume = data_isl['VOLUME']              # Square Arcseconds
    island_peakflux = data_isl['PEAK']              # mJy/As.Sq.
    island_sumflux = data_isl['SUM']                # mJy/As.Sq.

    # Island derived properties (From Steve code) -- Not distance independant!
    # n_islands = len(data_isl_der)
    # island_mass = data_isl_der['mass']              # Solar masses
    # island_radiusj = data_isl_der['rj']/206264.8    # AU
    # island_masspmj = data_isl_der['mass_over_mj']   # Unitless
    # island_conc = data_isl_der['conc']              # Unitless

    # Fragment properties (From Clumpfind)
    n_frags = len(data_frag)                        # Unitless
    frag_rapeak = data_frag['RA']                   # Degrees RA
    frag_decpeak = data_frag['DEC']                 # Degrees Dec
    frag_peakflux = data_frag['PEAK_FLUX']          # mJy/As.Sq.
    frag_sumflux = data_frag['SUM']                 # mJy/As.Sq.
    frag_volume = data_frag['VOLUME']               # Square Arcseconds

    # Read in Fragment protostar identifications, must be forced to dtype=U
    frag_class,frag_close = np.loadtxt('ClumpFind/'+region+'_frag_protos.txt',
                                        usecols=[2,3], dtype='S6').T
    frag_close = frag_close.astype('U6').astype(int)
    frag_class = frag_class.astype('U6')

    ############################################################################

    #############################
    # DERIVE FRAGMENT PROPERTIES:
    #############################
    print('Deriving Fragment properties...')

    # Convert fragment peak and sum flux to Jy:
    frag_peakflux_jy = frag_peakflux * (img_pixelres**2) / 1000.0
    frag_sumflux_jy = frag_sumflux * (img_pixelres**2) / 1000.0

    # Convert fragment peak flux to Jy/beam, known 229 sq.arcsec per beam.
    frag_peakflux_jybm = frag_peakflux * (229) / 1000.0
    noise_jybm = noise * (229) / 1000.0

    # Convert volume (area) into a "radius":
    frag_radius = distance * np.sqrt( frag_volume/np.pi ) / 206264.8 # parsecs

    # Calculate mass: see Mairs 2016
    frag_mass = 2.63 * ((distance/450)**2) * (frag_sumflux_jy) # Solar Masses

    # Calculate peak column number density assuming standard values, in cgs.
    frag_numcoldenspeak = 1.19E23 * frag_peakflux_jybm

    # Calculate MJeans: see Mairs 2016 (also includes T in units of 15K) -- bad way of doing it
    frag_massj = 2.9 * (frag_radius / 0.07) # Units of Msol

    # Compute the fragment mass density, in cgs
    frag_mass_density = (1.616E-23) * np.divide(frag_mass,np.power(frag_radius,3))

    # Calculate the Jeans Mass assuming the density is correct.
    # Calculate the density derived Jeans mass
    mH = 1.67E-27   # kg
    T = 15          # Kelvin
    mu = 2.33       # unitless
    kB = 1.381E-23  # SI
    cs = np.sqrt( (kB*T) / (mH*mu) )
    # Factor of 1000 accounts for CGS to SI.
    frag_massjdens = (2.691E-15) * (cs**3) * np.power(frag_mass_density*1000,-0.5)

    # Compute the jeans radius in parsecs
    frag_radiusj = (2.566E-11) *  np.divide(1,np.sqrt(frag_mass_density))

    # Compute the volumetric number density in cgs
    frag_numvoldenspeak = (2.570E23) * frag_mass_density

    # Calculate concentration
    frag_conc = 1 - np.divide( (1.13 * beamsize**2 * frag_sumflux_jy) ,
                    np.multiply( frag_volume , frag_peakflux_jybm ) )

    # Calculate mass in units of Jeans mass, use the density derived one.
    frag_masspmj = np.divide(frag_mass,frag_massjdens)

    # Calculate radius in units is Jeans radius
    frag_radiusprj = np.divide(frag_radius,frag_radiusj)

    # Find minimum fragment size, 7 pixels
    frag_minpix = 7 #pixels
    frag_minarea = frag_minpix * (img_pixelres)**2 #Square arcseconds
    frag_minradius = distance * np.sqrt( frag_minarea/np.pi ) / 206264.8 # Units of pc

    ############################################################################

    ###########################
    # DERIVE ISLAND PROPERTIES:
    ###########################
    print('Deriving Island properties...')

    # Convert island peak and sum flux to Jy:
    island_peakflux_jy = island_peakflux * (img_pixelres**2) / 1000.0
    island_sumflux_jy = island_sumflux * (img_pixelres**2) / 1000.0

    # Convert island peak flux to Jy/beam, known 229 sq.arcsec per beam.
    island_peakflux_jybm = island_peakflux * (229) / 1000.0

    # Convert volume (area) into a "radius":
    island_radius = distance * np.sqrt( island_volume/np.pi ) / 206264.8 # Units of pc

    # Calculate mass: see Mairs 2016
    island_mass = 2.63 * ((distance/450)**2) * (island_sumflux_jy) # Solar Masses

    # Compute the island jeans mass -- Uses bad Jeans mass.
    island_massj = 2.9 * (island_radius / 0.07) # Units of Msol

    # Compute the island mass density, in cgs
    island_mass_density = (1.616E-23) * np.divide(island_mass,np.power(island_radius,3))

    # Compute the Jeans Mass assuming the density is correct.
    island_massjdens = (2.691E-15) * (cs**3) * np.power(island_mass_density*1000,-0.5)

    # Compute the jeans radius in parsecs
    island_radiusj = (2.566E-11) *  np.divide(1,np.sqrt(island_mass_density))

    # Calculate concentration
    island_conc = 1 - np.divide( (1.13 * beamsize**2 * island_sumflux_jy) ,
                    np.multiply( island_volume , island_peakflux_jybm ) )

    # Recompute the mass in units of jeans mass using the
    island_masspmj = np.divide(island_mass, island_massjdens)

    # Calculate radius in units of Jeans radius.
    island_radiusprj = np.divide(island_radius,island_radiusj)

    # Compute the volumetric number density in cgs
    island_numvoldenspeak = (2.570E23) * island_mass_density

    # Find minimum island size, beamsize squared
    island_minarea = beamsize**2 #pixels
    island_minradius = distance * np.sqrt( island_minarea/np.pi ) / 206264.8 # Units of pc

    ############################################################################

    ######################
    # PROTOSTAR PROPERTIES
    ######################
    print('Getting YSO properties...')

    # Read in the Protostar catalog, either Dunham or Megeath depending on location.
    proto_catalog = 'Dependencies/Dunham_Catalog.FIT'
    if orion == True: proto_catalog = 'Dependencies/Megeath_Catalog.FIT'
    proto_data = fits.getdata(proto_catalog, 0)

    # Unpack data.
    n_proto = len(proto_data)
    proto_ra = proto_data['RA']
    proto_dec = proto_data['Dec']
    proto_class = proto_data['Class']

    # Use the image WCS objects to convert angular positions into pixel positions.
    proto_coords_deg = np.array([proto_ra,proto_dec,np.zeros(n_proto)]).T
    proto_xpix, proto_ypix = img_wcs.all_world2pix(proto_coords_deg, 1).T.astype(int)[:2]

    # Cut the catalog to only inlude protostars in the boundaries of the image.
    proto_cut_inimg = np.where((0 < proto_xpix) & (proto_xpix < img_hdr['NAXIS1']-5) &
                            (0 < proto_ypix) & (proto_ypix < img_hdr['NAXIS2']-5))[0]
    proto_xpix = proto_xpix[proto_cut_inimg]
    proto_ypix = proto_ypix[proto_cut_inimg]
    proto_class = proto_class[proto_cut_inimg]
    proto_ra = proto_ra[proto_cut_inimg]
    proto_dec = proto_dec[proto_cut_inimg]

    # Cut the catalog to only include protostars ontop of emission
    proto_cut_onflux = np.isfinite(img_hdu[0].data[0,proto_ypix,proto_xpix])
    proto_ypix = proto_ypix[proto_cut_onflux]
    proto_xpix = proto_xpix[proto_cut_onflux]
    proto_class = proto_class[proto_cut_onflux]
    proto_ra = proto_ra[proto_cut_onflux]
    proto_dec = proto_dec[proto_cut_onflux]
    jlastro.core.make_tab('ysos_ingbs.tab',proto_ra,proto_dec)

    # Define the flux density corresponding to each protostar.
    proto_flux = img_hdu[0].data[0,proto_ypix,proto_xpix] * (229) / 1000.0

    # Find the distance from each protostar to the nearest fragment peak.
    proto_skycoord = SkyCoord(ra=proto_ra, dec=proto_dec, unit=u.deg)
    frag_skycoord = SkyCoord(ra=frag_rapeak, dec=frag_decpeak, unit=u.deg)
    skyidx, skysep, skydist = proto_skycoord.match_to_catalog_sky(frag_skycoord)
    proto_fragdist = distance * skysep.deg * (np.pi/180)

    ############################################################################

    ###################
    # ISLAND PROPERTIES
    ###################
    print('Matching Fragments to Islands...')

    # Read in the Island boundary map.
    island_bounds_hdu = fits.open('ClumpFind/'+region+'_islands_out.fits')
    island_bounds_wcs = wcs.WCS(island_bounds_hdu[0].header)
    island_bounds_data = island_bounds_hdu[0].data

    # Convert fragment positions to pixels.
    island_coords_deg = np.array([island_rapeak,island_decpeak,np.zeros(n_islands)]).T
    frag_coords_deg = np.array([frag_rapeak,frag_decpeak,np.zeros(n_frags)]).T
    island_xpeak, island_ypeak = island_bounds_wcs.all_world2pix(island_coords_deg, 1).T[:2]
    frag_xpix, frag_ypix = island_bounds_wcs.all_world2pix(frag_coords_deg, 1).T.astype(int)[:2]

    # Check to see which islands are monolithic and which have no fragments.
    island_nfrags = np.zeros(n_islands, dtype=int)
    frag_islandnum = np.zeros(n_frags, dtype=int)
    frag_monolith = np.zeros(n_frags, dtype=int)

    # Loop over all fragments, determining which island they lie in.
    for i in range(n_frags):
        if np.isfinite(island_bounds_data[0,frag_ypix[i],frag_xpix[i]]):
            yind=int(frag_ypix[i])
            xind=int(frag_xpix[i])
            island_nfrags[int(island_bounds_data[0,yind,xind])-1] += 1
            frag_islandnum[i] = island_bounds_data[0,yind,xind]
        ##fi
    ###i

    #Check for monolithic fragments (belonging to complex or monolithic islands)
    for i in range(n_frags):
        if island_nfrags[frag_islandnum[i]-1] <= 1:
            frag_monolith[i]=1
        ##fi
    ###i

    # Cut islands that have no fragments.
    fragged_islands = np.where(island_nfrags > 0)[0]
    island_rapeak = island_rapeak[fragged_islands]
    island_decpeak = island_decpeak[fragged_islands]
    island_volume = island_volume[fragged_islands]
    island_xpeak = island_xpeak[fragged_islands]
    island_ypeak = island_ypeak[fragged_islands]
    island_nfrags = island_nfrags[fragged_islands]
    island_pident = island_id[fragged_islands]
    island_radius = island_radius[fragged_islands]
    island_radiusj = island_radiusj[fragged_islands]
    island_radiusprj = island_radiusprj[fragged_islands]
    island_mass = island_mass[fragged_islands]
    island_massj = island_massj[fragged_islands]
    island_massjdens = island_massjdens[fragged_islands]
    island_masspmj = island_masspmj[fragged_islands]
    island_conc = island_conc[fragged_islands]
    island_numvoldenspeak = island_numvoldenspeak[fragged_islands]
    n_islands = len(island_radiusj)
    island_monolith = np.zeros(n_islands, dtype=int)
    island_monolith[np.where(island_nfrags == 1)[0]] = 1
    island_skycoord = SkyCoord(ra=island_rapeak, dec=island_decpeak, unit=u.deg)

    ############################################################################

    ########################
    # ISLANDS AND PROTOSTARS
    ########################
    print('Matching Islands with YSOs...')

    island_class = np.empty(n_islands,dtype='U6')
    island_class[:] = 'None'
    island_close = np.zeros(n_islands,dtype=int)
    island_nproto = np.zeros(n_islands,dtype=int)

    # Look for the protostar associated with each island, only if there are
    # protostars in the map. Can crash in certain regions where there are no
    # protostars present.
    if len(proto_skycoord) > 0:
        for i in range(n_islands):
            # Distances to each protostar, find the closest one:
            ipidx, ipsep, ipdist = island_skycoord[i].match_to_catalog_sky(proto_skycoord)
            proto_close = ipidx
            proto_close_xpix = proto_xpix[proto_close]
            proto_close_ypix = proto_ypix[proto_close]

            # Find total number of protostars in this Island
            proto_in_isl = np.where(island_bounds_data[0,proto_ypix,proto_xpix])[0]
            island_nproto = len(proto_in_isl)

            # Find the closest YSO to the island peak
            proto_islnum = island_bounds_data[0,proto_close_ypix,proto_close_xpix]
            if proto_islnum == island_pident[i]:
                island_class[i] = proto_class[proto_close]
                if ipsep.arcsecond < 15:
                    island_close[i] = 1
                ##fi
            ##fi
        ###i
    ##fi

    ############################################################################

    ############
    # EXTINCTION
    ############
    print('Extinction and cumulative mass...')

    # Read in NICEST extinction data.
    ext_hdu = fits.open(extinction_map)
    ext_data = ext_hdu[0].data
    conf_ext_data = np.array(ext_data)
    ext_wcs = wcs.WCS(ext_hdu[0].header)

    # Get NICEST pixels in the SCUBA-2 footprint.
    where_in_s2 = np.isfinite(img_data)
    ext_evals = ext_data[ where_in_s2 ]

    # If testing is required then output some diagnostic plots.
    if extinction_test == True:
        test_data = np.zeros( ext_data.shape )
        test_data[:,:,:] = np.nan
        test_data[ where_in_s2 ] = ext_evals
        plt.clf()
        plt.imshow(test_data[0,::-1,:])
        plt.savefig('Test/'+distance_dir+'ext_data_select.png', dpi=300)
        plt.clf()
        plt.imshow(ext_data[0,::-1,:])
        plt.savefig('Test/'+distance_dir+'ext_data.png', dpi=300)
        plt.clf()
    ##fi

    # Get NICEST pixels in the Islands footprint.
    where_islands = np.isfinite( island_bounds_data.flatten() )
    where_sig_islands = np.in1d( island_bounds_data.flatten()[ where_islands ],
                                 fragged_islands+1 )
    #pdb.set_trace()
    where_islands = np.where( where_islands )[0]
    where_sig_islands = np.where( where_sig_islands )[0]
    where_good_islands = where_islands[ where_sig_islands ]
    s2_evals = ext_data.flatten()[ where_islands ][ where_sig_islands ]
    s2_fvals = img_data.flatten()[ where_islands ][ where_sig_islands ]
    # s2_evals = s2_evals.flatten()
    # s2_fvals = s2_fvals.flatten()

    # If testing is required then output some diagnostic plots.
    if extinction_test == True:
        test_s2_data = np.zeros( img_data.shape )
        test_s2_data[:,:,:] = np.nan
        test_s2_data = test_s2_data.flatten()
        test_s2_data[ where_good_islands ] = s2_fvals / 10
        test_s2_data = test_s2_data.reshape( img_data.shape )
        test_ext_data = np.zeros( ext_data.shape )
        test_ext_data[:,:,:] = np.nan
        test_ext_data = test_ext_data.flatten()
        test_ext_data[ where_good_islands ] = s2_evals
        test_ext_data = test_ext_data.reshape( ext_data.shape )
        masked_ext_data = np.ma.masked_where(test_ext_data==0, test_ext_data)
        #pdb.set_trace()
        plt.imshow(test_ext_data[0,::-1,:])
        plt.savefig('Test/'+distance_dir+'s2_ext_data_select.png', dpi=300)
        plt.clf()
        plt.imshow(test_s2_data[0,::-1,:])
        plt.savefig('Test/'+distance_dir+'s2_islands.png', dpi=300)
        plt.clf()
    ##fi

    # Get extinction of protostars. Using already culled catalog of protostars,
    # so only including protostars that lie on emission.

    proto_evals = ext_data[0, proto_ypix-1, proto_xpix-1]
    proto_evals = proto_evals.flatten()

    # Separate the protostar sample into the constituent classes
    if orion == True:
        where_p = np.where(proto_class == 'P')
        where_d = np.where(proto_class == 'D')
        where_fp = np.where(proto_class == 'FP')
        where_rp = np.where(proto_class == 'RP')

        p_evals = ext_data[0, proto_ypix[where_p]-1, proto_xpix[where_p]-1 ]
        d_evals = ext_data[0, proto_ypix[where_d]-1, proto_xpix[where_d]-1 ]
        fp_evals = ext_data[0, proto_ypix[where_fp]-1, proto_xpix[where_fp]-1 ]
        rp_evals = ext_data[0, proto_ypix[where_rp]-1, proto_xpix[where_rp]-1 ]

        p_evals = p_evals.flatten()
        d_evals = d_evals.flatten()
        fp_evals = fp_evals.flatten()
        rp_evals = rp_evals.flatten()
    ##fi
    if orion == False:
        where_c01 = np.where(proto_class == '0+1')
        where_f = np.where(proto_class == 'F')
        where_c2 = np.where(proto_class == '2')
        where_c3 = np.where(proto_class == '3')

        c01_evals = ext_data[0, proto_ypix[where_c01]-1, proto_xpix[where_c01]-1 ]
        f_evals = ext_data[0, proto_ypix[where_f]-1, proto_xpix[where_f]-1 ]
        c2_evals = ext_data[0, proto_ypix[where_c2]-1, proto_xpix[where_c2]-1 ]
        c3_evals = ext_data[0, proto_ypix[where_c3]-1, proto_xpix[where_c3]-1 ]

        c01_evals = c01_evals.flatten()
        f_evals = f_evals.flatten()
        c2_evals = c2_evals.flatten()
        c3_evals = c3_evals.flatten()
    ##fi

    # Make plots of cumulative mass as a function of extinction.

    # Convert all arrays to mass.
    # Create the pixel to square parsec conversion factor. Expect extinction
    # map resolution in arcminutes. Expect distance in parsecs.
    pcsqperpix = ( extinction_pixelres * distance / 206264.8 )**2

    ext_mvals = ext_evals * 183 * pcsqperpix               # Mass as in Mairs+16
    s2_mvals = s2_fvals * 0.00263 * ( img_pixelres ** 2 )  # Mass as in Mairs+16
    proto_mvals = np.ones_like( proto_evals ) * 0.5        # 0.5 Msol for each

    # Sort the arrays by extinction value.
    ext_argsort = np.argsort(ext_evals)
    s2_ext_argsort = np.argsort(s2_evals)
    proto_ext_argsort = np.argsort(proto_evals)

    ext_evals = ext_evals[ext_argsort][::-1]
    ext_mvals = ext_mvals[ext_argsort][::-1]
    s2_evals = s2_evals[s2_ext_argsort][::-1]
    s2_mvals = s2_mvals[s2_ext_argsort][::-1]
    proto_evals = proto_evals[proto_ext_argsort][::-1]
    proto_mvals = proto_mvals[proto_ext_argsort][::-1]

    # Ensure everything is positive, only S2 should be able to be negative,
    # and extinction could maybe be negative depending on calibrator:
    ext_evals = ext_evals[ ext_evals > 0 ]
    ext_mvals = ext_mvals[ ext_mvals > 0 ]
    s2_evals = s2_evals[ s2_mvals > 0 ]
    s2_mvals = s2_mvals[ s2_mvals > 0 ]


    # Define a non-failure of a cumulative sum.
    def jlcml(a):
        out = np.zeros(len(a))
        csum = 0
        for i,val in enumerate(a):
            csum += val
            out[i] = csum
        return out
    #def

    ext_cml = jlcml(ext_mvals)
    s2_cml = jlcml(s2_mvals)
    proto_cml = jlcml(proto_mvals)

    # Find the total mass.
    ext_norm_mass = np.sum(ext_mvals)
    s2_norm_mass = np.sum(s2_mvals)
    proto_norm_mass = np.sum(proto_mvals)

    # Find median column.
    ext_med_column = np.median(ext_evals) * 16.6
    s2_med_column = np.median(s2_evals) * 16.6
    proto_med_column = np.median(proto_evals) * 16.6

    # Normalize each.
    ext_cml_norm = ext_cml / ext_norm_mass
    s2_cml_norm = s2_cml / s2_norm_mass
    proto_cml_norm = proto_cml / proto_norm_mass

    # Make the whole extinction map
    where_whole_data = np.where( ext_data.flatten().astype(float) > 0 )
    whole_ext_data = ext_data.flatten()[ where_whole_data ]
    whole_evals = np.sort( whole_ext_data )[::-1]
    whole_ext_cml = np.cumsum( whole_evals )
    whole_ext_cml = jlcml( whole_evals )
    whole_ext_norm_mass = np.sum(whole_evals)
    whole_cml_norm = whole_ext_cml / whole_ext_norm_mass
    whole_med_column = np.median(whole_evals) * 16.6
    whole_ext_norm_mass *= (183 * pcsqperpix)

    ext_stats = open('Catalogs/'+distance_dir+'/stats_ext.txt', 'w')
    ext_stats.write( 'Total Mass: '+str(ext_norm_mass)+'\n' )
    ext_stats.write( 'Total Median Column [E21]: '+str(ext_med_column)+'\n' )
    ext_stats.write( 'Total Npix: '+str(len(ext_evals))+'\n\n' )
    ext_stats.write( 'Island Mass: '+str(s2_norm_mass)+'\n' )
    ext_stats.write( 'Island Median Column [E21]: '+str(s2_med_column)+'\n' )
    ext_stats.write( 'Island Npix: '+str(len(s2_evals))+'\n\n' )
    ext_stats.write( 'Protostar Mass: '+str(proto_norm_mass)+'\n' )
    ext_stats.write( 'Protostar Median Column [E21]: '+str(proto_med_column)+'\n' )
    ext_stats.write( 'Protostar Npix: '+str(len(proto_evals))+'\n\n' )
    ext_stats.write( 'No Mask Mass: '+str(whole_ext_norm_mass)+'\n' )
    ext_stats.write( 'No Mask Median Column [E21]: '+str(whole_med_column)+'\n' )
    ext_stats.write( 'No Mask Npix: '+str(len(whole_evals)) )
    ext_stats.write( 'Individual Protostars:\n\n' )

    if orion == True:
        p_mvals = np.ones_like( p_evals ) * 0.5
        d_mvals = np.ones_like( d_evals ) * 0.5
        fp_mvals = np.ones_like( fp_evals ) * 0.5
        rp_mvals = np.ones_like( rp_evals ) * 0.5

        p_ext_argsort = np.argsort(p_evals)
        d_ext_argsort = np.argsort(d_evals)
        fp_ext_argsort = np.argsort(fp_evals)
        rp_ext_argsort = np.argsort(rp_evals)

        p_evals = p_evals[p_ext_argsort][::-1]
        d_evals = d_evals[d_ext_argsort][::-1]
        fp_evals = fp_evals[fp_ext_argsort][::-1]
        rp_evals = rp_evals[rp_ext_argsort][::-1]

        # Take the absolute value to account for noise
        p_evals = np.absolute(p_evals)
        d_evals = np.absolute(d_evals)
        fp_evals = np.absolute(fp_evals)
        rp_evals = np.absolute(rp_evals)

        # Normalized cumulative.
        p_norm_mass = np.sum(p_mvals)
        d_norm_mass = np.sum(d_mvals)
        fp_norm_mass = np.sum(fp_mvals)
        rp_norm_mass = np.sum(rp_mvals)

        p_cml_norm = jlcml(p_mvals) / p_norm_mass
        d_cml_norm = jlcml(d_mvals) / d_norm_mass
        fp_cml_norm = jlcml(fp_mvals) / fp_norm_mass
        rp_cml_norm = jlcml(rp_mvals) / rp_norm_mass

        # Find median column.
        p_med_column = np.median(p_mvals) * 16.6
        d_med_column = np.median(d_mvals) * 16.6
        fp_med_column = np.median(fp_mvals) * 16.6
        rp_med_column = np.median(rp_mvals) * 16.6

        ext_stats.write( 'Protostar Mass: '+str(p_norm_mass)+'\n' )
        ext_stats.write( 'Protostar Median Column [E21]: '+str(p_med_column)+'\n' )
        ext_stats.write( 'Protostar Npix: '+str(len(p_evals)) )
        ext_stats.write( 'Disk Mass: '+str(d_norm_mass)+'\n' )
        ext_stats.write( 'Disk Median Column [E21]: '+str(d_med_column)+'\n' )
        ext_stats.write( 'Disk Npix: '+str(len(d_evals)) )
        ext_stats.write( 'Faint Protostar Mass: '+str(fp_norm_mass)+'\n' )
        ext_stats.write( 'Faint Protostar Median Column [E21]: '+str(fp_med_column)+'\n' )
        ext_stats.write( 'Faint Protostar Npix: '+str(len(fp_evals)) )
        ext_stats.write( 'Red Protostar Mass: '+str(rp_norm_mass)+'\n' )
        ext_stats.write( 'Red Protostar Median Column [E21]: '+str(rp_med_column)+'\n' )
        ext_stats.write( 'Red Protostar Npix: '+str(len(rp_evals)) )

    ##fi
    if orion == False:
        c01_mvals = np.ones_like( c01_evals ) * 0.5
        f_mvals = np.ones_like( f_evals ) * 0.5
        c2_mvals = np.ones_like( c2_evals ) * 0.5
        c3_mvals = np.ones_like( c3_evals ) * 0.5

        c01_ext_argsort = np.argsort(c01_evals)
        f_ext_argsort = np.argsort(f_evals)
        c2_ext_argsort = np.argsort(c2_evals)
        c3_ext_argsort = np.argsort(c3_evals)

        c01_evals = c01_evals[c01_ext_argsort][::-1]
        f_evals = f_evals[f_ext_argsort][::-1]
        c2_evals = c2_evals[c2_ext_argsort][::-1]
        c3_evals = c3_evals[c3_ext_argsort][::-1]

        # Take the absolute value to account for noise
        c01_evals = np.absolute(c01_evals)
        f_evals = np.absolute(f_evals)
        c2_evals = np.absolute(c2_evals)
        c3_evals = np.absolute(c3_evals)

        # Normalized cumulative.
        c01_norm_mass = np.sum(c01_mvals)
        f_norm_mass = np.sum(f_mvals)
        c2_norm_mass = np.sum(c2_mvals)
        c3_norm_mass = np.sum(c3_mvals)

        c01_cml_norm = jlcml(c01_mvals) / c01_norm_mass
        f_cml_norm = jlcml(f_mvals) / f_norm_mass
        c2_cml_norm = jlcml(c2_mvals) / c2_norm_mass
        c3_cml_norm = jlcml(c3_mvals) / c3_norm_mass

        # Find median column.
        c01_med_column = np.median(c01_evals) * 16.6
        f_med_column = np.median(f_evals) * 16.6
        c2_med_column = np.median(c2_evals) * 16.6
        c3_med_column = np.median(c3_evals) * 16.6

        ext_stats.write( 'Class 0+1 Mass: '+str(c01_norm_mass)+'\n' )
        ext_stats.write( 'Class 0+1 Median Column [E21]: '+str(c01_med_column)+'\n' )
        ext_stats.write( 'Class 0+1 Npix: '+str(len(c01_evals)) )
        ext_stats.write( 'Flat Mass: '+str(f_norm_mass)+'\n' )
        ext_stats.write( 'Flat Median Column [E21]: '+str(f_med_column)+'\n' )
        ext_stats.write( 'Flat Npix: '+str(len(f_evals)) )
        ext_stats.write( 'Class 2 Mass: '+str(c2_norm_mass)+'\n' )
        ext_stats.write( 'Class 2 Median Column [E21]: '+str(c2_med_column)+'\n' )
        ext_stats.write( 'Class 2 Npix: '+str(len(c2_evals)) )
        ext_stats.write( 'Class 3 Mass: '+str(c3_norm_mass)+'\n' )
        ext_stats.write( 'Class 3 Median Column [E21]: '+str(c3_med_column)+'\n' )
        ext_stats.write( 'Class 3 Npix: '+str(len(c3_evals)) )
    ##fi

    ext_stats.close()

    if extinction_test == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(ext_evals*16.6, ext_cml_norm, c='b', label='NICEST')
        ax.plot(s2_evals*16.6, s2_cml_norm, c='r', label='Islands')
        ax.plot(proto_evals*16.6, proto_cml_norm, c='m', label='Protostars')
        ax.plot(whole_evals*16.6, whole_cml_norm, '--b', label='NICEST No Mask')

        if orion == True:
            ax.plot(p_evals*16.6, p_cml_norm, c='g', label='P')
            ax.plot(d_evals*16.6, d_cml_norm, c='Olive', label='D')
            #ax.plot(fp_evals*16.6, fp_cml_norm, c='g', label='FP')
            #ax.plot(rp_evals*16.6, rp_cml_norm, '--b', label='RP')
        ##fi
        if orion == False:
            ax.plot(c01_evals*16.6, c01_cml_norm, c='g', label='Class 0+1')
            ax.plot(f_evals*16.6, f_cml_norm, c='Olive', label='Flat')
            #ax.plot(c2_evals*16.6, c2_cml_norm, c='g', label='Class 2')
            #ax.plot(c3_evals*16.6, c3_cml_norm, '--b', label='Class 3')
        ##fi

        ax.legend()
        ax.set_ylim(0,1.1)
        ax.set_xlim(0,100)
        ax.set_xlabel(r'Column [$10^{21}$ cm$^{-2}$]')
        ax.set_ylabel('Fractional Cumulative Mass')
        plt.savefig('Test/'+distance_dir+'test_cumulative.png', format='png', dpi=300)
        plt.clf()
    ##fi

    ############################################################################

    #####################
    # EXTINCTION CONTOURS
    #####################

    ext_img = aplpy.FITSFigure(img_name)
    ext_img.show_colorscale(cmap='OrRd')
    ext_img.show_contour(   extinction_map,
                            levels=[1/1.66,2.5/1.66,5/1.66],
                            colors=['SkyBlue','DodgerBlue','Blue'],
                            linewidths=2)
    ext_img.tick_labels.set_style('colons')
    ext_img.tick_labels.set_xformat('hh:mm:ss')
    ext_img.tick_labels.set_yformat('dd:mm:ss')
    ext_img.ticks.set_color('Black')
    ext_img.add_colorbar()
    ext_img.colorbar.set_axis_label_text('Flux (mJy/sq.arcsecond)')

    # Create legend handles.
    ext_contour_10 = mlines.Line2D([], [], color='SkyBlue',
                                    label=r'$10 \times 10^{22}$ cm$^{-2}$')
    ext_contour_25 = mlines.Line2D([], [], color='DodgerBlue',
                                    label=r'$25 \times 10^{22}$ cm$^{-2}$')
    ext_contour_50 = mlines.Line2D([], [], color='Blue',
                                    label=r'$50 \times 10^{22}$ cm$^{-2}$')
    if np.nanmax(ext_data) > 1/1.66:
        ext_handles = [ext_contour_10,]
    if np.nanmax(ext_data) > 2.5/1.66:
        ext_handles.append(ext_contour_25)
    if np.nanmax(ext_data) > 5/1.66:
        ext_handles.append(ext_contour_50)
    ##fi
    # Create the legend.
    if np.nanmax(ext_data) > 1/1.66:
        ext_img._ax1.legend(handles=ext_handles, loc='best')
    ##fi

    # Save the image.
    ext_img.save('Plots/'+distance_dir+'/'+region+'ext_img.pdf', format='pdf', dpi=300)

    ############################################################################

    #################
    # OUTPUT CATALOGS
    #################

    # Initialize the Island property table
    island_data_names = ('ID',             # ID number
                         'RA_PEAK',        # RA of peak
                         'DEC_PEAK',       # Dec of peak
                         'VOLUME',         # Area of island [Sq.As.]
                         'MASS',           # Mass in [Msol]
                         'MASS_JEANS',     # Jeans mass in [Msol]
                         'MASS_JEANS_DENS',# Alternate Jeans mass [Msol]
                         'MASSPMJ',        # Mass in Jeans masses
                         'RADIUS',         # Circular size [pc]
                         'RADIUS_JEANS',   # Jeans radius [pc]
                         'RADIUSPRJ',      # Circular size in Jeans radii
                         'CONCENTRATION',  # Concentration
                         'NUMVOLDENSPEAK', # Peak volumetric number density
                         'MONOLITH',       # Monolithic (one fragment?)
                         'CLASS',          # Class of closest YSO
                         'CLOSE')          # Closest YSO within 1 beam?
    island_data_dtypes = ('int','float','float','float','float','float','float',
                          'float','float','float','float','float','float',
                          'int','U6','int')

    island_data = table.Table( [island_pident,
                                island_rapeak,
                                island_decpeak,
                                island_volume,
                                island_mass,
                                island_massj,
                                island_massjdens,
                                island_masspmj,
                                island_radius,
                                island_radiusj,
                                island_radiusprj,
                                island_conc,
                                island_numvoldenspeak,
                                island_monolith,
                                island_class,
                                island_close],
                                names = island_data_names,
                                dtype = island_data_dtypes)

    # Initialize the Fragment property table
    fragment_data_names = ('RA_PEAK',           # RA of fragment peak
                           'DEC_PEAK',          # Dec of fragment peak
                           'PEAK_FLUX',         # Peak flux of fragment
                           'SUM_FLUX',          # Total flux of fragment
                           'VOLUME',            # Area of fragment [Sq.As.]
                           'MASS',              # Mass of fragment [Msol]
                           'MASS_JEANS',        # Jeans mass of fragment [Msol]
                           'MASS_JEANS_DENS',   # Alternate Jeans mass [Msol]
                           'MASSPMJ',           # Mass in Jeans masses
                           'RADIUS',            # Circular size [pc]
                           'RADIUS_JEANS',      # Jeans radius [pc]
                           'RADIUSPRJ',         # Circular size in Jeans radii
                           'CONCENTRATION',     # Concentration
                           'NUMVOLDENSPEAK',    # Peak volumetric number density
                           'MONOLITH',          # Monolithic (only fragment?)
                           'CLASS',             # Class of nearest YSO
                           'CLOSE')             # Nearest YSO within 1 beam?
    fragment_data_dtypes = ('float','float','float','float','float','float','float',
                            'float','float','float','float','float','float','float',
                            'int','U6','int')
    fragment_data = table.Table( [frag_rapeak,
                                  frag_decpeak,
                                  frag_peakflux,
                                  frag_sumflux,
                                  frag_volume,
                                  frag_mass,
                                  frag_massj,
                                  frag_massjdens,
                                  frag_masspmj,
                                  frag_radius,
                                  frag_radiusj,
                                  frag_radiusprj,
                                  frag_conc,
                                  frag_numvoldenspeak,
                                  frag_monolith,
                                  frag_class,
                                  frag_close],
                                  names = fragment_data_names,
                                  dtype = fragment_data_dtypes)
    # Initialize the protostar property table.
    protostar_data_names = ('CLASS',        # YSO class
                            'DISTANCE',     # Distance to nearest fragment peak
                            'FLUX')         # Flux density at protostar position
    protostar_data_dtypes = ('U6','float','float')
    protostar_data = table.Table( [proto_class,
                                   proto_fragdist,
                                   proto_flux],
                                   names = protostar_data_names,
                                   dtype = protostar_data_dtypes)

    # Initialize the cloud property table.
    if orion == True:
        cloud_data_names = ('FIELD EXT MASS',
                            'TOTAL EXT MASS',
                            'SCUBA2 MASS',
                            'YSO MASS',
                            'P MASS',
                            'D MASS',
                            'N ISLAND',
                            'ISLAND N MONOLITH',
                            'ISLAND N COMPLEX',
                            'ISLAND N YSO',
                            'N FRAGMENT',
                            'FRAGMENT N MONOLITH',
                            'FRAGMENT N COMPLEX',
                            'FRAGMENT N YSO',
                            )
        cloud_data_dtypes = (   'float','float','float','float','float','float',
                                'int','int','int','int','int','int','int','int')
        cloud_data_cols = [ ext_norm_mass,
                            whole_ext_norm_mass,
                            s2_norm_mass,
                            proto_norm_mass,
                            p_norm_mass,
                            d_norm_mass,
                            n_islands,
                            len(np.where(island_monolith == 1)[0]),
                            n_islands-len(np.where(island_monolith == 1)[0]),
                            n_islands-len(np.where(island_class == 'None')[0]),
                            n_frags,
                            len(np.where(frag_monolith == 1)[0]),
                            n_frags-len(np.where(frag_monolith == 1)[0]),
                            n_frags-len(np.where(frag_class == 'None')[0])
                            ]
        # Convert to numpy arrays, arghh!
        for i,val in enumerate(cloud_data_cols):
            cloud_data_cols[i] = np.array([val])
        #for

        cloud_data = table.table(   cloud_data_cols,
                                    names = cloud_data_names,
                                    dtype = cloud_data_dtypes)
    ##fi
    if orion == False:
        cloud_data_names = ('FIELD EXT MASS',
                            'TOTAL EXT MASS',
                            'SCUBA2 MASS',
                            'YSO MASS',
                            'C01 MASS',
                            'F MASS',
                            'C2 MASS',
                            'N ISLAND',
                            'ISLAND N MONOLITH',
                            'ISLAND N COMPLEX',
                            'ISLAND N YSO',
                            'N FRAGMENT',
                            'FRAGMENT N MONOLITH',
                            'FRAGMENT N COMPLEX',
                            'FRAGMENT N YSO',
                            )
        cloud_data_dtypes = (   'float','float','float','float','float','float',
                                'float','int','int','int','int','int','int',
                                'int','int')
        cloud_data_cols =  [ext_norm_mass,
                            whole_ext_norm_mass,
                            s2_norm_mass,
                            proto_norm_mass,
                            c01_norm_mass,
                            f_norm_mass,
                            c2_norm_mass,
                            n_islands,
                            len(np.where(island_monolith == 1)[0]),
                            n_islands-len(np.where(island_monolith == 1)[0]),
                            n_islands-len(np.where(island_class == 'None')[0]),
                            n_frags,
                            len(np.where(frag_monolith == 1)[0]),
                            n_frags-len(np.where(frag_monolith == 1)[0]),
                            n_frags-len(np.where(frag_class == 'None')[0])
                            ]
        # Convert to numpy arrays, arghh!
        for i,val in enumerate(cloud_data_cols):
            cloud_data_cols[i] = np.array([val])
        #for
        cloud_data = table.Table(   cloud_data_cols,
                                    names = cloud_data_names,
                                    dtype = cloud_data_dtypes)
    ##fi

    # Initialize the Extinction information binary:
    all_evals = [ext_evals, s2_evals, proto_evals, whole_evals]
    all_cml_norm = [ext_cml_norm, s2_cml_norm, proto_cml_norm, whole_cml_norm]
    if orion == True:
        yso_evals = [p_evals, d_evals, fp_evals, rp_evals]
        yso_cml_norm = [p_cml_norm, d_cml_norm, fp_cml_norm, rp_cml_norm]
    ##fi
    if orion == False:
        yso_evals = [c01_evals, f_evals, c2_evals, c3_evals]
        yso_cml_norm = [c01_cml_norm, f_cml_norm, c2_cml_norm, c3_cml_norm]
    ##fi

    # Write the data to file.
    island_data.write('Catalogs/'+distance_dir+'/'+region+'_island_properties.FIT',
                        format='fits', overwrite=True)
    fragment_data.write('Catalogs/'+distance_dir+'/'+region+'_fragment_properties.FIT',
                        format='fits', overwrite=True)
    protostar_data.write('Catalogs/'+distance_dir+'/'+region+'_protostar_properties.FIT',
                        format='fits', overwrite=True)
    cloud_data.write('Catalogs/'+distance_dir+'/'+region+'_cloud_propertes.FIT',
                        format='fits', overwrite=True)

    # Save extinction information in a numpy binary, then zip it.
    np.save('Catalogs/'+distance_dir+'/'+region+'_extinction.npy',
            np.array([all_evals, all_cml_norm, yso_evals, yso_cml_norm]))
    #os.system('gzip Catalogs/+'region+'_extinction.npy')

    # Write mass for each type of protostar and region.
    if orion == True:
        ext_mass_out = np.array([   ext_norm_mass,
                                    whole_ext_norm_mass,
                                    s2_norm_mass,
                                    proto_norm_mass,
                                    p_norm_mass,
                                    d_norm_mass])
        np.savetxt('Catalogs/'+distance_dir+'/'+region+'_extinction_mass.txt', ext_mass_out,
        header='mass: nicest gbs footprint, whole nicest map, s2 islands, yso, p, d')
    ##fi
    if orion == False:
        ext_mass_out = np.array([   ext_norm_mass,
                                    whole_ext_norm_mass,
                                    s2_norm_mass,
                                    proto_norm_mass,
                                    c01_norm_mass,
                                    f_norm_mass,
                                    c2_norm_mass])
        np.savetxt('Catalogs/'+distance_dir+'/'+region+'_extinction_mass.txt', ext_mass_out,
        header='mass: nicest gbs footprint, whole nicest map, s2 islands, yso, c0+1, flat, c2')
    ##fi

    ############################################################################
#def
