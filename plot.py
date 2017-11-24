# ----------------------------------------------------------------------------
#
# TITLE - plot.py
# AUTHOR - James Lane
# PROJECT - Archipelago
# CONTENTS:
#	1. IslandMass
#   2. IslandJeansMass
#   3. FragmentMass
#   4. FragmentJeansMass
#   5. PeakBrightnessProtostars
#   6. PeakDistanceProtostars
#   7. FragmentConcentration
#   8. IslandDensityRadius
#   9. FragmentDensityRadius
#   10. MonolithicJeansRadius
#   11. ComplexJeansRadius
#   12. ColumnMassFraction
#   13. PlottingProperties
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Plotting routines for the Archipelago analysis.'''
__author__ = "James Lane"

#Imports
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import patches as mpatches
from matplotlib.ticker import MaxNLocator
from astropy.io import fits
from astropy import units as u
from scipy.stats import norm
from scipy.optimize import curve_fit
import pdb
import os
import archipelago.tasks

################################################################################



def IslandMass( fig,
                island_data,
                outplots,
                region_properties,
                plot_individual = False,
              ):
    '''
    Island Mass:

    Plot a histogram of the logarithmic mass of the islands from the
    archipelago analysis.

    Args:
        fig (figure): matplotlib figure object.
        island_data (Table): Astropy table of island data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_individual (Boolean): Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data:
    beamsize = region_properties.beamsize
    region = region_properties.region
    noise = region_properties.noise
    distance = region_properties.distance
    island_mass = island_data['MASS'].data
    island_radius = island_data['RADIUS'].data

    # Power law tail fiducial
    cmf_m = np.array([10,1000])
    cmf_0_35 = 100*np.power(cmf_m,-0.35)

    # Detection limit. Island the size of beamsize squared
    island_mass_min = 1.299E-8 * (3*noise) * (distance**2) * (beamsize**2)

    # Incompleteness identification, 3-sigma clump at median clump size.
    med_island_radius = np.median(island_radius) * 206264.8 / distance
    island_mass_med = 1.299E-8 * (3*noise) * (distance**2) * (np.pi*med_island_radius**2)

    # Warn if large islands of interest.
    if max(island_mass) > 1000: print('Warning: hypermassive islands in '+region)

    # Make the histogram.
    ax = fig.add_subplot(111)
    num, bins, rects = ax.hist(np.log10(island_mass), bins=12,
                                range=[-3,3], color='b', alpha=0.75,
                                log=True)

    # Plot supporting information.
    ax.axvline(0, color='Black', linestyle='dashed')
    ax.axvline( np.log10(island_mass_min), color='Black', linestyle='dotted',
                linewidth=0.5, label='Min. detection limit')
    ax.axvline( np.log10(island_mass_med), color='Black', linestyle='dashdot',
                linewidth=0.5, label='Med. detection limit')
    ax.plot(np.log10(cmf_m), cmf_0_35, color='DodgerBlue', label=r'$\alpha=-0.35$')
    scale_relation = r'$M \propto d^{2}$'
    ax.annotate(scale_relation, (0.55,0.9), xycoords='axes fraction')

    # Adjust axes
    ax.set_ylabel('Number')
    ax.set_xlabel(r'Island Mass log(M$_{\odot}$)')
    ax.set_ylim(0.5, 200)
    ax.legend(loc='upper right', fontsize=8, labelspacing=0.5)
    ax.set_title(region+' Island Masses')

    # Save
    plt.savefig(outplots, format='pdf')

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('islandmass_'+region+'.png', dpi=300)
    ##fi

    plt.clf()
#def

def IslandJeansMass(fig,
                    island_data,
                    outplots,
                    region_properties,
                    plot_properties,
                    plot_individual = False,
                    ):
    '''
    Island Jeans Mass:

    Plot a histogram of the logarithmic mass, in units of Jeans mass of the
    islands from the archipelago analysis.

    Args:
        fig (figure): matplotlib figure object.
        island_data (Table): Astropy table of island data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        plot_individual (Boolean): Make an individual png of the image [False]

    Returns:
        None?
    '''

    # Set region properties and unpack data:
    beamsize = region_properties.beamsize
    region = region_properties.region
    distance = region_properties.distance
    noise = region_properties.noise
    island_masspmj = island_data['MASSPMJ'].data
    where_islandnyso = plot_properties.where_islandnyso
    islandyso = np.ones_like(island_masspmj)
    islandyso[where_islandnyso] = 0
    where_islandyso = np.where(islandyso == 1)[0]

    # Detection limit. Island the size of beamsize squared
    island_mass_min = 1.299E-8 * (3*noise) * (distance**2) * (beamsize**2)
    island_radius_min = (1/365594.8) * distance * np.sqrt( beamsize**2 )
    island_density_min = (1.616E-23) * np.divide(island_mass_min,
                                                np.power(island_radius_min,3))
    island_jeans_mass_min = (3.305E-8) * np.power(island_density_min*1000,-0.5)
    island_masspmj_min = island_mass_min / island_jeans_mass_min

    # Warn if islands of interest are present.
    if max(np.log10(island_masspmj)) > 3: print('Warning: hyperunstable islands in '+region)

    # Plot different regimes.
    ax = fig.add_subplot(111)
    # Stable:
    bin_range = [-3,0]
    snum, sbins, rects = ax.hist(np.log10(island_masspmj), bins=12, range=bin_range,
                                color='MediumBlue', alpha=1.0, log=True)
    _,_,_ = ax.hist(np.log10(island_masspmj[where_islandyso]), bins=12, range=bin_range,
                    color='ForestGreen', alpha=1.0, log=True)
    # Unstable:
    bin_range = [0,np.log10(4)]
    unum, ubins, rects = ax.hist(np.log10(island_masspmj), bins=3, range=bin_range,
                                color='MediumBlue', alpha=0.75, log=True)
    _,_,_ = ax.hist(np.log10(island_masspmj[where_islandyso]), bins=3, range=bin_range,
                    color='ForestGreen', alpha=1.0, log=True)
    # Significantly Unstable:
    bin_range = [np.log10(4),3]
    sunum, subins, rects = ax.hist(np.log10(island_masspmj), bins=9, range=bin_range,
                                color='MediumBlue', alpha=0.5, log=True)
    _,_,_ = ax.hist(np.log10(island_masspmj[where_islandyso]), bins=9, range=bin_range,
                    color='ForestGreen', alpha=1.0, log=True)

    # Plot supporting information.
    ax.axvline(0, 0, 1, color='k', linestyle='dashed')
    ax.axvline(np.log10(4), 0, 1, color='k', linestyle='dashed')
    ax.axvline(np.log10(island_masspmj_min), color='Black', linestyle='dotted',
                linewidth=0.5)
    scale_relation = r'$M/M_{Jeans} \propto d^{\frac{3}{2}}$'
    ax.annotate(scale_relation, (0.8,0.7), xycoords='axes fraction')

    # Adjust axes
    ax.set_ylabel('Number')
    ax.set_xlabel(r'log(M$_{island}$ / M$_{island,J}$)')
    ax.set_xlim(-2.5, 3)
    ax.set_ylim(0.5, 200)
    ax.set_title(region+' Island Stability')

    plt.savefig(outplots, format='pdf')

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('islandjeansmass_'+region+'.png', dpi=300)
    ##fi

    plt.clf()
#def

def FragmentMass(fig,
                 frag_data,
                 outplots,
                 region_properties,
                 plot_individual = False,
                 ):
    '''
    Fragment Mass:

    Plot a histogram of the logarithmic mass of the islands from the
    archipelago analysis.

    Args:
        fig (figure): matplotlib figure object.
        fragment_data (Table): Astropy table of fragmen data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_individual (Boolean): Make an individual png of the image [False]

    Returns:
        None?
    '''

    # Set region properties and unpack the data.
    beamsize = region_properties.beamsize
    region = region_properties.region
    distance = region_properties.distance
    img_pixelres = region_properties.pixel_resolution
    noise = region_properties.noise
    frag_mass = frag_data['MASS'].data
    frag_radius = frag_data['RADIUS'].data

    # Detection limit. A 3-sigma clump with size 7 pixels.
    frag_mass_min = 1.299E-8 * (3*noise) * (distance**2) * (7*img_pixelres**2)

    # Incompleteness identification, 3-sigma clump at median clump size.
    med_frag_radius = np.median(frag_radius) * 206264.8 / distance
    frag_mass_med = 1.299E-8 * (3*noise) * (distance**2) * (np.pi*med_frag_radius**2)

    # Power law tail fiducial
    cmf_m = np.array([10,1000])
    cmf_0_35 = 100*np.power(cmf_m,-0.35)

    # Warn about objects of interest
    if max(frag_mass) > 1000: print('Warning: hypermassive fragments in '+region)

    # Plot
    ax = fig.add_subplot(111)
    num, bins, rects = ax.hist(np.log10(frag_mass), bins=12,
                                range=[-3,3], color='b', alpha=0.75,
                                log=True)

    # Plot supporting information.
    ax.axvline(0, color='Black', linestyle='dashed')
    ax.axvline(np.log10(frag_mass_min), color='Black', linestyle='dotted',
                linewidth=0.5, label='Min. detection limit')
    ax.axvline(np.log10(frag_mass_med), color='Black', linestyle='dashdot',
                linewidth=0.5, label='Med. detection limit')
    ax.plot(np.log10(cmf_m), cmf_0_35, color='DodgerBlue', label=r'$\alpha=-0.35$')
    scale_relation = r'$M \propto d^{2}$'
    ax.annotate(scale_relation, (0.55,0.9), xycoords='axes fraction')

    # Adjust axes
    ax.set_ylabel('Number')
    ax.set_xlabel(r'Fragment Mass log(M$_{\odot}$)')
    ax.set_ylim(0.5, 200)
    ax.set_title(region+' Fragment Masses')
    ax.legend(loc='upper right', fontsize=8, labelspacing=0.5)

    # Save plot
    plt.savefig(outplots)

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('fragmentmass_'+region+'.png', dpi=300)

    plt.clf()
#def

def FragmentJeansMass(fig,
                      frag_data,
                      outplots,
                      region_properties,
                      plot_properties,
                      plot_individual = False,
                      ):
    '''
    Fragment Jeans Mass:

    Plot a histogram of the logarithmic mass, in units of Jeans mass of the
    fragments from the archipelago analysis.

    Args:
        fig (figure): matplotlib figure object.
        frag_data (Table): Astropy table of island data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_individual (Boolean):Make an individual png of the image [False]

    Returns:
        None?
    '''

    # Set region properties and unpack data:
    beamsize = region_properties.beamsize
    region = region_properties.region
    noise = region_properties.noise
    distance = region_properties.distance
    img_pixelres = region_properties.pixel_resolution
    frag_masspmj = frag_data['MASSPMJ'].data
    where_fragnyso = plot_properties.where_fragnyso
    fragyso = np.ones_like(frag_masspmj)
    fragyso[where_fragnyso] = 0
    where_fragyso = np.where(fragyso == 1)[0]

    # Detection limit. Fragment the size of 7 pixels
    frag_mass_min = 1.299E-8 * (3*noise) * (distance**2) * (7*img_pixelres**2)
    frag_radius_min = (1/365594.8) * distance * np.sqrt( 7*img_pixelres**2 )
    frag_density_min = (1.616E-23) * np.divide(frag_mass_min,
                                                np.power(frag_radius_min,3))
    frag_jeans_mass_min = (3.305E-8) * np.power(frag_density_min*1000,-0.5)
    frag_masspmj_min = frag_mass_min / frag_jeans_mass_min

    # Warn about objects of interest.
    if max(np.log10(frag_masspmj)) > 3: print('Warning: hyperunstable fragments in '+region)

    # Plot
    ax = fig.add_subplot(111)
    # Stable
    bin_range = [-3,0]
    snum, sbins, rects = ax.hist(np.log10(frag_masspmj), bins=12, range=bin_range,
                                color='MediumBlue', alpha=1.0, log=True)
    _,_,_ = ax.hist(np.log10(frag_masspmj[where_fragyso]), bins=12, range=bin_range,
                    color='ForestGreen', alpha=1.0, log=True)
    # Unstable
    bin_range = [0,np.log10(4)]
    unum, ubins, rects = ax.hist(np.log10(frag_masspmj), bins=3, range=bin_range,
                                color='MediumBlue', alpha=0.75, log=True)
    _,_,_ = ax.hist(np.log10(frag_masspmj[where_fragyso]), bins=3, range=bin_range,
                    color='ForestGreen', alpha=1.0, log=True)
    # Significantly Unstable
    bin_range = [np.log10(4),3]
    sunum, subins, rects = ax.hist(np.log10(frag_masspmj), bins=6, range=bin_range,
                                color='MediumBlue', alpha=0.5, log=True)
    _,_,_ = ax.hist(np.log10(frag_masspmj[where_fragyso]), bins=6, range=bin_range,
                    color='ForestGreen', alpha=1.0, log=True)

    # Plot supporting information
    ax.axvline(0, 0, 1, color='k', linestyle='dashed')
    ax.axvline(np.log10(4), 0, 1, color='k', linestyle='dashed')
    ax.axvline(np.log10(frag_masspmj_min), color='Black', linestyle='dotted',
                linewidth=0.5)
    scale_relation = r'$M/M_{Jeans} \propto d^{\frac{3}{2}}$'
    ax.annotate(scale_relation, (0.8,0.7), xycoords='axes fraction')

    ax.set_ylabel('Number')
    ax.set_xlabel(r'log(M$_{frag}$ / M$_{frag,J}$)')
    ax.set_xlim(-3, 3)
    ax.set_ylim(0.5, 200)
    ax.set_title(region+' Fragment Stability')

    # Save plot
    plt.savefig(outplots, format='pdf')

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('islandjeansmass_'+region+'.png', dpi=300)
    ##fi

    plt.clf()
#def

def PeakBrightnessProtostars(fig,
                            proto_data,
                            outplots,
                            region_properties,
                            plot_properties,
                            plot_individual = False,
                            ax = None,
                            save_ind = True,
                            spec_fig = False
                            ):
    '''
    PeakBrightnessProtostars:

    Plot a histogram of the peak brightness for each protostars in the region.

    Args:
        fig (figure): matplotlib figure object.
        proto_data (Table): Astropy table of protostar data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Hardcode some properties of the plot.
    binsize = 0.028 # 0.03 Jy/beam, about 3rms in Orion.
    binrange = [0,1.0]
    binnum = 34

    # Set region properties and unpack data.
    orion = region_properties.orion
    region = region_properties.region
    proto_flux = proto_data['FLUX'].data
    if orion == True:
        where_proto = plot_properties.where_proto
        where_disk = plot_properties.where_disk
        where_rp = plot_properties.where_rp
        where_fp = plot_properties.where_fp
    ##fi
    if orion == False:
        where_class01 = plot_properties.where_class01
        where_flat = plot_properties.where_flat
        where_class2 = plot_properties.where_class2
        where_class3 = plot_properties.where_class3
    ##fi

    # Force all protostars into the frame, protostars on negative pixels get
    # 0.001 Jy/beam flux, brighter than 1.0 Jy/beam get 0.999 Jy/beam.
    proto_flux[np.where(proto_flux < 0)[0]] = 0.001
    proto_flux[np.where(proto_flux > 0.999)[0]] = 0.999

    # Make an axis object.
    if ax == None:
        ax = fig.add_subplot(111)
    ##fi

    # Do plotting and based on whether the region is in Orion.
    if orion == True:
        #Disks
        dnums, dbins, drects = ax.hist(proto_flux[where_disk],
                                        bins=binnum, range=binrange,
                                        color='Olive', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Protostars
        pnums, pbins, prects = ax.hist(proto_flux[where_proto],
                                        bins=binnum, range=binrange,
                                        color='Green', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Faint protostars
        fpnums, fpbins, fprects = ax.hist(proto_flux[where_fp],
                                        bins=binnum, range=binrange,
                                        color='Blue', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Red protostars
        rpnums, rpbins, rprects = ax.hist(proto_flux[where_rp],
                                        bins=binnum, range=binrange,
                                        color='Red', alpha=0.75,
                                        histtype='step', linewidth=2)
        # Annotate what the out-of-frame maximum is
        ax.annotate(str(int(np.amax(dnums))),(0.05,0.95),
                    xycoords='axes fraction')
    ##fi

    if orion == False:
        #Class 01
        c01nums, c01bins, c01rects = ax.hist(proto_flux[where_class01],
                                        bins=binnum, range=binrange,
                                        color='Green', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Class 2
        c2nums, c2bins, c2rects = ax.hist(proto_flux[where_class2],
                                        bins=binnum, range=binrange,
                                        color='Olive', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Class 3
        c3nums, c3bins, c3rects = ax.hist(proto_flux[where_class3],
                                        bins=binnum, range=binrange,
                                        color='Purple', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Class Flat
        fnums, fbins, frects = ax.hist(proto_flux[where_flat],
                                        bins=binnum, range=binrange,
                                        color='Orange', alpha=0.75,
                                        histtype='step', linewidth=2)
        # Annotate what the out-of-frame maximum is.
        ax.annotate(str(int(np.amax(c2nums))),(0.05,0.95),
                    xycoords='axes fraction')
    ##fi

    ax.set_xlim(0,1)
    ax.set_ylim(0,120)
    ax.set_ylabel('N')
    ax.set_xlabel(r'Pixel Brightness (Jy Beam$^{-1}$)')
    ax.set_title(region+' Protostar Fluxes')

    # Do the legend depending on whether the region is in Orion.
    if orion == True:
        green_patch = mpatches.Patch(fc = 'Green', ec = 'k', label='P')
        olive_patch = mpatches.Patch(fc= 'Olive', ec = 'k', label='D')
        red_patch = mpatches.Patch(fc = 'Red', ec = 'k', label='RP')
        blue_patch = mpatches.Patch(fc = 'Blue', ec = 'k', label='FP')
        legend_handles = np.array([green_patch, olive_patch,
                                    blue_patch, red_patch],dtype='O')
    ##fi

    if orion == False:
        green_patch = mpatches.Patch(fc = 'Green', ec = 'k', label='0+1')
        orange_patch = mpatches.Patch(fc= 'Orange', ec = 'k', label='Flat')
        olive_patch = mpatches.Patch(fc = 'Olive', ec = 'k', label='2')
        purple_patch = mpatches.Patch(fc = 'Purple', ec = 'k', label='3')
        legend_handles = np.array([green_patch, orange_patch,
                                    olive_patch, purple_patch],dtype='O')
    ##fi

    # Make the legend and save the figure.
    plt.legend(handles=list(legend_handles),loc='upper right')
    scale_relation = r'$f_{peak}\ constant$'
    ax.annotate(scale_relation, (0.75,0.6), xycoords='axes fraction')
    if save_ind == True:
      plt.savefig(outplots, format='pdf')
    ##fi

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('protobrightness_'+region+'.png', dpi=300)

    plt.clf()
#def

def PeakDistanceProtostars(fig,
                           proto_data,
                           outplots,
                           region_properties,
                           plot_properties,
                           zoom = False,
                           plot_individual = False,
                           ax = None,
                           save_ind = True,
                           spec_fig = False
                           ):
    '''
    PeakDistanceProtostar:

    Plot a histogram of the fragment peak distance for each protostars in
    the region.

    Args:
        fig (figure): matplotlib figure object.
        proto_data (Table): Astropy table of protostar data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        zoom (Boolean): Do a zoomed in version of this plot.
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data.
    orion = region_properties.orion
    distance = region_properties.distance
    beamsize = region_properties.beamsize
    region = region_properties.region
    proto_fragdist = proto_data['DISTANCE'].data
    if orion == True:
        where_proto = plot_properties.where_proto
        where_disk = plot_properties.where_disk
        where_rp = plot_properties.where_rp
        where_fp = plot_properties.where_fp
    ##fi
    if orion == False:
        where_class01 = plot_properties.where_class01
        where_flat = plot_properties.where_flat
        where_class2 = plot_properties.where_class2
        where_class3 = plot_properties.where_class3
    ##fi

    # Hardcode some properties of the plot.
    binsize = distance * beamsize / (206264.8) # beamsize in pc
    binrange = [0,1.996]
    binnum = 61
    if zoom == True:
        binrange = [0,0.49]
        binnum = int(beamsize)
    ##fi

    # Force all protostars into the frame, set all distances larger than
    # 1.995 pc to be 1.995 pc.
    proto_fragdist[np.where(proto_fragdist > 1.995)[0]] = 1.995

    # Make an axis object.
    if ax == None:
        ax = fig.add_subplot(111)
    ##fi

    # Plot depending on whether the region is in Orion.
    if orion == True:
        #Disks
        dnums, dbins, drects = ax.hist(proto_fragdist[where_disk],
                                        bins=binnum, range=binrange,
                                        color='Olive', alpha=0.75,
                                        histtype='step', linewidth=2,
                                        zorder=3)
        #Protostars
        pnums, pbins, prects = ax.hist(proto_fragdist[where_proto],
                                        bins=binnum, range=binrange,
                                        color='g', alpha=0.75, histtype='step',
                                        linewidth=2)
        #Faint protostars
        fpnums, fpbins, fprects = ax.hist(proto_fragdist[where_fp],
                                        bins=binnum, range=binrange,
                                        color='b', alpha=0.75, histtype='step',
                                        linewidth=2)
        #Red protostars
        rpnums, rpbins, rprects = ax.hist(proto_fragdist[where_rp],
                                        bins=binnum, range=binrange,
                                        color='r', alpha=0.75, histtype='step',
                                        linewidth=2)
    ##fi

    if orion == False:
        #Flat
        fnums, fbins, frects = ax.hist(proto_fragdist[where_flat],
                                        bins=binnum, range=binrange,
                                        color='Orange', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Class 01
        c01nums, c01bins, c01rects = ax.hist(proto_fragdist[where_class01],
                                        bins=binnum, range=binrange,
                                        color='Green', alpha=0.75,
                                        histtype='step', linewidth=2)
        #Class 2
        c2nums, c2bins, c2rects = ax.hist(proto_fragdist[where_class2],
                                        bins=binnum, range=binrange,
                                        color='Olive', alpha=0.75,
                                        histtype='step', linewidth=3, zorder=2)
        #Class 3
        c3nums, c3bins, c3rects = ax.hist(proto_fragdist[where_class3],
                                        bins=binnum, range=binrange,
                                        color='Purple', alpha=0.75,
                                        histtype='step', linewidth=2)
    ##fi

    # If zooming, fit a skew-normal distribution.
    if zoom == True:
        # Define the skew-normal distribution.
        def skewnorm(x,a,b,sigma,gamma):
            arg = ( x - b ) / sigma
            return a * norm.pdf(arg) * norm.cdf(gamma * arg)
        #def

        # If in Orion fit and plot the disk distribution.
        if orion == True:
            try:
                dbin_cents = (dbins[:-1] + dbins[1:]) / 2
                skew_guess = [ np.amax(dnums), np.median(dbin_cents), 1, 1 ]
                skew_popt,_ = curve_fit(skewnorm, dbin_cents, dnums,
                                        p0=skew_guess, ftol=0.01, maxfev=50000)
                skew_x = np.arange(0,2,0.005)
                skew_y = skewnorm(skew_x,*skew_popt)
                ax.plot(skew_x, skew_y, color='Black', linewidth=1)
            except RuntimeError:
                pass
        ##fi

        # If not in Orion fit and plot the class 2 distribution.
        if orion == False:
            try:
                c2bin_cents = (c2bins[:-1] + c2bins[1:]) / 2
                skew_guess = [ np.amax(c2nums), np.median(c2bin_cents), 1, 1 ]
                skew_popt,_ = curve_fit(skewnorm, c2bin_cents, c2nums,
                                        p0=skew_guess, ftol=0.01, maxfev=10000)
                skew_x = np.arange(0,2,0.005)
                skew_y = skewnorm(skew_x, *skew_popt)
                ax.plot(skew_x, skew_y, color='Black', linewidth=1)
            except RuntimeError:
                pass
        ##fi
        ax.annotate(r'y = A pdf($x\prime$)cdf($x\prime \gamma$)',
                    (0.6,0.9), xycoords='axes fraction')
        ax.annotate(r'A = '+str(round(skew_popt[0],3)), (0.6,0.85),
                    xycoords='axes fraction')
        ax.annotate(r'b = '+str(round(skew_popt[1],3)), (0.6,0.8),
                    xycoords='axes fraction')
        ax.annotate(r'$\sigma$ = '+str(round(skew_popt[2],3)), (0.6,0.75),
                    xycoords='axes fraction')
        ax.annotate(r'$\gamma$ = '+str(round(skew_popt[3],3)), (0.6,0.7),
                    xycoords='axes fraction')
    ##fi

    # Set the Y axis range.
    bin_ymax = 120
    if orion == True and zoom == True:
        bin_ymax = 1.1*np.amax([dnums,pnums,fpnums,rpnums])
    if orion == False and zoom == True:
        bin_ymax = 1.1*np.amax([fnums,c01nums,c2nums,c3nums])

    # Set the X axis range.
    bin_xmax = 2
    if zoom == True:
        bin_xmax = 0.5

    ax.axvline(distance * 14.6 / (206264.8), color='Magenta', linewidth=0.25)
    ax.set_xlim(0,bin_xmax)
    ax.set_ylim(0,bin_ymax)
    ax.set_ylabel('N')
    ax.set_xlabel('Distance (pc)')
    ax.set_title(region+' Protostar Distances')
    scale_relation = r'$D_{proto} \propto d^{1}$'
    ax.annotate(scale_relation, (0.8,0.7), xycoords='axes fraction')
    if save_ind == True:
      plt.savefig(outplots, format='pdf')
    ##fi

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('protodistance_'+region+'.png', dpi=300)

    plt.clf()
#def

def FragmentConcentration(fig,
                          frag_data,
                          outplots,
                          region_properties,
                          plot_properties,
                          plot_individual = False,
                          ax = None,
                          save_ind = True,
                          spec_fig = False
                          ):
    '''
    FragmentConcentration:

    Plot fragment concentration as a function of mass in units of Jeans mass.

    Args:
        fig (figure): matplotlib figure object.
        frag_data (Table): Astropy table of protostar data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data.
    region = region_properties.region
    orion = region_properties.orion
    noise = region_properties.noise
    distance = region_properties.distance
    img_pixelres = region_properties.pixel_resolution
    frag_masspmj = frag_data['MASSPMJ'].data
    frag_conc = frag_data['CONCENTRATION'].data
    where_fragmono = plot_properties.where_fragmono
    where_fragcomp = plot_properties.where_fragcomp
    frag_fcolor = plot_properties.frag_fcolor
    frag_ecolor = plot_properties.frag_ecolor
    where_fragnyso = plot_properties.where_fragnyso
    if orion == True:
        where_fragpc = plot_properties.where_fragpc
        where_fragpnc = plot_properties.where_fragpnc
        where_fragd = plot_properties.where_fragd
    ##fi
    if orion == False:
        where_frag01c = plot_properties.where_frag01c
        where_frag01nc = plot_properties.where_frag01nc
        where_fragf = plot_properties.where_fragf
        where_frag2 = plot_properties.where_frag2
    ##fi

    # Detection limit. Fragment the size of 7 pixels
    frag_mass_min = 1.299E-8 * (3*noise) * (distance**2) * (7*img_pixelres**2)
    frag_radius_min = (1/365594.8) * distance * np.sqrt( 7*img_pixelres**2 )
    frag_density_min = (1.616E-23) * np.divide(frag_mass_min,
                                                np.power(frag_radius_min,3))
    frag_jeans_mass_min = (3.305E-8) * np.power(frag_density_min*1000,-0.5)
    frag_masspmj_min = frag_mass_min / frag_jeans_mass_min

    # Make the axes.
    if spec_fig == False:
        ax = plt.subplot2grid((4, 4), (1, 0), rowspan=3, colspan=3)
        xax = plt.subplot2grid((4, 4), (0, 0), colspan=3, sharex=ax)
        yax = plt.subplot2grid((4, 4), (1, 3), rowspan=3, sharey=ax)
    ##fi

    # Plot main scatter.
    ax.scatter(np.log10(frag_masspmj[where_fragmono]),
                        frag_conc[where_fragmono], s=25,
                        facecolor=frag_fcolor[where_fragmono],
                        marker='o', edgecolor=frag_ecolor[where_fragmono])
    ax.scatter(np.log10(frag_masspmj[where_fragcomp]),
                        frag_conc[where_fragcomp], s=25,
                        facecolor=frag_fcolor[where_fragcomp],
                        marker='D', edgecolor=frag_ecolor[where_fragcomp])

    nbins = 10
    setlog = False
    if spec_fig == False:
        if orion == True:
            # No YSO
            if len(where_fragnyso) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_fragnyso]), nbins,
                        histtype='step', color='Black', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_fragnyso], nbins,
                        orientation='horizontal', histtype='step', color='Black',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Close Protostars
            if len(where_fragpc) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_fragpc]), nbins,
                        histtype='step', color='Green', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_fragpc], nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Non-Close Protostars
            if len(where_fragpnc) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_fragpnc]), nbins,
                        histtype='step', color='Green', linestyle='dotted', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_fragpnc], nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        linestyle='dotted', log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Disks
            if len(where_fragd) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_fragd]), nbins,
                        histtype='step', color='Olive', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_fragd], nbins,
                        orientation='horizontal', histtype='step', color='Olive',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi
        ##fi
    ##fi
    if spec_fig == False:
        if orion == False:
            if len(where_fragnyso) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_fragnyso]), nbins,
                        histtype='step', color='Black', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_fragnyso], nbins,
                        orientation='horizontal', histtype='step', color='Black',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Close Class 0+1
            if len(where_frag01c) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_frag01c]), nbins,
                        histtype='step', color='Green', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_frag01c], nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Non-Close Class 0+1
            if len(where_frag01nc) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_frag01nc]), nbins,
                        histtype='step', color='Green', linestyle='dotted', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_frag01nc], nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        linestyle='dotted', log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Flat
            if len(where_fragf) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_fragf]), nbins,
                        histtype='step', color='Orange', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_fragf], nbins,
                        orientation='horizontal', histtype='step', color='Orange',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Disks
            if len(where_frag2) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_masspmj[where_frag2]), nbins,
                        histtype='step', color='Olive', log=setlog)
                _,_,poly_d = yax.hist(frag_conc[where_frag2], nbins,
                        orientation='horizontal', histtype='step', color='Olive',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi
        ##fi
    ##fi

    # Add extras and set axis properties.
    ax.axvline(0, linestyle='dashed', color='k', linewidth=0.5)
    ax.axvline(np.log10(4), linestyle='dashed', color='k', linewidth=0.5)
    ax.axvline(np.log10(frag_masspmj_min), color='Black', linestyle='dotted',
                linewidth=0.5)
    ax.axhline(0.5,linestyle='dashed',color='k', linewidth=0.5)
    ax.set_xlabel(r'log(M$_{frag}$ / M$_{frag,J}$)')
    ax.set_ylabel('Concentraction')
    ax.set_xlim(-3,3.0)
    ax.set_ylim(0.05,1.05)

    # X-axis histogram axis properties.
    if spec_fig == False:
        xax.set_xlim(-3,3.0)
        xax.tick_params(axis='both', which='both', direction='in')
        xax.tick_params(labelbottom='off')
        xax.set_ylim(0,1.1)
        xax.set_yticks([0.25,0.5,0.75,1.0])
        xax.set_yticklabels(['','0.5','','1'])
        xax.axvline(0, linestyle='dashed', color='k', linewidth=0.5)
        xax.axvline(np.log10(4), linestyle='dashed', color='k', linewidth=0.5)
        xax.set_ylabel(r'N/N$_{max}$')
    ##fi

    # Y-axis histogram axis properties.
    if spec_fig == False:
        yax.set_ylim(0.05,1.0)
        yax.tick_params(axis='both', which='both', direction='in')
        yax.tick_params(labelleft='off')
        yax.set_xlim(0,1.1)
        yax.set_xticks([0.25,0.5,0.75,1.0])
        yax.set_xticklabels(['','0.5','','1'])
        yax.axhline(0.5,linestyle='dashed',color='k', linewidth=0.5)
        yax.set_xlabel(r'N/N$_{max}$')
    ##fi

    # Global properties and save.
    if spec_fig == False:
        plt.subplots_adjust(wspace=0, hspace=0)
        fig.suptitle(region+' Fragment Concentrations')
        scale_relation = r'$M/M_{Jeans} \propto d^{\frac{3}{2}}$'
        ax.annotate(scale_relation, (0.1,0.9), xycoords='axes fraction')
    ##fi

    if save_ind == True:
        plt.savefig(outplots,format='pdf')
    ##fi

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('concentration_'+region+'.png', dpi=300)
    ##fi

    if spec_fig == True:
        return fig, ax
    ##fi

    plt.clf()
#def

def IslandDensityRadius(fig,
                        island_data,
                        outplots,
                        region_properties,
                        plot_properties,
                        plot_individual = False,
                        ax = None,
                        save_ind = True,
                        spec_fig = False
                        ):
    '''
    IslandDensityRadius:

    Plot island volumetric number density as a function of island radius in
    log-log.

    Args:
        fig (figure): matplotlib figure object.
        island_data (Table): Astropy table of protostar data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data.
    orion = region_properties.orion
    region = region_properties.region
    distance = region_properties.distance
    beamsize = region_properties.beamsize
    noise = region_properties.noise
    island_radius = island_data['RADIUS'].data
    island_numvoldenspeak = island_data['NUMVOLDENSPEAK'].data
    island_minarea = beamsize**2 #pixels
    island_minradius = distance * np.sqrt(island_minarea/np.pi) / 206264.8 # pc
    where_islandmono = plot_properties.where_islandmono
    where_islandcomp = plot_properties.where_islandcomp
    island_fcolor = plot_properties.island_fcolor
    island_ecolor = plot_properties.island_ecolor
    where_islandnyso = plot_properties.where_islandnyso
    if orion == True:
        where_islandpc = plot_properties.where_islandpc
        where_islandpnc = plot_properties.where_islandpnc
        where_islandd = plot_properties.where_islandd
    ##fi
    if orion == False:
        where_island01c = plot_properties.where_island01c
        where_island01nc = plot_properties.where_island01nc
        where_islandf = plot_properties.where_islandf
        where_island2 = plot_properties.where_island2
    ##fi

    # Make constant radius fiducials and the detection limit.
    constant_r = np.arange(0.001,1.0,0.001)
    constant_n1jeans = 169.04 * np.power(constant_r,-2)
    constant_n2jeans = 169.04 * np.power(constant_r,-2) * 4
    # detection_limit = 1.4593E9 * 3 * noise / ( constant_r )
    detection_limit = 7206 * 3 * noise / ( constant_r )

    # Label the filter scale: 600" diameter.
    filter_radius = distance * 300 / 206264.8

    # Make the axis objects.
    if spec_fig == False:
        ax = plt.subplot2grid((4, 4), (1, 0), rowspan=3, colspan=3)
        xax = plt.subplot2grid((4, 4), (0, 0), colspan=3, sharex=ax)
        yax = plt.subplot2grid((4, 4), (1, 3), rowspan=3, sharey=ax)
    ##fi

    # Plot main scatter.
    ax.scatter(np.log10(island_radius[where_islandmono]),
                np.log10(island_numvoldenspeak[where_islandmono]), s=50,
                facecolor=island_fcolor[where_islandmono],
                marker='o', edgecolor=island_ecolor[where_islandmono])
    ax.scatter(np.log10(island_radius[where_islandcomp]),
                np.log10(island_numvoldenspeak[where_islandcomp]), s=50,
                facecolor=island_fcolor[where_islandcomp],
                marker='D', edgecolor=island_ecolor[where_islandcomp])

    # Add fiducials
    ax.plot(np.log10(constant_r),np.log10(constant_n1jeans), color='Magenta',
            label=r'1R$_{J}$', linewidth=0.5)
    ax.plot(np.log10(constant_r),np.log10(constant_n2jeans), color='Blue',
            label=r'2R$_{J}$', linewidth=0.5)
    ax.plot(np.log10(constant_r),np.log10(detection_limit), color='Black',
            label='Detection limit', linewidth=0.5, linestyle='dotted')
    ax.axvline(np.log10(island_minradius), linestyle='dashed', color='k',
                linewidth=0.5, label='Minimum size')
    ax.axvline(np.log10(filter_radius), linestyle='dashdot', color='k',
                linewidth=0.5, label='Filter Scale')

    if spec_fig == False:
        ax.legend(loc='upper right', fontsize=10, markerscale=5)
    ##fi

    nbins = 10
    setlog = False
    if spec_fig == False:
        if orion == True:
            # No YSO
            if len(where_islandnyso) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_islandnyso]), nbins,
                        histtype='step', color='Black', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_islandnyso]),
                        nbins, orientation='horizontal', histtype='step', color='Black',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Close Protostars
            if len(where_islandpc) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_islandpc]), nbins,
                        histtype='step', color='Green', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_islandpc]),
                        nbins, orientation='horizontal', histtype='step', color='Green',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Non-Close Protostars
            if len(where_islandpnc) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_islandpnc]), nbins,
                        histtype='step', color='Green', linestyle='dotted', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_islandpnc]),
                        nbins, orientation='horizontal', histtype='step', color='Green',
                        linestyle='dotted', log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Disks
            if len(where_islandd) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_islandd]), nbins,
                        histtype='step', color='Olive', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_islandd]),
                        nbins, orientation='horizontal', histtype='step', color='Olive',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi
        ##fi
        if orion == False:
            # No YSOs
            if len(where_islandnyso) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_islandnyso]), nbins,
                        histtype='step', color='Black', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_islandnyso]),
                        nbins, orientation='horizontal', histtype='step', color='Black',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Close Class 0+1
            if len(where_island01c) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_island01c]), nbins,
                        histtype='step', color='Green', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_island01c]),
                        nbins, orientation='horizontal', histtype='step', color='Green',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Non-Close Class 0+1
            if len(where_island01nc) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_island01nc]), nbins,
                        histtype='step', color='Green', linestyle='dotted', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_island01nc]),
                        nbins, orientation='horizontal', histtype='step', color='Green',
                        linestyle='dotted', log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Flat
            if len(where_islandf) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_islandf]), nbins,
                        histtype='step', color='Orange', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_islandf]),
                        nbins, orientation='horizontal', histtype='step', color='Orange',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Disks
            if len(where_island2) > 0:
                _,_,poly_f = xax.hist(np.log10(island_radius[where_island2]), nbins,
                        histtype='step', color='Olive', log=setlog)
                _,_,poly_d = yax.hist(np.log10(island_numvoldenspeak[where_island2]),
                        nbins, orientation='horizontal', histtype='step', color='Olive',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi
        ##fi
    ##fi

    # Add general properties
    ax.set_xlabel(r'log [R$_{island}$ pc]')
    ax.set_ylabel(r'log [ n / cm$^{-3}$]')
    ax.set_xlim(-2.5,0.5)
    ax.set_ylim(3.3,6.3)

    # X-axis histogram axis properties.
    if spec_fig == False:
        xax.set_xlim(-2.5,0)
        xax.tick_params(axis='both', which='both', direction='in')
        xax.tick_params(labelbottom='off')
        xax.set_ylim(0,1.1)
        xax.set_yticks([0.25,0.5,0.75,1.0])
        xax.set_yticklabels(['','0.5','','1'])
        xax.axvline(island_minradius, linestyle='dashed', color='k',
                    linewidth=0.5, label='Minimum size')
        xax.set_ylabel(r'N/N$_{max}$')
    ##fi

    # Y-axis histogram axis properties.
    if spec_fig == False:
        yax.set_ylim(3.3,6.3)
        yax.tick_params(axis='both', which='both', direction='in')
        yax.tick_params(labelleft='off')
        yax.set_xlim(0,1.1)
        yax.set_xticks([0.25,0.5,0.75,1.0])
        yax.set_xticklabels(['','0.5','','1'])
        yax.set_xlabel(r'N/N$_{max}$')
    ##fi

    # Globals and save.
    if spec_fig == False:
        fig.suptitle(region+' Island Number Densities')
        scale_relation1 = r'$\rho_{N} \propto d^{-1}$'
        scale_relation2 = r'$R \propto d^{1}$'
        ax.annotate(scale_relation1, (0.3,0.2), xycoords='axes fraction')
        ax.annotate(scale_relation2, (0.3,0.1), xycoords='axes fraction')
        fig.subplots_adjust(wspace=0, hspace=0)
    ##fi

    if save_ind == True:
        plt.savefig(outplots,format='pdf')
    ##fi

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('concentration_'+region+'.png', dpi=300)
    ##fi

    if spec_fig == True:
        return fig, ax
    ##fi

    plt.clf()
#def

def FragmentDensityRadius(fig,
                          frag_data,
                          outplots,
                          region_properties,
                          plot_properties,
                          plot_individual = False,
                          ax = None,
                          save_ind = True,
                          spec_fig = False
                          ):

    '''
    FragmentDensityRadius:

    Plot fragment volumetric number density as a function of fragment radius in
    log-log.

    Args:
        fig (figure): matplotlib figure object.
        frag_data (Table): Astropy table of protostar data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data.
    orion = region_properties.orion
    region = region_properties.region
    distance = region_properties.distance
    beamsize = region_properties.beamsize
    noise = region_properties.noise
    img_pixelres = region_properties.pixel_resolution
    frag_radius = frag_data['RADIUS'].data
    frag_numvoldenspeak = frag_data['NUMVOLDENSPEAK'].data
    # Find minimum fragment size, 7 pixels
    frag_minpix = 7 #pixels
    frag_minarea = frag_minpix * (img_pixelres)**2 #Square arcseconds
    frag_minradius = distance * np.sqrt( frag_minarea/np.pi ) / 206264.8 # pc
    where_fragmono = plot_properties.where_fragmono
    where_fragcomp = plot_properties.where_fragcomp
    frag_fcolor = plot_properties.frag_fcolor
    frag_ecolor = plot_properties.frag_ecolor

    where_fragnyso = plot_properties.where_fragnyso
    if orion == True:
        where_fragpc = plot_properties.where_fragpc
        where_fragpnc = plot_properties.where_fragpnc
        where_fragd = plot_properties.where_fragd
    ##fi
    if orion == False:
        where_frag01c = plot_properties.where_frag01c
        where_frag01nc = plot_properties.where_frag01nc
        where_fragf = plot_properties.where_fragf
        where_frag2 = plot_properties.where_frag2
    ##fi

    # Make constant radius fiducials and the detection limit.
    constant_r = np.arange(0.001,1.0,0.001)
    constant_n1jeans = 169.04 * np.power(constant_r,-2)
    constant_n2jeans = 169.04 * np.power(constant_r,-2) * 4
    # detection_limit = 1.4593E9 * 3 * noise / ( distance**2 * constant_r )
    detection_limit = 7206 * 3 * noise / ( constant_r )

    # Label the filter scale: 600" diameter.
    filter_radius = distance * 300 / 206264.8

    # Make the axis objects.
    if spec_fig == False:
        ax = plt.subplot2grid((4, 4), (1, 0), rowspan=3, colspan=3)
        xax = plt.subplot2grid((4, 4), (0, 0), colspan=3, sharex=ax)
        yax = plt.subplot2grid((4, 4), (1, 3), rowspan=3, sharey=ax)
    ##fi

    ax.scatter(np.log10(frag_radius[where_fragmono]),
                        np.log10(frag_numvoldenspeak[where_fragmono]), s=50,
                        facecolor=frag_fcolor[where_fragmono],
                        marker='o', edgecolor=frag_ecolor[where_fragmono])
    ax.scatter(np.log10(frag_radius[where_fragcomp]),
                        np.log10(frag_numvoldenspeak[where_fragcomp]), s=50,
                        facecolor=frag_fcolor[where_fragcomp],
                        marker='D', edgecolor=frag_ecolor[where_fragcomp])
    ax.plot(np.log10(constant_r), np.log10(constant_n1jeans), color='Magenta',
            linewidth=0.5, label=r'1R$_{J}$')
    ax.plot(np.log10(constant_r), np.log10(constant_n2jeans), color='Blue',
            linewidth=0.5, label=r'2R$_{J}$')
    ax.plot(np.log10(constant_r), np.log10(detection_limit), color='Black',
            linewidth=0.5, linestyle='dotted', label='Detection Limit')
    ax.axvline(np.log10(frag_minradius), linestyle='dashed', color='Black',
                linewidth=0.5, label='Minimum Sizes')
    ax.axvline(np.log10(filter_radius), linestyle='dashdot', color='Black',
                linewidth=0.5, label='Filter Scale')

    if spec_fig == False:
        ax.legend(loc='upper right', fontsize=10, markerscale=5)
    ##fi

    nbins = 10
    setlog = False

    # Plot histograms, rescaled so each peaks at 1
    if spec_fig == False:
        if orion == True:
            # No YSO
            if len(where_fragnyso) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_fragnyso]),
                        bins=nbins, histtype='step', color='Black', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_fragnyso]),
                        bins=nbins, orientation='horizontal', histtype='step',
                        color='Black', log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##fi
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##fi
            ##fi

            # Close Protostars
            if len(where_fragpc) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_fragpc]), nbins,
                        histtype='step', color='Green', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_fragpc]), nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Non-Close Protostars
            if len(where_fragpnc) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_fragpnc]), nbins,
                        histtype='step', color='Green', linestyle='dotted', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_fragpnc]), nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        linestyle='dotted', log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Disks
            if len(where_fragd) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_fragd]), nbins,
                        histtype='step', color='Olive', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_fragd]), nbins,
                        orientation='horizontal', histtype='step', color='Olive',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi
        ##fi
        if orion == False:
            if len(where_fragnyso) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_fragnyso]), nbins,
                        histtype='step', color='Black', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_fragnyso]), nbins,
                        orientation='horizontal', histtype='step', color='Black',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Close Class 0+1
            if len(where_frag01c) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_frag01c]), nbins,
                        histtype='step', color='Green', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_frag01c]), nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Non-Close Class 0+1
            if len(where_frag01nc) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_frag01nc]), nbins,
                        histtype='step', color='Green', linestyle='dotted', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_frag01nc]), nbins,
                        orientation='horizontal', histtype='step', color='Green',
                        linestyle='dotted', log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Flat
            if len(where_fragf) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_fragf]), nbins,
                        histtype='step', color='Orange', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_fragf]), nbins,
                        orientation='horizontal', histtype='step', color='Orange',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi

            # Disks
            if len(where_frag2) > 0:
                _,_,poly_f = xax.hist(np.log10(frag_radius[where_frag2]), nbins,
                        histtype='step', color='Olive', log=setlog)
                _,_,poly_d = yax.hist(np.log10(frag_numvoldenspeak[where_frag2]), nbins,
                        orientation='horizontal', histtype='step', color='Olive',
                        log=setlog)
                for pf in poly_f:
                    pf.set_xy( pf.get_xy()*np.array([1, 1/max(pf.get_xy()[:,1])]) )
                ##pf
                for pd in poly_d:
                    pd.set_xy( pd.get_xy()*np.array([1/max(pd.get_xy()[:,0]), 1]) )
                ##pd
            ##fi
        ##fi
    ##fi

    # Plot generic
    ax.set_xlabel(r'log [R$_{frag}$ / pc]')
    ax.set_ylabel(r'log [ n / cm$^{-3}$]')
    ax.set_xlim(-3,0)
    ax.set_ylim(3.3,6.3)

    # X-axis histogram axis properties.
    #xax.set_xscale('log')
    if spec_fig == False:
        xax.set_xlim(-3,0)
        xax.set_ylim(0,1.1)
        xax.tick_params(axis='both', which='both', direction='in')
        xax.tick_params(labelbottom='off')
        xax.set_yticks([0.25,0.5,0.75,1.0])
        xax.set_yticklabels(['','0.5','','1'])
        xax.axvline(frag_minradius, linestyle='dashed', color='k',
                    linewidth=0.5, label='Minimum size')
        xax.set_ylabel(r'N/N$_{max}$')
    ##fi

    # Y-axis histogram axis properties.
    #yax.set_yscale('log')
    if spec_fig == False:
        yax.set_ylim(3.3,6.3)
        yax.set_xlim(0,1.1)
        yax.tick_params(axis='both', which='both', direction='in')
        yax.tick_params(labelleft='off')
        yax.set_xticks([0.25,0.5,0.75,1.0])
        yax.set_xticklabels(['','0.5','','1'])
        yax.set_xlabel(r'N/N$_{max}$')
    ##fi

    # Set global properties and save.
    if spec_fig == False:
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.suptitle(region+' Fragment Number Densities')
        scale_relation1 = r'$\rho_{N} \propto d^{-1}$'
        scale_relation2 = r'$R \propto d^{1}$'
        ax.annotate(scale_relation1, (0.2,0.2), xycoords='axes fraction')
        ax.annotate(scale_relation2, (0.2,0.1), xycoords='axes fraction')
    ##fi

    if save_ind == True:
        plt.savefig(outplots,format='pdf')
    ##fi

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('fragmentdensity_'+region+'.png', dpi=300)
    ##fi

    if spec_fig == True:
        return fig, ax
    ##fi

    plt.clf()
#def

def MonolithicJeansRadius(fig,
                          frag_data,
                          outplots,
                          region_properties,
                          plot_properties,
                          plot_individual = False,
                          ax = None,
                          save_ind = True,
                          spec_fig = False
                          ):
    '''
    MonolithicJeansRadius:

    Make fragment radius histogram in terms of Jeans radius for monolithic
    fragments.

    Args:
        fig (figure): matplotlib figure object.
        island_data (Table): Astropy table of protostar data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data.
    orion = region_properties.orion
    region = region_properties.region
    frag_radiusprj = frag_data['RADIUSPRJ'].data
    frag_monolith = frag_data['MONOLITH'].data
    frag_radiusprj_mono = frag_radiusprj[np.where(frag_monolith == 1)[0]]
    if orion == True:
        where_fragpc = plot_properties.where_fragpc
        where_fragmonopc = np.where(frag_monolith[where_fragpc] == 1)[0]
        frag_radiusprj_monop = frag_radiusprj[where_fragpc][where_fragmonopc]
    ##fi
    if orion == False:
        where_frag01c = plot_properties.where_frag01c
        where_fragmono01c = np.where(frag_monolith[where_frag01c] == 1)[0]
        frag_radiusprj_monop = frag_radiusprj[where_frag01c][where_fragmono01c]
    ##fi

    nf_mono = len(frag_radiusprj_mono)
    if ax == None:
        ax = fig.add_subplot(111)
    ##fi
    binsize = 0.25
    binrange = [0,3]
    binnum = int((binrange[1]-binrange[0])/binsize)
    nums, bins, rects = ax.hist(frag_radiusprj_mono, bins=binnum, range=binrange,
                                    color='y')
    ax.set_xlim(0,3)
    ax.set_ylim(0,np.amax(nums)*1.25)
    nums, bins, rects = ax.hist(frag_radiusprj_monop, bins=binnum, range=binrange,
                                    color='olive')
    ax.annotate('Total N: '+str(nf_mono),xy=(0.7,0.7),xycoords='axes fraction')
    ax.set_xlabel(r"R$_{frag}$ / R$_{frag,J}$")
    ax.set_ylabel('N')
    ax.set_title(region+' Monolithic Fragment Radius Histogram')
    scale_relation = r'$R/R_{J} \propto d^{\frac{1}{2}}$'
    ax.annotate(scale_relation, (0.75,0.8), xycoords='axes fraction')
    if save_ind == True:
      plt.savefig(outplots, format='pdf')
    ##fi


    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('monolithic_'+region+'.png', dpi=300)

    plt.clf()
#def

def ComplexJeansRadius(fig,
                       frag_data,
                       outplots,
                       region_properties,
                       plot_properties,
                       plot_individual = False,
                       ax = None,
                       save_ind = True,
                       spec_fig = False
                       ):
    '''
    ComplexJeansRadius:

    Make fragment radius histogram in terms of Jeans radius for monolithic
    fragments.

    Args:
        fig (figure): matplotlib figure object.
        island_data (Table): Astropy table of protostar data.
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data.
    orion = region_properties.orion
    region = region_properties.region
    frag_radiusprj = frag_data['RADIUSPRJ'].data
    frag_monolith = frag_data['MONOLITH'].data
    frag_radiusprj_comp = frag_radiusprj[np.where(frag_monolith == 0)[0]]

    if orion == True:
        where_fragpc = plot_properties.where_fragpc
        where_fragmonopc = np.where(frag_monolith[where_fragpc] == 0)[0]
        frag_radiusprj_compp = frag_radiusprj[where_fragpc][where_fragmonopc]
    ##fi
    if orion == False:
        where_frag01c = plot_properties.where_frag01c
        where_fragmono01c = np.where(frag_monolith[where_frag01c] == 0)[0]
        frag_radiusprj_compp = frag_radiusprj[where_frag01c][where_fragmono01c]
    ##fi

    nf_comp = len(frag_radiusprj_comp)
    if ax == None:
        ax = fig.add_subplot(111)
    ##fi
    binsize = 0.25
    binrange = [0,3]
    binnum = int((binrange[1]-binrange[0])/binsize)
    nums, bins, rects = ax.hist(frag_radiusprj_comp, bins=binnum, range=binrange,
                                    color='y')
    ax.set_xlim(0,4)
    ax.set_ylim(0,np.amax(nums)*1.25)
    nums, bins, rects = ax.hist(frag_radiusprj_compp, bins=binnum, range=binrange,
                                    color='olive')
    ax.annotate('Total N: '+str(nf_comp),xy=(0.7,0.7),xycoords='axes fraction')
    ax.set_xlabel(r"R$_{frag}$ / R$_{frag,J}$")
    ax.set_ylabel('N')
    ax.set_title(region+' Complex Fragment Radius Histogram')
    scale_relation = r'$R/R_{J} \propto d^{\frac{1}{2}}$'
    ax.annotate(scale_relation, (0.75,0.8), xycoords='axes fraction')
    if save_ind == True:
      plt.savefig(outplots, format='pdf')
    ##fi

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('complex_'+region+'.png', dpi=300)

    plt.clf()
#def

def ColumnMassFraction( fig,
                        extinction_data,
                        outplots,
                        region_properties,
                        plot_properties,
                        plot_individual = False,
                        ax = None,
                        save_ind = True,
                        spec_fig = False
                        ):
    '''
    ColumnMassFraction:

    Make fragment radius histogram in terms of Jeans radius for monolithic
    fragments.

    Args:
        fig (figure): matplotlib figure object.
        extinction_data (Table): Numpy
        outplots (PdfPages): The PdfPages object on which to write the plot.
        region_properties (RegionProperties): The RegionProperties object.
        plot_properties (PlotProperties): The PlotProperties object.
        plot_individual (Boolean):Make an individual png of the image [False]
        ax (Object): The axis object to use in the plotting routine
        save_ind (Boolean): Save each individual figure to the outplot
        spec_fig (Boolean): The program is being called for making paper Figure.

    Returns:
        None?
    '''

    # Set region properties and unpack data
    orion = region_properties.orion
    region = region_properties.region
    all_evals, all_cml_evals, yso_evals, yso_cml_norm = extinction_data
    ext_evals, s2_evals, proto_evals, whole_evals = all_evals
    ext_cml_norm, s2_cml_norm, proto_cml_norm, whole_cml_norm = all_cml_evals
    if orion == True:
        p_evals, d_evals, fp_evals, rp_evals = yso_evals
        p_cml_norm, d_cml_norm, fp_cml_norm, rp_cml_norm = yso_cml_norm
    ##fi
    if orion == False:
        c01_evals, f_evals, c2_evals, c3_evals = yso_evals
        c01_cml_norm, f_cml_norm, c2_cml_norm, c3_cml_norm = yso_cml_norm
    ##fi

    if ax == None:
        ax = fig.add_subplot(111)
    ##fi
    ax.plot(ext_evals*16.6, ext_cml_norm, c='DodgerBlue', label='NICEST')
    ax.plot([-1,ext_evals[-1]*16.6], [1,1], c='DodgerBlue')
    ax.plot(whole_evals*16.6, whole_cml_norm, c='DodgerBlue', linestyle='dotted',
            label='NICEST No Mask')
    ax.plot([-1,whole_evals[-1]*16.6], [1,1], c='DodgerBlue', linestyle='dashed')
    ax.plot(s2_evals*16.6, s2_cml_norm, c='Red', label='Islands')
    ax.plot([-1,s2_evals[-1]*16.6], [1,1], c='Red')
    #ax.plot(proto_evals*16.6, proto_cml_norm, c='m', label='Protostars')

    if orion == True:
        ax.plot(p_evals*16.6, p_cml_norm, c='DarkGreen', linestyle='dashed',
                alpha=0.5, label='P')
        if len(p_evals) > 0:
            ax.plot([-1,p_evals[-1]*16.6], [1,1], c='DarkGreen', linestyle='dashed',
                    alpha=0.5)
        ##fi
        ax.plot(d_evals*16.6, d_cml_norm, c='Olive', linestyle='dashed',
                alpha=0.5, label='D')
        if len(d_evals) > 0:
            ax.plot([-1,d_evals[-1]*16.6], [1,1], c='Olive', linestyle='dashed',
                    alpha=0.5)
        ##fi
        #ax.plot(fp_evals*16.6, fp_cml_norm, c='g', label='FP')
        #ax.plot(rp_evals*16.6, rp_cml_norm, '--b', label='RP')
    ##fi
    if orion == False:
        ax.plot(c01_evals*16.6, c01_cml_norm, c='DarkGreen', linestyle='dashed',
                alpha=0.5, label='Class 0+1')
        if len(c01_evals) > 0:
            ax.plot([-1,c01_evals[-1]*16.6], [1,1], c='DarkGreen', linestyle='dashed',
                    alpha=0.5)
        ##fi
        ax.plot(f_evals*16.6, f_cml_norm, c='Orange', linestyle='dashed',
                alpha=0.5, label='Flat')
        if len(f_evals) > 0:
            ax.plot([-1,f_evals[-1]*16.6], [1,1], c='Orange', linestyle='dashed',
                    alpha=0.5)
        ##fi
        ax.plot(c2_evals*16.6, c2_cml_norm, c='Olive', linestyle='dashed',
                alpha=0.5, label='Class 2')
        if len(c2_evals) > 0:
            ax.plot([-1,c2_evals[-1]*16.6], [1,1], c='Olive', linestyle='dashed',
                    alpha=0.5)
        ##fi
        #ax.plot(c3_evals*16.6, c3_cml_norm, '--b', label='Class 3')
    ##fi

    if os.path.isfile('Catalogs/'+region+'_extinction_mass.txt'):
        ext_mass = np.loadtxt('Catalogs/'+region+'_extinction_mass.txt')
        frac_isl_tot = ext_mass[2] / ext_mass[0]
        ax.annotate(r'$M_{tot.}$ = '+r'{:.0f}'.format(ext_mass[0])+r'$M_{\odot}$',
                    (0.35,0.9), xycoords='axes fraction')
        ax.annotate(r'$M_{Isl.}/M_{tot.}$ = '+'{:.1f}%'.format(frac_isl_tot*100),
                    (0.35,0.85), xycoords='axes fraction')

        # Different for Orion and not.
        if orion == True:
            frac_proto_tot = ext_mass[4] / ext_mass[0]
            frac_proto_isl = ext_mass[4] / ext_mass[2]
            frac_proto_disk = ext_mass[4] / ext_mass[5]
            ax.annotate(r'$M_{P}/M_{tot.}$ = '+'{:.1f}%'.format(frac_proto_tot*100),
                        (0.35,0.8), xycoords='axes fraction')
            ax.annotate(r'$M_{P}/M_{Isl.}$ = '+'{:.1f}%'.format(frac_proto_isl*100),
                        (0.35,0.75), xycoords='axes fraction')
            ax.annotate(r'$N_{P}/N_{D}$ = '+'{:.1f}%'.format(frac_proto_disk*100),
                        (0.35,0.7), xycoords='axes fraction')

        ##fi
        if orion == False:
            frac_c01_tot = ext_mass[4] / ext_mass[0]
            frac_c01_isl = ext_mass[4] / ext_mass[2]
            frac_c01_c2 = ext_mass[4] / ext_mass[6]
            ax.annotate(r'$M_{0+1}/M_{tot}$ = '+'{:.1f}%'.format(frac_c01_tot*100),
                        (0.35,0.8), xycoords='axes fraction')
            ax.annotate(r'$M_{0+1}/M_{Isl.}$ = '+'{:.1f}%'.format(frac_c01_isl*100),
                        (0.35,0.75), xycoords='axes fraction')
            ax.annotate(r'$N_{0+1}/N_{2}$ = '+'{:.1f}%'.format(frac_c01_c2*100),
                        (0.35,0.7), xycoords='axes fraction')
        ##fi
    ##fi

    ax.legend()
    ax.set_ylim(0,1.5)
    ax.set_xlim(0,100)
    ax.set_xlabel(r'Column [$10^{21}$ cm$^{-2}$]')
    ax.set_ylabel('Fractional Cumulative Mass')
    ax.set_title(region+' Column')
    scale_relation = r'$\Sigma \propto d^{2}$'
    ax.annotate(scale_relation, (0.4,0.95), xycoords='axes fraction')
    if save_ind == True:
      plt.savefig(outplots, format='pdf')
    ##fi

    # If the individual plot keyword is set then plot a .png as well:
    if plot_individual == True:
        plt.savefig('column_'+region+'.png', dpi=300)

    plt.clf()
#def


class PlotProperties:
    '''
    PlottingProperties:

    Contains information required for plotting for each set of fragment, island
    and protostar data.

    Attributes:
        n_frags (int): The number of fragments
        n_island (int): The number of islands
        n_proto (int): The number of protostars

        All following properties are arrays of varying sizes.

        Fragment/Island protostar masks for fragment/island data:

        [Orion]
        where_(frag/island)p: protostars
        where_...pc: close protostars (within beam)
        where_...pnc: not close protostars
        where_...d: disks
        where_...fp: faint protostars
        where_...rp: red protostars
        where_...nyso: No YSO
        [Not Orion]
        where_...01: Class 0/1
        where_...01c: close Class 0/1 (within beam)
        where_...01nc: not close Class 0/1
        where_...f: Flat spectrum
        where_...2: Class 2
        where_...3: Class 3
        where_...nyso: No YSO

        Fragment/Island monolithic masks for fragment/island data:

        where_...mono: Monolithic island or fragment in island
        where_...comp: Complex island or fragment in island

        Fragment/Island coloring masks for fragment/island data:
        ..._ecolor: Marker edge color
        ..._fcolor: Marker face color

        Protostar masks for protostar data:

        [Orion]
        where_proto: Protostars
        where_disk: Disks
        where_rp: Red protostars
        where_fp: Faint protostars
        [Not Orion]
        where_class01: Class 0/1
        where_flat: Flat spectrum
        where_class2: Class 2
        where_class3: Class 3

    Args:
        region_properties (RegionProperties): Describes the region of interest
        frag_data (Table): Astropy table of fragment data
        island_data (Table): Astropy table of island data
        proto_data (Table): Astropy table of protostar data

    Returns:
        None

    Raises:
        None
    '''


    def __init__(   self,
                    region_properties,
                    frag_data,
                    island_data,
                    proto_data
                ):

        #####################
        # FRAGMENT PROPERTIES
        #####################

        # Set region properties.
        orion = region_properties.orion

        # Fragment, Island and protostar numbers.
        self.n_frags = len(frag_data)
        self.n_islands = len(island_data)
        self.n_protos = len(proto_data)

        # Make where's for fragments and their associated protostars.
        frag_class = frag_data['CLASS'].data
        frag_close = frag_data['CLOSE'].data

        # Get masks for Megeath data.
        if orion == True:
            self.where_fragp = np.where(frag_class == 'P')[0]
            self.where_fragpc = np.where( (frag_class == 'P') &
                                          (frag_close == 1))[0]
            self.where_fragpnc = np.where( (frag_class == 'P') &
                                           (frag_close == 0))[0]
            self.where_fragd = np.where(frag_class == 'D')[0]
            self.where_fragfp = np.where(frag_class == 'FP')[0]
            self.where_fragrp = np.where(frag_class == 'RP')[0]
            self.where_fragnyso = np.where(frag_class == 'None')[0]
        ##fi

        # Get masks for Dunham Data.
        if orion == False:
            self.where_frag01 = np.where(frag_class == '0+1')[0]
            self.where_frag01c = np.where( (frag_class == '0+1') &
                                      (frag_close == 1))[0]
            self.where_frag01nc = np.where( (frag_class == '0+1') &
                                       (frag_close == 0))[0]
            self.where_fragf = np.where(frag_class == 'F')[0]
            self.where_frag2 = np.where(frag_class == '2')[0]
            self.where_frag3 = np.where(frag_class == '3')[0]
            self.where_fragnyso = np.where(frag_class == 'None')[0]
        ##fi

        # Make monolithic and complex masks.
        frag_monolith = frag_data['MONOLITH'].data
        self.where_fragmono = np.where(frag_monolith == 1)[0]
        self.where_fragcomp = np.where(frag_monolith == 0)[0]

        # Make an array of colors, representing classes.
        frag_ecolor = np.empty(self.n_frags, dtype='U10')

        # Set colors for Megeath data.
        if orion == True:
            frag_ecolor[self.where_fragnyso] = 'Black'
            frag_ecolor[self.where_fragd] = 'Olive'
            frag_ecolor[self.where_fragrp] = 'Red'
            frag_ecolor[self.where_fragfp] = 'Blue'
            frag_ecolor[self.where_fragpc] = 'Green'
            frag_ecolor[self.where_fragpnc] = 'Green'
        ##fi

        # Set colors for Dunham data.
        if orion == False:
            frag_ecolor[self.where_fragnyso] = 'Black'
            frag_ecolor[self.where_frag01c] = 'Green'
            frag_ecolor[self.where_frag01nc] = 'Green'
            frag_ecolor[self.where_fragf] = 'Orange'
            frag_ecolor[self.where_frag2] = 'Olive'
            frag_ecolor[self.where_frag3] = 'Purple'
        ##fi

        self.frag_ecolor = frag_ecolor

        # Set the facecolor as green for close protostars.
        frag_fcolor = np.empty(self.n_frags, dtype='U10')
        frag_fcolor[:] = 'none'
        if orion == True:
            frag_fcolor[self.where_fragpc] = 'Green'
        ##fi
        if orion == False:
            frag_fcolor[self.where_frag01c] = 'Green'
        ##fi

        self.frag_fcolor = frag_fcolor

        ###################
        # ISLAND PROPERTIES
        ###################

        # Make where's for islands and ther associated protostars.
        island_class = island_data['CLASS'].data
        island_close = island_data['CLOSE'].data

        # Get masks for Megeath data.
        if orion == True:
            self.where_islandp = np.where(island_class == 'P')[0]
            self.where_islandpc = np.where((island_class == 'P') &
                                        (island_close == 1))[0]
            self.where_islandpnc = np.where((island_class == 'P') &
                                         (island_close == 0))[0]
            self.where_islandd = np.where(island_class == 'D')[0]
            self.where_islandfp = np.where(island_class == 'FP')[0]
            self.where_islandrp = np.where(island_class == 'RP')[0]
            self.where_islandnyso = np.where(island_class == 'None')[0]
        ##fi

        # Get masks for Dunham data.
        if orion == False:
            self.where_island01 = np.where(island_class == '0+1')[0]
            self.where_island01c = np.where((island_class == '0+1') &
                                         (island_close == 1))[0]
            self.where_island01nc = np.where((island_class == '0+1') &
                                          (island_close == 0))[0]
            self.where_islandf = np.where(island_class == 'F')[0]
            self.where_island2 = np.where(island_class == '2')[0]
            self.where_island3 = np.where(island_class == '3')[0]
            self.where_islandnyso = np.where(island_class == 'None')[0]
        ##fi

        # Make monolithic and complex masks.
        island_monolith = island_data['MONOLITH'].data
        self.where_islandmono = np.where(island_monolith == 1)[0]
        self.where_islandcomp = np.where(island_monolith == 0)[0]

        # Make an array of colors, representing classes.
        island_ecolor = np.empty(self.n_islands, dtype='U10')

        # Set colors for Megeath data.
        if orion == True:
            island_ecolor[self.where_islandnyso] = 'Black'
            island_ecolor[self.where_islandd] = 'Olive'
            island_ecolor[self.where_islandrp] = 'Red'
            island_ecolor[self.where_islandfp] = 'Blue'
            island_ecolor[self.where_islandpc] = 'Green'
            island_ecolor[self.where_islandpnc] = 'Green'
        ##fi

        # Set colors for Dunham data.
        if orion == False:
            island_ecolor[self.where_islandnyso] = 'Black'
            island_ecolor[self.where_island01c] = 'Green'
            island_ecolor[self.where_island01nc] = 'Green'
            island_ecolor[self.where_islandf] = 'Orange'
            island_ecolor[self.where_island2] = 'Olive'
            island_ecolor[self.where_island3] = 'Purple'
        ##fi

        self.island_ecolor = island_ecolor

        # Set facecolors as green for close protostars.
        island_fcolor = np.empty(self.n_islands, dtype='U10')
        island_fcolor[:] = 'none'
        if orion == True:
            island_fcolor[self.where_islandpc] = 'Green'
        ##fi
        if orion == False:
            island_fcolor[self.where_island01c] = 'Green'
        ##fi

        self.island_fcolor = island_fcolor

        ######################
        # PROTOSTAR PROPERTIES
        ######################

        # Also do protostar data for both catalogs:
        proto_class = proto_data['CLASS'].data

        # Megeath:
        if orion == True:
            self.where_proto = np.where(proto_class == 'P')[0]
            self.where_disk = np.where(proto_class == 'D')[0]
            self.where_rp = np.where(proto_class == 'RP')[0]
            self.where_fp = np.where(proto_class == 'FP')[0]
        ##fi

        # Dunham:
        if orion == False:
            self.where_class01 = np.where(proto_class == '0+1')[0]
            self.where_flat = np.where(proto_class == 'F')[0]
            self.where_class2 = np.where(proto_class == '2')[0]
            self.where_class3 = np.where(proto_class == '3')[0]
        ##fi
    #def
#cls
