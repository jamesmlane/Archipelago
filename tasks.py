# ----------------------------------------------------------------------------
#
# TITLE - tasks.py
# AUTHOR - James Lane edited from Steve Mairs
# PROJECT - archipelago
# CONTENTS:
#	1. KillNoisyEdges
#   2. GetInOrion
#   3. GetDistance
#   4. GetAnalysisFiles
#   5. ZipAnalysisFiles
#   6. StackCatalogs
#   7. NanToZero
#   8. PopulateDirectories
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''
General tasks in the Archipelago project.
'''

__author__ = "James Lane"

#Imports
import numpy as np
import subprocess
import os
import sys
import pdb
from astropy.io import fits
from astropy import table

def KillNoisyEdges(filename, min_var):
    '''
    KillNoisyEdges:

    Cut the noisy edges off a map. Default output is filename with '_clipped'
    appended.

    Args:
        filename (str): Name of the (fits) file.
        min_var (float): Variance level above which all data will be trimmed.

    Returns:
        None
    '''

    # Get the data.
    filebase = filename.split('.', 1)[0]
    img = fits.open(filebase+'.fits')

    # Perform the cut
    bad_ind = np.where( img[1].data > min_var )
    img[0].data[ bad_ind ] = np.nan

    # Write to file.
    img.writeto(filebase+'_clipped.fits', clobber=True)

#def

def GetInOrion(region):
    '''
    GetInOrion:

    Determine whether or not the region is in the Orion field, and therefore
    the Megeath protostar catalog should be used.

    Args:
        region (str): Name of the region

    Returns:
        in_orion (Boolean): Is the region in Orion (ie: use Megeath protostars)

    Raises:
        RuntumeError: The supplied region name is not identified.
    '''

    orion_fields = ['OrionA',
                    'OrionAS',
                    'OrionAN',
                    'OrionB',
                    'OrionB_L1622',
                    'OrionB_N2023',
                    'OrionB_N2068']

    non_orion_fields = ['Aquila',
                        'Auriga',
                        'CepheusL1228',
                        'CepheusL1251',
                        'CepheusSouth',
                        'CrA',
                        'IC5146',
                        'OphScoMain',
                        'OphScoN2',
                        'OphScoN3',
                        'OphScoN6',
                        'PerseusIC348',
                        'PerseusWest',
                        'PipeB59',
                        'PipeE1',
                        'Serpens',
                        'SerpensAquila',
                        'SerpensE',
                        'SerpensMWC297',
                        'SerpensMain',
                        'SerpensN',
                        'Taurus',
                        'TaurusL1495',
                        'TaurusSouth',
                        'TaurusTMC']

    if region in orion_fields: in_orion = True
    elif region in non_orion_fields: in_orion = False
    else: raise RuntimeError(region+', region not known!')

    return in_orion

def GetDistance(reg):
    '''
    GetDistance

    Get the canonical distance for each region

    Args:
        reg (string) - The name of the region

    Returns:
        distance (float) - The distance in pc
    '''

    regions = np.array(['Aquila',
                        'Auriga',
                        'CepheusL1228',
                        'CepheusL1251',
                        'CepheusSouth',
                        'CrA',
                        'IC5146',
                        'OphScoMain',
                        'OphScoN2',
                        'OphScoN3',
                        'OphScoN6',
                        'OrionA',
                        'OrionAN',
                        'OrionAS',
                        'OrionB',
                        'OrionB_L1622',
                        'OrionB_N2023',
                        'OrionB_N2068',
                        'PerseusIC348',
                        'PerseusWest',
                        'PipeB59',
                        'PipeE1',
                        'Serpens',
                        'SerpensAquila',
                        'SerpensE',
                        'SerpensMain',
                        'SerpensMWC297',
                        'SerpensN',
                        'Taurus',
                        'TaurusL1495',
                        'TaurusSouth',
                        'TaurusTMC'])
    distances = np.array([436,      # Aquila
                          450,      # Auriga
                          0,        # CepheusL1228
                          0,        # CepheusL1251
                          0,        # CepheusSouth
                          0,        # CrA
                          1000,     # IC5146
                          137,      # OphScoMain
                          0,        # OphN2
                          0,        # OphN3
                          0,        # OphN6
                          0,        # OrionA
                          388,      # OrionAN
                          420,      # OrionAS
                          388,      # OrionB
                          388,      # OrionB_L1622
                          388,      # OrionB_N2023
                          388,      # OrionB_N2068
                          303,      # PerseusIC348
                          271,      # PerseusWest
                          0,        # PipeB59
                          0,        # PipeE1
                          436,      # Serpens
                          436,      # SerpensAquila
                          436,      # SerpensE
                          436,      # SerpensMain
                          436,      # SerpensMWC297
                          436,      # SerpensN
                          140,      # Taurus
                          140,      # TaurusL1495
                          140,      # TaurusSouth
                          140       # TaurusTMC
                          ])

    if not reg in regions:
        raise RuntimeError(reg+', region not known!')
    where_reg = np.where(regions == reg)[0][0]
    dist = distances[where_reg]
    return dist

def GetPrintName(reg):
    '''
    GetPrintName

    Get a printable name for each region

    Args:
        reg (string) - The name of the region

    Returns:
        name (str) - The printable name of the region
    '''

    regions = np.array(['Aquila',
                        'Auriga',
                        'CepheusL1228',
                        'CepheusL1251',
                        'CepheusSouth',
                        'CrA',
                        'IC5146',
                        'OphScoMain',
                        'OphScoN2',
                        'OphScoN3',
                        'OphScoN6',
                        'OrionA',
                        'OrionAN',
                        'OrionAS',
                        'OrionB',
                        'OrionB_L1622',
                        'OrionB_N2023',
                        'OrionB_N2068',
                        'PerseusIC348',
                        'PerseusWest',
                        'PipeB59',
                        'PipeE1',
                        'Serpens',
                        'SerpensAquila',
                        'SerpensE',
                        'SerpensMain',
                        'SerpensMWC297',
                        'SerpensN',
                        'Taurus',
                        'TaurusL1495',
                        'TaurusSouth',
                        'TaurusTMC'])
    printable_regions = np.array([  436,      # Aquila
                                    'Aur',    # Auriga
                                    0,        # CepheusL1228
                                    0,        # CepheusL1251
                                    0,        # CepheusSouth
                                    0,        # CrA
                                    'IC5146', # IC5146
                                    'OphSco', # OphScoMain
                                    0,        # OphN2
                                    0,        # OphN3
                                    0,        # OphN6
                                    0,        # OrionA
                                    388,      # OrionAN
                                    'OriAS',  # OrionAS
                                    'OriB',   # OrionB
                                    388,      # OrionB_L1622
                                    388,      # OrionB_N2023
                                    388,      # OrionB_N2068
                                    'IC348',  # PerseusIC348
                                    'PersW',  # PerseusWest
                                    0,        # PipeB59
                                    0,        # PipeE1
                                    436,      # Serpens
                                    'Ser',    # SerpensAquila
                                    436,      # SerpensE
                                    436,      # SerpensMain
                                    436,      # SerpensMWC297
                                    436,      # SerpensN
                                    'Tau',    # Taurus
                                    'L1495',  # TaurusL1495
                                    140,      # TaurusSouth
                                    140       # TaurusTMC
                                    ])

    if not reg in regions:
        raise RuntimeError(reg+', region not known!')
    where_reg = np.where(regions == reg)[0][0]
    printreg = printable_regions[where_reg]
    return printreg


def GetAnalysisFiles(region):
    '''
    GetAnalysisFiles:

    Get the files required for archipelago analysis.

    Args:
        region (str): Name of the region (no .fits extension)

    Returns:
        None
    '''

    file1 = 'Transfer/'+region+'_analysis.tar.gz'

    scp_command = 'scp -C jlane93@turbot.phys.uvic.ca:'
    scp_command = scp_command+'~/Projects/Archipelago/Islands/Data/IR3/'+region+'/'
    scp_command = scp_command+file1+' ./'
    os.system(scp_command)

    unzip_command = 'tar -zxvf '+region+'_analysis.tar.gz'
    os.system(unzip_command)
    os.system('rm -f '+region+'_analysis.tar.gz')

def ZipAnalysisFiles(region):
    '''
    ZipAnalysisFiles:

    Zip the files required for archipelago analyss

    Args:
        region (str): Name of the region (no .fits extension)

    Returns:
        None
    '''

    file1 = 'island_info/'+region+'_islands_out.sdf'
    file2 = 'island_info/'+region+'_islands_outcat.FIT'
    file3 = 'fragment_info/'+region+'_peak-cat.fits'
    file4 = 'fragment_info/'+region+'_frag_protos.txt'

    os.system('mkdir Transfer')
    os.system('mkdir ClumpFind')
    os.system('cp -f '+file1+' ClumpFind')
    os.system('cp -f '+file2+' ClumpFind')
    os.system('cp -f '+file3+' ClumpFind')
    os.system('cp -f '+file4+' ClumpFind')
    zip_name = 'Transfer/'+region+'_analysis.tar.gz'
    zip_command = 'tar -zcvf '+zip_name+' ClumpFind'

    os.system(zip_command)

def StackCatalogs(fileout,*args):
    '''
    StackCatalogs:

    Stack multiple archipelago catalogs of the same type.

    Args:
        fileout (string) - Name of the output file (with .FIT ending)
        *args (strings) - The names of the catalogs to be stacked.

    Returns:
        None

    Outputs:

    '''

    # Empty list for catalogs
    catalogs = []
    for i,filename in enumerate(args[0]):
        catalog = fits.getdata(filename,1)
        # Turn into Table object
        catalogs.append( table.Table(catalog) )
    ####

    # stack into the new catalog
    try:
       input = raw_input
    except NameError:
       pass
    while fileout[-4:] != '.FIT':
        print('Error! Filename ending is not .FIT')
        fileout = input('Enter a new name...')

    # Make the
    new_catalog = table.vstack(catalogs)
    new_catalog.write(fileout, format='fits', overwrite=True)
#def

def StackExtinction(fileout,*args):
    '''
    StackExtinction:

    Stack multiple archipelago extinction-mass binaries.

    Args:
        *args (strings) - The names of the binaries to be stacked.

    Returns:
        None

    Outputs:

    '''

    # Empty list for catalogs
    catalogs = []

    for i,filename in enumerate(args):
        catalog = np.load(filename)

    ####

    # stack into the new catalog
    new_catalog = table.vstack(catalogs)
#def

def NanToZero(filename,hdu_index=0):
    '''
    NanToZero:

    Convert all the NaNs in a fits image to zeros.

    Args:
        filename (string) - The name of the fits file.
        hdu_index (int) - Use an HDU index other than 0

    Returns:
        filename+'_zeros' - The output file
    '''

    hdu = fits.open(filename)
    data = hdu[hdu_index].data
    hdr = hdu[hdu_index].header

    where_nan = np.isnan(data)
    data[where_nan] = 0

    fileout = os.path.splitext(filename)[0]+'_zeros.fits'
    fits.writeto(fileout,data,header=hdr,clobber=True)
#def

def PopulateDirectories():
    '''
    PopulateDirectores:

    Create the directories necessary to hold all of the data for the
    Archipelago data.

    Args:
        None

    Returns:
        None
    '''
    os.system('mkdir {ClumpFind,Catalogs,Dependencies,Plots,Test}')
#def
