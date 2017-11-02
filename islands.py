# ----------------------------------------------------------------------------
#
# TITLE - islands.py
# AUTHOR - James Lane
# PROJECT - archipelago
# CONTENTS:
#	1. FindIslands
#   2. GetClumpCatalog
#   3. Beam2Pix
#   4. Pix2Beam
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''
Functions for running the island and fragment finding routines for the
Archipelago project. Created by Steve Mairs; assembled and presented by
James Lane
'''

__author__ = "James Lane"

# Imports
import numpy as np
import os
from astropy.io import fits
import subprocess

def FindIslands(region,noise):
    '''
    FindIslands:

    Args:
        region (string) - Region for finding sources in.
        noise (float) - Noise in the map, in mJy/arcsecond sq.

    Returns:
        None

    Outputs:
        Both JSA_CATALOGUE and Findclumps output many files, see:
        http://www.starlink.ac.uk/docs/sun255.htx/sun255ss5.html#Q1-11-37
        http://www.starlink.ac.uk/docs/sun265.htx/sun265ss9.html

    History:
        November 2016: Created from Steve Mairs code - James Lane
        July 2017:  Heavily edited, prepared for Python3, Removed excess catalog
                    creation - James Lane
    '''

    # Keywords
    beamwidth = 15.0 # SCUBA-2 850um beam in arcseconds.
    method='ClumpFind' # CUPID clumpfinding method.
    sdffile=region+'.sdf'
    fitsfile=region+'.fits'

    # Determine the smallest sized island to keep, the beam FWHM in pixels.
    hdr = fits.open(fitsfile)[0].header
    if 'CDELT2' in hdr.keys():
        pix_size_degrees=hdr['CDELT2']
        pix_size_arcsecs=206264.806*pix_size_degrees*(np.pi/180.0)
        beam_fwhm=beamwidth
        beam_fwhm_pix=beam_fwhm/pix_size_arcsecs
    else:
        pix_size_arcsecs=3.0
        beam_fwhm=beamwidth
        beam_fwhm_pix=beam_fwhm/pix_size_arcsecs
    ##fi

    ### Begin
    print('')
    print('')
    print('######')
    print('######')
    print('')
    print('Now finding ISLANDS...')
    print('')
    print('Beamwidth   = '+str(beamwidth)+ ' arcseconds')
    print('Beamsize    = '+str(round(beam_fwhm_pix,1))+' pixels')
    print('')
    print('######')
    print('######')

    # Clumpfinding algorithm parameters.
    ClumpFind_AllowEdge=0
    ClumpFind_DeltaT=3.0*noise#20.0*noise
    ClumpFind_FwhmBeam=0.0#beam_fwhm_pix/2.0 <-- For detecting Islands
    ClumpFind_IDLAlg=0
    ClumpFind_Level1=3.0*noise
    ClumpFind_Level2=3.0*noise
    ClumpFind_MaxBad=1.0
    ClumpFind_MinPix=(beam_fwhm_pix)**2.0 #Beam squared. 225 as.sq. vs 229 as.sq. actual beam.
    #ClumpFind_Naxis= Just run with the default for now
    ClumpFind_RMS=noise
    ClumpFind_Tlow=3.0*noise
    #ClumpFind_VeloRes= JUST USE DEFUALT FOR NOW

    # Construct the dictionary to hold the parameters.
    paramdict={ 'ClumpFind_AllowEdge':ClumpFind_AllowEdge,
                'ClumpFind_FwhmBeam':ClumpFind_FwhmBeam,
                'ClumpFind_IDLAlg':ClumpFind_IDLAlg,
                'ClumpFind_MaxBad':ClumpFind_MaxBad,
                'ClumpFind_MinPix':ClumpFind_MinPix,
                'ClumpFind_RMS':ClumpFind_RMS,
                'ClumpFind_Level1':ClumpFind_Level1,
                'ClumpFind_Level2':ClumpFind_Level2}

    #Create the configuration file, change underscores to periods -- Why??
    configparams=open("config.txt","w")
    dummy=0
    for i in paramdict.keys():
        dummy+=1
        i_dot=i.replace('_','.')
        if dummy == len(paramdict.keys()):
            configparams.write(i_dot+"="+str(round(paramdict[i],6)))
        else:
            configparams.write(i_dot+"="+str(round(paramdict[i],6))+',')
        ##fi
    ###i
    configparams.close()

    # Name the output files.
    outputfile=region+'_islands_out'
    outcatfile=region+'_islands_outcat'
    logfilename=region+'_islands.log'
    configfile="config.txt"

    #Construct the shell script to run the clumpfind command
    shellscript=open("clumpfinding.sh","w")
    shellscript.write('#!/bin/sh\n')
    shellscript.write("$STARLINK_DIR/bin/cupid/findclumps '"+sdffile+"' '"+outputfile+"' '"+outcatfile+"' "+method+' CONFIG=^'+configfile+' LOGFILE='+logfilename+' REPCONF=TRUE RMS='+str(noise)+' DECONV=FALSE SHAPE=Polygon WCSPAR=TRUE BACKOFF=TRUE')
    shellscript.close()

    #Run the clumpfind command
    subprocess.call("sh ./clumpfinding.sh", shell=True)

    #Clean up all the files produced and involved in running ClumpFind
    os.system('rm config.txt')
    os.system('rm clumpfinding.sh')
    os.system('mkdir island_info/')
    os.system('mv '+region+'*_islands* island_info/')

    ### Now find fragments.
    print('')
    print('')
    print('######')
    print('######')
    print('')
    print('Now finding FRAGMENTS...')
    print('')
    print('######')
    print('######')
    print('')
    print('')

    # Create the shellscript for running JSA_CATALOGUE
    shellscript_frag=open("fragfinding.sh","w")
    shellscript_frag.write('#!/bin/sh\n')
    shellscript_frag.write("${ORAC_DIR}/etc/picard_start.sh -nodisplay -log hs -verbose JSA_CATALOGUE "+sdffile)
    shellscript_frag.close()

    #Run the JSA_CATALOGUE command
    subprocess.call("sh ./fragfinding.sh", shell=True)

    #Clean up all the files produced
    os.system('rm fragfinding.sh')
    os.system('mkdir fragment_sdfs/')
    os.system('mv *peak.sdf fragment_sdfs/')
    os.system('rm *clump.sdf')
    os.system('rm *extent*png')
    os.system('rm disp.dat')
    os.system('rm log.group')
    os.system('rm rules.badobs')
    os.system('mkdir fragment_info/')
    os.system('mv *mask.sdf fragment_info/')
    os.system('mv *extent-cat.fits fragment_info/')
    os.system('mv *peak-cat.fits fragment_info/')
    os.system('mv *snr*.sdf fragment_info/')
    os.system('rm *temp2dmap.sdf')

    print('')
    print('You now have islands and fragments for '+region+'! :D')
    print('')
    print('')
#def








# #Get all the FIT files generated by clumpfind
# clumpinfofiles=os.listdir('.') #lists all files in the current directory
# rid=[]
# for i in range(len(clumpinfofiles)): #make sure we are only dealing with .FIT files
#   if clumpinfofiles[i][-1] != 'T':
#     rid.append(i)
# if len(rid)>0:
#   for j in rid:
#     clumpinfofiles[j]=0
#   for k in range(len(rid)):
#     k=0
#     clumpinfofiles.remove(0)
# for clumpinfofile in clumpinfofiles:
#   GetClumpCatalog(clumpinfofile,fitsfile,CATALOG='yes',DISTANCE=distance,TEMP=temp,BEAMSIZE=beamwidth,OPACITY=opacity,PIXEL_ARCSEC=pix_size_arcsecs,UNITS=units)

# def Beam2Pix(   value,
#                 pix_length,
#                 beamwidth,
#                 PIXTYPE='arcsecs',
#                 BEAMTYPE='arcsecs',
#                 DISTANCE=450):
#     '''
#     Beam2Pix:
#
#     Convert the pixel length and beamwidth to arcseconds.
#
#     Args:
#         values
#         pix_length
#         beamwidth
#         PIXTYPE
#         BEAMTYPE
#         DISTANCE
#
#     Returns:
#         jyperpixel
#     '''
#     #Convert pix_length and beamwidth to arcseconds
#     if PIXTYPE=='degrees':
#         pix_length=206264.806*pix_length*(np.pi/180.0)
#     if PIXTYPE=='radians':
#         pix_length=206264.806*pix_length
#     if PIXTYPE=='parsecs':
#         pix_length=(pix_length/DISTANCE)*206264.806
#     if BEAMTYPE=='degrees':
#         beamwidth=206264.806*beamwidth*(np.pi/180.0)
#     if BEAMTYPE=='radians':
#         beamwidth=206264.806*beamwidth
#     #TotalFlux=Peak*2*pi*int(x*exp(-x^2/(2*sigma^2))) -- int is evaluated from 0 -> inf so, by the substitution z=x/sigma (then m=z^2), our integral evaluates to 1 and we have:  TotalFlux=Peak*2*Pi*sigma^2=Peak*beam area
#     #Total/area = Peak flux of a pixel. Duh.
#     sigma=np.sqrt((-(beamwidth/2.0)**2.0)/(2.0*np.log(0.5)))
#     beam_area=2.0*np.pi*sigma**2.0
#     #Find the number of pixels in a beam
#     pixels_in_beam = beam_area/pix_length**2
#     #convert from jy/beam to jy/pixel
#     jyperpixel=value/pixels_in_beam
#     return jyperpixel
# #def
#
#     def pix2beam(   value,
#                     pix_length,
#                     beamwidth,
#                     PIXTYPE='arcsecs',
#                     BEAMTYPE='arcsecs',
#                     DISTANCE=450):
#     '''
#     Pix2Beam:
#
#     Converts the
#     '''
#     #Convert pix_length and beamwidth to arcseconds
#     if PIXTYPE=='degrees':
#      pix_length=206264.806*pix_length*(np.pi/180.0)
#     if PIXTYPE=='radians':
#      pix_length=206264.806*pix_length
#     if PIXTYPE=='parsecs':
#      pix_length=(pix_length/DISTANCE)*206264.806
#     if BEAMTYPE=='degrees':
#      beamwidth=206264.806*beamwidth*(np.pi/180.0)
#     if BEAMTYPE=='radians':
#      beamwidth=206264.806*beamwidth
#     #TotalFlux=Mean*2*pi*int(x*exp(-x^2/(2*sigma^2))) -- int is evaluated from 0 -> inf so, by the substitution z=x/sigma, our integral evaluates to 1 and we have:  TotalFlux=Mean*2*Pi*sigma^2=Mean*beam area
#     #Total/area = mean flux of a pixel. Duh.
#     sigma=np.sqrt((-(beamwidth/2.0)**2.0)/(2.0*np.log(0.5)))
#     beam_area=2.0*np.pi*sigma**2.0
#     #Find the number of pixels in a beam
#     pixels_in_beam = beam_area/pix_length**2
#     #convert from jy/pixel to jy/beam
#     jyperbeam=value/(1.0/pixels_in_beam)
#     return jyperbeam
# #def

# def GetClumpCatalog(clumpinfofile,
#                     fitsfile,
#                     CATALOG='yes',
#                     DISTANCE=450.0,
#                     TEMP=15.0,
#                     BEAMSIZE=14.6,
#                     OPACITY=0.012,
#                     PIXEL_ARCSEC=3.0,
#                     UNITS='MJYA2'):
#     '''
#     GetClumpCatalog:
#
#
#     '''
#
#     fitsheader=pyfits.getheader(fitsfile)
#
#     #Import all the information from each FIT file:
#     clumptable=atpy.Table(clumpinfofile)
#     clumpnumber=clumptable['PIDENT']
#     xcent=clumptable['Cen1']
#     ycent=clumptable['Cen2']
#     xpeak=clumptable['Peak1']
#     ypeak=clumptable['Peak2']
#     area=clumptable['Volume']
#     totalflux=clumptable['Sum']
#     size1=clumptable['Size1']
#     size2=clumptable['Size2']
#     peakflux=clumptable['Peak']
#     clumpinfo={}
#     catalogue=open(clumpinfofile.split('.FIT')[0]+'_derived.cat',"w")
#     catalogue.write('Clump Number        x_cent   y_cent     x_peak_cent    y_peak_cent       Radius (AU)       Mass (M_sun)      Mass/M_J     Concentration\n')
#     #Find the masses and stabilities like we did with Mairs et al. 2014.
#     if 'CDELT2' in fitsheader.keys():
#     cdelt2=float(fitsheader['CDELT2'])
#     pixel_arcsec=206264.806*float(fitsheader['CDELT2'])*(np.pi/180.0) #Arcseconds per pixel
#     else:
#     pixel_arcsec=PIXEL_ARCSEC
#
#
#     fitscolnames=['ID','x_cent','y_cent','x_peak','y_peak','rj','mass','mass_over_mj','conc']
#     fitscolformats=['D','D','D','D','D','D','D','D','D']
#
#     derived_dict={}
#
#     for i in range(len(fitscolnames)):
#       derived_dict[fitscolnames[i]]=[]
#
#     for i in range(len(area)):
#     #SARAH'S WAY OF CALCULATING M_J - IN SOLAR MASSES. The R_j she uses is the clump's effective radius - not the correct value, but all she had to go on.
#     r_j_pc=DISTANCE*(np.sqrt((area[i]/pixel_arcsec**2.0)/np.pi)*pixel_arcsec/206264.806) #450.0*(sqrt(area/pi) * pixel_arcsec /206264.806) #S = RTHETA, area[i]/pixel_arcsec**2.0 - turns AREA from arcsec^2 to pixels^2
#     r_j_au=r_j_pc*206264.806 #CONVERSION FACTOR IS AU IN A PARSEC DUH
#     m_j=1.9*(TEMP/10.0)*(r_j_pc/0.07)
#     #;;;;;;;;;;;;;;;;;
#     #SUM(T)=Peak*2*pi*int(x*exp(-x^2/(2*sigma^2))) -> int is evaluated from 0 -> inf so, by the substitution z=x/sigma, our integral evaluates to 1 and we have:  SumT=Peak*2*Pi*sigma^2 - the calculation of sigma is performed, based on a 20 arcsec beam, then the beam area is calculated so we can find the number of pixels in a beam
#     #since 0.5=exp(-(fwhm/2)^2/(2*sigma^2)):
#     #SO WE ARE JUST GETTING
#     if UNITS =='MJYA2':
#       coreflux=totalflux[i]*pixel_arcsec**2.0/1000.0
#       peakflux_converted=peakflux[i]*pixel_arcsec**2.0/1000.0
#     if UNITS == 'JYBM':
#       coreflux=beam2pix(float(totalflux[i]),pixel_arcsec,BEAMSIZE,DISTANCE=DISTANCE,PIXTYPE='arcsecs') #NOW IN JANSKY'S PER PIXEL
#       peakflux_converted=beam2pix(float(peakflux[i]),pixel_arcsec,BEAMSIZE,DISTANCE=DISTANCE,PIXTYPE='arcsecs') #NOW IN JANSKY'S PER PIXEL
#     if UNITS == 'JYPIX':
#       coreflux=float(totalflux[i]) #SUM(T) FROM THE LOG
#       peakflux_converted=float(peakflux[i])
#     coremass=0.074*coreflux*(DISTANCE/100.0)**2.0*(OPACITY/0.01)**(-1.0)*(np.exp(17.0/TEMP)-1) #IN SOLAR MASSES
#
#     concentration=1-(1.13*BEAMSIZE**2.0*coreflux/(area[i]*pix2beam(float(peakflux_converted[i]),pixel_arcsec,BEAMSIZE,DISTANCE=DISTANCE,PIXTYPE='arcsecs'))) #Area is in arcseconds squared
#
#     #GET PIXEL COORDINATES
#     xrefvalue=fitsheader['CRVAL1']
#     yrefvalue=fitsheader['CRVAL2']
#     xrefpix=fitsheader['CRPIX1']
#     yrefpix=fitsheader['CRPIX2']
#
#     xcentcoord=xrefpix-((xcent[i]-xrefvalue)*(np.pi/180.0)*206264.806)/pixel_arcsec
#     ycentcoord=yrefpix+((ycent[i]-yrefvalue)*(np.pi/180.0)*206264.806)/pixel_arcsec
#     xpeakcoord=xrefpix-((xpeak[i]-xrefvalue)*(np.pi/180.0)*206264.806)/pixel_arcsec
#     ypeakcoord=yrefpix+((ypeak[i]-yrefvalue)*(np.pi/180.0)*206264.806)/pixel_arcsec
#
#
#     derived_dict['ID'].append(clumpnumber[i])
#     derived_dict['x_cent'].append(xcentcoord)
#     derived_dict['y_cent'].append(ycentcoord)
#     derived_dict['x_peak'].append(xpeakcoord)
#     derived_dict['y_peak'].append(ypeakcoord)
#     derived_dict['rj'].append(r_j_au)
#     derived_dict['mass'].append(coremass)
#     derived_dict['mass_over_mj'].append(coremass/m_j)
#     derived_dict['conc'].append(concentration)
#
#
#     catalogue.write(str(clumpnumber[i])+'                    '+str(round(xcentcoord,2))+'   '+str(round(ycentcoord,2))+'       '+str(round(xpeakcoord,2))+'           '+str(round(ypeakcoord,2))+'         '+str(r_j_au)+'     '+str(coremass)+'   '+str(coremass/m_j)+'     '+str(round(concentration,2))+'\n')
#     clumpinfo[i]={'logfile':CATALOG,'index':clumpnumber,'peak1_coord':xpeakcoord,'peak2_coord':ypeakcoord,'cen1':xcentcoord,'cen2':ycentcoord,'size1':size1[i],'size2':size2[i],'sum':totalflux[i],'peak':peakflux[i],'area':area[i]}
#     catalogue.close()
#
#     #Make a fits file of the derived parameters like rj and mass
#
#     new_cols=[]
#     for eachcolumn in range(len(fitscolnames)):
#       colname=fitscolnames[eachcolumn]
#       coldata=derived_dict[fitscolnames[eachcolumn]]
#       colformat=fitscolformats[eachcolumn]
#       new_cols.append(pyfits.ColDefs([pyfits.Column(name=colname,format=colformat,array=coldata)]))
#     updatecolumns=new_cols[0]
#     for eachnewcol in np.array(range(len(new_cols)-1))+1:
#       updatecolumns=updatecolumns+new_cols[eachnewcol]
#     hdu = pyfits.BinTableHDU.from_columns(updatecolumns)
#     hdu.writeto(clumpinfofile.split('.FIT')[0]+'_derived.fits')
#
#     return clumpinfo
