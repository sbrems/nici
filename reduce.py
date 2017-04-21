import os
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import register_translation
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.interpolation import shift
from nici.params import *
from nici import misc,fix_pix,data,darkflat,cut_star,badpixel

import ipdb

def do(wavebands,target=None,date=None,pardir = None, datadir =None,darkdir=None,flatdir=None,
       ncpu=10,keepfrac=0.8,fkeepfrac = 0.9,minhalfsize=None):
    '''A Pipeline to reduce NICI data. Provide darks in darkdir, flats in flatdir.
    Assuming dithering, so taking the neighbouring two images with a different dither
    position and interpolating between the pixels to to get the flat.
    BP: Using flats for BP detection: Values changing differently than the others are
    being replaced by median of good neighbours. Taken from Jordans pipeline.'''
    known_wavebands = ['LP','CH4H1S','H20L','BLOCK','K','H']#red,blue
    dic_fitmode={'CH4H1S':'pos',
                 'LP':'double',
                 'K':'double',
                 'H':'neg',
                 'H20L':'pos',
    }
    dic_minhalfsize = {'CH4H1S':240,
                       'LP':244,
                       'K':0,
                       'H':250,
                       'H20L':243,
    }
    wavebands = [i.upper() for i in wavebands]
    wavefolder = ''.join(wavebands)
    dic_hexpt = {'red':'ITIME_R','blue':'ITIME_B'}
    exptime = 0.38 #only use images with this exptime
    wb2cam = {wavebands[0]:'red',wavebands[1]:'blue'}#always usses two wavebands
    print('Starting reduction. Assuming band %s was taken with camera %s and %s with camera %s'\
          %(wavebands[0],wb2cam[wavebands[0]],wavebands[1],wb2cam[wavebands[1]]))
    if wavebands[1] == 'BLOCK':
        wavebands = [wavebands[0]]
    for waveband in wavebands[::-1]:
        if minhalfsize == None:
            minhalfsize = dic_minhalfsize[waveband]
        cam = wb2cam[waveband]
        if cam == 'blue':
            print('Going to flip x-axis as recommended for cam {}'.format(cam))
            flipx = True
        else: 
            flipx = False
        hexpt = dic_hexpt[cam]
        print('\n Processing waveband %s from camera %s'%(waveband,cam))
        if waveband not in known_wavebands: raise ValueError('Your waveband %s not known. Choose from %s'\
                                                   %(waveband,known_wavebands))   
        if target == None: target = 'HD97048'
        if date   == None: date ='20120401'
        if pardir == None: pardir  = '/disk1/brems/NICI/'
        if datadir== None: datadir = os.path.join(pardir,date,target,'raw',wavefolder)
#        if darkdir== None: darkdir = os.path.join(pardir,date,'darks')
        if flatdir== None: flatdir = os.path.join(datadir,'flats')
        outdir = os.path.join(pardir,date,target,'reduced',waveband)
        intermdir = os.path.join(outdir,'interm')
        finaldir =  outdir
        dewarpdir = os.path.join(intermdir,'dewarped')
        for direct in [outdir,intermdir,finaldir,dewarpdir]:
            if not os.path.exists(direct):
                os.makedirs(direct)
        fnbpm = 'bpm_'+cam+'.fits'
    
        #apply header corrections and distortion from #Hayward+2014 before!
        print('\nCAUTION! ASSUMING YOU ALREADY FIXED THE HEADERS FOLLOWING E.G. HAYWARD+2014!\n')
        #getting and sorting darks and flats
#        ddata,dTexp = darkflat.sort(darkdir,hexpt,cam=cam) #no darks needed
        fnmflat = os.path.join(intermdir,'masterflat_'+cam+'.fits')
        fnbpm   = os.path.join(intermdir,'bpm_'+cam+'.fits')
        if os.path.exists(fnmflat) and os.path.exists(fnbpm):
            print('Using found masterflat and bpm in {} and {}'.format(fnmflat,fnbpm))
            masterflat = fits.getdata(fnmflat)
            bpm        = fits.getdata(fnbpm)
        else:
            masterflat,bpm = darkflat.masterflat_lamp(flatdir,intermdir,cam,hexpt,wavefolder,
                                                      exptime=exptime,sigma=3)
#        fdata = darkflat.subtract_dark(fdata,fTexp,ddata,dTexp)

        #finding the star positions and time of observation
        filetable = data.return_quad(datadir,hexpt,bpm,cam)
        #determine which ones are the fluxframes and separate them
        filetable = data.ident_fluxframes(filetable)
        #get the right observation for background subtraction based on positions
        filetable = data.determine_neighbours(filetable)
        #getting science data. Flatfielding, bkgrnd and fixpixing it and saving in intermdir
        data.flatfield_bkgrnd(masterflat,filetable,bpm,intermdir,cam=cam)
        ##############################################
        #dewarp the images
        data.dewarp(cam,pardir,indir=intermdir,outdir=dewarpdir)
        #cut out the star, align them and get the paralactic angles. First for fluxframes,
        #then for sat frames
        fitmode = dic_fitmode[waveband]
        fluxtable = cut_star.align_stars(dewarpdir,finaldir,fflux=1,ncpu=ncpu,keepfrac=fkeepfrac,
                                         maxhalfsize=130,flipx=flipx,fitmode='pos')
        cut_star.align_stars(dewarpdir,finaldir,fluxtable=fluxtable,fflux=0,minhalfsize=minhalfsize,
                             ncpu=ncpu,keepfrac=keepfrac,flipx=flipx,fitmode =fitmode)
        
    print('DONE WITH ALL FOR {}'.format(waveband))
