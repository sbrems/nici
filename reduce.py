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

def do(waveband,target=None,date=None,pardir = None, datadir =None,darkdir=None,flatdir=None,
       ncpu=10,keepfrac=0.7):
    '''A Pipeline to reduce NICI data. Provide darks in darkdir, flats in flatdir.
    Assuming dithering, so taking the neighbouring two images with a different dither
    position and interpolating between the pixels to to get the flat.
    BP: Using flats for BP detection: Values changing differently than the others are
    being replaced by median of good neighbours. Taken from Jordans pipeline.'''
    wavebands = ['KH','LCH4H1S']
    waveband = waveband.upper()
    if waveband not in wavebands: raise ValueError('Your waveband %s not known. Choose from %s'\
                                                   %(waveband,wavebands))
    if target == None: target = 'HD97048'
    if date   == None: date ='20120401'
    if pardir == None: pardir  = '/disk1/brems/NICI/'
    if datadir== None: datadir = pardir+date+'/'+target+'/'+waveband+'/'
    if darkdir== None: darkdir = pardir+date+'/darks/'+waveband+'/'
    if flatdir== None: flatdir = pardir+date+'/flats/'+waveband+'/'
    intermdir = datadir+'interm/'
    finaldir =  datadir+'results/'
    for direct in [intermdir,finaldir]:
        if not os.path.exists(direct):
            os.mkdir(direct)
    fnbpm = 'bpm_'+waveband+'.fits'

    #getting and sorting darks and flats
    ddata,dTexp = darkflat.sort(darkdir,hexpt)
    fdata,fTexp = darkflat.sort(flatdir,hexpt)
    fdata = darkflat.subtract_dark(fdata,fTexp,ddata,dTexp)
    
    #making the bpm

    if os.path.exists(pardir+fnbpm):
        print('Found bpm in %s'%(pardir+fnbpm))
        bpm = fits.getdata(pardir+fnbpm)
    else:
        bpm = badpixel.make_bpm(flatdir)
        fits.writeto(pardir+fnbpm,bpm)
    
    #finding the star positions and time of observation
    filetable = data.return_quad(datadir,hexpt,ddata,dTexp,fdata,fTexp,bpm)
    #get the right observation for background subtraction based on positions
    filetable = data.determine_neighbours(filetable)
    #getting science data. Flatfielding, bkgrnd and fixpixing it and saving in intermdir
    data.flatfield_bkgrnd(fdata,fTexp,filetable,bpm,intermdir)
##############################################
    
    #cut out the star, align them and get the paralactic angles
    cut_star.align_stars(indir=intermdir,outdir=finaldir,ncpu=ncpu,keepfrac=keepfrac)
