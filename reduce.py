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


def do(waveband,target=None,date=None,pardir = None, datadir =None,darkdir=None,ncpu=10,
       keepfrac=0.7):
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
    if pardir == None: pardir  = '/home/sbrems/NICI/'
    if datadir== None: datadir = pardir+'/'+date+'/'+target+'/'+waveband+'/'
    if darkdir== None: darkdir = pardir+'/'+date+'/'+darks+'/'+waveband+'/'
    if flatdir== None: flatdir = pardir+'/'+date+'/'+flats+'/'+waveband+'/'
    intermdir = datadir+'interm/'
    finaldir =  datadir+'results/'

    #getting and sorting darks and flats
    ddata,dTexp = darkflat.sort(darkdir,hexpt)
    fdata,fTexp = darkflat.sort(flatdir,hexpt)
    fdata = darkflat.subtract_dark(fdata,fTexp,ddata,dTexp)
    
    #making the bpm
    if os.exists(flatdir+fnbpm):
        bpm = fits.getdata(flatdir+fnbpm)
    else:
        bpm = badpixel.make_bpm(fdata)
        fits.writeto(flatdir+fnbpm)
    
    #finding the star positions and time of observation
    filetable = data.return_quad(datadir,hexpt,ddata,dTexp,fdata,fTexp,bpm)
    #get the right observation for background subtraction based on positions
    filetable = data.determine_neighbours(filetable)
    #getting science data. Flatfielding, bkgrnd and fixpixing it and saving in intermdir
    data.flatfield_bkgrnd(fdata,fTexp,filetable,bpm,intermdir)
##############################################
    
    #cut out the star, align them and get the paralactic angles
    cut_star.align_stars(indir=intermdir,outdir=finaldir,ncpu=ncpu,keepfrac=keepfrac)
