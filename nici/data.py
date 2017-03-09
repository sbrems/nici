import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from misc import read_fits
from parameters import *



def find_quad(frame,bpm):
    ylen,xlen = frame.shape
    threshold = 0.9*np.median(frame) + 3*np.std(frame)
    test = np.logical_and(image > threshold, bpm == 0)
    y_high_counts, x_high_counts = test.nonzero()
    
    if np.median(y_high_counts) > np.floor(ylen/2.):
        side = 't'
    elif np.median(y_high_counts)  < np.floor(ylen/2.):
        side = 'b'
    else: #just at center
        side = 'c'
    if np. median(x_high_counts) > np.floor(xlen/2.):
        side += 'r'
    elif np.median(x_high_counts)> np.floor(xlen/2.):
        side += 'l'
    else: #just at center
        side += 'c'

    try:
        quad = side2quad(side)
    except:
        quad = 5

    return quad


def return_quad(datadir,hexpt,ddata,dTexp,fdata,fTexp,bpm):
    '''Give datadir. Returns astropy table containing filenames and the position of the
    star (0:topleft,1:bottomleft,2:bottomright,3:topright). 5 if not found.'''
    fnames,data,headers = read_fits(datadir,only_first=True)
    nframes = data.shape[-2]
    data_med = np.median(data,axis=0)
    #checking for same exptime. Small sanity check and getting obstime
    exptimes = []
    obstimes = []
    for expt,obst in headers[hexpt],headers['DATE-OBS']+headers['TIME-OBS']:
        exptimes.append(expt)
        obstimes.append(Time(obst))
    exptime = np.unique(exptimes)
    if len(exptime) > 1:
        raise ValueError('Different Exposure times: (%s) in Science frames encountered.\n'+\
                         'This is not included!'%np.unique(exptimes))
    else:
        exptime = exptime[0]

    quad = []
    roughy = []
    roughx = []
    for iframe in nframes:
        fnames[iframe] = datadir+fnames[iframe]
        data[iframe,:,:] -= data_med
        iquad,iy,ix = find_quad(data[iframe,:,:],bpm)
        quad.append(itquad)
        roughy.append(iy)
        roughx.append(ix)
        

    table = Table([fnames,quad,obstimes,roughx,roughy],\
                  names =('fname','quad','obstime','roughx','roughy'))

    return table.sort('obstime')

def determine_neighbours(table):
    '''Get the table with the quads. For each one, give two values of the
    neighbouring observations that should be used for background 
    subtraction. Returns the table'''
    obsbef = []
    obsaft = []
    for iobs in range(len(table)):
        tbefore = table[0:iobs]
        try:
            obsbefore = tbefore[ tbefore != table[iobs]['quad'] ][-1]['quad']
        except:
            obsbefore= None
        tafter = table[iobs:]
        try:
            obsafter = tafter[ tafter != table[iobs]['quad'] ][0]['quad']
        except:
            obsafter = obsbefore
        if obsbefore == None:
            obsbefore = obsafter
            if obsafter == None:
                raise ValueError('Cant sort observations for bkgrnd subtraction'\+
                                 'Check obsnr %s!'%iobs)
        obsbef.append(obsbefore)
        obsaft.append(obsafter)
    
    table['obsbefore'] = obsbef
    table['obsafter']  = obsafter
    
    return table

def flatfield_bkgrnd(fdata,fTexp,filetable,bpm,intermdir):
    print('Flatfielding')
    allsdata = []
    allsheader=[]
    allsmed = []
    fninterm = []
    nfiles = len(filetable)
    print('Flatfieldiing the %d files'%nfiles)
    for ifile in range(nfiles):
        sdata,sheader = fits.getdata(filetable['fname'][ifile],header=True)
        allsdata.append(divide_flat(sdata,sheader[hexpt],fdata,fTexp))
        allsheader.append(sheader)
        allsmed.append(np.median(sdata,axis=0))
    print('Creating fixpix map')
    bad_and_neighbors = fix_pix.find_neighbors(bpm,n=4)
    print('Removing background and fixing pixels for %d files'%nfiles)
    for ifile in range(nfiles):
        print('Doing bgrnd+fixpix for %d / %d'%(ifile,nfiles))
        for iim in range(allsdata[ifile].shape[0]):
            fn = filetable['fname'][ifile].split('/')[-1]
            sdata =  allsdata[ifile][iim,:,:] -\
                     (allsmed[filetable['obsbefore']] +\
                      allsmed[filetable['obsafter']]) /2.
            fix_pix.correct_with_precomputed_neighbors(sdata,bad_and_neighbors)
            fninterm.append(intermdir+fn)
    filetable['fninterm'] = fninterm

    ascii.write(filetable,output=intermdir+'filetable_bkgrnd.csv',delimiter=',')
    print('Done Removing background. Stored results in %s'%intermdir)
#    return filetable
    
