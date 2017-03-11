import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii
from misc import read_fits
import darkflat
import fix_pix
from params import *



def find_quad(frame,bpm):
    ylen,xlen = frame.shape
    threshold = 0.9*np.median(frame) + 3*np.std(frame)
    test = np.logical_and(frame > threshold, bpm == 0)
    y_high_counts, x_high_counts = test.nonzero()
    ypos =np.median(y_high_counts)
    xpos = np.median(x_high_counts)
    if ypos > np.floor(ylen/2.):
        side = 't'
    elif ypos  < np.floor(ylen/2.):
        side = 'b'
    else: #just at center
        side = 'c'
    if xpos > np.floor(xlen/2.):
        side += 'r'
    elif xpos < np.floor(xlen/2.):
        side += 'l'
    else: #just at center
        side += 'c'

    try:
        quad = side2quad[side]
    except:
        quad = 5#at center in one direct

    return quad,ypos,xpos


def return_quad(datadir,hexpt,ddata,dTexp,fdata,fTexp,bpm):
    '''Give datadir. Returns astropy table containing filenames and the position of the
    star (0:topleft,1:bottomleft,2:bottomright,3:topright). 5 if not found.'''
    fnames,data,headers = read_fits(datadir,only_first=True)
    nframes = data.shape[0]
    data_med = np.median(data,axis=0)
    #checking for same exptime. Small sanity check
    exptimes = []
#    obstimes = []
    for head in headers:
        exptimes.append(head[hexpt])
#        obstimes.append(Time(head['DATE-OBS']+head['TIME-OBS']))
    exptime = np.unique(exptimes)
    if len(exptime) > 1:
        raise ValueError('Different Exposure times: (%s) in Science frames encountered.\n'+\
                         'This is not included!'%np.unique(exptimes))
    else:
        exptime = exptime[0]

    quad = []
    roughy = []
    roughx = []
    for iframe in range(nframes):
        fnames[iframe] = datadir+fnames[iframe]
        data[iframe,:,:] -= data_med
        iquad,iy,ix = find_quad(data[iframe,:,:],bpm)
        quad.append(iquad)
        roughy.append(iy)
        roughx.append(ix)
        

    table = Table([fnames,quad,roughx,roughy],\
                  names =('fname','quad','roughx','roughy'))
    table.sort('fname')
    table['fnumber'] = range(len(table))
    return table

def determine_neighbours(table):
    '''Get the table with the quads. For each one, give two values of the
    neighbouring observations that should be used for background 
    subtraction. Returns the table'''
    obsbef = []
    obsaft = []
    for iobs in range(len(table)):
        tbefore = table[0:iobs]
        try:
            obsbefore = tbefore[ tbefore['quad'] != table[iobs]['quad'] ]['fnumber'][-1]
        except:
            obsbefore= None
        tafter = table[iobs:]
        try:
            obsafter = tafter[ tafter['quad'] != table[iobs]['quad'] ]['fnumber'][0]
        except:
            obsafter = obsbefore
        if obsbefore == None:
            obsbefore = obsafter
            if obsafter == None:
                raise ValueError('Cant sort observations for bkgrnd subtraction'+\
                                 'Check obsnr %s!'%iobs)
        obsbef.append(obsbefore)
        obsaft.append(obsafter)
    
    table['obsbefore'] = obsbef
    table['obsafter']  = obsafter
    
    return table

def flatfield_bkgrnd(fdata,fTexp,filetable,bpm,intermdir):
    allsdata = []
    allsheader=[]
    allsmed = []
    fninterm = []
    nfiles = len(filetable)
    print('Flatfieldiing the %d files'%nfiles)
    for ifile in range(nfiles):
        fn,sdata,sheader = read_fits(filetable['fname'][ifile],verbose=False)
        sdata = sdata[0]
        sheader= sheader[0]
        allsdata.append(darkflat.divide_flat(sdata,[sheader[hexpt]],fdata,fTexp))
        allsheader.append(sheader)
#        allsmed.append(np.median(sdata,axis=0)) images already med
    print('Creating fixpix map')
    bad_and_neighbors = fix_pix.find_neighbors(bpm,n=4)
    sdata_red = []
    print('Removing background and fixing pixels for %d files'%nfiles)
    for ifile in range(nfiles):
        if ifile%10==0:print('Doing bgrnd+fixpix for %d / %d files'%(ifile,nfiles))
        if len(allsdata[ifile].shape) ==3:
            sdata = []
            for iim in range(allsdata[ifile].shape[0]):
                fn = filetable['fname'][ifile].split('/')[-1]
                sdata.append(allsdata[ifile][iim,:,:] -\
                         (allsdata[filetable['obsbefore'][ifile]] +\
                          allsdata[filetable['obsafter'][ifile]]) /2.)
                fix_pix.correct_with_precomputed_neighbors(sdata,bad_and_neighbors)#inplace
            sdata = np.array(sdata)
        elif len(allsdata[ifile].shape) ==2:
            fn = filetable['fname'][ifile].split('/')[-1]
            sdata =  allsdata[ifile] - \
                     (allsdata[filetable['obsbefore'][ifile]] +\
                     allsdata[filetable['obsafter'][ifile]]) /2.
            fix_pix.correct_with_precomputed_neighbors(sdata,bad_and_neighbors)#inplace
        fninterm.append(intermdir+fn)
        fits.writeto(intermdir+fn,sdata,header=allsheader[ifile],clobber=True)
    filetable['fninterm'] = fninterm
    ascii.write(filetable,output=intermdir+'filetable_bkgrnd.csv',delimiter=',')
    print('Done Removing background. Stored results in %s'%intermdir)
#    return filetable
    
