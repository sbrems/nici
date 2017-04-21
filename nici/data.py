import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii
import os
from .misc import read_fits
from .darkflat import divide_flat
from .fix_pix import find_neighbors,correct_with_precomputed_neighbors
from .params import *



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


def return_quad(datadir,hexpt,bpm,cam):
    '''Give datadir. Returns astropy table containing filenames and the position of the
    star (0:topleft,1:bottomleft,2:bottomright,3:topright). 5 if not found.'''
    fnames,data,headers = read_fits(datadir,only_first=True,cam=cam)
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
    ignore_iims =[]
    for iframe in range(nframes):
        fnames[iframe] = os.path.join(datadir,fnames[iframe])
        data[iframe,:,:] -= data_med
        iquad,iy,ix = find_quad(data[iframe,:,:],bpm)
        if not np.isfinite(ix):
            ignore_iims.append(iframe)
        else:
            quad.append(iquad)
            roughy.append(iy)
            roughx.append(ix)

    print('Not using the following {} images, as no star was found:'.format(len(ignore_iims)))
    for iim in sorted(ignore_iims,reverse=True):
        print(fnames[iim])
        del fnames[iim]

    table = Table([fnames,quad,roughx,roughy],\
                  names =('fname','quad','roughx','roughy'))
    table.sort('fname')
    table['fnumber'] = range(len(table))
    return table

def ident_fluxframes(filetable):
    '''Find the fluxframes and also group them if other frames were taken inbetween'''
    #group the observations if change from flux to sat and back
    filetable['obsblock'] = -1
    filetable['flux'] = -1
    nflux = 0
    nsat  = 0
    obsblock = 0
    wasflux = None
    for ifile,fn in enumerate(filetable['fname']):
        head = fits.getheader(fn)
        if (head['NDFW'].upper() == 'OPEN') and \
           (head['LPSTATE'].upper() !='IDLE'):
            nsat += 1
            filetable['flux'][ifile] = 0
            if wasflux ==False or wasflux ==None:
                wasflux = False
            else:
                obsblock += 1
                wasflux = False
        else:
            nflux += 1
            filetable['flux'][ifile] = 1
            if wasflux ==True or wasflux == None:
                wasflux = True
            else:
                obsblock += 1
                wasflux = True
        filetable['obsblock'][ifile] = obsblock
    print('Found {} fluxframes and {} saturated in {} blocks'.format(nflux,nsat,obsblock))

    #sanitycheck
    if np.min(filetable['obsblock']) < 0 or np.min(filetable['flux']) < 0:
        raise ValueError('obsblock or fluxframe/satframe could not be determined correctly')

    return filetable


def determine_neighbours(filetable):
    '''Get the table with the quads. For each one, give two values of the
    neighbouring observations that should be used for background 
    subtraction. Returns the table.'''
    filetable = filetable.group_by('obsblock')
    obsbef = []
    obsaft = []
    for tablegr in filetable.groups:
        for iobs in range(len(tablegr)):
            tbefore = tablegr[0:iobs]
            try:
                obsbefore = tbefore[ tbefore['quad'] != tablegr[iobs]['quad'] ]['fnumber'][-1]
            except:
                obsbefore= None
            tafter = tablegr[iobs:]
            try:
                obsafter = tafter[ tafter['quad'] != tablegr[iobs]['quad'] ]['fnumber'][0]
            except:
                obsafter = obsbefore
            if obsbefore == None:
                obsbefore = obsafter
                if obsafter == None:
                    raise ValueError('Cant sort observations for bkgrnd subtraction'+\
                                     'Check obsnr %s!'%iobs)
            obsbef.append(obsbefore)
            obsaft.append(obsafter)
        
    filetable['obsbefore'] = obsbef
    filetable['obsafter']  = obsaft
    
    return filetable

def flatfield_bkgrnd(masterflat,filetable,bpm,intermdir,cam=None,verbose=True):
    allsdata = []
    allsheader=[]
    allsmed = []
    fninterm = []
    nfiles = len(filetable)
    print('Flatfielding the %d files'%nfiles)
    for ifile in range(nfiles):
        fn,sdata,sheader = read_fits(filetable['fname'][ifile],verbose=False,cam=cam)
        sdata = sdata[0]
        sheader= sheader[0]
        allsdata.append(sdata/masterflat)
        allsheader.append(sheader)
        verbose=False
#        allsmed.append(np.median(sdata,axis=0)) images already med
    print('Creating fixpix map')
    bad_and_neighbors = find_neighbors(bpm,n=4)
    sdata_red = []
    print('Removing background and fixing pixels for %d files'%nfiles)
    for ifile in range(nfiles):
        if ifile%10==0:print('Doing bgrnd+fixpix for %d / %d files'%(ifile,nfiles))
        if len(allsdata[ifile].shape) ==3:
            sdata = []
            for iim in range(allsdata[ifile].shape[0]):
                fn = os.path.split(filetable['fname'][ifile])[-1]
                sdata.append(allsdata[ifile][iim,:,:] -\
                         (allsdata[filetable['obsbefore'][ifile]] +\
                          allsdata[filetable['obsafter'][ifile]]) /2.)
                correct_with_precomputed_neighbors(sdata,bad_and_neighbors)#inplace
            sdata = np.array(sdata)
        elif len(allsdata[ifile].shape) ==2:
            fn = os.path.split(filetable['fname'][ifile])[-1]
            sdata =  allsdata[ifile] - \
                     (allsdata[filetable['obsbefore'][ifile]] +\
                     allsdata[filetable['obsafter'][ifile]]) /2.
            correct_with_precomputed_neighbors(sdata,bad_and_neighbors)#inplace
        fninterm.append(os.path.join(intermdir,fn))
        fits.writeto(os.path.join(intermdir,fn),sdata,header=allsheader[ifile],overwrite=True)
    filetable['fninterm'] = fninterm
    ascii.write(filetable,output=os.path.join(intermdir,'filetable_bkgrnd.csv'),delimiter=',',
                overwrite=True)
    print('Done Removing background. Stored results in %s'%intermdir)
#    return filetable
    
def dewarp(cam,pardir,indir=None,outdir=None):
    '''Do the dewarping following the NICI campaign paper 
    ( http://iopscience.iop.org/article/10.1086/679508/pdf ). This uses IDL,
    so make sure it is installed. Wrapper is pidly'''
    import pidly
    
    fnwarp = os.path.join(pardir,'niciwarp_'+cam+'.sav')
    print('Dewarping using IDL. Assuming dewarp file is {}'.format(fnwarp))
    filetable = ascii.read(os.path.join(indir,'filetable_bkgrnd.csv'),delimiter=',')
    filetable['fndewarped'] = 'error_in_dewarping_the_filename'*3 #make sure string is long enough
    nfiles = len(filetable)
    print('Dewarping the {} files'.format(nfiles))
    idl = pidly.IDL()
    idl('restore,"'+fnwarp+'"')
    for ifile in range(nfiles):
        if (ifile % 10) ==0: 
            #somehow IDL crashes after ~25 images. restart it more often
            print('Dewarping file {} / {} files'.format(ifile,nfiles))
            #if (ifile != 0) and (ifile != nfiles-1): idl.close()   
            #idl = pidly.IDL()
            #idl('restore,"'+fnwarp+'"')
        fn = filetable['fninterm'][ifile]
        fnout = os.path.join(outdir,os.path.split(fn)[-1])
        filetable['fndewarped'][ifile] = fnout
        idl('fits_open,"'+fn+'", fcb')
        idl('fits_read, fcb, im, hdr')
        idl('fits_close,fcb')
        idl('im_w = poly_2d(im, kx, ky, 2, cubic=-0.5)')
        idl('fits_open,"'+fnout+'",fout,/WRITE')
        idl('fits_write,fout,im_w,hdr')
        idl('fits_close,fout')
        
    idl.close()
    ascii.write(filetable,output=os.path.join(outdir,'filetable_bkgrnd.csv'),delimiter=',',
                overwrite=True)
