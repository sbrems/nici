import numpy as np
import os.path
from .misc import read_fits
from astropy.stats import sigma_clip
from astropy.io import fits
from .badpixel import make_masterflat_bpm

def masterflat_lamp(flatdir,intermdir,cam,hexpt,wavefolder,exptime=0.38,sigma = 3):
    print('Making BPM and Masterflat for cam {}'.format(cam))
    if wavefolder not in ['H20LBLOCK']:
        print('Assuming Images from {} use lampflats'.format(wavefolder))
        
        fns,data,header = read_fits(flatdir,cam=cam)
        ilampoff = []
        ilampon  = []
        if wavefolder == 'LPCH4H1S':
            #Note: the blue flats will be overwritten later, as all frames are dark there :(
            prefix = 'S20120401S'
            dic_fnrs = {'red' :{'bright':['0019','0020','0021','0022'],
                                'dark'  :['0023','0024','0025','0026']},
                        'blue':{'dark'  :['0283','0284','0285','0286'],
                                'bright':['0287','0288','0289','0290']},
                    }
            
            lampon = []
            lampoff= []
            fns,data,header = read_fits(flatdir,cam=cam)
            for fn,dat in zip(fns,data):
                if fn in [prefix+fnr+'.fits' for fnr in dic_fnrs[cam]['bright']]:
                    lampon.append(dat)
                elif fn in [prefix+fnr+'.fits' for fnr in dic_fnrs[cam]['dark']]:
                    lampoff.append(dat)
            medlampon = np.median(np.array(lampon),axis=0)
            medlampoff= np.median(np.array(lampoff),axis=0)
            fits.writeto(os.path.join(intermdir,'medlampon_'+cam+'.fits'),medlampon,overwrite=True)
            fits.writeto(os.path.join(intermdir,'medlampoff_'+cam+'.fits'),medlampoff,overwrite=True)
        else:
            for ii,head in enumerate(header):
                if head[hexpt] == exptime:
                    gcalshut =  head['GCALSHUT'].strip().upper()
                    if 'OPEN' == gcalshut:
                        ilampon.append(ii)
                    elif 'CLOSED' == gcalshut:
                        ilampoff.append(ii)
                    else:
                        print('WARNING! UNKNOWN PROPERTY {} IN GCALSHUT! IS LAMP \
ON OR OFF? CHECK!'.format(gcalshut))
            medlampoff = np.median(data[ilampoff,:,:],axis=0)
            medlampon  = np.median(data[ilampon,:,:],axis=0)
        
        gradmap = medlampon - medlampoff
        bpm = np.full(medlampon.shape,0,dtype=np.int64)
        bpm[sigma_clip(gradmap,sigma=sigma).mask] = 1
         
        masterflat = gradmap / np.median(gradmap[bpm ==0])
        masterflat[bpm ==1] = 1
    if wavefolder == 'H20LBLOCK':
        print('Processing band {}. Assuming skyflats')
        masterflat,bpm = make_masterflat_bpm(flatdir,cam=cam,sigma=sigma)
    if cam =='blue' and wavefolder =='LPCH4H1S':
         print('\nTHERE WERE NO BRIGHT FLATS FOR CAM {}. \
 ONLY MAKING ROUGHT BPM, NO FLAT!'.format(cam))
         fullmedian = np.median(np.concatenate((np.array(lampoff),np.array(lampon)),axis=0),axis=0)
         bpm = np.zeros(medlampon.shape,dtype=np.int64)
         masterflat = np.full(medlampon.shape,1.,dtype=np.int64)
         bpm[np.abs(fullmedian) >= 1000.] = 1   
         fits.writeto(os.path.join(intermdir,'medlampon_'+cam+'.fits'),medlampon,overwrite=True)
         fits.writeto(os.path.join(intermdir,'medlampoff_'+cam+'.fits'),medlampoff,overwrite=True)
        
        
                

    fits.writeto(os.path.join(intermdir,'bpm_'+cam+'.fits'),bpm,overwrite=True)
    fits.writeto(os.path.join(intermdir,'masterflat_'+cam+'.fits'),masterflat,
                 overwrite=True)

    
    return masterflat,bpm
    

def sort(direct,hexpt,cam=None):
    '''returns the median of all same exposures and a list with the exposure times (sorted)'''
    fns,data,header = read_fits(direct,cam=cam)
    exptimes = []
    for head in header:
        exptimes.append(head[hexpt])
    uexpts = sorted(np.unique(exptimes))
    darks = np.full(np.hstack((len(uexpts),data[0,:,:].shape)),np.nan)
    for ii,uexpt in enumerate(uexpts):
        darks[ii,:,:] = np.nanmedian( data[np.where(exptimes==uexpt)[0],:,:] , axis=0)
    print('Found exposure times {} in {}'.format(uexpts,direct))
    import ipdb;ipdb.set_trace()
    return darks,uexpts


def subtract_dark(data,dataTexp,ddata,dTexp,verbose=True):
    '''Subtracts darks using the nearest dark exposure time'''
    for ii,dtime in enumerate(dataTexp):
        dtuse = min(dTexp, key = lambda x:abs(x-dtime))
        idtuse = np.where(dTexp == dtuse)[0][0]
        if (dtuse/dtime < 10) and (dtuse/dtime > 0.1):
            if len(data.shape) == 3:
                for iframe in range(data.shape[0]):
                    data[iframe,:,:] /= fdata[iftuse,:,:]
            elif len(data.shape)==2:
                data /= fdata[iftuse,:,:]
            else:
                raise ValueError('Data has unknown shape')
            data[ii,:,:] -= ddata[idtuse,:,:]
        else:
            if verbose: print('Not removing dark as no good exposure time was found:\
Darkexpts: {}, Dataexpts: {}'.format(dTexp,dataTexp))

    return data

def divide_flat(data,dataTexp,fdata,fTexp,verbose=True):
    for ii, time in enumerate(dataTexp):
        ftuse = min(fTexp, key= lambda x:abs(x-time))
        iftuse  = np.where(np.array(fTexp) == ftuse)[0][0]
        if (ftuse/time < 10) and (ftuse/time > 0.1):
            if len(data.shape) == 3:
                for iframe in range(data.shape[0]):
                    data[iframe,:,:] /= fdata[iftuse,:,:]
            elif len(data.shape)==2:
                data /= fdata[iftuse,:,:]
            else:
                raise ValueError('Data has unknown shape')
        else:
            if verbose: print('Not flatfielding as no good exposure time was found:\
Flatexpts: {}, Dataexpts: {}'.format(fTexp,dataTexp))
    return data
