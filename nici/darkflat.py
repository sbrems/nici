import numpy as np
from misc import read_fits


def sort(direct,hexpt):
    '''returns the median of all same exposures and a list with the exposure times (sorted)'''
    fns,data,header = read_fits(direct)
    exptimes = []
    for ihead in header:
        exptimes.append(ihead[hexpt])
    uexpts = sorted(np.unique(exptimes))
    darks = np.full(np.hstack((len(uexpts),data[0,:,:].shape)),np.nan)
    for ii,uexpt in enumerate(uexpts):
        darks[ii,:,:] = np.nanmedian( data[np.where(exptimes==uexpt)[0],:,:] , axis=0)
        
    return darks,uexpts


def subtract_dark(data,dataTexp,ddata,dTexp):
    '''Subtracts darks using the nearest dark exposure time'''
    for ii,dtime in enumerate(dataTexp):
        dtuse = min(dTexp, key = lambda x:abs(x-dtime))
        idtuse = np.where(dTexp == dtuse)[0][0]
        data[ii,:,:] -= ddata[idtuse,:,:]
    
    return data

def divide_flat(data,dataTexp,fdata,fTexp):
    for ii, time in enumerate(dataTexp):
        ftuse = min(dataTexp, key= lambda x:abs(x-time))
        iftuse  = np.where(np.array(fTexp) == ftuse)[0][0]
        if len(data.shape) == 3:
            for iframe in range(data.shape[0]):
                data[iframe,:,:] /= fdata[iftuse,:,:]
        elif len(data.shape)==2:
            data /= fdata[iftuse,:,:]
        else:
            raise ValueError('Data has unknown shape')
    
    return data
