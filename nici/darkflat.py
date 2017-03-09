import numpy as np
from misc import read_fits


def sort(darkdir,hexpt):
    '''returns the median of all same exposures and a list with the exposure times (sorted)'''
    fns,data,header = read_fits(darkdir)
    exptimes = []
    for expt in header[hexpt]:
        exptimes.append(expt)
    uexpts = sorted(np.unique(expt))
    darks = np.full(np.hstack(len(uexpts),data[0,:,:].shape),np.nan)
    for ii,uexpt in enumerate(uexpts):
        darks[ii,:,:] = np.median( data[np.where(exptimes==uexpt),:,:] , axis=0)
        
    return darks,uexpts


def subtract_dark(data,dataTexp,ddata,dTexp):
    '''Subtracts darks using the nearest dark exposure time'''
    for ii,dtime in enumerate(dataTexp):
        dtuse = min(dTexp, key = lambda x:abs(x-time))
        idtuse = np.where(dTexp == dtuse)[0]
        data[ii,:,:] -= ddata[idtuse,:,:]
    
    return data

def divide_flat(data,dataTexp,fdata,fTexp):
    for ii, time in enumerate(dataTexp):
        ftuse = min(dTexp, key= lambda x:abs(x-time))
        iftuse  = np.where(fTexp == ftuse)[0]
        if len(data.shape) == 3:
            for iframe in range(data.shape[0]):
                data[iframe,:,:] /= fdata[idtuse,:,:]
        elif len(data.shape)==2:
            data /= fdata[idtuse,:,:]
        else:
            raise ValueError('Data has unknown shape')
    
    return data
