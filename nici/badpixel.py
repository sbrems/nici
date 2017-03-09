import numpy as np
from astropy.stats import sigma_clip


def make_bpm(flats,sigma=3):
    '''Assuming the flats are sorted. E.g. in increasing or decreasing intensity. Then clip the
    Values which are of. Using for loops. So slow. Returns a matrix which is 1 whre a bp is found.'''
    bpm = np.full(flats.shape[0,:,:],0).astype(int)
    gradmap = np.full(flats.shape[0,:,:],np.nan)
    nflats = flats.shape[0]
    iflats = np.linspace(0,nflats-1,nflats).astype(int)
    for yy in flats.shape[1]:
        for xx in flats.shape[2]:
            values = sigma_clip(flats[:,yy,xx],sig=3) #filter cosmic rays
            gradmap[yy,xx] = np.polyfit(iflats[~values.mask],\
                                        values[~values.mask],deg=1)[1]

    return bpm[~sigma_clip(gradmap,sig=sigma).mask] = 1
    
    
