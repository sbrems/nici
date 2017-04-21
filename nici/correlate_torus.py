import numpy as np
from scipy.signal import fftconvolve
from .misc import find_max_star

def correlate_torus(im,rin=2,rout=7):
    '''Correlate a torus to an image'''
    if len(im.shape) !=2:
        raise ValueError('Image has unknown shape!')
    if rin >= rout:
        raise ValueError('rin: {} >= rout: {}'.format(rin,rout))
    #make a torus map
    ny = im.shape[0]
    nx = im.shape[1]
    ceny = np.int(np.floor(ny/2.))
    cenx = np.int(np.floor(nx/2.))
    
    y,x = np.ogrid[-ceny:ny-ceny, -cenx:nx-cenx]
    mask = x*x + y*y <= rout*rout
    mask[x*x + y*y <= rin*rin] = False
    
    #CC the torus
    cor = fftconvolve(im,mask,mode='same')
    return cor
