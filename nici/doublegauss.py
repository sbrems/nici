import numpy as np

def shared_center(xy,x0,y0,sigposx,sigposy,signegx,signegy,amppos,ampneg,offset):
    '''fit a positive and a negative gauss with the same center'''
    x,y =xy
    x0 = float(x0)
    y0 = float(y0)    
    gpos = amppos * np.exp(-( (x-x0)**2/2./sigposx**2 + (y-y0)**2/2./sigposy**2 ))
    gneg = ampneg * np.exp(-( (x-x0)**2/2./signegx**2 + (y-y0)**2/2./signegy**2 ))
    g    = offset + gpos + gneg
    
    return g.ravel()

