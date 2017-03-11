import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from multiprocessing import Pool
from scipy.ndimage.interpolation import shift
from scipy.signal import fftconvolve
from gaussfitter import gaussfit
from params import *



class subreg:
    def __init__(self,reference):
        self.reference=reference
    def __call__(self,im):
        kernel=self.reference[::-1,::-1]
        cor = fftconvolve(im,kernel,mode='same')
        y,x=find_max_star(cor)
        g=gaussfit(cor[max(0, y-40):min(y+40, cor.shape[0]),
                       max(0, x-40):min(x+40, cor.shape[1])])
        shiftx=np.rint(cor.shape[1]/2.) - max(0,x-40)-g[2]
        shifty=np.rint(cor.shape[0]/2.) - max(0,y-40)-g[3]
        #shifts.append((shifty,shiftx))
        return (shifty,shiftx)
def bin_median(arr,smaller_by_factor=1,returnStd=False):
    '''bin an array arr by creating super pixels of size
    smaller_by_factor*smaller_by_factor and taking the median. 
    Can optionally also return the standard deviation within
    the super-pixels.
    INPUTS:
    arr: 2d array
    smaller_by_factor: integer
    returnStd: bool, default False
    RETURNS binned_array OR binned_array, bin_std'''
    sub_arrs0=[]
    for i in xrange(smaller_by_factor):
        for j in xrange(smaller_by_factor):
            sub_arrs0.append(arr[i::smaller_by_factor,j::smaller_by_factor])
    sub_arrs=[s[:sub_arrs0[-1].shape[0]-1,:sub_arrs0[-1].shape[1]-1] for s in sub_arrs0]
    if returnStd:
        #zip truncates each sub-arr at the length of the minimum length subarr.
        #this ensures every bin has the same number of datapoints, but throws
        #away data if the last bin doesn't have a full share of datapoints.
        return np.median(sub_arrs,axis=0),np.std(sub_arrs,axis=0)
    else:
        return np.median(sub_arrs,axis=0)

def find_max_star(image):
    '''Median smooth an image and find the max pixel.
    The median smoothing helps filter hot pixels and 
    cosmic rays. The median is taken by using bin_median
    with a smaller_by_factor=16'''
    image[np.isnan(image)]=np.median(image[~np.isnan(image)])
    binned=bin_median(image,smaller_by_factor=16)
    y,x=np.transpose((binned==binned.max()).nonzero())[0]
    y*=16
    x*=16
    while True:
        x0=max(x-15,0)
        y0=max(y-15,0)
        patch = image[y0:min(y+15,image.shape[0]),x0:min(x+15,image.shape[1])]
        dy,dx=np.transpose(np.where(patch==patch.max()))[0]
        dy-=(y-y0)
        dx-=(x-x0)
        y=min( max(y+dy,0), image.shape[0]-1)
        x=min( max(x+dx,0), image.shape[1]-1)
        if (dx==0) and (dy==0):
            break
    return y,x


def get_center(image):
    indices=np.indices(image.shape)
    test=np.where((indices[0]-int(np.ceil(image.shape[0]/2.)))**2+\
                  (indices[1]-int(np.ceil(image.shape[1]/2.)))**2 < 10**2)
    m=np.zeros_like(image)
    m[test]=1
    smoothed=fftconvolve(image,m,'same')
    y,x=np.where(smoothed==np.max(smoothed))
    return x[0], y[0]        


def align_stars(indir,outdir=None,ncpu=1,keepfrac=0.7):
    if outdir == None:
        outdir=indir
    pxhalf = 100 #px right and left of image. final size 2*pxhalf+1
    filetable = ascii.read(indir+'filetable_bkgrnd.csv',delimiter=',')
    nimages = len(filetable)
    #read in the images
    PAs = []
    images = []
    for ii,fn in enumerate(filetable['fninterm']):
        data,head = fits.getdata(fn,header=True)
        if not len(data.shape) == 2:
            raise ValueError('Unknown data fromat at file %s'%fn)
        starx = int(filetable[ii]['roughx'])
        stary = int(filetable[ii]['roughy'])
#        nims = data.shape[0]
        datacut = np.full([2*pxhalf+1,2*pxhalf+1],np.nan)
        data = data[max(0, stary-pxhalf) : min(stary+pxhalf+1, data.shape[1]),
                    max(0, starx-pxhalf) : min(starx+pxhalf+1, data.shape[0]),
                ]
        datacut[0:data.shape[0],0:data.shape[1]] = data
        PAs.append(head[hpa])
        images.append(datacut)
#        for hh in range(nims):
#            images.append( data[hh,:,:])
    
    
    
    #/////////////////////////////////////////////////////////
    #median combine and first xreg
    #/////////////////////////////////////////////////////////

    print('register number getting medianed: ',len(images))
    print(images[0].size)
    print(images[0].shape)
    print(images[0].dtype)
    first_median=np.median(images, axis=0)


    pool=Pool(ncpu)
    get_shifts=subreg(first_median)
    shifts=pool.map(get_shifts,images)
    first_shifts = shifts
    pool.close()

    for hh in range(len(images)): 
        images[hh]=shift(images[hh], shifts[hh])

 
    #/////////////////////////////////////////////////////////
    #keep only the best of images
    #/////////////////////////////////////////////////////////
    cross_reg=[]
    for im in images:
        cross_reg.append(np.sum((im-first_median)**2.))
        
    sorted_cross_reg=np.argsort(cross_reg)
    selected_cross_reg=sorted_cross_reg[0:int(keepfrac*len(images))]
    n_selected=len(selected_cross_reg)

    
    #/////////////////////////////////////////////////////////
    #median combine and second xreg
    #/////////////////////////////////////////////////////////

    images=np.array(images)[selected_cross_reg,:,:]
    second_median=np.median(images,axis=0)

    print 'second subreg'
    pool=Pool(ncpu)
    get_shifts=subreg(second_median)
    shifts=pool.map(get_shifts,images)
    second_shifts = shifts
    pool.close()

    #get center for images. Move to pxhalf
    xycen = []
    for h in range(n_selected):
#        xycen.append( get_center(images[h,:,:]) )
        xycen.append( gaussfit(images[h,pxhalf-5:pxhalf+5,
                                        pxhalf-5:pxhalf+5],\
                               params=(0.,-200.,6.,6.,3.,3.,0.))[2:4] - [5.,5.])
    yxcenshift = np.median(xycen,axis=0)[::-1] #- [pxhalf,pxhalf]
    print('General offset of images: %s' %yxcenshift)
    for h in range(n_selected): 
        shifts[h] += yxcenshift
        images[h,:,:]=shift(images[h,:,:], shifts[h])
    
    #/////////////////////////////////////////////////////////
    #save
    #/////////////////////////////////////////////////////////
    images = np.stack(images,axis=0)
    PAs = np.array(PAs)
    PAs_sel = PAs[selected_cross_reg]
    filet_sel = filetable[selected_cross_reg]
    filet_sel['orig_nr'] = selected_cross_reg
    fits.writeto(outdir+'center_im_sat.fits',images,clobber=True)
    fits.writeto(outdir+'rotnth.fits',PAs,clobber=True)

    print('DONE WITH ALL. SAVED DEROTATED IMAGES AND ANGLES IN %s. Their shape is %s'%(outdir,images.shape))
