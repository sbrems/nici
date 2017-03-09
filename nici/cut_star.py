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

def get_center(image):
    indices=np.indices(image.shape)
    test=np.where((indices[0]-150)**2+(indices[1]-150)**2<10**2)
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
        if not len(data.shape) == 3:
            raise ValueError('Unknown data fromat at file %s'%fn)
        starx = filetable[ii]['roughx']
        stary = filetable[ii]['roughy']
        nims = data.shape[0]
        datacut = np.full([nims,2*pxhalf+1,2*pxhalf+1],np.nan)
        data = data[:,
                       min(0, stary-pxhalf) : max(stary+pxhalf+1, data.shape[1]),
                       min(0, starx-pxhalf) : max(starx+pxhalf+1, data.shape[0]),
                ]
        datacut[:,0:data.shape[1],0:data.shape[2]]
        PAs.append(nims*[list(head[hpa])])
        for hh in range(nims):
            images.append = data[hh,:,:]
    
    
    
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
    selected_cross_reg=sorted_cross_reg[0:int(keepFrac*len(images))]
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

    #get center for images. Should be pxhalf
    xycen = []
    for h in range(n_selected):
        xycen.append( get_center(images[h,:,:]) )
    yxcenshift = np.median(xycen,axis=0)[::-1] - [pxhalf,pxhalf]
            
    for h in range(n_selected): 
        shifts[h] += xycenshift
        images[h,:,:]=shift(images[h,:,:], shifts[h])
    
    #/////////////////////////////////////////////////////////
    #save
    #/////////////////////////////////////////////////////////
    images = np.stack(images,axis=0)
    PAs_sel = PAs[selected_cross_reg]
    filet_sel = filetable[selected_cross_reg]
    filet_sel['orig_nr'] = selected_cross_reg
    fits.writeto(outdir+'center_im_sat.fits',images)
    fits.writeto(outdir+'rotnth.fits')

    print('DONE WITH ALL. SAVED DEROTATED IMAGES AND ANGLES IN %s. Their shape is %'%(outdir,images.shape))
