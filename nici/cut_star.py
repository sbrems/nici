import numpy as np
import os
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table,vstack
from multiprocessing import Pool
from scipy.signal import fftconvolve
from .parallactic_angle import get_pa
from .gaussfitter import gaussfit
from .params import *
from .misc import find_max_star
from .correlate_torus import correlate_torus
import scipy.optimize as opt
from .doublegauss import shared_center



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
    test=np.where((indices[0]-int(np.ceil(image.shape[0]/2.)))**2+\
                  (indices[1]-int(np.ceil(image.shape[1]/2.)))**2 < 10**2)
    m=np.zeros_like(image)
    m[test]=1
    smoothed=fftconvolve(image,m,'same')
    y,x=np.where(smoothed==np.max(smoothed))
    return x[0], y[0]        

def get_rough_center(im,n=3,neg = 'neg'):
    '''Return the rough center of the most negative/positive points.
    If neg: Use a torus-correlated map to find the max'''
    if neg not in ['neg','pos','double']:
        raise ValueError('Value {} for neg unknown. Choose pos or neg'.format(neg))
    if neg in ['neg','double']:
        #yy,xx = np.unravel_index(im.flatten().argsort()[0:n],im.shape)
        im = correlate_torus(im,rin=2,rout=10)
    imsize = im.shape[-1]
    fitres = gaussfit(im,
                      params = (np.min(im), np.max(im)-np.min(im),
                                int(np.floor(imsize/2.)),int(np.floor(imsize/2.)),
                                imsize/4.,imsize/4.,0.))
                                
    x,y = fitres[2:4]
    x = int(round(x))
    y = int(round(y))
    #yy,xx = np.unravel_index(im.flatten().argsort()[-n::],im.shape)
    #y = int(round(np.median(yy)))
    #x = int(round(np.median(xx)))
    return y,x

def fit_gaussian(fitim,ceny,cenx,neg='pos'):
    '''Fit a negative gaussian to an array of 2xhalsize width. Use the median of the
    positions of the n most negative points'''
    if neg == 'neg':
        offset = np.max(fitim)
        amp  = -np.max(fitim)
    elif neg =='pos':
        offset = np.min(fitim)
        amp    = np.max(fitim)
    else:
        raise ValueError('Value {} for neg unknown. Choose pos,neg here'.format(neg))
    fitres = gaussfit(fitim,
                      params = (offset,amp,5,5,
                                3.,3.,0.))
    if neg and fitres[1] >= 0:
        UserWarning('WARNING! SHOULD BE A NEG GAUSSIAN BUT {} WAS FOUND AS\
AMPLITUDE!'.format(fitres[1]))
    if not neg and fitres[1] <= 0:
        UserWarning('WARNING! SHOULD BE A POS GAUSSIAN BUT {} WAS FOUND AS\
AMPLITUDE!'.format(fitres[1]))
        
    yxcen = fitres[2:4][::-1]
    return yxcen


def align_stars(indir,outdir,fluxtable=None,fflux=None,ncpu=6,keepfrac=0.7,
                maxhalfsize=np.inf,minhalfsize=0,flipx=False,fitmode='neg'):
    '''align the stars and cut the maximum possible square around them'''
    from scipy.ndimage.interpolation import shift
    if outdir == None:
        outdir=indir
    dic_fflux = {0:'sat',1:'flux'}
    print('\n****Cutting out {}-frames****\n'.format(dic_fflux[fflux]))
    pxhalf = 30 #px right and left of image. final size 2*pxhalf+1, planet at ~242px
    filetable = ascii.read(os.path.join(indir,'filetable_bkgrnd.csv'),delimiter=',')
    filetable = filetable[filetable['flux'] ==fflux]
    #only continue if there is at least one file
    if len(filetable) == 0:
        print('No files found for type {}'.format(dic_fflux[fflux]))
    else:
        filetable['PA'] = np.nan
        nimages = len(filetable)
        imdim = fits.getdata(filetable['fndewarped'][0]).shape[-1]
        #read in the images
        images = []
        full_images = []
        remove_ims = [] #which are too close to corner
        for ii,fn in enumerate(filetable['fndewarped']):
            data,head = fits.getdata(fn,header=True)
            full_images.append(data)
            if not len(data.shape) == 2:
                raise ValueError('Unknown data fromat at file %s'%fn)
            starx = int(filetable[ii]['roughx'])
            stary = int(filetable[ii]['roughy'])
    #        nims = data.shape[0]
            datacut = np.full([2*pxhalf+1,2*pxhalf+1],np.nan)
            data = data[max(0, stary-pxhalf) : min(stary+pxhalf+1, data.shape[1]),
                        max(0, starx-pxhalf) : min(starx+pxhalf+1, data.shape[0]),
                    ]
            filetable['PA'][ii] = get_pa(head,verbose=False)
            try:
                datacut[0:datacut.shape[0],0:datacut.shape[1]] = data
            except:
                print('Star too close to center. Ignoring image %s'%fn)
                remove_ims.append(ii)
                continue
            filetable['PA'][ii] = get_pa(head,verbose=False)
            images.append(datacut)
        full_images = np.array(full_images)
        filetable.remove_rows(remove_ims)
    #    filetable.write('filetable_removed_close_borders.csv',delimiter=',',overwrite=True)
        print('Ignored {} images due to wrongly placed star.'.format(len(remove_ims)))
        nimages = len(filetable)
        
        #/////////////////////////////////////////////////////////
        #median combine and first xreg. Keep orig images and shift
        #them only once at the end!
        #/////////////////////////////////////////////////////////
    
        print('register number getting medianed: ',len(images))
        print(images[0].size, images[0].shape, images[0].dtype)
    
        first_median=np.median(images, axis=0)
    
        pool=Pool(ncpu)
        get_shifts=subreg(first_median)
        first_shifts=pool.map(get_shifts,images)
        pool.close()
        shifted_images = []
        for hh in range(len(images)): 
            shifted_images.append(shift(images[hh], first_shifts[hh]))
    
     
        #/////////////////////////////////////////////////////////
        #keep only the best of images
        #/////////////////////////////////////////////////////////
        cross_reg=[]
        for im in shifted_images:
            cross_reg.append(np.sum((im-first_median)**2.))
            
        sorted_cross_reg=np.argsort(cross_reg)
        selected_cross_reg=sorted(sorted_cross_reg[0:int(keepfrac*len(images))])
        n_selected=len(selected_cross_reg)
        
        #/////////////////////////////////////////////////////////
        #median combine and second xreg
        #/////////////////////////////////////////////////////////
    
        images=np.array(images)[selected_cross_reg,:,:]
        shifted_images = np.array(shifted_images)[selected_cross_reg,:,:]
        second_median=np.median(shifted_images,axis=0)
    
        print('second subreg')
        pool=Pool(ncpu)
        get_shifts=subreg(second_median)
        #second shifts is the shift to get all at the same center. need to absolutely
        #center them later
        second_shifts=pool.map(get_shifts,images)
        pool.close()
        
        shifted_images =[]
        for hh in range(n_selected): 
            shifted_images.append(shift(images[hh], second_shifts[hh]))
        shifted_images = np.array(shifted_images)
        
        #get center for images. Move to pxhalf
        print('Finding the absolute center')
        yxoff = []
        smallsizehalf = 20
        fitsizehalf = 6
        nfitworked = 0
        for hh in range(n_selected):
    #        xycen.append( get_center(images[h,:,:]) )
            small_im = shifted_images[hh,pxhalf-smallsizehalf:pxhalf+smallsizehalf+1,
                                         pxhalf-smallsizehalf:pxhalf+smallsizehalf+1]
            yxroughcen = np.array(get_rough_center(small_im,n=3,neg = fitmode))
            yxroughcen += pxhalf -smallsizehalf
            fit_im = shifted_images[hh,yxroughcen[0]-fitsizehalf: 
                                       yxroughcen[0]+fitsizehalf+1,
                                       yxroughcen[1]-fitsizehalf:
                                       yxroughcen[1]+fitsizehalf+1]
            if fitmode == 'double':
                xx = range(fit_im.shape[1])
                yy = range(fit_im.shape[0])
                xx, yy = np.meshgrid(xx,yy)
                amp = np.max(fit_im) - np.min(fit_im)
                initial_guess = (fitsizehalf,fitsizehalf,#center all in x,y
                                 fitsizehalf,fitsizehalf,#sigpos
                                 fitsizehalf/2.,fitsizehalf/2.,#signeg
                                 amp, -amp/2.,0.,) #amppos,ampneg,offset
                try:
                    popt, pcov = opt.curve_fit(shared_center, (xx,yy), fit_im.flatten(),
                                               p0=initial_guess)
                    yxcen = popt[0:2][::-1]
                    nfitworked +=1
                    yxoff.append(yxroughcen + (yxcen-fitsizehalf) - pxhalf)
                except:
                    continue
                    
            elif fitmode in ['neg','pos']:
                yxcen = np.array(fit_gaussian(fit_im,fitsizehalf,fitsizehalf,
                                              neg=fitmode))
            else:
                raise ValueError('Fitmode {} unknown!'.format(fitmode))
            if fitmode != 'double':
                yxoff.append(yxroughcen + (yxcen -fitsizehalf) - pxhalf)
#        del shifted_images
#        del images
        if fitmode == 'double': 
            print('Fitting precise double gaussian used for absolute centering \
worked {} of {} times'.format(nfitworked,n_selected))
        yxcenshift = np.median(-np.array(yxoff),axis=0)
        print('General offset of images: %s' %yxcenshift)
        all_shifts  =second_shifts + yxcenshift
        
        #store the precise centers
        filetable['precisey'] = np.nan
        filetable['precisex'] = np.nan
        filetable['pxtoborder']=0
        for istar,tstar,thisshift in zip(selected_cross_reg,
                                     filetable[selected_cross_reg],
                                     all_shifts):
            filetable['precisey'][istar]  = tstar['roughy'] - thisshift[0]
            filetable['precisex'][istar]  = tstar['roughx'] - thisshift[1]
            filetable['pxtoborder'][istar]= int(np.floor(np.min(\
                                            (filetable['precisex'][istar],
                                             filetable['precisey'][istar],
                                             imdim - filetable['precisex'][istar],
                                             imdim - filetable['precisey'][istar]))-0.5))
        #////////////////////////////////////////////////////////
        #after we got the precise centers cut the maximum area
        #////////////////////////////////////////////////////////
        #determine the minimum distance to any border
        min_dist = np.min((filetable['pxtoborder'][selected_cross_reg]))
        if min_dist >= minhalfsize:
            if min_dist <= maxhalfsize:
                maxsize = min_dist
                print('Cutting images to size {}, which is the maximum without discarding \
images'.format(2*maxsize+1))
            else:
                maxsize = maxhalfsize
                print('Cutting images to selected size {}'.format(2*maxsize+1))
            selected_final = selected_cross_reg
            n_final = n_selected
        else:
            maxsize = np.min(filetable['pxtoborder'][filetable['pxtoborder'] >=minhalfsize])
            selected_border = np.where(filetable['pxtoborder'] > maxsize)[0]
            selected_final = [ii for ii in selected_cross_reg if ii in selected_border]
            n_final = len(selected_final)
            print('Cutting images to chosen size {}, which means {} more images are \
being discarded.'.format(2*maxsize+1,n_selected-n_final))
        final_ims = np.full((n_final,2*maxsize+1,2*maxsize+1),np.nan)
    
        #1. round the center 2.shift the image by residuum 3.cut floor of min dist
        for ii,tstar in zip(range(n_final),filetable[selected_final]):
            inty = int(np.round(tstar['precisey']))
            intx = int(np.round(tstar['precisex']))
            tempim = shift(full_images[selected_final[ii],:,:],
                           (inty - tstar['precisey'],intx-tstar['precisex']))
            final_ims[ii,:,:]=tempim[inty -maxsize: inty +maxsize+1,
                                     intx -maxsize: intx +maxsize+1]
        #/////////////////////////////////////////////////////////
        #save
        #/////////////////////////////////////////////////////////
        filetable['selected_final'] = 0
        filetable['selected_final'][selected_final] = 1
        if flipx:
            final_ims = final_ims[:,:,::-1]
            filetable['PA'] = -filetable['PA']
        if fflux == 1:
            fits.writeto(os.path.join(outdir,'median_unsat.fits'),np.median(final_ims,axis=0),
                         overwrite=True)
            fits.writeto(os.path.join(outdir,'cube_unsat.fits'),final_ims,overwrite=True)
            print('Saved median images to {}'.format(outdir))
            return filetable
        else:
            if fluxtable != None:
                filetable = vstack([fluxtable,filetable])
    
            fits.writeto(os.path.join(outdir,'center_im.fits'),final_ims,overwrite=True)
            fits.writeto(os.path.join(outdir,'rotnth.fits'),filetable['PA'][selected_final],
                         overwrite=True)
            ascii.write(filetable,output=os.path.join(indir,'filetable_stars.csv'),delimiter=',',
                        overwrite=True)
    
        print('Cut stars. SAVED CENTERED IMAGES AND ANGLES IN {}.\
Their shape is {}'.format(outdir,final_ims.shape))
