import os
from astropy.io import fits
import numpy as np
from .params import *

def read_fits(directory,verbose=True,only_first=False,cam=None):
    '''This routine reads all fits files data into a big data cube and all header files
    into a big header cube. The order is the same and is alphabetically.Filenames and headers
    are multiple if it was an imagecube. So all have the same length. So it returns:
    (fits)filenames,datacube,headercube.
    If you give the path of a fits-file,only this is read out.
    Also compatible with NICI data (2 cubes, 3headers)'''
    #to avoid buffer overflows we need the number images first, which is not the same as
    #filenumber, as some images are cubes and some are not
    n_images = 0
    form = []
    if directory.endswith('.fits'):
        files= [os.path.split(directory)[-1]]
        directory = os.path.dirname(directory)+'/'
    else:
        files = sorted(os.listdir(directory))
    for fl in files:
        if fl.endswith('.fits'):
            form = fits.getdata(os.path.join(directory,fl)).shape
            if len(form) == 3 : #image cube
                if only_first:
                    n_images += 1
                else:
                    n_images += form[0]
            elif len(form) == 2: #one image
                n_images += 1
            else:
                raise ValueError('Fits file has unknown format!')
    if verbose: print('Found ',n_images,' frames in ',directory)
    #now make the array
    filenames = []
    headers = []
    all_data = np.full((n_images,form[-2],form[-1]),np.nan,dtype=np.float64) #float16 for memory
    n = 0
    for fl in files:
        if fl.endswith('.fits'):
            hdulist = fits.open(os.path.join(directory,fl))
            #check hdu shape
            if cam != None:
                if len(hdulist) !=3:#nici special                
                    raise ValueError('Unknown file format for NICI camera!file:%s'%fl)
                header = hdulist[0].header+hdulist[cam2head[cam]].header
                data = hdulist[cam2head[cam]].data
            else:
                if len(hdulist) !=1: raise ValueError('Unknown file format at file %s'%fl)
                header = hdulist[0].header
                data = hdulist[0].data
            #check data shape
            if (len(data.shape) == 3) & (only_first):
                data = data[0,:,:]
            if len(data.shape) == 3 :#image cube
                all_data[n:n+data.shape[0],:,:] = data
                headers.extend(header for ii in range(data.shape[0]))
                filenames.extend(fl for ii in range(data.shape[0]))
                n += data.shape[0]
            elif len(data.shape) == 2: #one image
                all_data[n,:,:] = data
                headers.append(header)
                filenames.append(fl)
                n += 1
            else:
                raise ValueError('Fits file has unknown format!')

    return filenames, all_data,headers

def find_max_star(image,smaller_by_factor = 16):
    '''Median smooth an image and find the max pixel.
    The median smoothing helps filter hot pixels and 
    cosmic rays. The median is taken by using bin_median
    with a smaller_by_factor=16'''
    image[np.isnan(image)]=np.median(image[~np.isnan(image)])
    binned=bin_median(image,smaller_by_factor=smaller_by_factor)
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
    for i in range(smaller_by_factor):
        for j in range(smaller_by_factor):
            sub_arrs0.append(arr[i::smaller_by_factor,j::smaller_by_factor])
    sub_arrs=[s[:sub_arrs0[-1].shape[0]-1,:sub_arrs0[-1].shape[1]-1] for s in sub_arrs0]
    if returnStd:
        #zip truncates each sub-arr at the length of the minimum length subarr.
        #this ensures every bin has the same number of datapoints, but throws
        #away data if the last bin doesn't have a full share of datapoints.
        return np.median(sub_arrs,axis=0),np.std(sub_arrs,axis=0)
    else:
        return np.median(sub_arrs,axis=0)


