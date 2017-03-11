import os
from astropy.io import fits
import numpy as np

def read_fits(directory,verbose=True,only_first=False):
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
    for file in files:
        if file.endswith('.fits'):
            form = fits.getdata(directory+file).shape
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
    for file in files:
        if file.endswith('.fits'):
            hdulist = fits.open(directory+file)
            if len(hdulist) ==3:
                header = hdulist[0].header+hdulist[1].header
            elif len(hdulist) ==2:#standard
                header = hdulist[0].header
            data = hdulist[1].data
            if (len(data.shape) == 3) & (only_first):
                data = data[0,:,:]
            if len(data.shape) == 3 :#image cube
                all_data[n:n+data.shape[0],:,:] = data
                headers.extend(header for ii in xrange(data.shape[0]))
                filenames.extend(file for ii in xrange(data.shape[0]))
                n += data.shape[0]
            elif len(data.shape) == 2: #one image
                all_data[n,:,:] = data
                headers.append(header)
                filenames.append(file)
                n += 1
            else:
                raise ValueError('Fits file has unknown format!')

    return filenames, all_data,headers
