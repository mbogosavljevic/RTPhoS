#!/usr/bin/python
# M. Bogosavljevic, AOB, March 2015
# working on automation  of offsets.py by Z. Ioannou 
# to be used in the sequential mode when fits files are coming in one by one
# *** work in progress ****
"""
--- Zach:
This routine gets a list of FITS files that belong to a timeseries photometry run
and calculates the pixel offsets of each file in relation to the first file.

The offsets are calculated by cross correlating a blured masked of the reference
frame (usually the first image in a sequence) with a blured mask of the following
image. The blured masks are necesseary to make the stars appear as prominent as
possible without any background noise, hot pixels or cosmic rays.

"""
def zach_offsets(dataref,data2):

  import pyfits
  from scipy import signal, ndimage
  import numpy

  xshift = 0
  yshift = 0

  # Crop the image by 10 pixels on each side
  xsize1 = dataref.shape[0]
  ysize1 = dataref.shape[1]

  xstart = 9
  xend   = xsize1-10	
  ystart = 9
  yend   = ysize1-10

  croped1 = dataref[xstart:xend,ystart:yend]
  median1 = pyfits.np.median(croped1)
  
  # Create image 1 mask and blur it using a Gausian filter of
  # FWHM of 1 pixel. Then any pixel with a value less than 100
  # is set to zero.

  # WARNING - This is completely arbitrary but usually
  # pixels that correspond to actual stars will have values
  # a lot greater than 100.

  mask1   = croped1
  mask1[mask1 < 1.5*median1] = 0.0
  blured1 = ndimage.gaussian_filter(croped1, sigma=1)
  mask1   = blured1
  mask1[mask1 < 100.0] = 0.0
  
  # Get the size of the croped masked arrays
  xsize = mask1.shape[0]
  ysize = mask1.shape[1]
  
  # Create collapsed image arrays for reference image dataref
  xvals1 = mask1.sum(axis=0)
  yvals1 = mask1.sum(axis=1)

  # Crop the data2 image by 10 pixes on each side
  xsize2 = data2.shape[0]
  ysize2 = data2.shape[1]
  xstart = 9
  xend   = xsize2-10
  ystart = 9
  yend   = ysize2-10
  croped2 = data2[xstart:xend,ystart:yend]
  median2 = pyfits.np.median(croped2)
    	
  # Create image 2 mask and blur it
  mask2   = croped2
  mask2[mask2 < 1.5*median2] = 0.0
  blured2 = ndimage.gaussian_filter(croped2, sigma=1)
  mask2   = blured2
  mask2[mask2 < 100.0] = 0.0
    
  # Create collapsed image arrays for image 2
  xvals2=mask2.sum(axis=0)
  yvals2=mask2.sum(axis=1)
  
  # Calculate the x and y shift of the image in pixels using cross correlation
  xshift = (numpy.argmax(signal.correlate(xvals1,xvals2)))-(ysize-1)
  yshift = (numpy.argmax(signal.correlate(yvals1,yvals2)))-(xsize-1)

  return (xshift, yshift)

if __name__ == "__main__":

  import sys

  # check if path to watch ends with '/'
  path_to_watch = sys.argv[1]
  if not(path_to_watch.endswith('/')):
         path_to_watch = path_to_watch + '/'
