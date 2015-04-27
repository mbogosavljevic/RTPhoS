#!/usr/bin/python
####################################
## M. Bogosavljevic, SQU, Jan 21 2015
## Plot a radial profile and return fwhm
## Intended to be used as a DS9 analysis task
## Make a get_fwhm.ans file with the following 4 lines which are between ---- :
## ----
## Radial Profile
## *
## bind r
## get_fwhm.py $xpa_method | $text
## ----
## This will bind the r key to launch this code from DS9. 
## The code expects that there is a selected circle region
## It will center it, and plot its radial profile, using the region size as radius
## It will also display a plot, with FWHM in title, and wait for the plot to be closed 
## This later part is a "feature" of matplotlib objects, and can be changed only
## if one includes fork process and such.

import pyregion
import pyfits
import ds9
import numpy as np
import matplotlib.pyplot as plt
from   scipy.optimize import curve_fit

##############################################################################

# Define model function to be used to fit to the data:
# in this case its a Gauss with a  center Xc, sigma, amplitude A and base level B
def gauss(x, *p):
    A, xc, sigma, B  = p
    return A*np.exp(-(x-xc)**2/(2.*sigma**2)) + B

##############################################################################

def fwhm_radial_plot(image):
#   Calculate the azimuthal radial profile.
#   Input image should be a small stamp from the 2D image (already cropped)

    # Calculate the indices from the image
    y, x = np.indices(image.shape)
    center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])
    r = np.hypot(x - center[0], y - center[1])
    
    # select only elements inclosed in a circle of radius r, not a box of side 2r
    mr = np.amax(r)/ 2.0**0.5
    take = np.where(r <= mr)
    cr  = r[take]
    im2 = image[take]

    # Get sorted radii
    ind = np.argsort(cr.flat)
    r_sorted = cr.flat[ind]
    i_sorted = im2.flat[ind]

    # mirror the curve 
    neg_r = -1 * r_sorted[::-1]
    neg_i = i_sorted[::-1]
    fully = np.concatenate([neg_i,i_sorted])
    fullx = np.concatenate([neg_r,r_sorted])

    # Fit a function 
    # p0 is the initial guess for the fitting coefficients
    # (A, xc, sigma, B in above in gauss function)
    # so make some educated guesses:
    n = len(fullx)
    sigma = 1.
    base = np.mean(i_sorted)
    peak = np.amax(fully) - base

    p0 = [peak, 0., sigma, base]
    coeff, var_matrix = curve_fit(gauss, fullx, fully, p0=p0)

    # convert sigma to fwhm:
    fwhm = 2.3548 * coeff[2] 
    # format it as 3 decimapl place string
    fwhms = "{:.3f}".format(fwhm)

    # create a more dense fake array for plotting the fit
    r2 = max(r_sorted)
    fakex = np.arange(0,r2,0.1)
    fit = gauss(fakex, *coeff)
  
    # plot setup
    plt.ylabel('Intensity')
    plt.xlabel('Pixel')
    plt.title('Radial profile, fwhm='+fwhms)

    # plot fit in red dashed
    plt.plot(fakex, fit, 'r--')
    # plot data in blue circles
    plt.plot(r_sorted, i_sorted, 'bo')
    plt.show()
    # must wait for window to close, since window is interactive

    return fwhm, coeff

#############################################################################

def get_fwhm(xpapoint):
# given a xpapoint handle
# find the selected region (should be only one)
# centroid it, and calculate the fwhm of the object inside
# the circular radius is set by the region object itself

    # create a ds9 object linked with an XPA point
    win = ds9.ds9(xpapoint)

    # load the image from ds9 (not from disk!)
    hdu_link = win.get_pyfits() 
    image = hdu_link[0].data

    # initally, target should be the only selected region
    # ceontroid 50 times (DS9 algo for centroiding is not great)
    for x in range(0, 49):
        win.set("regions centroid")
    targetrows = win.get("regions selected") 
    target = pyregion.parse(targetrows)
    xt = target[0].coord_list[0]
    yt = target[0].coord_list[1]
    r  = target[0].coord_list[2]


    print "Xcen Ycen", xt, yt

    # check if target is length of one
    if len(target)!=1:
       print(" Select one region (only)! ")
       raise Exception(" Select one region (only)! ")
    #else:
    #   print('%-15s %7.2f %7.2f' % ("Target (X,Y):", xt, yt) )

    x1 = xt - r
    x2 = xt + r
    y1 = yt - r
    y2 = yt + r
   
    # x and y are inverted??? found this for HTCas example image
    crop_image = image[y1:y2,x1:x2]
    #pyfits.writeto("temp.fits", crop_image, clobber=1)

    fwhm, coefs = fwhm_radial_plot(crop_image)
    ## these output valeues can be used later on for something
    print('%-6s %7.3f' % ('FWHM: ', fwhm))

if __name__ == "__main__":

    import sys
    xpaname = sys.argv[1]
    print " Calculating centroid and FWHM "
    get_fwhm(xpaname)

# Earlier version of this code had been taking the radius as input
# now it determines the radius from the size of the 
# circle region which is currently selected in DS9
#    radius = float(sys.argv[2])
#    get_fwhm(xpaname,radius)
