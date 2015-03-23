#!/usr/bin/python
# M. Bogosavljevic, AOB, March 2015
### WORK IN PROGRESS!!!! ###
### Contains 
### gauss : just a gaussian function definision for PSF fitting
### fwhm_from_star : fit the gauss and returns FWM (in pixels)
### get_comps_fwhm : do the above fit for a number of comparison stars,
###                  returns average (in pixels)
### run_rtphos : connects to ds9 via xpapoint, gets regions
###              figures out comparisons and targets
### seekfits:
# run_rtphos.ans needs some improvement so that the text output goes to terminal in real time
#RUN RTPHOS
#*
#menu
#run_rtphos.py $xpa_method

import pyregion
import pyfits
import ds9
import numpy as np
import os, time
from   scipy.optimize import curve_fit

##############################################################################
def namesplit(s, leader, trailer):
    end_of_leader = s.index(leader) + len(leader)
    start_of_trailer = s.index(trailer, end_of_leader)
    out = s[end_of_leader:start_of_trailer]
    return out

# Define model function to be used for PSF fit to the stars:
# in this case its a Gauss with a center Xc, sigma, amplitude A and base level B
def gauss(x, *p):
    A, xc, sigma, B  = p
    return A*np.exp(-(x-xc)**2/(2.*sigma**2)) + B

##############################################################################

def fwhm_from_star(image):
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
    fwhm = int(2.3548 * coeff[2] * 100) / 100.0
    return fwhm

#############################################################################

def get_comps_fwhm(comparisons, xpapoint):

    # create a ds9 object linked with an XPA point
    win = ds9.ds9(xpapoint)
    # load the image from ds9 (not from disk!)
    hdu_link = win.get_pyfits() 
    image = hdu_link[0].data

    nc = len(comparisons)
    tot_fwhm = 0.
    # comparisons is alist of tuples, in each tuple first element has x,y,r
    print "I will use", nc, " stars for FWHM"
    for i in range(0,nc):
        xt = comparisons[i][0][0]
        yt = comparisons[i][0][1]
        r  = comparisons[i][0][2]
        x1 = xt - r
        x2 = xt + r
        y1 = yt - r
        y2 = yt + r   
        # caution - I don't know why order of x and y is inverted here
        crop_image = image[y1:y2,x1:x2]
        fwhm = fwhm_from_star(crop_image)
        sname = (comparisons[i][1])
        cname = namesplit(sname,'{','}')
        print i+1, cname, fwhm
        tot_fwhm = tot_fwhm + fwhm
        
    mean_fwhm = tot_fwhm / nc
    return mean_fwhm

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

def seekfits(dataref,imdir,tsleep, comparisons, targets, psf_fwhm):

    before = dict ([(f, None) for f in os.listdir (imdir)])
    
    # counter needed just to print out progress
    count = 0.
    
    try:
      while 1:
           count = count + 1
           print int(count), "*LISTENING*", imdir, time.strftime('%X %x %Z')
           after = dict ([(f, None) for f in os.listdir (imdir)])
           added = [f for f in after if not f in before]

           if added: 
             print "Added: ", ", ".join (added)
             for filein in enumerate(added):                  
               # check if it is a fits file
               filename = imdir + '/' + filein[1] 
               print filename
               # WARNING - hardcoded the '.fits' or '.fit' extensions
               if (filename.endswith('.fits') or filename.endswith('.fit')):
                 # Can load both data and header with this trick
                 data2, hdr = pyfits.getdata(filename, header=True)              
    
                 # Get Observation Header for image and display it
                 dateobs     = hdr['DATE-OBS']
                 timeobs     = hdr['TIME-OBS']
                 exposure    = hdr['EXPTIME']
                 CCDfilter   = hdr['FILTER']
                 print dateobs, timeobs, exposure, CCDfilter
                 ################################################
                 # now initiate the calibration, offsets and photometry

                 # test if the image has been calibrated previously
                 ### test_for_calib.py ...
                 ### if not processed, do calib.py ...

                 # find offsets from dataref
                 thisoffset = zach_offsets(dataref,data2)
                 #print "Offsets ", thisoffset

                 # create optphot init files
                 ### RUN OPTPHOT ......
                 #!!!!!!
                 

                 ################################################
               else:
                 print "Whatsaaap"

           before = after
           time.sleep (tsleep)   # Wait for tsleep seconds before repeating
    except KeyboardInterrupt:
     pass

def run_rtphos(xpapoint):

    # create a ds9 object linked with an XPA point
    win = ds9.ds9(xpapoint)
    # load the image which is displayed in ds9 (from disk!)
    ref_filename = win.get("file")
    dataref, hdr = pyfits.getdata(ref_filename, header=True)     
    
    # Check if the file has been processed (bias,dark,flat)
    # (will use standard IRAF keywords for this)
    ## ### TBD !!!! ####
    ### test_for_calib.py ...
    ### if not processed, do calib.py ...
    ### return the processed filename
    
    # extract path to where other images are expected to be coming in
    imdir = os.path.dirname(os.path.abspath(ref_filename))
    
    # select and centroid all created regions
    # make lists of xc,yc of comparison stars (name must start with 'C-'
    # get xt,yt of target 
    
    win.set("regions select all")

    # I do it twice with different radii on purpose
    win.set("regions centroid radius 20")
    win.set("regions centroid iteration 20")
    win.set("regions centroid")
    win.set("regions centroid radius 5")
    win.set("regions centroid iteration 5")
    win.set("regions centroid")
    # put back the default
    win.set("regions centroid radius 20")
    win.set("regions centroid iteration 20")

    win.set("regions format ds9")
    win.set("regions save "+ ref_filename + ".reg")
    win.set("saveimage jpeg " + ref_filename + ".jpg 75 ")
    
    sourcelist = win.get("regions selected") 
    sources    = pyregion.parse(sourcelist)
    n = len(sources)
    # find out which are the comparison stars (they must have "C-" in name)
 
    comparisons  = [(s.coord_list,s.comment) for s in sources if "C-" in s.comment]
    nc = len(comparisons)
    targets =  [(s.coord_list,s.comment) for s in sources if not("C-" in s.comment)]
    nt = len(targets)
    
    if nc == 0:
        print("Must have one comparison star labeled as 'C-<name>' ")
        raise Exception("Must have one comparison star labeled as 'C-<name>' ")
    
    if nt == 0:
        print("No target selected!")
        raise Exception("No target selected!")
    
    # get FWHM of stellar PSF using comparison stars
    psf_fwhm = get_comps_fwhm(comparisons, xpapoint)
    print "Found PSF FWHM:",  ("%.2f" % psf_fwhm)

    # WARNING - HARDCODING - deal with optional parameter later
    # take 3 seconds default sleep time
    tsleep = 3 
    seekfits(dataref, imdir, tsleep, comparisons, targets, psf_fwhm)
    
if  __name__ == "__main__":

    global xpapoint 
    import sys
    xpapoint       = sys.argv[1]
    run_rtphos(xpapoint)