#!/usr/bin/python
# M. Bogosavljevic, AOB, May 2015
### WORK IN PROGRESS!!!! ###
### Contains 
### gauss : just a gaussian function definition for PSF fitting
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
import astropy.io.fits as pyfits
import pyds9
import numpy as np
import os, time
import ccdcalib
from   scipy.optimize import curve_fit
from   subprocess import Popen, PIPE

##############################################################################
def dict_of_floats(list_of_strings, num_items):
    dict_of_floats={}
    
    for i in range(num_items):
        for j in range(num_items):
            dummy = [float(x) for x in list_of_strings[j].split()]
            dict_of_floats[j]=dummy[1:3]
            seeing = dummy[3]

    result = (dict_of_floats, seeing)
    return (result)

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
# requires gauss
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
# requires fwhm_from_star, namesplit

    # create a ds9 object linked with an XPA point
    win = pyds9.DS9(xpapoint)
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

#############################################################################

def zach_offsets(dataref,data2red):

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
    xsize2 = data2red.shape[0]
    ysize2 = data2red.shape[1]
    xstart = 9
    xend   = xsize2-10
    ystart = 9
    yend   = ysize2-10
    croped2 = data2red[xstart:xend,ystart:yend]
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

#############################################################################

def write_optphot_init(imdir, comparisons, targets):

    text_file1 = open(imdir+"psf.cat", "w")
    text_file2 = open(imdir+"stars.cat", "w")
    text_file1.write("!\n!\n!\n")
    text_file2.write("!\n!\n!\n")
   
    nt = len(targets)
    print("Targets nt:", nt)
    for i in range(0,nt):
        x = targets[i][0][0]
        y = targets[i][0][1]
        # write the target as the first source in stars.cat
        text_file2.write('%-5s %8.1f %8.1f \n' % ("1", x, y) )

    nc = len(comparisons)
    print("Comparisons nc:", nc)
    for k in range(0,nc):
        x = comparisons[k][0][0]
        y = comparisons[k][0][1]
        text_file1.write('%-5i %8.1f %8.1f \n' % (k+1, x, y) )
        text_file2.write('%-5i %8.1f %8.1f \n' % (k+2, x, y) )
        #print('%-5i %8.1f %8.1f' % (k+1, x, y) )

    text_file1.close()
    text_file2.close()
    print ("Wrote:", imdir+"psf.cat")
    print ("Wrote:", imdir+"stars.cat")
    return (1)


#############################################################################
def run_photometry(dirs, inputfile):

    os.chdir(dirs['reduced'])  # Move to the reduced image directory
    p = Popen(["optimal"], stdin=PIPE, stdout=PIPE)

    filename   = inputfile
    psfpos     = " psf.cat "
    starpos    = "stars.cat "
    verbose    = "N"
    badskyskew = "-1 "
    badskychi  = "-1 "
    fwhm       = "2.0 "
    clip       = "5.0 "
    aprad      = "10.0 "
    iopt       = "1 "
    searchrad  = "5.0 "
    adu        = "2.0"

    input_txt=[]
    input_txt.append(filename+psfpos+starpos+verbose)
    input_txt.append(badskyskew+badskychi+fwhm+clip+aprad+iopt+searchrad+adu)
    print input_txt[0]
    print input_txt[1]

    data_out = p.communicate(input_txt[0]+"\n"
                             +input_txt[1]+"\n")[0]

    results=data_out.split("\n")
    total_records = len(results)-1
    total_stars=total_records/2

    optimal_data = results[0:total_stars]
    aperture_data = results[total_stars:total_records]

    print optimal_data
    print
    print aperture_data
    print

    # Gets a float dictionary from a list of results (optimal or aperture)
    optimal_res=dict_of_floats(optimal_data, total_stars)
    optimal_stars=optimal_res[0]
    seeing = optimal_res[1]

    aperture_res=dict_of_floats(aperture_data, total_stars)
    aperture_stars=aperture_res[0]
    seeing = aperture_res[1]

    print
    print optimal_stars
    print seeing

    print
    print aperture_stars
    print seeing

    os.chdir(dirs['data'])  # Move back to the raw data directory
    return


#############################################################################

def seekfits(dataref, dirs, tsleep, comparisons, targets, psf_fwhm):
# requires zach_offsets, write_optphot_init

    print("IMDIR: "+dirs['data'])

    before = dict ([(f, None) for f in os.listdir(dirs['data'])])
    
    # counter needed just to print out progress
    count = 0.
    
    try:
      while 1:
           count = count + 1
           print int(count), "*LISTENING*", dirs['data'], time.strftime('%X %x %Z')
           after = dict ([(f, None) for f in os.listdir(dirs['data'])])
           added = [f for f in after if not f in before]

           if added: 
             print "Added files: ", ", ".join (added)
             for filein in enumerate(added):                  
               # check if it is a fits file
               filename = dirs['data']+'/'+filein[1] 
               print filename
               # WARNING - hardcoded the '.fits' or '.fit' extensions
               if (filename.endswith('.fits') or filename.endswith('.fit')):
                 # Can load both data and header with this trick
                 data2, hdr = pyfits.getdata(filename, header=True)              
                 print("I read image: "+filename)

                 # Get Observation Header for image and display it
                 #dateobs     = hdr['DATE-OBS']
                 #timeobs     = hdr['TIME-OBS']
                 #exposure    = hdr['EXPTIME']
                 #CCDfilter   = hdr['FILTER']
                 #print dateobs, timeobs, exposure, CCDfilter
                 ################################################

                 # Now initiate the calibration, offsets and photometry.
                 # ccdcalib will either calibrate the image and place the
                 # calibrated image file in the '/reduced/' directory or
                 # if the image did not require calibration just copy the image
                 # to the '/reduced/' directory. In either case the image will
                 # have a 'c_' prefix to indicate that ccdcalib has seen it.
                 calib_data=ccdcalib.calib(dirs, filename, data2, hdr)
                 data2 = calib_data[0]
                 hdr = calib_data[1]
                 calib_fname = calib_data[2]

                 # find offsets from dataref
                 thisoffset = zach_offsets(dataref,data2)
                 print "Offsets ", thisoffset

                 # create optphot init files
                 print("Targets here", targets)
                 t = write_optphot_init(dirs['reduced'], comparisons, targets)
                 print "Wrote Opphot init files"
                 # call optimal and do the photometry.
                 run_photometry(dirs, calib_fname)

           before = after
           time.sleep (tsleep)   # Wait for tsleep seconds before repeating
    except KeyboardInterrupt:
     pass

def run_rtphos(xpapoint):
# requires get_comps_fwhm, seekfits

    # create a ds9 object linked with an XPA point
    win = pyds9.DS9(xpapoint)
    # load the image which is displayed in ds9 (from disk!)
    ref_filename = win.get("file")

    # Set up input and output directories
    path, filename = os.path.split(ref_filename)    
    data_dir = path # Set data directory as the directory of the DS9 image.
    os.chdir("..")  # Move up a directory
    current_dir = os.path.abspath(os.curdir)
    bias_dir = current_dir+"/bias/"    # Ultimately these dirs need to be set by
    dark_dir = current_dir+"/dark/"    # the user using a DS9 input window.
    flat_dir = current_dir+"/flat/"    #  -""-
    reduced_dir = data_dir+"/reduced/" #  -""-
    if not os.path.exists(reduced_dir): os.makedirs(reduced_dir)
    os.chdir(data_dir) # Move back to the data directory
    # Make a dictionary with all the required directories
    dirs = {'current':current_dir, 'bias':bias_dir, 'dark':dark_dir, 'flat':flat_dir, 'data':data_dir, 'reduced':reduced_dir}

    dataref, hdr = pyfits.getdata(ref_filename, header=True)     
    print("... Working ...")

    # Check image calibration and calibrate if required.   
    result = ccdcalib.calib(dirs, ref_filename, dataref, hdr)
    dataref = result[0]
    hdr = result[1]
    calib_fname = result[2]

    ### return the processed filename - This is calib_fname but why return it here???
    
    # select and centroid all created regions
    # make lists of xc,yc of comparison stars (name must start with 'C-'
    # get xt,yt of target 
    
    win.set("regions select all")

    # I do it twice with different radii on purpose
    win.set("regions centroid radius 20")
    win.set("regions centroid iteration 20")
    # not happy at all with DS9 centering so repeating it 50 times
    for x in range(0, 49):
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
    print("I see "+ str(n) + " sources")
    # find out which are the comparison stars (they must have "C-" in name)
 
    comparisons  = [(s.coord_list,s.comment) for s in sources if "C-" in s.comment]
    nc = len(comparisons)
    targets =  [(s.coord_list,s.comment) for s in sources if not("C-" in s.comment)]
    nt = len(targets)
    #print("Comparisons:")
    #print(comparisons)
    #print("Targets:")
    #print(targets)

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
    seekfits(dataref, dirs, tsleep, comparisons, targets, psf_fwhm)
    
if  __name__ == "__main__":

    import sys
    xpapoint       = sys.argv[1]
    run_rtphos(xpapoint)
