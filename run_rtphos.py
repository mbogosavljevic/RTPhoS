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
#run_rtphos.py $xpa_method rtphos.defaults

import pyregion
import astropy.io.fits as pyfits
from   astropy.time import Time
from   datetime import datetime
import pyds9
import numpy as np
from numpy import inf
import os, time, math
import ccdcalib
from   scipy.optimize import curve_fit
from   scipy import signal, ndimage
from   subprocess import call, Popen, PIPE
import f2n
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.patches as mpatches

# NOTE: where else can we put this?
#sys.path.append("~/pythoncode/f2n/f2n") # The directory that contains f2n.py and f2n_fonts !

##############################################################################
def dict_of_floats(list_of_strings, num_items):
    dict_of_floats={}
    xypos={}
    
    for i in range(num_items):
        for j in range(num_items):
            #dummy = [float(x) for x in list_of_strings[j].split()]
            # Keep values as strings since star name is also included
            dummy = [x for x in list_of_strings[j].split()]
            #print ("Dummy", dummy)
            dict_of_floats[j]=dummy[1:3]
            seeing = dummy[3]
            xypos[j] = dummy[4:6]
                        
    result = (dict_of_floats, seeing, xypos)
    return (result)

##############################################################################
# Define model function to be used for PSF fit to the stars:
# in this case its a Gauss with a center Xc, sigma, amplitude A and base level B
def gauss(x, *p):
    A, xc, sigma, B  = p
    return A*np.exp(-(x-xc)**2/(2.*sigma**2)) + B

##############################################################################

def make_png(dirs, ref_filename, data, rbin):
# requires f2n installed

    frame_name = os.path.splitext(os.path.basename(ref_filename))[0]
    png_image_name = dirs['png'] + frame_name + ".png"
    png_image = f2n.f2nimage(numpyarray=data)  # give the image data to f2n class
    png_image.setzscale("flat","flat")  # works best to my liking
    png_image.rebin(rbin)
    png_image.makepilimage("lin")       # linear image scaling.
    # we can play with marking the star and comparisons in the image, like so
    # png_image.drawcircle(112, 101, r=15) 
    # TBD later, for now just label the frame
    png_image.writetitle(frame_name, colour=(200, 200, 0))
    png_image.tonet(png_image_name)     # write the png.

    
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
    print ("I will use", nc, " stars for FWHM")
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
        cname = (comparisons[i][1])
        print (i+1, cname, fwhm)
        tot_fwhm = tot_fwhm + fwhm
        
    mean_fwhm = tot_fwhm / nc
    return mean_fwhm

#############################################################################
def stripdate(longjd):

    # reduce a Julian Date type date format to 2 significant figures before the decimal.
    twosig = int(longjd/100.0)*100.0
    shortjd =  jd - twosig

    return (shortjd, twosig)

#############################################################################
def barytime(checklist):

    # Deconstruct the Date string from the DATE Keyword
    date  = checklist['DATE'].split('-')
    date = " ".join(date)
    # Deconstruct the Time string from the TIME Keyword
    time  = checklist['TIME'].split(':')
    time = " ".join(time)
    # Deconstruct the RA string from the RA Keyword
    if (checklist['RA'] != "Invalid"):
       if ':' in checklist['RA']: ra  = checklist['RA'].split(':')
       if ' ' in checklist['RA']: ra  = checklist['RA'].split(' ')
    else:
       ra = "NaN"
    ra = " ".join(ra)
    # Deconstruct the Dec string from the DEC Keyword
    if (checklist['DEC'] != "Invalid"):
       if ':' in checklist['DEC']: dec  = checklist['DEC'].split(':')
       if ' ' in checklist['DEC']: dec  = checklist['DEC'].split(' ')
    else:
       dec = "NaN"
    dec = " ".join(dec)
    # Get the Exposure time
    exposure = str(checklist['EXPOSURE'])
    if exposure=="Invalid": exposure = "0.0"

    equin  = '2000.0'
    utcorr = '0.0'

#   For Later:
#   Get Equinox form Header
#   Get Time Correction from Header

#   The input to barycor should be the following:
#   equin, rah, ram, ras, decd, decm, decs
#   year, month, day, hrs, mins, secs, utcorr, exposure
    datetimeinfo = []
    inputline1 = equin+" "+" "+ra+  " "+dec
    inputline2 = date +" "+" "+time+" "+" "+utcorr+" "+" "+exposure

#    print inputline1
#    print inputline2

    p = Popen(["barycor"], stdin=PIPE, stdout=PIPE)

    time_BDJD = p.communicate(inputline1+"\n"
                             +inputline2)[0]

    return time_BDJD


#############################################################################
def zach_offsets(dataref,data2red):
    
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
    median1 = np.median(croped1)
    
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
    median2 = np.median(croped2)
      	
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
    xshift = (np.argmax(signal.correlate(xvals1,xvals2)))-(ysize-1)
    yshift = (np.argmax(signal.correlate(yvals1,yvals2)))-(xsize-1)
    
    return (xshift, yshift)

#############################################################################

def write_optphot_init(rtdefs, imdir, comparisons, targets, thisoffset):

    text_file1 = open(imdir+rtdefs['psfs'], "w")
    text_file2 = open(imdir+rtdefs['stars'], "w")
    text_file1.write("!\n!\n!\n")
    text_file2.write("!\n!\n!\n")
   
    xoff = float(thisoffset[0])
    yoff = float(thisoffset[1])

    nt = len(targets)
    for i in range(0,nt):
        x = targets[i][0][0]
        y = targets[i][0][1]
        name = targets[i][1]

        # write the targets as first source(s) in stars.cat
        text_file2.write('%-5i %8.1f %8.1f %-15s \n' % (i+1, x-xoff, y-yoff, name) )
        #text_file2.write('%-5i %8.1f %8.1f \n' % (i+1, x-xoff, y-yoff))

    nc = len(comparisons)
    for k in range(0,nc):
        x = comparisons[k][0][0]
        y = comparisons[k][0][1]
        name = comparisons[k][1]

        # write comparisons to PSF coordinates file.
        text_file1.write('%-5i %8.1f %8.1f %-15s \n' % (k+1, x-xoff, y-yoff, name) )
        #text_file1.write('%-5i %8.1f %8.1f \n' % (k+1, x-xoff, y-yoff))

        # add comparisons as additional sources to stars coordinates file.
        text_file2.write('%-5i %8.1f %8.1f %-15s \n' % (k+nt+1, x-xoff, y-yoff, name) )
        #text_file2.write('%-5i %8.1f %8.1f \n' % (k+nt+1, x-xoff, y-yoff))

    text_file1.close()
    text_file2.close()
    print ("Wrote " + str(nc) + " in " + imdir + rtdefs['psfs'])
    print ("Wrote " + str(nc+nt) + " in "  + imdir + rtdefs['stars'])
    return (1)


#############################################################################
def run_photometry(rtdefs, dirs, inputfile, psf_fwhm):

    filename   = inputfile+" "
    psfpos     = rtdefs['psfs']+" "    
    starpos    = rtdefs['stars']+" "   
    badskyskew = rtdefs['skyskew']+" "
    badskychi  = rtdefs['skyfit']+" "
    clip       = rtdefs['cradius']+" "
    aprad      = rtdefs['aradius']+" "
    iopt       = rtdefs['starnumber']+" "
    searchrad  = rtdefs['sradius']+" "
    adu        = rtdefs['gain']
    verbose    = rtdefs['verbose']
    fwhm       = str(psf_fwhm)+" "
    
    input_txt=[]
    input_txt.append(filename+psfpos+starpos+verbose)
    input_txt.append(badskyskew+badskychi+fwhm+clip+aprad+iopt+searchrad+adu)

    #print filename+psfpos+starpos+verbose
    #print badskyskew+badskychi+fwhm+clip+aprad+iopt+searchrad+adu

    os.chdir(dirs['reduced'])  # Move to the reduced image directory
    p = Popen(["optimal"], stdin=PIPE, stdout=PIPE)
    data_out = p.communicate(input_txt[0]+"\n"
                             +input_txt[1]+"\n")[0]
 
    if verbose=='Y':
        print " ### OPTHOT OUTPUT START ###"
        print data_out
        print " ### OPTPHOT OUTPUT END ###"

    results=data_out.split("\n")
    total_records = len(results)-1
    total_stars=total_records/2

    optimal_data = results[0:total_stars]
    aperture_data = results[total_stars:total_records]

    # Debug
    #print ("Optimal data", optimal_data)
    #print ("Aperture data", aperture_data)

    # Gets a float dictionary from a list of results (optimal or aperture)
    optimal_res=dict_of_floats(optimal_data, total_stars)
    optimal_stars=optimal_res[0]
    seeing = optimal_res[1]
    xypos = optimal_res[2]
    
    aperture_res=dict_of_floats(aperture_data, total_stars)
    aperture_stars=aperture_res[0]
    seeing = aperture_res[1]
    xypos = aperture_res[2]

    photometry_result = (optimal_stars, aperture_stars, seeing, xypos)

    os.chdir(dirs['data'])  # Move back to the raw data directory
    return photometry_result

#############################################################################
def outputfiles(alltargets, optimalist, aperatlist, seeing, \
                frame_time, frame_timerr, pdatetime, filename, count):

    # Loop will write out files with the following outputformat:
    # seq. number, frame_time, frame_timerr, flux, flux error, seeing
    for i in range(0,len(alltargets)):
        # First write the optimal photometry data
        with open(alltargets[i][1]+".opt", "a") as outfile:
            outdata = (count, pdatetime, frame_time, frame_timerr*86400.0, optimalist[i][0], \
                       optimalist[i][1], float(seeing), filename)
            fmtstring = '%5i %20s %15.6f %6.2f %12s %9s %6.2f %s \n'
            outfile.write(fmtstring % outdata)
#            outfile.write(str(count)+" "+pdatetime+" "+str(frame_time)+\
#            " "+ str(frame_timerr)+" "+str(optimalist[i][0])+" "+\
#            str(optimalist[i][1])+" "+str(seeing)+" "+filename+" \n")
        # Now write the aperture photometry data
        with open(alltargets[i][1]+".dat", "a") as outfile:
            outdata = (count, pdatetime, frame_time, frame_timerr*86400.0, aperatlist[i][0], \
                       aperatlist[i][1], float(seeing), filename)
            fmtstring = '%5i %20s %15.6f %6.2f %12s %9s %6.2f %s \n'
            outfile.write(fmtstring % outdata)
#            outfile.write(str(count)+" "+pdatetime+" "+str(frame_time)+\
#            " "+ str(frame_timerr)+" "+str(aperatlist[i][0])+" "+\
#            str(aperatlist[i][1])+" "+str(seeing)+" "+filename+" \n")

    ### Files are always appended. This might be a problem when running a
    ### new round of rtphos of the same data. Perhaps we should make rtphos
    ### to first check is *.opt and *.dat files exist in the output directory
    ### delete them and then proceed with the reduction(?)
    ### Should we place the output files in a new directory (e.g. output)?

    return


#############################################################################
def seekfits(rtdefs, dataref, dirs, tsleep, comparisons, targets, psf_fwhm, ax1, ax2, ax3, ax4, ax5):
# requires zach_offsets, write_optphot_init
    
    # Reduced data lists
    xdata=[]         # X-axis data (time)
    yseeing=[]       # Seeing data
    yrawtarget=[]    # Raw target counts
    yrawtargeterr=[] # Raw target error bars
    yrawcomp=[]      # Raw comparison counts
    yrawcomperr=[]   # Raw comparison error bars
    ydflux=[]        # Differential photometry counts
    ydfluxerr=[]     # Differential photometry error bars

    before = sorted(os.listdir(dirs['data']))
    # a switch to first check if files found on startup
    # need to be reduced or not (if the reduced output exists)
    reducebefore = True

    # counter needed just to print out progress
    count = 0

    try:
       while 1:
           print "*LISTENING*", dirs['data'], time.strftime('%X %x %Z')
           after = sorted(os.listdir(dirs['data']))
           added = [f for f in after if not f in before]

           if reducebefore:
               added = added + before
               reducebefore = False

           if added:   
               for filein in added:                  
                   # check if it is a fits file
                   filename = dirs['data']+'/'+filein
                   
                   # WARNING - hardcoded the '.fits' or '.fit' extensions
                   if (filename.endswith('.fits') or filename.endswith('.fit')):
                       # Can load both data and header with this trick
                       data2, hdr = pyfits.getdata(filename, header=True) 
                       print("I read image: "+filename)

                       # Data array used for plotting the current image. 
                       dataplt = data2  

                       # First check that all the required header keywords are in
                       # the FITS file and then get the time stamp for this frame.
                       # If the RA and DEC of the image are in the headers then
                       # time will be in Barycentric Dynamical Julian Date. If not 
                       # then time will be in plain simple Julian Date.
                       checklist = ccdcalib.makechecklist(hdr)

                       count = count + 1 
                       if checklist['RA']=="Invalid" or checklist['DEC']=="Invalid":  
                           mdatetime = checklist['DATE']+" "+checklist['TIME']
                           pdatetime = checklist['DATE']+"|"+checklist['TIME']
                           #print ("HERE",mdatetime)
                           # this is needed to get plot-able UTC time
                           sdatetime = datetime.strptime(mdatetime,  "%Y-%m-%d %H:%M:%S")                        
                           mdatetime = Time(mdatetime, format='iso', scale='utc')
                           time_frame  = mdatetime.jd # in days
                           exposure = float(checklist['EXPOSURE']) # in seconds
                           midexp = (exposure / 2.0) / 86400. # in days
                           # Make frame time the middle of the exposure
                           # these are now float variables, not string
                           frame_time = float(time_frame) + midexp
                           frame_timerr = midexp
                       else:
                           time_BDJD = barytime(checklist)
                           frame_time = float(time_BDJD[1])
                           #frame_timerr  = format(float(time_BDJD[2]), '.15g')
                           frame_timerr = midexp

                       # Strip JD and reduced it to 2 significant figures.
                       stripjd = int(float(frame_time)/100.0)*100.0
                       twosig_time =  float(frame_time) - stripjd

                       # Now initiate the calibration, offsets and photometry.
                       # ccdcalib will either calibrate the image and place the
                       # calibrated image file in the '/reduced/' directory or
                       # if the image did not require calibration just copy the image
                       # to the '/reduced/' directory. In either case the image will
                       # have a 'c_' prefix to indicate that ccdcalib has seen it.
                       calib_data = ccdcalib.calib(rtdefs, dirs, filename, data2, hdr)
                       (data2, hdr, calib_fname) = calib_data

                       # find offsets from dataref
                       thisoffset = zach_offsets(dataref,data2)
                       print "Offsets: (x,y) ", thisoffset

                       # create opphot init files
                       t = write_optphot_init(rtdefs, dirs['reduced'], comparisons, targets, thisoffset)
                       print "Wrote Opphot init files"
                       # call optimal and do the photometry.
                       frame_photometry = run_photometry(rtdefs, dirs, calib_fname, psf_fwhm)
  
                       # Deconstruct the photometry results from optimal.f90
                       (optimaldict, aperatdict, seeing, xypos) = frame_photometry
                       optimalist = optimaldict.values()
                       aperatlist = aperatdict.values()
                       xyposlist  = xypos.values()

                       # File output
                       junk, sfilename = os.path.split(filename)
                       alltargets = targets + comparisons
                       outputfiles(alltargets, optimalist, aperatlist, seeing, \
                                   frame_time, frame_timerr, pdatetime, sfilename, count)

                       # Fill the data lists
                       xdata.append(twosig_time)
                       yseeing.append(seeing)
                       yrawtarget.append(float(optimalist[0][0]))
                       yrawtargeterr.append(float(optimalist[0][1]))
                       yrawcomp.append(float(optimalist[1][0]))
                       yrawcomperr.append(float(optimalist[1][1]))
                       # Do the differential photometry calculations
                       tcounts    = float(optimalist[0][0])
                       terror     = float(optimalist[1][1])
                       ccounts    = float(optimalist[1][0])
                       cerror     = float(optimalist[1][1])
                       ydfluxs    = (tcounts/ccounts)
                       ydfluxerrs = ydfluxs*math.sqrt( ((terror/tcounts)**2.0) + \
                                                     ((cerror/ccounts)**2.0) )
                       ydflux.append(ydfluxs)                       
                       ydfluxerr.append(ydfluxerrs)

                       # Begin graphics output calculations
                       dataplt[dataplt == -inf] = 0.0             # Remove inf values
                       # Crop a 100px square around the target
                       targetx = float(xyposlist[0][0])
                       targety = float(xyposlist[0][1])
                       target_crop = dataplt[targety-50:targety+50,targetx-50:targetx+50]
                       medianintens = np.median(target_crop)
                       target_crop[target_crop==0] = medianintens # Remove zero values
                       target_crop = np.log(target_crop)          # Use for log scale plotting
                       cropmin = np.amin(target_crop)
                       cropmax = np.amax(target_crop)
                       # Attempt for a reasonable intensity scale 
                       maxintens = ((cropmax-cropmin)/2.0)+cropmin
                       # Crop a 100px square around the first comparison star
                       compx = float(xyposlist[1][0])
                       compy = float(xyposlist[1][1])
                       comp_crop = dataplt[compy-50:compy+50,compx-50:compx+50]
                       comp_crop[comp_crop==0] = medianintens     # Remove zero values
                       comp_crop = np.log(comp_crop)              # Use for log scale plotting

                       # Target thumbnail
                       ax1.plot([50,50],[0,100],'r:')             # Plot cross-hairs
                       ax1.plot([0,100],[50,50],'r:')             #      -""-
                       ax1.imshow(target_crop, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))
                       ax1.text(1, 5, targets[0][1], color='yellow', fontsize=12)
                       #ax1.text(100, 5, "Current frame: "+filein, color='yellow', fontsize=12)
                       #ax1.text(100, 25,"Frames processed: "+str(count), color='yellow', fontsize=12)
                       ax1.set_xlabel("Physical X: "+xyposlist[0][0])
                       ax1.set_ylabel("Physical Y: "+xyposlist[0][1])
                       ax1.set_xticklabels([])
                       ax1.set_yticklabels([])

                       # Comparison thumbnail
                       ax2.plot([50,50],[0,100],'r:')             # Plot cross-hairs
                       ax2.plot([0,100],[50,50],'r:')             #      -""-
                       ax2.imshow(comp_crop, cmap='gray', norm=LogNorm(vmin=cropmin, vmax=maxintens))
                       ax2.text(1, 5, comparisons[0][1], color='yellow', fontsize=12)
                       ax2.set_xlabel("Physical X: "+xyposlist[0][0])
                       ax2.set_ylabel("Physical Y: "+xyposlist[0][1])
                       ax2.set_xticklabels([])
                       ax2.set_yticklabels([])

                       # Target Raw counts
                       ax3.errorbar(xdata, yrawtarget, yerr=yrawtargeterr, label='Test', fmt='go')
                       ax3.errorbar(xdata, yrawcomp,   yerr=yrawcomperr, fmt='ro')                      
                       green = mpatches.Patch(color='green', label=targets[0][1])
                       red = mpatches.Patch(color='red', label=comparisons[0][1])
                       ax3.legend(handles=[green, red])                       
                       ax3.grid(True)
                       ax3.set_ylabel('Raw Counts')

#                       # Comparison Raw counts
#                       plt.subplot(6,1,4)
#                       plt.errorbar(xdata, yrawcomp, yerr=yrawcomperr, fmt='ro')
#                       #plt.title(targets[0][1]+' flux')
#                       plt.xticks([])
#                       plt.grid(True)
#                       plt.ylabel('Raw Counts')
#                       plt.tight_layout()
                       # Target Relative Flux
                       ax4.errorbar(xdata, ydflux, yerr=ydfluxerr, fmt='bo')
                       blue = mpatches.Patch(color='blue', label=targets[0][1])
                       ax4.legend(handles=[blue])                       
                       ax4.grid(True)
                       ax4.set_ylabel('Relative Flux')

                       # Seeing plot
                       ax5.scatter(xdata,yseeing)
                       blue = mpatches.Patch(color='blue', label='Seeing')
                       ax5.legend(handles=[blue])                       
                       ax5.grid(True)
                       ax5.set_xlabel('Time JD+ '+str(stripjd))
                       ax5.set_ylabel('Pixels')

                       plt.pause(1) # Small time delay to allow for matplotlib to plot the graphs

                       # Text output
                       print "============================================"
                       print "FILENAME ", filename
                       print "FRAME_TIME ", frame_time, frame_timerr
                       #print "Optimal Photometry Results:"
                       for i in range(0,len(optimalist)):
                           print alltargets[i][1], optimalist[i][0], optimalist[i][1], \
                                 seeing, xyposlist[i][0], xyposlist[i][1]
                       print "------------------------------"
                       #print "Aperature Photometry Results:"
                       for i in range(0,len(aperatlist)):
                           print alltargets[i][1], aperatlist[i][0], aperatlist[i][1], \
                                 seeing, xyposlist[i][0], xyposlist[i][1]
                       print

           before = after
           time.sleep(tsleep)   # Wait for tsleep seconds before repeating

    except KeyboardInterrupt:
       pass

    return

#############################################################################
def run_rtphos(rtphosdir, xpapoint, pathdefs):
# requires get_comps_fwhm, seekfits
  
    print " ==== RTPhoS START ==== " + time.strftime('%X %x %Z')

    ######## SETTING UP DIRECTORIES AND LINKS
    # create a ds9 object linked with an XPA point
    win = pyds9.DS9(xpapoint)
    # image name which is displayed in ds9 
    ref_filename = win.get("file")

    # Set up input and output directories
    path, filename = os.path.split(ref_filename) 

    # Read in default values from file
    with open (pathdefs, "r") as defsfile:
        defs = defsfile.readlines()
        print " Read defaults from: "+ pathdefs
        print " "
        for l in range(0,len(defs)):
            print ("Def ",l," ",defs[l][:-2]) # without the \n
        print " ========================================== "

    ## If you would want to grep the defaults file and find match
    ## this would be one way to do it
    #Set data directory from defaults file
    #findthis = "# Raw data frames"
    #matching = [s for s in defs if findthis in s]
    #data_dir = matching[0].split()[0]

    # Will use positional assignement for defaults
    # These are all string types:    
    data_dir   = defs[1].split()[0]
    bias_dir   = defs[2].split()[0]
    dark_dir   = defs[3].split()[0]
    flat_dir   = defs[4].split()[0] 
    biaswc     = defs[5].split()[0]
    darkwc     = defs[6].split()[0]
    flatwc     = defs[7].split()[0]
    mbias      = defs[8].split()[0]
    mdark      = defs[9].split()[0]
    mflat      = defs[10].split()[0]
    psfs       = defs[11].split()[0]
    stars      = defs[12].split()[0]
    cprefix    = defs[13].split()[0]
    sradius    = defs[14].split()[0]
    aradius    = defs[15].split()[0]
    cradius    = defs[16].split()[0]
    starnumber = defs[17].split()[0]
    skyskew    = defs[18].split()[0]
    skyfit     = defs[19].split()[0]
    gain       = defs[20].split()[0]
    # These are numbers:
    verbose    =   int(defs[21].split()[0])
    tsleep     =   int(defs[22].split()[0])

    reduced_dir  = data_dir+"/reduced/"
    png_dir      = data_dir+"/png/"

    # Current root working dir
    current_dir = os.path.abspath(os.path.join(data_dir, os.pardir))

    if not os.path.exists(reduced_dir): 
        os.makedirs(reduced_dir)
        print "Created output reduced dir: " + reduced_dir
    if not os.path.exists(png_dir): 
        os.makedirs(png_dir)
        print "Created output png dir: " + png_dir
    
    # Create the symbolic links required for running barycor.f90
    os.chdir(reduced_dir)
    call(['ln', '-s', rtphosdir+'/Timing/jpleph.dat', 'JPLEPH'])
    call(['ln', '-s', rtphosdir+'/Timing/leap.dat', 'leapdat'])
    os.chdir(data_dir) # Move back to the data directory

    # Convert verbose switch value to a string
    if verbose==1:
       verbose='Y'
    else:
       verbose='N'

    # Make a dictionary with all the required directories.
    dirs = {'current':current_dir, 'bias':bias_dir, 'dark':dark_dir, \
            'flat':flat_dir, 'data':data_dir, 'reduced':reduced_dir, 'png':png_dir}

    # Make a dictionary with all the required parameters.
    rtdefs = {'biaswc':biaswc, 'darkwc':darkwc, 'flatwc':flatwc, 'mbias':mbias,\
               'mdark':mdark,   'mflat':mflat, 'psfs':psfs, 'stars':stars, 'cprefix':cprefix, \
              'sradius':sradius, 'aradius':aradius, 'cradius':cradius, 'starnumber':starnumber, \
              'skyskew':skyskew, 'skyfit':skyfit, 'gain':gain, 'verbose':verbose, 'tsleep':tsleep}

    ######## DS9 SOURCE IDENTIFICATION AND FWHM ESTIMATE ###########
    # Ceontrid regions in DS9
    win.set("regions select all")

    # I do it twice with different radii on purpose
    win.set("regions centroid radius 20")
    win.set("regions centroid iteration 20")
    # not happy at all with DS9 centering so repeating it 20 times
    for x in range(0, 19):
        win.set("regions centroid")
    win.set("regions centroid radius 5")
    win.set("regions centroid iteration 5")
    win.set("regions centroid")
    # put back the default
    win.set("regions centroid radius 20")
    win.set("regions centroid iteration 20")

    # save regions file for later
    win.set("regions format ds9")
    win.set("regions save "+ ref_filename + ".reg")

    # Get source (target, comparison) lists from regions selected
    sourcelist = win.get("regions selected") 
    sources    = pyregion.parse(sourcelist)
    n = len(sources)
    print(" In DS9 I see "+ str(n) + " sources labeled")

    for l in range(0,n):
       if (sources[l].comment is not None):
            sources[l].comment = sources[l].comment[sources[l].comment.find("{")+1:\
            sources[l].comment.find("}")]
            print (l+1, sources[l].comment, sources[l].coord_list)
       else:
            print "You have some unlabeled sources!"
            print sources[l].coord_list
            print "====  Exiting. ==== "
            return

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
    
    print "Targets:"
    print targets
    print "Comparisons:"
    print comparisons

    # get FWHM of stellar PSF using comparison stars
    psf_fwhm = get_comps_fwhm(comparisons, xpapoint)
    print "Found PSF FWHM:",  ("%.2f" % psf_fwhm)

    ##### START PIPELINE #############################################
    # This is where the pipiline looks at the data for the first time!
    print " DS9: Starting with file " + ref_filename
    dataref, hdr = pyfits.getdata(ref_filename, header=True)     

    # Check first DS9 image calibration and calibrate if required.   
    result = ccdcalib.calib(rtdefs, dirs, ref_filename, dataref, hdr)
    (dataref, hdr, calib_fname) = result

    # Do the photometry
    plt.ion()
    plt.figure(figsize=(12,12))
    ax1  = plt.subplot(421)
    ax2  = plt.subplot(422)
    ax3  = plt.subplot(412)
    ax4  = plt.subplot(413, sharex=ax3)
    ax5  = plt.subplot(414, sharex=ax3)

    #fig1 = plt.figure(figsize=(8,8))
    #ax1  = fig1.add_subplot(221)
    #ax2  = fig1.add_subplot(222)
    #fig2 = plt.figure(figsize=(12,12))

    seekfits(rtdefs, dataref, dirs, tsleep, comparisons, targets, psf_fwhm, ax1, ax2, ax3, ax4, ax5)
    
if  __name__ == "__main__":

    import sys
    rtphosdir, dummy = os.path.split(sys.argv[0]) 
    xpapoint       = sys.argv[1]
    pathdefs       = sys.argv[2]
    run_rtphos(rtphosdir, xpapoint, pathdefs)
