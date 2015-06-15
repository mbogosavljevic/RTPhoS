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
from   astropy.time import Time
import pyds9
import numpy as np
import os, time
import ccdcalib
from   scipy.optimize import curve_fit
from   scipy import signal, ndimage
from   subprocess import call, Popen, PIPE
import matplotlib.pyplot as plt
import f2n
import sys

# NOTE: where else can we put this?
#sys.path.append("~/pythoncode/f2n/f2n") # The directory that contains f2n.py and f2n_fonts !

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

def write_optphot_init(imdir, comparisons, targets, thisoffset):

    text_file1 = open(imdir+"psf.cat", "w")
    text_file2 = open(imdir+"stars.cat", "w")
    text_file1.write("!\n!\n!\n")
    text_file2.write("!\n!\n!\n")
   
    xoff = float(thisoffset[0])
    yoff = float(thisoffset[1])

    nt = len(targets)
    for i in range(0,nt):
        x = targets[i][0][0]
        y = targets[i][0][1]
        name = targets[i][1]
        name2 = name[name.find("{")+1:name.find("}")]
        # write the targets as first source(s) in stars.cat
        text_file2.write('%-5i %8.1f %8.1f %-15s \n' % (i, x-xoff, y-yoff, name2) )

    nc = len(comparisons)
    for k in range(0,nc):
        x = comparisons[k][0][0]
        y = comparisons[k][0][1]
        name = comparisons[k][1]
        name2 = name[name.find("{")+1:name.find("}")]
        # write comparisons to psf.cat
        text_file1.write('%-5i %8.1f %8.1f %-15s \n' % (k+1, x-xoff, y-yoff, name2) )
        # add comparisons as additional sources to stars.cat
        text_file2.write('%-5i %8.1f %8.1f %-15s \n' % (k+nt+1, x-xoff, y-yoff, name2) )

    text_file1.close()
    text_file2.close()
    print ("Wrote " + str(nc) + " in " + imdir + "psf.cat")
    print ("Wrote " + str(nc+nt) + " in "  + imdir + "stars.cat")
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

    data_out = p.communicate(input_txt[0]+"\n"
                             +input_txt[1]+"\n")[0]

#    print input_txt[0]
#    print input_txt[1]

    results=data_out.split("\n")
    total_records = len(results)-1
    total_stars=total_records/2

    optimal_data = results[0:total_stars]
    aperture_data = results[total_stars:total_records]

    # Gets a float dictionary from a list of results (optimal or aperture)
    optimal_res=dict_of_floats(optimal_data, total_stars)
    optimal_stars=optimal_res[0]
    seeing = optimal_res[1]

    aperture_res=dict_of_floats(aperture_data, total_stars)
    aperture_stars=aperture_res[0]
    seeing = aperture_res[1]

    photometry_result = (optimal_stars, aperture_stars, seeing)

    os.chdir(dirs['data'])  # Move back to the raw data directory
    return photometry_result


#############################################################################

class seekfits():
# requires zach_offsets, write_optphot_init
    
    #Suppose we know the x range
    #min_x = 0
    #max_x = 10

    def on_launch(self, dirs):
        #Set up plot
        self.figure, self.ax = plt.subplots()
        self.lines, = self.ax.plot([],[], 'o')
#        self.ax.errorbar([],[],yerr=[])
#        self.lines, = self.ax.errorbar([],[],[],[],fmt='o')
#        self.line, (down, up), verts = self.ax.errorbar([],[],yerr=[],ftm='o')
        #Autoscale on unknown axis and known lims on the other
        self.ax.set_autoscaley_on(True)
        self.ax.set_autoscalex_on(True)
        #self.ax.set_xlim(self.min_x, self.max_x)
        #Other stuff
        self.ax.grid()
        self.ax.set_xlabel('Time')
        self.ax.set_ylabel('Counts')
        print("IMDIR: "+dirs['data'])

    def on_running(self, xdata, ydata, yerror):
        #Update data (with the new _and_ the old points)
#        self.lines, = self.ax.plot(xdata, ydata, 'o')
#        self.ax.errorbar(xdata, ydata, yerr=yerror)
#        self.lines, = self.ax.plot(xdata,ydata, 'o')
        self.lines.set_xdata(xdata)
        self.lines.set_ydata(ydata)
#        down.set_ydata()
#    bottoms.set_ydata(y - yerr)
#    tops.set_ydata(y + yerr)

        #Need both of these in order to rescale
        self.ax.relim()
        self.ax.autoscale_view()
        #We need to draw *and* flush
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()


    def __call__(self, dataref, dirs, tsleep, comparisons, targets, psf_fwhm ):
        self.on_launch(dirs)
        xdata = []
        ydata = []
        xerror = []
        yerror = []

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
                  count = count + 1             
                  for filein in added:                  
                     # check if it is a fits file
                     filename = dirs['data']+'/'+filein
                     print filename
                     # WARNING - hardcoded the '.fits' or '.fit' extensions
                     if (filename.endswith('.fits') or filename.endswith('.fit')):

                        # Can load both data and header with this trick
                        data2, hdr = pyfits.getdata(filename, header=True) 
                        print("I read image: "+filename)

                        # First check that all the required header keywords are in
                        # the FITS file and then get the time stamp for this frame.
                        # If the RA and DEC of the image are in the headers then
                        # time will be in Barycentric Dynamical Julian Date. If not 
                        # then time will be in plain simple Julian Date.
                        checklist = ccdcalib.makechecklist(hdr)

                        if checklist['RA']=="Invalid" or checklist['DEC']=="Invalid":  
                           datetime = checklist['DATE']+" "+checklist['TIME']
                           datetime = Time(datetime, format='iso', scale='utc')
                           time_frame  = datetime.jd
                           exp = float(checklist['EXPOSURE'])
                           midexp = exp/2.0
                           # Make frame time the middle of the exposure
                           frame_time = format(time_frame+midexp, '.15g')
                           frame_timerr = format(midexp/86400.0, '.15g')   
                        else:
                           time_BDJD = barytime(checklist)
                           frame_time = format(float(time_BDJD[1]), '.15g')
                           frame_timerr  = format(float(time_BDJD[2]), '.15g')

                        # Strip JD and reduced it to 2 significant figures.
                        # stripjd = int(float(frame_time)/100.0)*100.0
                        # twosig_time =  float(frame_time) - stripjd

                        # Now initiate the calibration, offsets and photometry.
                        # ccdcalib will either calibrate the image and place the
                        # calibrated image file in the '/reduced/' directory or
                        # if the image did not require calibration just copy the image
                        # to the '/reduced/' directory. In either case the image will
                        # have a 'c_' prefix to indicate that ccdcalib has seen it.
                        calib_data = ccdcalib.calib(dirs, filename, data2, hdr)
                        (data2, hdr, calib_fname) = calib_data

                        # find offsets from dataref
                        thisoffset = zach_offsets(dataref,data2)
                        print "Offsets: (x,y) ", thisoffset

                        # create optphot init files
                        t = write_optphot_init(dirs['reduced'], comparisons, targets, thisoffset)
                        print "Wrote Opphot init files"
                        # call optimal and do the photometry.
                        frame_photometry = run_photometry(dirs, calib_fname)

                        # Deconstruct the photometry results from optimal.f90
                        (optimaldict, aperatdict, seeing) = frame_photometry
                        optimalist  = optimaldict.values()
                        aperatlist  = aperatdict.values()

                        # Screen output
                        print "============================================"
                        print "FILENAME ", filename
                        print "FRAME_TIME ", frame_time, frame_timerr, count, len(optimalist)
                        #print "Optimal Photometry Results:"
                        for i in range(0,len(optimalist)):
                            print "OPT_STAR ", i+1, optimalist[i][0], optimalist[i][1], seeing
                        print "------------------------------"
                        #print "Aperature Photometry Results:"
                        for i in range(0,len(aperatlist)):
                            print "APER_STAR ", i+1, aperatlist[i][0], aperatlist[i][1], seeing
                        print

                        xdata.append(frame_time)
                        ydata.append(optimalist[1][0])

                        xer = 1
                        xerror.append(xer)
                        yerror.append(optimalist[1][1])

                        self.on_running(xdata, ydata, yerror)

                before = after
                time.sleep(tsleep)   # Wait for tsleep seconds before repeating

        except KeyboardInterrupt:
           pass

        return # xdata, ydata

#############################################################################
def run_rtphos(xpapoint):
# requires get_comps_fwhm, seekfits

    # create a ds9 object linked with an XPA point
    win = pyds9.DS9(xpapoint)
    # load the image which is displayed in ds9 (from disk!)
    ref_filename = win.get("file")

    # Set up input and output directories
    path, filename = os.path.split(ref_filename)    
    data_dir = path # Set data directory as the directory of the DS9 image.
    current_dir = os.path.abspath(os.path.join(data_dir, os.pardir))
    bias_dir = current_dir+"/bias/"    # Ultimately these dirs need to be set by
    dark_dir = current_dir+"/dark/"    # the user using a DS9 input window.
    flat_dir = current_dir+"/flat/"    #  -""-
    reduced_dir  = data_dir+"/reduced/" #  -""-
    png_dir      = data_dir+"/png/"

    if not os.path.exists(reduced_dir): os.makedirs(reduced_dir)
    if not os.path.exists(png_dir): os.makedirs(png_dir)

    # Create the symbolic links required for running barycor.f90
    # WARNING: The links are hardwired into the code. This should be made
    # such that the code checks to see if the links are true. If they are
    # proceed as normal otherwise return just Julian Days and not BDJD.
    os.chdir(reduced_dir)
    call(['ln', '-s', '/opt/star-kapuahi/etc/jpleph.dat', 'JPLEPH'])
    call(['ln', '-s', '/home/zac/Software/Ark/data/leap.dat', 'leapdat'])
    os.chdir(data_dir) # Move back to the data directory

    # Make a dictionary with all the required directories
    dirs = {'current':current_dir, 'bias':bias_dir, 'dark':dark_dir, \
            'flat':flat_dir, 'data':data_dir, 'reduced':reduced_dir, 'png':png_dir}

    # This is where the pipiline looks at the data for the first time!
    dataref, hdr = pyfits.getdata(ref_filename, header=True)     
    print " ==== RTPhoS START ==== " + time.strftime('%X %x %Z')
    print " DS9: Starting with file " + ref_filename

    # Check image calibration and calibrate if required.   
    result = ccdcalib.calib(dirs, ref_filename, dataref, hdr)
    (dataref, hdr, calib_fname) = result
    
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

    #win.set("saveimage jpeg " + ref_filename + ".jpg 75 ")
    # Instead of using win.set, which saves the current screen in DS9 
    # and has a number of limitations and flaws, I will use f2n
    # which I wrapped with some settings in make_png, 
    # incorporated for all frames in ccdcalib.py
   
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

    # Do the photometry
    tsleep = 3                    # Arbitrary time delay
    photometry = seekfits()
    photometry(dataref, dirs, tsleep, comparisons, targets, psf_fwhm)
    
if  __name__ == "__main__":

    import sys
    plt.ion()
    xpapoint       = sys.argv[1]
    run_rtphos(xpapoint)
