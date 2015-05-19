#!/usr/bin/python
#
"""
Currently the code below will check to see if a data frame has been calibrated
by looking to see if the IRAF FITS Headers BIASCOR, DARKCOR and FLATCOR exist.
If they exist the code will assume that the data is calibrated and will do nothing.
If it hasn't been calibrated it will attempt to calibrate the frame
by looking for bias, dark and flat frames in the same directory. It will then
create the appropriate calibration files and apply them to the data image.
The calibrated image is saved with the prefix 'c_' in the same directory.

Things that need to be done to further improve the calibration process:

- If flats of the same filter are found but they have different exposures calculate
each flat seperately and then combine them.

- Find out if there are any pixels with zero value in the flats and flag them as
bad pixels.

- Include the cosmic ray detection script cosmics to clean the image.

- Put in various other checks such as observation epoch.

"""
import astropy.io.fits as pyfits
import glob
import numpy as np
import os
import datetime
#import sys

#===============================================================================
# Make a dictionary of selected header keywords and their values.
#===============================================================================
def makechecklist(hdr):

# ------------------------------------------------------------------------------
# Inputs are:
# hdr       - Type: Class      - Header class created by astropy.io.fits

# Outputs are:
# checklist - Type: Dictionary - Dictionary of selected keywords

# Modules required:
# - astropy.io.fits
# ------------------------------------------------------------------------------       

    checklist = {}
    keylist = hdr.keys()

    # Get the image size header values
    checklist['NAXIS1'] = hdr['NAXIS1']
    checklist['NAXIS2'] = hdr['NAXIS2']

    # Get the filter type from the header
    if "FILTER" in keylist: checklist['FILTER'] = hdr['FILTER']

    # Get the integration time from the header
    if "EXPOSURE" in keylist:  
       checklist['EXPOSURE'] = hdr['EXPOSURE']
    elif "EXPTIME" in keylist: 
       checklist['EXPOSURE'] = hdr['EXPTIME']
    else:
       checklist['EXPOSURE'] = 0
       print "WARNING: Exposure keyword not found!"
       print "         This might cause dark subtraction and flat fielding errors!"         

    # Get the Date from the header
    if "DATE-OBS" in keylist: checklist['DATE'] = hdr['DATE-OBS']

    return checklist



#===============================================================================
# Write a FITS file with updated headers.
#===============================================================================
def writefits(data, hdr, filename):

# ------------------------------------------------------------------------------
# Inputs are:
# data      - Type: Array   - 2-D image data to be written to a FITS file.
# hdr       - Type: Class   - Header class created by astropy.io.fits
# filename  - Type: String  - Output filename.

# Outputs are:
# Returns nothing. Writes a file to disk.

# Modules required:
# - astropy.io.fits
# ------------------------------------------------------------------------------       

    # Make the FITS header for the output file.
    hdr_out = hdr.copy(strip=True)  
    hdu_out = pyfits.PrimaryHDU(data, hdr_out)
#   hdu_out.scale(type='float32', bzero=32768, bscale=1)    # For compatibility use 32 bits per data pixel. This does not seem to work as it rounds up and precision is lost.
    hdu_out.writeto(filename, output_verify='warn', clobber=True)
    del data, hdr, hdr_out                  # Free some memory
    return


#===============================================================================
# Median combine images read from a file list.
#===============================================================================
def mediancomb(filenamesin):

# ------------------------------------------------------------------------------
# Inputs are:
# filenamesin   - Type: List    - List of filenames for data input

# Outputs are:
# median_image  - Type: Array   - A 2-D image array

# Modules required:
# - astropy.io.fits
# - numpy
# ------------------------------------------------------------------------------

    # Count the number of filenames
    filenum = len(filenamesin)

	# Determine the size of the image based on the size of the first file.
    nx = pyfits.getval(filenamesin[0], 'NAXIS1')
    ny = pyfits.getval(filenamesin[0], 'NAXIS2')

	# Set up the array that will hold all the images.
    images = np.ndarray((filenum,ny,nx),dtype=float)

    # Go round this loop and fill the images array.
    for i in range(filenum):
        images[i] = pyfits.getdata(filenamesin[i])

    # If only one file is found skip the median calculation and set itself as
    # the output file.
    if filenum == 1:
       median_image = pyfits.getdata(filenamesin[0])
    else:
       # Median combine all the images.
       median_image = np.median(images, axis=0)

    # Memory Freedom!  
    del images

    return median_image


#===============================================================================
# Select FITS files based on a filename wildcard.
#===============================================================================
def makefilelist(wildcard, directory="."):

# ------------------------------------------------------------------------------
# Input arguments are:
# wildcard   - Type: String    - A wildcard string to search for files.
# directory  - Type: String    - Directory of the filenames (optional).

# Output arguments are:
# filelist   - Type: List      - A list of filenames

# Modules required:
# - glob
# - os
# ------------------------------------------------------------------------------

    # Move to the specified directory
    os.chdir(directory)

    # Search for all files in the current directory that satisfy the wildcard provided
    filelist = sorted(glob.glob(wildcard))

    # Discard the files without a .fit or .fits extension and count what is left.
    filelist = [f for f in filelist if f.endswith('.fits') or f.endswith('.fit')]
    filenum = len(filelist)     
    print "Found ", filenum, wildcard, "frames!"

    return filelist


#===============================================================================
# Groups FITS files according to their image array size and compares them to
# a given size. Returns a file list with file names that correspond to the files
# that match the required size.
#===============================================================================
def checksize(filelistin, dsize, calltxt):

# ------------------------------------------------------------------------------
# Inputs:
# filelistin    - Type: List     - The list with input filenames
# dsize         - Type: Tuple    - The size of the original image data array
# calltxt       - Type: String   - A flag showing where the routine was called from

# Outputs:
# required_filelist - Type: List - A list containing the matching size filenames 

# This routine depends on:
# mediancomb - Median combine images read from a list of files.

# Modules required:
# - astropy.io.fits
# - numpy
# ------------------------------------------------------------------------------

    # First make a list of all the frame sizes that are available
    sizes=[]
    filenum = len(filelistin)
    for i in range(filenum):
        dum = pyfits.getdata(filelistin[i])
        dummy = np.shape(dum)
        sizes.append(dummy)

    # Convert size tuples to strings so that np.unique will see them as pairs.
    strsizes = map(str, sizes)

    # Next find how many unique size groups there are.
    uniquevals = np.unique(strsizes)
    sizevers = np.size(uniquevals)

    # Convert numpy array type to a list and get the matching size index in the list.
    # If no matching size is found then return an empty file list.
    required_filelist = []
    listvals = np.array(uniquevals).tolist()
    try:
        pos = listvals.index(str(dsize))
    except ValueError:
        print "WARNING: No matching ", calltxt, " frames found!" 
        print "         Data will not be ",calltxt," calibrated!"
        return required_filelist

    # If more than 1 size available inform the user.
    if sizevers>1:
       print "Found ",calltxt," frames of ",sizevers," size(s): ",uniquevals

    # Build a file dictionary. It will be of the form:
    # ['size1': [filenames array1], 'size2': [filenames array2]]
    filedict = {}
    for i in uniquevals:
        filedict[i] = []
        for j in filelistin:
            dum = pyfits.getdata(j)
            dummy = np.shape(dum)
            if str(dummy)==i:
                #print "Found ", str(dummy), " in ", i, j
                filedict[i].append(j)

    del dum, dummy # Memory freedom!

    # Make lists of the filenames for every available size 
    filelists = []
    filelists = filedict.values()
 
    # Select the filelist that matches the required size
    required_filelist = filelists[pos]

    # For now, create calibration files based on matching only the size of the image.
    # In the future we need to modify the code so that if there are any calibration frames found
    # that are larger than the incoming image the code will look for header keywords
    # that will allow for cropping the calibration images to match the data image.
    # The code commented out below can help in this process by creating
    # several master calibration frames for every different calibration frame size found.

    # Make a list of output filenames.
#    filesout = []
#    if biasvers>1:
#        for i in range(biasvers):
#            filesout.append('masterbias_%i.fits' %(i+1)) 
#    else:
#        filesout.append('masterbias.fits')

    # Construct the COMMENT header text to be inserted to the output files.
#    biastxt = "Bias frame created by RTPhoS on "

    # Now go round this loop and create the master bias frames for every
    # available size.
#    for i in range(biasvers):
#        biasfiles = filelists[i]
#        outfilename = filesout[i]
#        mediancomb(biasfiles, outfilename, biastxt)

    return required_filelist



#===============================================================================
# Create a Masterbias frame from available bias frames.
#===============================================================================
def makebias(dsize, dirs):

# ------------------------------------------------------------------------------
# Inputs:
# dsize       - Type: Tuple       - The size of the original image data array
# dirs        - Type: Dictionary  - Input and output directories

# Outputs:
# bias     - Type: Tuple - Tuple containing a boolean and the bias image array 

# This routine depends on:
# mediancomb - Median combine images read from a list of files.
# writefits  - Write the masterbias to a file.

# Modules required:
# - astropy.io.fits
# - numpy
# - datetime
# - os
# ------------------------------------------------------------------------------

    # Check to see that bias directory exits. If it does not exit bias calibration.
    if not os.path.exists(dirs['bias']):
       print "WARNING: Bias frames directory does not exist!"
       print "         Bias calibration aborted!"
       biascheck = False
       masterbias = 0
       return (biascheck, masterbias)

    # Move to the bias files directory
    prev_dir = os.path.abspath(os.curdir)
    os.chdir(dirs['bias'])

    # Get all the available bias files into a file list.
    biasfiles = makefilelist('*bias*')
    biasnum = len(biasfiles)

    if biasnum==0:
       biascheck = False
       masterbias = 0
       print "WARNING: No bias frames found data will not be bias calibrated!"
       os.chdir(prev_dir)
       return (biascheck, masterbias)

    # Check to see that all bias frames are of the same size and make a list
    # of all frames of the same size.
    goodbiasfiles = checksize(biasfiles, dsize, 'bias')

    # Check to see if the list returned has anything in it.
    biasnum = len(goodbiasfiles)
    if biasnum>0:
       print "Match Found! Creating masterbias using ",biasnum," available frames!"
    else:
       biascheck = False
       masterbias = 0
       os.chdir(prev_dir)
       return (biascheck, masterbias)

    # Create the median image.
    masterbias = mediancomb(goodbiasfiles)
    print "Bias Median: ", np.median(masterbias)

    # Construct a COMMENT keyword and add it to the output image header.
    # First, get the header of the first file to use as a header for the output file.
    # Output filename and header text construct.
    filenameout = 'masterbias.fits'
    hdr_in = pyfits.getheader(goodbiasfiles[0])
    hdr_out = hdr_in.copy(strip=True)

    # Append the supplied text string with the current date and time
    biastxt = "Bias frame created by RTPhos on "+ datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Update the output header
    hdr_out['COMMENT'] = (biastxt)          

    # Output the image to a FITS file using the output filename supplied.
    writefits(masterbias, hdr_out, dirs['reduced']+filenameout)
    biascheck = True
    bias = (biascheck, masterbias)

    os.chdir(prev_dir)
    return bias


#===============================================================================
# Create a Masterdark frame from available dark frames.
#===============================================================================
def makedark(dsize, exposure, dirs):
    
# ------------------------------------------------------------------------------
# Inputs:
# dsize    - Type: Tuple  - The size of the original image data array
# exposure - Type: String - The exposure of the image data
# dirs     - Type: Dictionary  - Input and output directories

# Outputs:
# dark     - Type: Tuple  - Tuple containing a boolean and the dark image array

# This routine depends on:
# mediancomb - Median combine images read from a list of files.
# masterbias - Routine to create a masterbias frame
# writefits  - Write the masterdark to a file.

# Modules required:
# - astropy.io.fits
# - numpy
# - datetime
# - os
# ------------------------------------------------------------------------------

    # Check to see that dark frames directory exits. If it does not exit dark calibration.
    if not os.path.exists(dirs['dark']):
       print "WARNING: Dark frames directory does not exist!"
       print "         Dark calibration aborted!"
       darkcheck = False
       masterdark = 0
       return (darkcheck, masterdark)

    # Move to the dark files directory
    prev_dir = os.path.abspath(os.curdir)
    os.chdir(dirs['dark'])

    # Get all the available dark files into a file list.
    darkfiles = makefilelist('*dark*')
    darknum = len(darkfiles)
    if darknum==0:
       darkcheck = False
       masterdark = 0
       print "WARNING: No dark frames found data will not be dark calibrated!"
       os.chdir(prev_dir)
       return (darckcheck, masterdark)

    # Check to see that all dark frames are of the same size and make a list
    # of all frames of the same size.
    posdarkfiles = checksize(darkfiles, dsize, 'dark')

    # Check to see if the list returned has anything in it.
    darknum = len(posdarkfiles)
    if darknum<1:
       darkcheck = False
       masterdark = 0
       os.chdir(prev_dir)
       return (darkcheck, masterdark)  

    # Check that all the dark files in the file list have the same exposure.
    # If they don't create a dictionary of filenames and exposures.
    # First make a list of all the exposures that are available
    exposures=[]
    for i in range(darknum):
        dum = pyfits.getheader(posdarkfiles[i])
        keylist = dum.keys()
        if "EXPOSURE" in keylist:  
           dummy = dum['EXPOSURE']
           keyword = 'EXPOSURE'
        elif "EXPTIME" in keylist: 
           dummy = dum['EXPTIME']
           keyword = 'EXPTIME'         
        exposures.append(dummy)

    del keylist, dum, dummy   # Memory freedom!

    # Next find how many unique exposures there are.
    uniquevals = np.unique(exposures)
    expvers = np.size(uniquevals)

    # If more than 1 exposure available inform the user.
    if expvers>1: print "Found files with ",expvers," different exposures: ",uniquevals    

    # Build a file dictionary. It will be of the form:
    # ['Exposure 1': [filenames list1], 'Exposure 2': [filenames list2]]
    filedict = {}
    for i in uniquevals:
        filedict[i] = []
        for j in posdarkfiles:
            dum = pyfits.getheader(j)
            dummy = dum[keyword]
            if dummy==i:
                #print "Found ", str(dummy), " in ", i, j
                filedict[i].append(j)

    #exposure = ' 7'     # For testing different exposures (This is a string)
    # Find if any of the dark exposures match the data exposure and make the
    # appropriate file list.
    darkexps = [x for x in uniquevals]
    darkexps_fl = [float(x) for x in uniquevals] # Convert string to float
    exposure_fl = float(exposure)                # Convert string to float
    print "Data image exposure is:",exposure
    try:
        pos = darkexps.index(exposure)
        # pos = int(pos)
        print "Exposure match found!"
        gooddarkfiles = filedict[exposure]
        nomatch = False
    except ValueError:
        # The darks do not have the same exposure as the data.
        # Will select the closest one and scale the darks accordingly.
        expdiff = [abs(x-exposure_fl) for x in darkexps_fl]
        mindiff = min(expdiff)
        # This gets the index of the minimum value of the expdiff list.
        pos = [k for k, l in enumerate(expdiff) if l == mindiff ]
        pos = int(pos[0])
        bestexp = darkexps[pos]
        gooddarkfiles = filedict[bestexp]
        # Calculate the scaling factor.
        factor = exposure_fl/float(bestexp)
        nomatch = True
        print "No similar exposure found!"
        print "Closest Dark frame exposure found is:",bestexp
        print "Dark will be adjusted by a factor of", factor

    #print gooddarkfiles

    # Create the median image.
    masterdark = mediancomb(gooddarkfiles)
    
	# Check to see if a masterbias frame exists and if it does subract it from the dark.
    if os.path.isfile(dirs['reduced']+'masterbias.fits'):	
       masterbias = pyfits.getdata(dirs['reduced']+'masterbias.fits')
       # Check to see if the bias is of the same size as the dark.
       bias_size = np.shape(masterbias)
       if bias_size == np.shape(masterdark):
           masterdark = masterdark - masterbias
       else:
           bias = makebias(dsize, dirs)
           biascheck = bias[0]
           masterbias = bias[1]
           if biascheck: masterdark = masterdark - masterbias
    else:
       bias = makebias(dsize, dirs)
       biascheck = bias[0]
       masterbias = bias[1]
       if biascheck: masterdark = masterdark - masterbias

    # If the dark needs to be scaled use the scaling factor calculated above.
    if nomatch: masterdark = masterdark*factor 

    # Output filename and header text construct.
    outfilename = 'masterdark.fits'

	# Construct a COMMENT keyword and add it to the output dark header
	# First, get the header of the first file to use as a header for the output file.
    hdr_dark = pyfits.getheader(darkfiles[0])
    hdr_out = hdr_dark.copy(strip=True)
    # Make a text string with the current date and time
    darktxt = "Dark frame created by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Update the output header
    hdr_out['COMMENT'] = (darktxt)		

    print "Dark Median: ", np.median(masterdark)
    # Output the masterdark frame.
    # Filename hardwired as masterdark.fits
    writefits(masterdark, hdr_out, dirs['reduced']+outfilename)
    darkcheck = True
    dark = (darkcheck, masterdark)

    os.chdir(prev_dir)
    return dark

#===============================================================================
# Create a Masterflat frame from available flat frames.
#===============================================================================
def makeflat(dsize, obsfilter, dirs):
  
# ------------------------------------------------------------------------------
# Inputs:
# dsize     - Type: Tuple  - The size of the original image data array
# obsfilter - Type: String - The observation filter
# dirs      - Type: Dictionary  - Input and output directories

# Outputs:
# flat      - Type: Tuple  - Tuple containing a boolean and the flat image array

# This routine depends on:
# mediancomb - Median combine images read from a list of files.
# masterbias - Routine to create a masterbias frame.
# masterdark - Routine to create a masterdark frame.
# writefits  - Write the masterdark to a file.

# Modules required:
# - astropy.io.fits
# - numpy
# - datetime
# - os
# ------------------------------------------------------------------------------

    # Check to see that flats directory exits. If it does not exit flat calibration.
    if not os.path.exists(dirs['flat']):
       print "WARNING: Flat frames directory does not exist!"
       print "         Flat fielding aborted!"
       flatcheck = False
       masterflat = 0
       return (flatcheck, masterflat)

    # Move to the flat files directory
    prev_dir = os.path.abspath(os.curdir)
    os.chdir(dirs['flat'])
 
    # Get all the available flat files into a file list.
    flatfiles = makefilelist('*flat*')
    flatnum = len(flatfiles)
    if flatnum==0:
       flatcheck = False
       print "WARNING: No flat field frames found data will not be flat fielded!"
       return masterflat

    # Check to see that all flat frames are of the same size and make a list
    # of all frames of the same size.
    posflatfiles = checksize(flatfiles, dsize, 'flat')

    # Check to see if the list returned has anything in it.
    flatnum = len(posflatfiles)
    if flatnum<1:
       flatcheck = False
       masterflat = 0
       os.chdir(prev_dir)
       return (flatcheck, masterflat)    

    # Check that all the flat field files in the file list have the same filter.
    # If they don't create a dictionary of filenames and filters.
    # First make a list of all the filters that are available.
    filters=[]
    for i in range(flatnum):
        dum = pyfits.getheader(posflatfiles[i])
        keylist = dum.keys()
        if "FILTER" in keylist: 
           dummy = dum['FILTER']
           filters.append(dummy)
        else:
           print "No FILTER header keyword for file", posflatfiles[i]

    del keylist, dum, dummy   # Memory freedom!

    # Next find how many unique filter values there are.
    uniquevals = np.unique(filters)
    filtervers = np.size(uniquevals)
    filtervals = np.array(uniquevals).tolist() # Convert array to list

    # If more than 1 filter available inform the user.
    if filtervers>1: print "Found files corresponding to",filtervers,"different filters:",uniquevals  
  
    # Build a file dictionary. It will be of the form:
    # ['Filter 1': [filenames list1], 'Filter 2': [filenames list2]]
    filedict = {}
    for i in uniquevals:
        filedict[i] = []
        for j in posflatfiles:
            dum = pyfits.getheader(j)
            dummy = dum['FILTER']
            if dummy==i:
                #print "Found ", str(dummy), " in ", i, j
                filedict[i].append(j)

    # Make lists of the filenames for every filter
    filelists = []
    filelists = filedict.values()

    # Convert numpy array type to a list and get the matching size index in the list.
    # If no matching size is found then skip flat fielding.
    #obsfilter = 'R'     # For testing different filters (This is a string)
    try:
        pos = filtervals.index(obsfilter)
    except ValueError:
        print "WARNING: No matching flat fielding frames found!" 
        print "         Data will not be flat fielded!"
        flatcheck = False
        masterflat = 0
        os.chdir(prev_dir)
        return (flatcheck, masterflat)
 
    # Select the filelist that matches the required size
    goodflatfiles = filelists[pos]
     
    # Now check the good filter files that they all have the same exposure time!
    flatnum = len(goodflatfiles) 
    exposures=[]
    for i in range(flatnum):
        dum = pyfits.getheader(goodflatfiles[i])
        keylist = dum.keys()
        if "EXPOSURE" in keylist:  
           dummy = dum['EXPOSURE']
           keyword = 'EXPOSURE'
        elif "EXPTIME" in keylist: 
           dummy = dum['EXPTIME']
           keyword = 'EXPTIME'         
        exposures.append(dummy)

    flatexp = exposures[0]
    del keylist, dum, dummy   # Memory freedom!

    # Next find how many unique exposures there are.
    uniquevals = np.unique(exposures)
    expvers = np.size(uniquevals)
    
    # If more than 1 exposure for the same filter. Quit it is getting too
    # complicated for now. Will deal with this later.
    if expvers>1: 
       print "WARNING: Found Flat frames of the same filter but with different"
       print "         exposures. I can't handle this. No flat fielding will be done!"
       flatcheck = False
       masterflat = 0
       os.chdir(prev_dir)
       return (flatcheck, masterflat)

    #print goodflatfiles

    # Create the median image.
    masterflat = mediancomb(goodflatfiles)

	# Check to see if a masterbias frame exists and if it does subract it from the flat.
    if os.path.isfile(dirs['reduced']+'masterbias.fits'):	
       masterbias = pyfits.getdata(dirs['reduced']+'masterbias.fits')
       # Check to see if the bias is of the same size as the dark.
       bias_size = np.shape(masterbias)
       if bias_size == np.shape(masterflat):
           masterflat = masterflat - masterbias
       else:
           bias = makebias(dsize, dirs)
           biascheck = bias[0]
           masterbias = bias[1]
           if biascheck: masterflat = masterflat - masterbias
    else:
       bias = makebias(dsize)
       biascheck = bias[0]
       masterbias = bias[1]
       if biascheck: masterflat = masterflat - masterbias

	# Check to see if a masterdark frame exists and if it does subract it from the flat.
    if os.path.isfile(dirs['reduced']+'masterdark.fits'):	
       masterdark, darkhdr = pyfits.getdata(dirs['reduced']+'masterdark.fits', header=True)
       # Check to see if the dark is of the same size as the flat.
       dark_size = np.shape(masterdark)
       if dark_size == np.shape(masterflat):
          # Check to see if the Flat exposure time is the same as the dark.
          keylist = darkhdr.keys()
          if "EXPOSURE" in keylist:  
            darkexp = darkhdr['EXPOSURE']
          elif "EXPTIME" in keylist: 
            darkexp = darkhdr['EXPTIME']
          if darkexp == flatexp:
            masterflat = masterflat - masterdark
          else:
            # Adjust the dark frame to match the flat exposure time.
            darkexp_fl=float(darkexp)
            flatexp_fl=float(flatexp)            
            factor = flatexp_fl/darkexp_fl
            masterdark = masterdark * factor
            masterflat = masterflat - masterdark
       else:
          dark = makedark(dsize, flatexp, dirs)
          darkcheck = dark[0]
          masterdark = dark[1]
          if darkcheck: masterflat = masterflat - masterdark
    else:
       dark = makedark(dsize, flatexp, dirs)
       darkcheck = dark[0]
       masterdark = dark[1]
       if darkcheck: masterflat = masterflat - masterdark

    # Now find the median of the flat frame and divide by it to create a 
    # a coefficient frame.
    flatmedian = np.median(masterflat)
    print "Flat Median: ", np.median(masterflat)
    masterflat = masterflat/flatmedian
    print "Nomalized Flat Median: ", np.median(masterflat)

    # Output filename and header text construct.
    outfilename = 'masterflat.fits'

    # Construct a COMMENT keyword and add it to the output flat field header
    # First, get the header of the first file to use as a header for the output file.
    hdr_flat = pyfits.getheader(flatfiles[0])
    hdr_out = hdr_flat.copy(strip=True)
    # Make a text string with the current date and time
    flattxt = "Flat field frame created by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Update the output header
    hdr_out['COMMENT'] = (flattxt)          

    # Output the masterflat frame.
    # Filename hardwired as masterflat.fits
    writefits(masterflat, hdr_out, dirs['reduced']+outfilename)
    flatcheck = True
    flat = (flatcheck, masterflat)

    os.chdir(prev_dir)
    return flat


#===============================================================================
# ----------------------CALIBRATION INITIALIZATION------------------------------
#===============================================================================
def calib(dirs, ref_filename, dataref, hdr_data):

# ------------------------------------------------------------------------------
# Inputs:
# dirs         - Type: Dictionary - Holds all the relevant directory paths
# ref_filename - Type: String - The image data filename
# dataref      - Type: Array  - The image 2-D data array
# hdr_data     - Type: Class  - The image header

# Outputs:
# Calibrates an image are creates calibration master frames and writes the
# calibrated image to disk.

# This routine depends on:
# makechecklist - Routine that makes a dictionary of selected keywords and their values.
# masterbias - Routine to create a masterbias frame.
# masterdark - Routine to create a masterdark frame.
# masterflat - Routine to create a masterflat frame.
# writefits  - Write the masterdark to a file.

# Modules required:
# - astropy.io.fits
# - numpy
# - datetime
# ------------------------------------------------------------------------------

    print "Checking image calibration..."

    # Strip the filename of its path
    ref_filename = os.path.basename(ref_filename)

    # Check to see that this is a 2-D image if not stop.
    if hdr_data['NAXIS'] != 2:
       print "WARNING: Not a 2D image file! Proceeding to next frame..."
       return    

    # Get the size of the image.
    dsize = np.shape(dataref)

    # Create the list of headers to check for file compatibility.
    datacheck = makechecklist(hdr_data)

    # Construct a list of all the header keywords.
    keywordlist = hdr_data.keys()

    # Get the Exposure and Filter values from the data header 
    exposure  = datacheck['EXPOSURE']  
    obsfilter = datacheck['FILTER']

    # Set calibration flags. Default setting should be False.
    # Change to True if you want to exclude some parts of the code for testing.
    biascor = False
    darkcor = False
    flatcor = False

    # Look for these IRAF keywords, if they exist assume that the image has 
    # been calibrated.
    if "ZEROCOR" in keywordlist:
       biascor = True
    if "DARKCOR" in keywordlist:
       darkcor = True
    if "FLATCOR" in keywordlist:
       flatcor = True

    # If all keywords exist assume that the frame is calibrated and return to 
    # the main program.
    if biascor and darkcor and flatcor == True:
       return

    # Do the calibration according to which part is missing. The calibration assumes
    # that the calibration frames masterbias, masterdark and masterflat have already
    # been created and are in the same directory. Their names are hardwired in the
    # code. If no calibration frames are found then they are made using any available
    # calibration files.
    biascheck = False
    if not biascor:
       print "Frame is not Bias calibrated"
       print "Proceeding with removing the bias..."
       if os.path.isfile(dirs['reduced']+'masterbias.fits'):   
          masterbias = pyfits.getdata(dirs['reduced']+'masterbias.fits')
          if np.shape(masterbias) != dsize:
             biascheck = False
             print "WARNING: Masterbias frame is of different size than data frames!"
             print "         Attempting to make new masterbias frame..."
             bias = makebias(dsize, dirs)
             biascheck = bias[0]
             masterbias = bias[1]
          else:
             biascheck = True
       else:
          bias = makebias(dsize, dirs)
          biascheck = bias[0]
          masterbias = bias[1]
       if biascheck:
          dataref = dataref - masterbias
          biastxt = "Bias was subtracted by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
          print "Bias calibration successfull!"

    darkcheck = False
    if not darkcor:
       print "Frame is not Dark calibrated"
       print "Proceeding with removing the dark..."
       if os.path.isfile(dirs['reduced']+'masterdark.fits'):   
          masterdark = pyfits.getdata(dirs['reduced']+'masterdark.fits')
          if np.shape(masterdark) != dsize:
             darkcheck = False
             print "WARNING: Masterdark frame is of different size than data frames!"
             print "         Attempting to make new masterdark frame..."
             dark = makedark(dsize, exposure, dirs)
             darkcheck = dark[0]
             masterdark = dark[1]
          else:
             darkcheck = True
       else:
          dark = makedark(dsize, exposure, dirs)
          darkcheck = dark[0]
          masterdark = dark[1]
       if darkcheck:
          dataref = dataref - masterdark
          darktxt = "Dark was subtracted by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
          print "Dark calibration successfull!"

    flatcheck = False
    if not flatcor:
       print "Frame is not flat fielded!"
       print "Proceeding with flat fielding..."
       if os.path.isfile(dirs['reduced']+'masterflat.fits'):   
          masterdark = pyfits.getdata(dirs['reduced']+'masterflat.fits')
          if np.shape(masterflat) != dsize:
             flatcheck = False
             print "WARNING: Masterflat frame is of different size than data frames!"
             print "         Attempting to make new masterflat frame..."
             flat = makeflat(dsize, obsfilter, dirs)
             flatcheck = flat[0]
             masterflat = flat[1]
          else:
             flatcheck = True
       else:
          flat = makeflat(dsize, obsfilter, dirs)
          flatcheck = flat[0]
          masterflat = flat[1]          
       if flatcheck:
          dataref = dataref/masterflat
          flattxt = "Frame was flat fielded by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
          print "Flat fielding successfull!"
 
    # Update the headers. 
    hdr_out = hdr_data.copy(strip=True)
    if biascheck: hdr_out['BIASCOR'] = (biastxt)          
    if darkcheck: hdr_out['DARKCOR'] = (darktxt)
    if flatcheck: hdr_out['FLATCOR'] = (flattxt)

    # Write the output to disk giving the prefix of 'c_' to the calibrated frame.    
    if biascheck or darkcheck or flatcheck: writefits(dataref, hdr_out, dirs['reduced']+'c_'+ref_filename)

    return   

if __name__ == "__main__":

#  For testing    

   # Set up input and output directories
   current_dir = os.path.abspath(os.curdir)
   bias_dir = current_dir+"/bias/"
   dark_dir = current_dir+"/dark/"
   flat_dir = current_dir+"/flat/"
   data_dir = current_dir+"/data/"
   reduced_dir = data_dir+"/reduced/"
   if not os.path.exists(reduced_dir): os.makedirs(reduced_dir)

   # Make a dictionary with all the directories
   dirs = {'current':current_dir, 'bias':bias_dir, 'dark':dark_dir, 'flat':flat_dir, 'data':data_dir, 'reduced':reduced_dir}

   ref_filename = dirs['data']+'J1753.fits'
   dataref, hdr_data = pyfits.getdata(ref_filename, header=True)     
   calib(dirs, ref_filename, dataref, hdr_data)
   print "Calibration Done!"

