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

- Check to see that all calibration files are of the same shape (Currently this 
is done only for the bias frames)

- Check to see if the darks correspond to the same exposure time as the data (and
the flats). If they do not correspond to the same exposure then scale the dark
accordingly before applying it. 

- If directory contains darks of different exposure times then create masterdarks
for each exposure time, determine which one is more appropriate for the data frame and
either apply it directly or scale it and apply it.

- Check to see that the flats correspond to the same filter as the data. If the
directory contains flats for more than one filter then create masterflats for each
filter, determine which one corresponds to the data and apply it.

- Find out if there are any pixels with zero value in the flats and flag them as
bad pixels.

- Include the cosmic ray detection script cosmics to clean the image.

"""
import astropy.io.fits as pyfits
import glob
import numpy as np
import os
import datetime


# Make a dictionary of selected header keywords and their values.
def makechecklist(hdr):

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

    # Get the Date from the header
    if "DATE-OBS" in keylist: checklist['DATE'] = hdr['DATE-OBS']

    return checklist



# Given a keyword this function will return the value from a FITS header
def gethdrval(header,key):

    value=header[key]

    return value

# Make a list of fits filenames and makes a list of header values
def makelist(filenames, key):

    listout=[]
    for i in range(len(filenames)):
        dum=pyfits.getheader(filenames[i])
        dummy=gethdrval(dum,key)
        listout.append(dummy)

    del dum, dummy

    return listout


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

    hdr_out = hdr.copy(strip=True)  
    hdu_out = pyfits.PrimaryHDU(data, hdr_out)
#   hdu_out.scale(type='float32', bzero=32768, bscale=1)    # For compatibility use 32 bits per data pixel. This does not seem to work as it rounds up and precision is lost.
    hdu_out.writeto(filename, output_verify='warn', clobber=True)
    del data, hdr, hdr_out                  # Free some memory
    return


#===============================================================================
# Median combine images read from a file list and writes them to an output file.
#===============================================================================
def mediancomb(filenamesin, filenameout, hdrtext):

# ------------------------------------------------------------------------------
# Inputs are:
# filenamesin   - Type: List    - List of filenames for data input
# filenamesout  - Type: List    - List of filenames for writing output files
# hdrtext       - Type: String  - String to construct a COMMENT keyword in the output files.

# Outputs are:
# median_image  - Type: Array   - A 2-D image array

# This routine depends on:
# writefits     - to write output files.

# Modules required:
# - astropy.io.fits
# - numpy
# - datetime
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

    # Construct a COMMENT keyword and add it to the output image header.
    # First, get the header of the first file to use as a header for the output file.
    hdr_in = pyfits.getheader(filenamesin[0])
    hdr_out = hdr_in.copy(strip=True)

    # Append the supplied text string with the current date and time
    hdrtext = hdrtext + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Update the output header
    hdr_out['COMMENT'] = (hdrtext)          

    # Output the image to a FITS file using the output filename supplied.
    writefits(median_image, hdr_out, filenameout)

    # Memory Freedom!  
    del images

    return median_image


#===============================================================================
# Select FITS files based on a filename wildcard.
#===============================================================================
def makefilelist(wildcard):

# ------------------------------------------------------------------------------
# Input arguments are:
# wildcard   - Type: String    - A wildcard string to search for files.

# Output arguments are:
# filelist   - Type: List      - A list of filenames

# Modules required:
# - glob
# ------------------------------------------------------------------------------

    # Search for all files in the current directory that satisfy the wildcard provided
    filelist = sorted(glob.glob(wildcard))

    # Discard the files without a .fit or .fits extension and count what is left.
    filelist = [f for f in filelist if f.endswith('.fits') or f.endswith('.fit')]
    filenum = len(filelist)     
    print "Found ", filenum, wildcard, "frames!"

    return filelist


#===============================================================================
# Create a Masterbias frame from available bias frames.
#===============================================================================
def makebias(dsize):

# ------------------------------------------------------------------------------
# Inputs:
# dsize      - Type: Tuple    - The size of the original image data array

# Outputs:
# masterbias - Type: Array   - A 2-D image array 

# This routine depends on:
# mediancomb - Median combine images read from a list of files.

# Modules required:
# - glob
# - astropy.io.fits
# - numpy
# ------------------------------------------------------------------------------

    # Get all the available bias files into a file list.
    biasfiles = makefilelist('*bias*')
    biasnum = len(biasfiles)

    if biasnum==0:
       masterbias=0.0
       print "WARNING: No bias frames found data will not be bias calibrated!"
       return masterbias

    # Check to see that all the bias frames are of the same size.
    # First make a list of all the bias frame sizes that are available
    bsizes=[]
    testdict = {}
    for i in range(biasnum):
        dum = pyfits.getdata(biasfiles[i])
        dummy = np.shape(dum)
        bsizes.append(dummy)

    del dum, dummy # Memory freedom!

    # Convert size tuples to strings so that np.unique will see them as pairs.
    strbsizes = map(str, bsizes)

    # Next find how many unique size groups there are.
    uniquevals = np.unique(strbsizes)
    biasvers = np.size(uniquevals)

    # Convert numpy array type to a list and get the matching size index in the list.
    # If no matching size is found then exit and return a zero bias frame.
    listvals = np.array(uniquevals).tolist()
    try:
        pos = listvals.index(str(dsize))
    except ValueError:
        masterbias=0.0
        print "WARNING: No matching bias frames found! Data will not be bias calibrated!"
        return masterbias

    # If more than 1 size available inform the user.
    if biasvers>1:
       print "Found bias frames of ",biasvers," size(s): ",uniquevals

    # Build the bias file dictionary. It will be of the form:
    # ['size1': [filenames array1], 'size2': [filenames array2]]
    biasdict = {}
    for i in uniquevals:
        biasdict[i] = []
        for j in biasfiles:
            dum = pyfits.getdata(j)
            dummy = np.shape(dum)
            if str(dummy)==i:
                #print "Found ", str(dummy), " in ", i, j
                biasdict[i].append(j)

    del dum, dummy # Memory freedom!

    # Make lists of the filenames for every available size 
    filelists = []
    filelists = biasdict.values()
 
    # Input and output filenames and header text construct.
    biasfiles = filelists[pos]
    outfilename = 'masterbias.fits'
    biastxt = "Bias frame created by RTPhoS on "

    print "Match Found! Creating masterbias using ",len(biasfiles)," available frames!"

    # Create the median image.
    masterbias = mediancomb(biasfiles, outfilename, biastxt)

    # For now, create a masterbias based on matching only the size of the image.
    # In the future we need to modify the code so that if the bias frames found
    # are larger than the incoming image the code will look for header keywords
    # that will allow for cropping the bias image to match the data image.
    # The code commented out below can help in this process by creating
    # several masterbias frames for every different frame size found.

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

    return masterbias


#===============================================================================
# Create a Masterdark frame from available dark frames.
#===============================================================================
def makedark(dsize):
    
    # Get all the available dark files into a file list.
    darkfiles = makefilelist('*dark*')
    darknum = len(darkfiles)
    if darknum==0:
       masterdark=0.0
       print "WARNING: No dark frames found data will not be bias calibrated!"
       return masterdark

    # Determine the size of the image based on the size of the first dark file.
    nx = pyfits.getval(darkfiles[0], 'NAXIS1')
    ny = pyfits.getval(darkfiles[0], 'NAXIS2')

    # Set up the array that will hold all the dark images.
    dark_images = np.ndarray((darknum,nx,ny),dtype=float)

    # Go round this loop and fill the bias images array.
    for i in range(0,darknum):
#       print 'Reading file',i
        dark_images[i] = pyfits.getdata(darkfiles[i])

    # Median combine the dark frame
    masterdark = np.median(dark_images, axis=0)

    # Check to see if a bias exists and if it does subract it from the dark.
    if os.path.isfile('masterbias.fits'):   
       biasimg = pyfits.open('masterbias.fits')
       masterbias = biasimg[0].data
       masterdark = masterdark - masterbias
    # If the bias does not exist, create it.
    else:
       masterbias = makebias(dataref)
       masterdark = masterdark - masterbias

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
    writefits(masterdark, hdr_out, 'masterdark.fits')
  
    # Free the memory!  
    del dark_images, darkfiles

    return masterdark

   #testsize = np.unique(naxis1)
    #print testsize, np.size(testsize)

    #biasnx1 = dict(zip(biasfiles,naxis1))

    #print biasnx1


# Create a Masteflat frame from available flat field frames.
#-------------------------------------------------------------
def makeflat(dataref, dsize):
    
    # Search for all flat field files in the current directory
    # Warning: will search for all files with the text 'flat' in their filename.
    flatfiles = sorted(glob.glob('*flat*'))
    flatnum = len(flatfiles)                
    print "Found ",flatnum, "flat frames"
    if flatnum==0:
       masterflat=1.0
       print "WARNING: No flat field frames found data will not be flat fielded!"
       return masterflat

    # Determine the size of the image based on the size of the first flat field file.
    nx = pyfits.getval(flatfiles[0], 'NAXIS1')
    ny = pyfits.getval(flatfiles[0], 'NAXIS2')

    # Set up the array that will hold all the flat field images.
    flat_images = np.ndarray((flatnum,nx,ny),dtype=float)

    # Go round this loop and fill the flat field images array.
    for i in range(0,flatnum):
#       print 'Reading file',i
        flat_images[i] = pyfits.getdata(flatfiles[i])

    # Median combine the flat field frames
    masterflat = np.median(flat_images, axis=0)

    # Check to see if a masterbias exists and if it does subract it from the flat field.
    if os.path.isfile('masterbias.fits'):   
       biasimg = pyfits.open('masterbias.fits')
       masterbias = biasimg[0].data
       masterflat = masterflat - masterbias
    # If the bias does not exist, create it.
    else:
       masterbias = makebias(dataref)
       masterflat = masterflat - masterbias

    # Check to see if a masterdark exists and if it does subract it from the flat field.
    if os.path.isfile('masterdark.fits'):   
       darkimg = pyfits.open('masterdark.fits')
       masterdark = darkimg[0].data
       masterflat = masterflat - masterdark
    # If the dark does not exist, create it.
    else:
       masterdark = makedark(dataref)
       masterflat = masterflat - masterdark

    flatmedian = np.median(masterflat)
    print "Flat Median: ", np.median(masterflat)
    masterflat = masterflat/flatmedian
    print "Nomalized Flat Median: ", np.median(masterflat)

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
    writefits(masterflat, hdr_out, 'masterflat.fits')

    # Free the memory!  
    del flat_images, flatfiles

    return masterflat



# ---------------CALIBRATION INITIALIZATION-------------------------------------
def calib(ref_filename, dataref, hdr_data):

    print "Checking image calibration..."

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

    # Set calibration flags
    biascor = False 
    darkcor = True
    flatcor = True

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
    biascal = False
    if not biascor:
       print "Frame is not Bias calibrated"
       print "Proceeding with removing the bias..."
       if os.path.isfile('masterbias.fits'):   
          biasimg = pyfits.open('masterbias.fits')
          masterbias = biasimg[0].data
          if np.shape(masterbias) != dsize:
             masterbias=0.0
             print "WARNING: Masterbias frame is of different size than data frames!"
             print "         Bias calibration has not be performed!"
          else:
             dataref = dataref - masterbias
             print "Bias calibration performed!"
             biascal = True
       else:
          masterbias = makebias(dsize)
          dataref = dataref - masterbias
       if masterbias.all() !=0.0:
             biascal = True
             biastxt = "Bias was subtracted by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    darkcal = False
    if not darkcor:
       print "Frame is not Dark calibrated"
       print "Proceeding with removing the dark..."
       if os.path.isfile('masterdark.fits'):   
          darkimg = pyfits.open('masterdark.fits')
          masterdark = darkimg[0].data
          if np.shape(masterdark) != dsize:
             masterdark=0.0
             print "WARNING: Dark frames are of different size than data frames!"
             print "         Dark calibration has not be performed!"
          else:
             dataref = dataref - masterdark
             print "Bias calibration performed!"
             darkcal = True
       else:
          masterdark = makedark(dataref, dsize)
          dataref = dataref - masterdark
          print "Dark calibration performed!"
       if masterdark.all() !=0.0:
          darkcal = True
          darktxt = "Dark was subtracted by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    flatcal = False
    if not flatcor:
       print "Frame is not Flat fielded"
       print "Proceeding with Flat fielding..."
       if os.path.isfile('masterflat.fits'):   
          flatimg = pyfits.open('masterflat.fits')
          masterflat = flatimg[0].data
          if np.shape(masterflat) != dsize:
             masterflat=1.0
             print "WARNING: Masterflat frame is of different size than data frames!"
             print "         Flat fielding has not be performed!"
          else:
             dataref = dataref/masterbias
             print "Flat fielding performed!"
             flatcal = True
       else:
          masterflat = makeflat(dataref, dsize)
          dataref = dataref/masterflat
          print "Flat fielding performed!"
       if masterflat.all() !=1.0:
          flatcal = True
          flattxt = "Flatfield applied by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Update the headers. 
    hdr_out = hdr_data.copy(strip=True)
    if biascal: hdr_out['BIASCOR'] = (biastxt)          
    if darkcal: hdr_out['DARKCOR'] = (darktxt)
    if flatcal: hdr_out['FLATCOR'] = (flattxt)
    # Write the output to disk giving the prefix of 'c_' to the calibrated frame.
    writefits(dataref, hdr_out, 'c_'+ref_filename)

    return   

if __name__ == "__main__":

#  For testing
   ref_filename = 'J1753.fits'
   dataref, hdr_data = pyfits.getdata(ref_filename, header=True)     
   calib(ref_filename, dataref, hdr_data)
   print "Calibration Done!"

