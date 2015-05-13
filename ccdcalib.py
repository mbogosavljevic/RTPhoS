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
import sys
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


# Write a FITS file with updated headers.
def writefits(data, hdr, filename):
        
    hdr_out = hdr.copy(strip=True)  
    hdu_out = pyfits.PrimaryHDU(data, hdr_out)
#   hdu_out.scale(type='float32', bzero=32768, bscale=1)    # For compatibility use 32 bits per data pixel. This does not seem to work as it rounds up and precision is lost.
    hdu_out.writeto(filename, output_verify='warn', clobber=True)
    del data, hdr, hdr_out                  # Free some memory
    return



# Create a Masterbias frame from available bias frames.
#-------------------------------------------------------------
def makebias(dataref, dsize):

    # Search for all bias files in the current directory
    # Warning: will search for all files with the text 'bias' in their filename.
    biasfiles = sorted(glob.glob('*bias*'))

    # Check to see that all selected files have a fits file name extension
    biasfiles = [f for f in biasfiles if f.endswith('.fits') or f.endswith('.fit')]
    biasnum = len(biasfiles)     
    print "Found ", biasnum, "bias frames"
    if biasnum==0:
       masterbias=0.0
       print "WARNING: No bias frames found data will not be bias calibrated!"
       return masterbias

    # Check to see that all the bias frames are of the same size.
    # First make a list of all the bias frame sizes that are available
    bsizes=[]
    goodbiasfiles=[]
    for i in range(biasnum):
        dum = pyfits.getdata(biasfiles[i])
        dummy = np.shape(dum)
        bsizes.append(dummy)
        if dummy == dsize:
           goodbiasfiles.append(biasfiles[i]) 
    
    del dum, dummy # Memory freedom!

    # Convert size tuples to strings so that np.unique will see them as pairs.
    bsizes = map(str, bsizes)
    #print bsizes

    # Next find how many unique size groups there are.
    uniquevals = np.unique(bsizes)
    biasvers = np.size(uniquevals)
    
    # If more than 1 size available inform the user.
    if biasvers>1:
       print "Found bias frames of ",biasvers," size(s): ",uniquevals

    # Count how many frames found with the same size as the data
    matches = len(goodbiasfiles)
    print "Number of matches: ", matches

    # If no bias frames are found with the same size as the data then exit
    # This should be made to check if the bias frames found are larger and
    # if they are the program should proceed with making a masterbias frame
    # and then crop it (if the header keywords exist)
    if matches == 0:
       masterbias=0.0
       print "WARNING: Bias frames are of different shape than data frames!"
       print "         Bias calibration has not be performed!"
       return masterbias

    # Make lists of all NAXIS1 and NAXIS2 header values.
    # nx1 = makelist(biasfiles,'NAXIS1')
    # nx2 = makelist(biasfiles,'NAXIS2')


    #testsize = np.unique(naxis1)
    #print testsize, np.size(testsize)

    #biasnx1 = dict(zip(biasfiles,naxis1))

    #print biasnx1

    # Set up the array that will hold all the bias images.
    bias_images = np.ndarray((biasnum,dsize[0],dsize[1]),dtype=float)


    # Go round this loop and fill the bias images array.
    for i in range(0,matches):
#       print 'Reading bias file',i
        bias_images[i] = pyfits.getdata(goodbiasfiles[i])

    # If only one file is found skip the median routine and set it as the
    # masterbias files.
    if matches == 1:
       masterbias = pyfits.getdata(goodbiasfiles[0])
    else:
       # Median combine the bias frame
       masterbias = np.median(bias_images, axis=0)

       print "Bias Median: ", np.median(masterbias)

    # Construct a COMMENT keyword and add it to the output bias header
    # First, get the header of the first file to use as a header for the output file.
    hdr_bias = pyfits.getheader(biasfiles[0])
    hdr_out = hdr_bias.copy(strip=True)
    # Make a text string with the current date and time
    biastxt = "Bias frame created by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Update the output header
    hdr_out['COMMENT'] = (biastxt)          

    # Output the masterbias frame.
    # Filename hardwired as masterbias.fits
    writefits(masterbias, hdr_out, 'masterbias.fits')

    # Memory Freedom!  
    del bias_images, biasfiles

    return masterbias



# Create a Masterdark frame from available dark frames.
#-------------------------------------------------------------
def makedark(dataref, dsize):
    
    # Search for all dark files in the current directory
    # Warning: will search for all files with the text 'dark' in their filename.
    darkfiles = sorted(glob.glob('*dark*'))
    darknum = len(darkfiles)                
    print "Found ",darknum, "dark frames"
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
             print "WARNING: Bias frames are of different size than data frames!"
             print "         Bias calibration has not be performed!"
          else:
             dataref = dataref - masterbias
             print "Bias calibration performed!"
             biascal = True
       else:
          masterbias = makebias(dataref, dsize)
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

