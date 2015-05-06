#!/usr/bin/python
#

import pyfits
import glob
import numpy as np
import sys
import os

# Create a Masterbias frame from available bias frames.
#-------------------------------------------------------------
def makebias(dataref):
	# Initialize biasframes dictionary
	biasframes = {}
    
	# Search for all bias files in the current directory
	# Warning: will search for all files with the text 'bias' in their filename.
	biasfiles = sorted(glob.glob('*bias*'))
	biasnum = len(biasfiles)	
	print "Found ", biasnum, "bias frames"
	if biasnum==0:
		masterbias=0.0
		print "Warning: No bias frames found data will not be bias calibrated!"
		return masterbias

	# Determine the size of the image based on the size of the first bias file.
	nx = pyfits.getval(biasfiles[0], 'NAXIS1')
	ny = pyfits.getval(biasfiles[0], 'NAXIS2')

	# Set up the array that will hold all the bias images.
	bias_images = np.ndarray((biasnum,nx,ny),dtype=float)

	# Go round this loop and fill the bias images array.
	for i in range(0,biasnum):
   		print 'Reading file',i
    		bias_images[i] = pyfits.getdata(biasfiles[i])

	# Median combine the bias frame
	masterbias = np.median(bias_images, axis=0)
	print "Bias Median: ", np.median(masterbias)

	# Output the masterbias frame.
	# Filename hardwired as masterbias.fits
	hdu=pyfits.PrimaryHDU(masterbias)
	hdu.writeto('masterbias.fits')

	# Free the memory!  
	del bias_images

	return masterbias

# Create a Masterdark frame from available bias frames.
#-------------------------------------------------------------
def makedark(dataref):
	# Initialize darkframes dictionary
	darkframes = {}
    
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
		print 'Reading file',i
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

	print "Dark Median: ", np.median(masterdark)
	# Output the masterdark frame.
	# Filename hardwired as masterdark.fits
	hdu=pyfits.PrimaryHDU(masterdark)
	hdu.writeto('masterdark.fits')
  
	# Free the memory!  
	del dark_images

	return masterdark

#ref_filename = win.get("file")
ref_filename = 'J1753.fits'
dataref, hdr = pyfits.getdata(ref_filename, header=True)     
print("... Working ...")

# Check if the file has been processed (bias,dark,flat)
# (will use standard IRAF keywords for this)
keywordlist = hdr.keys()
biascor = darkcor = flatcor = False  # set calib flags to False

# Look for these IRAF keywords, if they exist assume that the image has 
# been calibrated.
if "ZEROCOR" in keywordlist:
	biascor = True
if "zach" in keywordlist:
	darkcor = True
if "FLATCOR" in keywordlist:
	flatcor = True

# Do the calibration according to which part is missing. The calibration assumes
# that the calibration frames masterbias, masterdark and masterflat have already
# been created and are in the same directory. Their names are hardwired in the
# code. Zach currently working on creating the the bias, darks and flats automatically.
if not biascor:
	print "Frame is not Bias calibrated"
	print "Proceeding with removing the bias..."
	if os.path.isfile('masterbias.fits'):	
		biasimg = pyfits.open('masterbias.fits')
		masterbias = biasimg[0].data
		dataref = dataref - masterbias
		print "Bias calibration performed"
	else:
		masterbias = makebias(dataref)
		dataref = dataref - masterbias
		print "Bias calibration performed"

if not darkcor:
	print "Frame is not Dark calibrated"
	print "Proceeding with removing the dark..."
	if os.path.isfile('masterdark.fits'):	
		darkimg = pyfits.open('masterdark.fits')
		masterdark = darkimg[0].data
		dataref = dataref - masterdark
		print "Dark calibration performed"
	else:
		masterdark = makedark(dataref)
		dataref = dataref - masterdark
		print "Dark calibration performed"

"""
if not flatcor:
	flatimg = pyfits.open('masterflat.fits')
	flatref = flatimg[0].data
	dataref = dataref / flatref
	print "Uncalibrated image: Removed Flat"
"""
# Write calibrated file to disk
hdu=pyfits.PrimaryHDU(dataref)
hdu.writeto('testcalib.fits')
print "All Done!"



"""
# Create a Masterdark frame from available dark frames.
-------------------------------------------------------------

# Create a Masterflat frame from available flat field frames.
-------------------------------------------------------------
# Initialize flat frames dictionary
flatframes = {}
    
# Search for all flat field files in the current directory
# Warning: will search for all files with the text 'flat' in their filename.
flatfiles = sorted(glob.glob('*flat*'))
print "Found ",len(flatfiles), "flat field frames"

# Determine the size of the image based on the size of the first flat field file.
nx = pyfits.getval(flatfiles[0], 'NAXIS1')
ny = pyfits.getval(flatfiles[0], 'NAXIS2')

# Set up the array that will hold all the flat field images.
flatnum = len(flatfiles)
flat_images = np.ndarray((flatnum,nx,ny),dtype=float)

# Go round this loop and fill the flat field images array.
for i in range(0,flatnum):
    print 'Reading file',i
    flat_images[i] = pyfits.getdata(flatfiles[i])

# Median combine the flat field frames
masterflat_cts = np.median(flat_images, axis=0)

# Find the median value of the masterflat_cts frame and divide by it to creat the
# masterflat image.
median = masterflat_cts.median()
print median
masterflat = masterflat_cts/median

# Output the masterflat frame.
# Filename hardwired as masterflat.fits
hdu=pyfits.PrimaryHDU(masterflat)
hdu.writeto('masterflat.fits')
  
# Free the memory!  
del flat_images
"""

