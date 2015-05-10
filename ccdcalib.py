#!/usr/bin/python
#

import astropy.io.fits as pyfits
import glob
import numpy as np
import sys
import os
import datetime

# Write a FITS file with updated headers.
def writefits(data, hdr, filename):
	
	hdr_out = hdr.copy(strip=True)	
	hdu_out = pyfits.PrimaryHDU(data, hdr_out)
	hdu_out.scale('float32', 'old')	# For compatibility use 32 bits per data pixel
	hdu_out.writeto(filename, output_verify='warn', clobber=True)
	del data, hdr, hdr_out 			# Free some memory
	return

# Create a Masterbias frame from available bias frames.
#-------------------------------------------------------------
def makebias(dataref):

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
#   		print 'Reading bias file',i
    		bias_images[i] = pyfits.getdata(biasfiles[i])

	# Median combine the bias frame
	masterbias = np.median(bias_images, axis=0)
	print "Bias Median: ", np.median(masterbias)

	# Constract the BIASCOR keyword and add it to the output bias header
	# First, get the header of the first file to use as a header for the output file.
	hdr_bias = pyfits.getheader(biasfiles[0])
	hdr_out = hdr_bias.copy(strip=True)
	# Make a text string with the current date and time
	biastxt = "Bias frame created by RTPhoS on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	# Constuct the keyword
	hdr_out['BIASCOR'] = (biastxt)		

	# Output the masterbias frame.
	# Filename hardwired as masterbias.fits
	writefits(masterbias, hdr_out, 'masterbias.fits')

	# Memory Freedom!  
	del bias_images, biasfiles

	return masterbias

# Create a Masterdark frame from available dark frames.
#-------------------------------------------------------------
def makedark(dataref):
    
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
	writefits(masterdark, 'masterdark.fits')
  
	# Free the memory!  
	del dark_images

	return masterdark

# Create a Masteflat frame from available flat field frames.
#-------------------------------------------------------------
def makeflat(dataref):
    
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
		print 'Reading file',i
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

	# Output the normalized masterflat frame.
	# Filename hardwired as masterflat.fits
	writefits(masterflat, 'masterflat.fits')
  
	# Free the memory!  
	del flat_images

	return masterflat


# --------------- MAIN ---------------------------------------
#ref_filename = win.get("file")
ref_filename = 'J1753.fits'
dataref, hdr_data = pyfits.getdata(ref_filename, header=True)     
print("... Working ...")

# Check if the file has been processed (bias,dark,flat)
# (will use standard IRAF keywords for this)
keywordlist = hdr_data.keys()
biascor = False  # set calibration flags
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
		print "Bias calibration performed!"
	else:
		masterbias = makebias(dataref)
		dataref = dataref - masterbias
		print "Bias calibration performed!"

if not darkcor:
	print "Frame is not Dark calibrated"
	print "Proceeding with removing the dark..."
	if os.path.isfile('masterdark.fits'):	
		darkimg = pyfits.open('masterdark.fits')
		masterdark = darkimg[0].data
		dataref = dataref - masterdark
		print "Dark calibration performed!"
	else:
		masterdark = makedark(dataref)
		dataref = dataref - masterdark
		print "Dark calibration performed!"

if not flatcor:
	print "Frame is not Flat fielded"
	print "Proceeding with Flat fielding..."
	if os.path.isfile('masterflat.fits'):	
		flatimg = pyfits.open('masterflat.fits')
		masterflat = flatimg[0].data
		dataref = dataref/masterflat
		print "Flat fielding performed!"
	else:
		masterflat = makeflat(dataref)
		dataref = dataref/masterflat
		print "Flat fielding performed!"

# Write calibrated file to disk
writefits(dataref, hdr_data, 'testcalib.fits')
print "All Done!"




