"""
This routine gets a list of FITS files that belong to a timeseries photometry run
and calculates the pixel offsets of each file in relation to the first file.

The offsets are calculated by cross correlating a blured masked of the reference
frame (usually the first image in a sequence) with a blured mask of the following
image. The blured masks are necesseary to make the stars appear as prominent as
possible without any background noise, hot pixels or cosmic rays.

"""
import os
import pyfits
import matplotlib.pyplot as plt
from scipy import signal, fftpack, ndimage
import numpy

# Ask user for the file containing the list of FITS images to process
filein = raw_input('Enter file with the list of images to be processed:   ')
print
print "Reading filelist from:", filein
print

# Read the list of files
file = open(filein,'r')
filelist = file.read().splitlines()
files = sum(1 for line in open(filein))
i=0
xshift=0
yshift=0

# Open the output file and write the header information and the first image info
fileout=open('offsets.dat','w')
fileout.write("! File created by offsets.py File format is:"+'\n')
fileout.write("! Filename"+'\n')
fileout.write("! Image#, X-Shift, Y-Shift"+'\n')
fileout.write(filelist[i]+'\n')
fileout.write("1 0 0"+'\n')

# Difine arrays to store the offsets
xoff=numpy.arange(0,files,1)
yoff=numpy.arange(0,files,1)

i=0
for n in range(files):

	if i==0:
		# Open the first image and read its data array
		image1 = pyfits.open(filelist[0])
		data1=image1[0].data
		#image1.info()

		# Crop the image by 10 pixels on each side
		xsize1=data1.shape[0]
		ysize1=data1.shape[1]
		xstart=9
		xend=xsize1-10	
		ystart=9
		yend=ysize1-10
		croped1=data1[xstart:xend,ystart:yend]
		median1 = pyfits.np.median(croped1)
		print "Cropped Image shape: ", croped1.shape
		print "Cropped Image Median: ", median1

		# Create image 1 mask and blur it using a Gausian filter of
		# FWHM of 1 pixel. Then any pixel with a value less than 100
		# is set to zero. This is completely arbitrary but usually
		# pixels that correspond to actual stars will have values
		# a lot greater than 100.
		mask1=croped1
		mask1[mask1 < 1.5*median1] = 0.0
		blured1=ndimage.gaussian_filter(croped1, sigma=1)
		mask1=blured1
		mask1[mask1 < 100.0] = 0.0

		# Get the size of the croped masked arrays
		xsize=mask1.shape[0]
		ysize=mask1.shape[1]

		# Create collapsed image arrays for image 1
		xvals1=mask1.sum(axis=0)
		yvals1=mask1.sum(axis=1)
		i=i+1
		xoff[i]=0
		yoff[i]=0

		# 	Write images to file
		#	hdu=pyfits.PrimaryHDU(mask1)
		#	hdu.writeto("mask1.fits")

	else:	

		# For the rest of the files load the data array
		image2 = pyfits.open(filelist[i])
		data2=image2[0].data
		#image2.info()
	
		# Crop the image by 10 pixes on each side
		xsize2=data2.shape[0]
		ysize2=data2.shape[1]
		xstart=9
		xend=xsize2-10
		ystart=9
		yend=ysize2-10
		croped2=data2[xstart:xend,ystart:yend]
		median2 = pyfits.np.median(croped2)

		#	print "Cropped Image shape: ", croped2.shape
		#	print "Cropped Image Median: ", median2
		#	print
	
		# Create image 2 mask and blur it
		mask2=croped2
		mask2[mask2 < 1.5*median2] = 0.0
		blured2=ndimage.gaussian_filter(croped2, sigma=1)
		mask2=blured2
		mask2[mask2 < 100.0] = 0.0

		#	hdu=pyfits.PrimaryHDU(mask2)
		#	hdu.writeto("mask2.fits")

		# Create collapsed image arrays for image 2
		xvals2=mask2.sum(axis=0)
		yvals2=mask2.sum(axis=1)

		# Plot sum of y-values of image 1
		#testx=numpy.arange(0,xsize,1)
		#plt.plot(testx,yvals1)
		#plt.show()

		# Plot sum of y-values of image 2
		#testx=numpy.arange(0,xsize,1)
		#plt.plot(testx,yvals2)
		#plt.show()

		# Plot the cross correlation function for the x-values
		#testxwidth=ysize*2
		#test=signal.correlate(xvals1,xvals2)
		#testx=numpy.arange(0,testxwidth-1,1)
		#plt.plot(testx,test)
		#plt.show()

		# Calculate the x and y shift of the image in pixels using cross correlation
		xshift = (numpy.argmax(signal.correlate(xvals1,xvals2)))-(ysize-1)
		yshift = (numpy.argmax(signal.correlate(yvals1,yvals2)))-(xsize-1)

		print "Shifts for image ", filelist[i]
		print "X-shift is: ", xshift
		print "Y-shift is: ", yshift

		# Store the offsets in arrays
		xoff[i]=xshift
		yoff[i]=yshift

#		Write the output to file
		fileout=open('offsets.dat','a')
		fileout.write(filelist[i]+'\n')
#		print >> [i], xshift, yshift
		fileout.write("%s %s %s" % (i+1, xshift, yshift)+'\n')

		i=i+1

fileout.close()

# Plot the X and Y offsets
print "Plot shows X-offsets"
print
plotx=numpy.arange(0,files,1)
plt.plot(plotx,xoff)
plt.plot(plotx,yoff)
plt.legend(['X-Shifts', 'Y-Shifts'], loc='upper left')
plt.show()

print "Processed ",i,"file names"
print

# Things for later:
# Figure out the x-y size discrepancy for the xshift and yshift calculation. Currently,
# to calculate the xshift i need to subract the ysize of the image rather than the xsize.
# This of course affects only non-square images.

# Properly integrate the whole routine into the reduction pipeline.

# Write the X and Y shifts into an ascii file.

# Make an X Y plot showing the offset of the center of the image.



